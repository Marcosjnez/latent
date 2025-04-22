#' @title
#' Latent Class Analysis.
#' @description
#'
#' Estimate latent class models with categorical and continuous data.
#'
#' @usage
#'
#' lca(data, model = rep("multinomial", ncol(data)), nclasses = 2L,
#' control = list(opt = "lbfgs", rstarts = 30L, cores = 1L))
#'
#' @param data data.frame or matrix of response.
#' @param model Character vector with the model for each item.
#' @param nclasses Number of latent classes.
#' @param control List of control parameters for the optimization algorithm. See 'details' for more information.
#'
#' @details \code{lca} estimates models with categorical and continuous data.
#'
#' @return List with the following objects:
#' \item{parameters}{The model for the probabilities of the classes.}
#' \item{f}{Logarithm likelihood at the maximum.}
#'
#' @references
#'
#' None yet.
#'
#' @export
lca <- function(data, nclasses = 2L, model = rep("multinomial", ncol(data)),
                fullmodel = NULL, control = NULL, do.fit = TRUE) {

  # Check that data is either a data.frame or a matrix:
  if(!is.data.frame(data) & !is.matrix(data)) {
    stop("data must be a matrix or data.frame")
  }

  nclasses <- as.integer(nclasses) # Ensure that nclasses is an integer
  # Check that nclasses is a positive integer:
  if(nclasses < 1L) {
    stop("nclasses must be a positive integer")
  }

  # Create defaults for the control of the optimizer:
  control <- lca_control(control)

  # Process the data:
  dt <- data.table::as.data.table(data) # Convert data to a data.table object
  # Create a matrix with the unique response patterns:
  counts_dt <- dt[, .(index = .I[1], count = .N), by = names(dt)]
  Y <- as.matrix(counts_dt[, names(dt), with = FALSE])
  npatterns <- nrow(counts_dt) # Number of unique response patterns
  n <- counts_dt$count # Number of repetitions of each response pattern
  # Get the position of the unique patterns in the original data:
  uniq_indices <- counts_dt$index
  # Indices to map the matrix of unique patterns to the original data:
  map2full <- match(do.call(paste, dt), do.call(paste, counts_dt[, -c("index", "count"), with = FALSE]))

  nobs <- nrow(data) # Sample size

  # Check that the model is a character vector with a string for each column of data:
  if(!is.character(model) || length(model) != ncol(data)) {

    stop("model must be a character vector with as many elements as columns in data")

  }

  if(is.null(fullmodel)) {

    # Get the full model specification with labels for each parameter:
    # Fix the first parameter of a constrained vector:
    free <- FALSE
    if(control$opt == "em") free <- TRUE # Free all the parameters if EM is used
    fullmodel <- getmodel(data = data, model = model, nclasses = nclasses, free = free)

  }

  classes <- fullmodel$classes           # Get the model for the classes
  conditionals <- fullmodel$conditionals # Get the model for the item parameters

  opt <- control$opt                     # Optimizer

  # Generate all the necessary objects to evaluate the estimator:
  data_list <- vector("list")
  data_list$Y <- Y
  data_list$n <- n
  data_list$uniq_indices <- uniq_indices
  data_list$map2full <- map2full

  args <- getargs_lca(data_list = data_list, model = model, classes = classes,
                      conditionals = conditionals, control = control)
  estimator_setup <- args$estimator_setup

  # Generate all the necessary objects to project onto the manifold:
  # All the parameters are fixed to belong to the Euclidean manifold:
  indices <- unique(unlist(lapply(estimator_setup, function(x) x$indices)))
  manifolds <- c("euclidean")
  arguments <- list()
  arguments[[1]] <- list(indices = indices)
  manifold_setup <- setup_all_manifolds(manifolds, arguments)

  # Initial parameters estimates and posterior probabilities for the optimizer:
  parameters <- args$parameters # For derivative-based algorithms
  posterior <- args$posterior   # For the EM algorithm

  # Give additional information to the optimizer:
  control$parameters <- parameters
  control$posterior <- posterior
  control$S <- estimator_setup[[1]]$S # Number of unique patterns
  control$nlatent <- length(estimator_setup[[1]]$classes) # Number of latent variables

  result <- vector("list") # Initialize the object to be returned

  # Fit the model or just get the model specification:
  if(!do.fit) {

    # If do.fit is FALSE, just return the model setup for the optimization:

    opt <- list(estimator_setup = estimator_setup,
                manifold_setup = manifold_setup,
                control = control,
                # init_param = parameters,  # Already in control
                # init_post = posterior,    # Already in control
                data = data)

    modelInfo <- list(model_vector = model,
                      model = fullmodel,
                      nobs = nobs,
                      nparam = length(parameters),
                      npatterns = npatterns,
                      df = npatterns - length(parameters))

    result$modelInfo <- modelInfo
    result$opt <- opt

    return(result) # Return information without fitting the model

  }

  if(control$opt == "em-lbfgs") {

    # Perform the optimization with EM (fit the model):
    control$opt <- "em"
    x <- lca(data, nclasses = nclasses, model = model,
             fullmodel = NULL, control = control, do.fit = TRUE)

    # Perform the optimization with L-BFGS (fit the model):
    control$opt <- "lbfgs"
    control$parameters <- list(x$opt$parameters) # THESE PARAMETERS SHOULD
    control$posterior <- list(x$opt$posterior)
    control$rstarts <- 1L
    control$cores <- 1L
    # which(fullmodel)
    x <- lca(data, nclasses = nclasses, model = model,
             fullmodel = fullmodel, control = control, do.fit = TRUE)

  } else {

    # Perform the optimization (fit the model):
    x <- optimizer(control_estimator = estimator_setup,
                   control_manifold = manifold_setup,
                   control_optimizer = control)

  }

  # Collect all the information about the optimization:
  opt <- list(parameters = x$parameters,
              posterior = x$posterior,
              estimator_setup = estimator_setup,
              manifold_setup = manifold_setup,
              control = control,
              init_param = parameters,
              init_post = posterior,
              iterations = x$iterations,
              convergence = x$convergence,
              ng = x$ng,
              data = data)

  loglik <- -x$f # Logarithm likelihood
  classes <- x$vectors[[1]][[1]] # Class probabilities
  names(classes) <- paste("Class", 1:nclasses, sep = "")

  nitems <- ncol(data) # Number of items in the data

  # Initialize the outputs:
  ClassConditional <- RespConditional <- vector("list", length = nitems)
  names(ClassConditional) <- names(RespConditional) <- paste("Item", 1:nitems)
  type <- list()

  gauss <- "gaussian" %in% model     # Check the existence of gaussian items
  multin <- "multinomial" %in% model # Check the existence of multinomial items
  rm(posterior) # Remove the posterior from the environment REMOVE?

  # Fill the objects that were initialized depending on the nature of the items:
  k <- 1

  if(gauss & multin) {

    # If there are both gaussian and multinomial items...

    # S1 <- estimator_setup[[1]]$S # Number of patterns in gaussian data
    # S2 <- estimator_setup[[2]]$S # Number of patterns in multinomial data
    # gauss_map2full <- estimator_setup[[1]]$map2full
    # categ_map2full <- estimator_setup[[2]]$map2full

    # Log likelihood of gaussian data conditioning on the latent parameters:
    gauss_latentloglik <- matrix(x$matrices[[1]][[2]], nrow = npatterns, ncol = nclasses)
    # Log likelihood of multinomial data conditioning on the latent parameters:
    categ_latentloglik <- matrix(x$matrices[[2]][[2]], nrow = npatterns, ncol = nclasses)

    # gauss_latentloglik <- gauss_latentloglik[gauss_map2full, ]
    # categ_latentloglik <- categ_latentloglik[categ_map2full, ]

    # Log likelihood of data conditioning on the latent parameters:
    latentloglik <- gauss_latentloglik + categ_latentloglik

    # Joint lkelihood of data and latent parameters:
    jointp <- t(exp(t(latentloglik) + log(classes)))
    # Marginal likelihood of latent parameters:
    marginalp <- rowSums(jointp)
    # Posterior probabilities of the latent parameters:
    posterior <- jointp / marginalp
    # Logarithm likelihood by pattern:
    loglik_case <- log(marginalp)

    # Most probable class for each pattern:
    state <- apply(posterior, MARGIN = 1, which.max)

    colnames(posterior) <- paste("Class", 1:nclasses, sep = "")

  }

  if(gauss) {

    type$gaussian <- list()
    # For gaussian items:
    indices <- which(model == "gaussian")
    gaussian_conditionals <- lapply(x$list_matrices[[1]][[1]], FUN = \(x) {
      # P(y|X):
      result <- matrix(x, nrow = 2, ncol = nclasses)
      if(control$opt == "em") {
        result[2, ] <- sqrt(result[2, ]) # Variance to Standard deviation
      } else {
        result[2, ] <- exp(result[2, ]) # Standard deviation
      }
      colnames(result) <- paste("Class", 1:nclasses, sep = "")
      rownames(result) <- c("Mean", "Std")
      return(result)
    })
    ClassConditional[indices] <- gaussian_conditionals

    uniq_indices <- estimator_setup[[1]]$uniq_indices
    observed <- estimator_setup[[1]]$n
    pattern <- estimator_setup[[1]]$Y
    logliks <- -x$vectors[[1]][[2]]
    liks <- exp(logliks)
    estimated <- round(liks/sum(liks) * sum(observed), 2)
    map2full <- estimator_setup[[1]]$map2full
    if(!exists("posterior")) {
      post <-  matrix(unlist(x$matrices[[1]][1]), nrow = nrow(pattern), ncol = nclasses)
      colnames(post) <- paste("Class", 1:nclasses, sep = "")
      stat <- apply(post, MARGIN = 1, which.max)
      posterior <- post
      state <- stat
      loglik_case <- logliks[map2full]
    } else {
      post <- posterior[uniq_indices, ]
      stat <- state[uniq_indices]
    }
    Freqs <- cbind(pattern,
                   observed = observed,
                   estimated = estimated,
                   posterior = post,
                   state = stat,
                   loglik = logliks)
    Freqs <- as.data.frame(Freqs)
    type$gaussian$loglik_case <- logliks[map2full]
    Freqs <- Freqs[do.call(order, as.data.frame(Freqs)), ]
    rownames(Freqs) <- paste("pattern", 1:nrow(Freqs), sep = "")
    type$gaussian$Pattern <- Freqs
    k <- k+1

  }

  if(multin) {

    type$multinomial <- list()
    # For multinomial items:
    indices <- which(model == "multinomial")
    multinomial_conditionals <- lapply(x$list_matrices[[k]][[1]], FUN = \(x) {
      # Extract P(y|X):
      ncategories <- length(x)/nclasses
      result <- matrix(x, nrow = ncategories, ncol = nclasses)
      colnames(result) <- paste("Class", 1:nclasses, sep = "")
      rownames(result) <- paste("Category", 1:ncategories, sep = "")
      return(result)
    })
    nitemscat <- length(multinomial_conditionals)
    names(multinomial_conditionals) <- paste("Item", 1:nitemscat)
    type$multinomial$ClassConditional <- multinomial_conditionals

    ClassConditional[indices] <- multinomial_conditionals

    probCat <- lapply(multinomial_conditionals, FUN = \(mat) {
      # Calculate P(y|X)*P(X), the joint probability:
      jointp <- t(mat) * classes
      # Calculate P(y), the denominator of the posterior:
      probCat <- colSums(jointp)
      return(probCat)
    })

    type$multinomial$probCat <- probCat

    RespConditional[indices] <- lapply(multinomial_conditionals, FUN = \(mat) {
      # Calculate P(y|X)*P(X), the joint probability:
      jointp <- t(mat) * classes
      # Calculate P(y), the denominator of the posterior:
      probCat <- colSums(jointp)
      # Calculate P(X|y) = P(y|X)*P(X)/P(y), the posterior:
      posterior <- t(jointp) / probCat
      return(posterior)
    })

    type$multinomial$RespConditional <- RespConditional[indices]
    names(type$multinomial$RespConditional) <- paste("Item", 1:nitemscat)

    uniq_indices <- estimator_setup[[k]]$uniq_indices
    observed <- estimator_setup[[k]]$n
    pattern <- estimator_setup[[k]]$Y + 1
    logliks <- -x$vectors[[k]][[2]]
    estimated <- round(exp(logliks) * sum(observed), 0)
    map2full <- estimator_setup[[k]]$map2full
    if(!exists("posterior")) {
      post <-  matrix(unlist(x$matrices[[k]][1]), nrow = nrow(pattern), ncol = nclasses)
      colnames(post) <- paste("Class", 1:nclasses, sep = "")
      stat <- apply(post, MARGIN = 1, which.max)
      posterior <- post
      state <- stat
      loglik_case <- logliks[map2full]
    } else {
      post <- posterior[uniq_indices, ]
      stat <- state[uniq_indices]
    }
    Freqs <- cbind(pattern,
                   observed = observed,
                   estimated = estimated,
                   posterior = post,
                   state = stat,
                   loglik = logliks)
    Freqs <- as.data.frame(Freqs)
    # type$multinomial$posterior <- post[map2full, ]
    type$multinomial$loglik_case <- logliks[map2full]
    Freqs <- Freqs[do.call(order, as.data.frame(Freqs)), ]
    rownames(Freqs) <- paste("pattern", 1:nrow(Freqs), sep = "")
    type$multinomial$Pattern <- Freqs

  }

  result$parameters$classes <- classes
  result$parameters$items <- ClassConditional
  result$parameters <- replace_near_zero(result$parameters)
  # result$classes <- classes
  # result$ClassConditional <- ClassConditional
  result$posterior <- posterior
  result$state <- state
  result$loglik <- loglik
  result$loglik_case <- loglik_case
  result$details <- type

  modelInfo <- list(model_vector = model,
                    model = fullmodel,
                    nobs = nobs,
                    nparam = length(x$parameters),
                    npatterns = npatterns,
                    df = npatterns - length(x$parameters))
  result$modelInfo <- modelInfo

  result$opt <- opt
  result$elapsed <- x$elapsed
  # table(x$classification)/nrow(data)
  # x$classes
  # table(dat$trueclass)/nrow(data)

  class(result) <- "lca"

  return(result)

}


# Andres Chull, fit indices
# input output lavaan
# missing data CFA EFA LCA
# 1:nc

## Estimated values

