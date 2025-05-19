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
#' @param itemmodel Character vector with the model for each item.
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
lca <- function(data, nclasses = 2L, item_model = rep("multinomial", ncol(data)),
                model = NULL, control = NULL, do.fit = TRUE, constraints = TRUE) {

  # library(data.table, quietly = TRUE)
  # library(lavaan, quietly = TRUE)

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
  if(!is.character(item_model) || length(item_model) != ncol(data)) {

    stop("item_model must be a character vector with as many elements as columns in data")

  }

  # Get the full model specification with labels for each parameter:
  # Fix the first parameter of a constrained vector:
  model <- getmodel(data = data, item_model = item_model, nclasses = nclasses,
                    model = model, constraints = constraints)

  # Generate all the necessary objects to evaluate the estimator:
  data_list <- vector("list")
  data_list$Y <- Y
  data_list$n <- n
  data_list$uniq_indices <- uniq_indices
  data_list$map2full <- map2full
  data_list$dt <- dt

  args <- getargs_lca(data_list = data_list, item_model = item_model, model = model,
                      control = control)
  manifold_setup <- args$manifold_setup
  transform_setup <- args$transform_setup
  estimator_setup <- args$estimator_setup

  # Initial parameters estimates and posterior probabilities for the optimizer:
  parameters <- args$parameters # For derivative-based algorithms
  posterior <- args$posterior   # For the EM algorithm
  transparameters <- args$transparameters

  # Give additional information to the optimizer:
  control$parameters <- parameters
  control$transparameters <- transparameters
  control$posterior <- posterior
  control$S <- estimator_setup[[1]]$S # Number of unique patterns
  control$nlatent <- length(estimator_setup[[1]]$classes) # Number of latent variables

  result <- vector("list") # Initialize the object to be returned

  # Fit the model or just get the model specification:
  if(!do.fit) {

    # If do.fit is FALSE, just return the model setup for the optimization:

    opt <- list(manifold_setup = manifold_setup,
                transform_setup = transform_setup,
                estimator_setup = estimator_setup,
                control = control,
                # init_param = parameters,  # Already in control
                # init_post = posterior,    # Already in control
                data = data)

    modelInfo <- list(item_model = item_model,
                      model = model$log_model,
                      nobs = nobs,
                      nparam = length(parameters[[1]]),
                      ntransparam = length(transparameters[[1]]),
                      npatterns = npatterns,
                      df = npatterns - length(parameters[[1]]))

    result$modelInfo <- modelInfo
    result$opt <- opt

    return(result) # Return information without fitting the model

  }

  if(control$opt == "em-lbfgs") {

    # Perform the optimization with EM (fit the model):
    control$opt <- "em"
    x <- lca(data, nclasses = nclasses, item_model = item_model,
             model = NULL, control = control, do.fit = TRUE)

    # Perform the optimization with L-BFGS (fit the model):
    control$opt <- "lbfgs"
    control$parameters <- list(x$opt$parameters) # THESE PARAMETERS SHOULD
    control$posterior <- list(x$opt$posterior)
    control$rstarts <- 1L
    control$cores <- 1L
    x <- lca(data, nclasses = nclasses, item_model = item_model,
             model = model$log_model, control = control, do.fit = TRUE)

  } else {

    # Perform the optimization (fit the model):
    x <- optimizer(control_transform = transform_setup,
                   control_estimator = estimator_setup,
                   control_manifold = manifold_setup,
                   control_optimizer = control)

  }

  if(control$opt != "em") {

    parameters <- x$parameters
    transparameters <- c(x$transparameters)

  } else {

    copy <- model$log_model$conditionals
    k <- 0

    if("gaussian" %in% item_model) {
      k <- k+1
      indices_gauss <- which(item_model == "gaussian")
      copy[indices_gauss] <- x$list_matrices[[k]][[1]]
    }

    if("multinomial" %in% item_model) {
      k <- k+1
      indices_mult <- which(item_model == "multinomial")
      copy[indices_mult] <- x$list_matrices[[k]][[1]]
    }

    parameters <- NULL
    transparameters <- c(x$vectors[[1]][[1]], unlist(copy))

  }

  # Collect all the information about the optimization:
  opt <- list(parameters = parameters,
              transparameters = transparameters,
              posterior = x$posterior,
              manifold_setup = manifold_setup,
              transform_setup = transform_setup,
              estimator_setup = estimator_setup,
              control = control,
              init_param = parameters,
              init_transparam = transparameters,
              init_post = posterior,
              iterations = x$iterations,
              convergence = x$convergence,
              ng = x$ng,
              data = data,
              outputs = x$outputs)

  loglik <- -x$f # Logarithm likelihood
  logliks <- x$outputs$estimators$vectors[[1]][[3]] # Loglik of response pattern
  loglik_case <- n * logliks # Sum of logliks by response pattern
  posterior <- x$posterior
  colnames(posterior) <- paste("P(", "Class", 1:nclasses, "|Y)", sep = "")
  state <- apply(posterior, MARGIN = 1, FUN = which.max)

  # Create a summary table with information for each response pattern:
  estimated <- round(exp(logliks) * sum(n), 0)
  state <- apply(x$posterior, MARGIN = 1, FUN = which.max)
  Freqs <- cbind(Pattern = Y + 1,
                 Observed = n,
                 Estimated = estimated,
                 Posterior = posterior,
                 State = state,
                 loglik = logliks,
                 total_loglik = loglik_case)
  Freqs <- as.data.frame(Freqs)
  Freqs <- Freqs[do.call(order, as.data.frame(Freqs)), ]
  rownames(Freqs) <- paste("pattern", 1:nrow(Freqs), sep = "")

  vec <- vector(length = args$ntransparam)
  if(!is.null(parameters)) {
    vec[args$indices_full_param_vector] <- parameters
    vec[args$indices_full_fixed_vector] <- args$full_fixed_vector
    result$parameters <- fill_list_with_vector(model$log_model, vec)
    result$parameters <- convert_all_to_numeric(result$parameters)
  }
  # result$parameters <- replace_near_zero(result$parameters)
  result$transformed_parameters <- fill_list_with_vector(model$prob_model, x$transparameters)
  result$transformed_parameters <- convert_all_to_numeric(result$transformed_parameters)
  # result$transformed_parameters <- replace_near_zero(result$transformed_parameters)

  result$posterior <- posterior[map2full, ]
  result$state <- state[map2full]
  result$loglik <- loglik
  result$loglik_case <- logliks[map2full]
  result$summary_table <- Freqs

  rownames(result$posterior) <- rownames(data)

  gauss <- "gaussian" %in% item_model     # Check the existence of gaussian items
  multin <- "multinomial" %in% item_model # Check the existence of multinomial items

  # Outputs only for multinomial models:
  if(multin & !gauss) {

    classes <- transparameters[1:nclasses]

    # For multinomial items:
    indices <- which(item_model == "multinomial")
    # Initialize the outputs:
    nitems <- length(indices) # Number of items in the data
    ClassConditional <- RespConditional <- vector("list", length = nitems)
    names(ClassConditional) <- names(RespConditional) <- paste("Item", 1:nitems)
    multinomial_conditionals <- lapply(x$outputs$estimators$list_matrices[[1]][[1]], FUN = \(x) {
      # Extract P(y|X):
      ncategories <- length(x)/nclasses
      result <- matrix(x, nrow = ncategories, ncol = nclasses)
      colnames(result) <- paste("Class", 1:nclasses, sep = "")
      rownames(result) <- paste("Category", 1:ncategories, sep = "")
      return(result)
    })
    names(multinomial_conditionals) <- paste("Item", 1:nitems)

    ClassConditional[indices] <- multinomial_conditionals

    probCat <- lapply(multinomial_conditionals, FUN = \(mat) {
      # Calculate P(y|X)*P(X), the joint probability:
      jointp <- t(mat) * classes
      # Calculate P(y), the denominator of the posterior:
      probCat <- colSums(jointp)
      return(probCat)
    })


    RespConditional[indices] <- lapply(multinomial_conditionals, FUN = \(mat) {
      # Calculate P(y|X)*P(X), the joint probability:
      jointp <- t(mat) * classes
      # Calculate P(y), the denominator of the posterior:
      probCat <- colSums(jointp)
      # Calculate P(X|y) = P(y|X)*P(X)/P(y), the posterior:
      posterior <- t(jointp) / probCat
      return(posterior)
    })

    result$ClassConditional <- multinomial_conditionals
    result$RespConditional <- RespConditional[indices]
    names(result$RespConditional) <- paste("Item", 1:nitems)
    result$probCat <- probCat

  }

  nparam <- length(parameters)
  modelInfo <- list(item_model = item_model,
                    model = model$log_model,
                    prob_model = model$prob_model,
                    nobs = nobs,
                    nparam = nparam,
                    npatterns = npatterns,
                    df = npatterns - nparam)
  result$modelInfo <- modelInfo

  result$opt <- opt
  result$elapsed <- x$elapsed

  class(result) <- "lca"

  return(result)

}


# Andres Chull, fit indices
# input output lavaan
# missing data CFA EFA LCA
# 1:nclasses
# Estimated values

