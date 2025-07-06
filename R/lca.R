# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 22/05/2025
#'
#' @title
#' Latent Class Analysis.
#' @description
#'
#' Estimate latent class models with gaussian and multinomial item models.
#'
#' @usage
#'
#' lca(data, item = rep("gaussian", ncol(data)), nclasses = 2L,
#' model = NULL, control = list(opt = "lbfgs", rstarts = 30L, cores = 1L),
#' do.fit = TRUE, constraints = TRUE)
#'
#' @param data data frame or matrix.
#' @param nclasses Number of latent classes.
#' @param item Character vector with the model for each item (i.e., "gaussian" or "multinomial"). Defaults to "gaussian" for all the items.
#' @param model List of parameter labels. See 'details' for more information.
#' @param control List of control parameters for the optimization algorithm. See 'details' for more information.
#' @param do.fit TRUE to fit the model and FALSE to return only the model setup. Defaults to TRUE.
#' @param constraints Should the model be checked for identification? Defaults to TRUE.
#'
#' @details \code{lca} estimates models with categorical and continuous data.
#'
#' @return List with the following objects:
#' \item{parameters}{The model for the logarithm probabilities of the classes.}
#' \item{f}{Logarithm likelihood at the maximum.}
#'
#' @references
#'
#' None yet.
#'
#' @export
lca <- function(data, nclasses = 2L, item = rep("gaussian", ncol(data)),
                model = NULL, control = NULL, do.fit = TRUE, constraints = TRUE) {

  # Check that data is either a data.frame or a matrix:
  if(!is.data.frame(data) & !is.matrix(data)) {
    stop("data must be a matrix or data.frame")
  }

  nclasses <- as.integer(nclasses) # Ensure that nclasses is an integer
  # Check that nclasses is a positive integer:
  if(nclasses < 1L) {
    stop("nclasses must be a positive integer")
  }

  ## store original call
  mc  <- match.call()

  # Create defaults for the control of the optimizer:
  control <- lca_control(control)

  ## Process the data ##
  # Sample size:
  nobs <- nrow(data)
  # Convert data to a data.table object:
  dt <- data.table::as.data.table(data)
  ## Collect some information from the data ##
  counts_dt <- dt[, .(index = .I[1], count = .N), by = names(dt)]
  # Data matrix with the unique response patterns:
  Y <- as.matrix(counts_dt[, names(dt), with = FALSE])
  # Number of unique response patterns:
  npatterns <- nrow(counts_dt)
  # Counts of each response pattern:
  n <- counts_dt$count
  # Indices to map the original data to the matrix of unique patterns:
  uniq_indices <- counts_dt$index
  # Indices to map the matrix of unique patterns to the original data:
  map2full <- match(do.call(paste, dt), do.call(paste, counts_dt[, -c("index", "count"), with = FALSE]))

  # Check that the itemmodel is a character vector with a string for each column of data:
  if(!is.character(item) || length(item) != ncol(data)) {

    stop("item must be a character vector with as many elements as columns in data")

  }

  # Get the full model specification (in logarithm and probability scale) with
  # labels for each parameter:
  model <- getmodel(data = data, item = item, nclasses = nclasses,
                    model = model, constraints = constraints)

  # Put in a list the objects generated form the data:
  data_list <- vector("list")
  data_list$dt <- dt
  data_list$Y <- Y
  data_list$n <- n
  data_list$uniq_indices <- uniq_indices
  data_list$map2full <- map2full

  # Generate all the necessary objects to evaluate the estimators:
  args <- getargs_lca(data_list = data_list, item = item,
                      model = model, control = control)
  manifold_setup <- args$manifold_setup
  transform_setup <- args$transform_setup
  estimator_setup <- args$estimator_setup
  control_setup <- args$control_setup

  # Extract the initial parameters estimates and posterior probabilities for the
  # optimizer:
  parameters <- args$parameters # For derivative-based algorithms
  posterior <- args$posterior   # For the EM algorithm
  transparameters <- args$transparameters
  nparam <- length(parameters[[1]])
  ntransparam <- length(transparameters[[1]])

  #result <- vector("list") # Initialize the object to be returned

  # Fit the model or just get the model specification:
  if(!do.fit) {

    # If do.fit is FALSE, just return the model setup:

    # Model information:
    modelInfo <- list(item = item,
                      model = model$log_model,
                      prob_model = model$prob_model,
                      nobs = nobs,
                      nparam = nparam,
                      ntransparam = ntransparam,
                      npatterns = npatterns,
                      df = npatterns - nparam)

    #result$modelInfo <- modelInfo

    # Data for the optimization algorithms:
    opt <- list(data = data,
                data_list = data_list,
                manifold_setup = manifold_setup,
                transform_setup = transform_setup,
                estimator_setup = estimator_setup,
                control_setup = control_setup,
                args = args)

    #result$opt <- opt
    llca <- new("llca",
                   version            = as.character( packageVersion('latent') ),
                   call               = mc, # matched call
                   timing             = NULL, # timing information
                   modelInfo          = modelInfo, # modelInfo
                   Optim              = opt, # opt
                   parameters         = NULL,
                   transformed_pars   = NULL,
                   posterior          = NULL,
                   state              = NULL,
                   loglik             = NULL, # loglik values and info
                   loglik_case        = NULL,
                   summary_table      = NULL,
                   ClassConditional   = NULL,
                   RespConditional    = NULL,
                   probCat            = NULL
                   )

    return(llca) # Return information without fitting the model

  }

  # Perform the optimization (fit the model):
  x <- optimizer(control_manifold = manifold_setup,
                 control_transform = transform_setup,
                 control_estimator = estimator_setup,
                 control_optimizer = control_setup)

  if(control_setup$opt != "em") {

    parameters <- c(x$parameters)
    names(parameters) <- args$parameter_vector
    transparameters <- c(x$transparameters)
    names(transparameters) <- args$transparameter_vector

  } else {

    # IF the EM algorithm was used...

    copy <- model$log_model$conditionals
    k <- 0

    if("gaussian" %in% item) {
      k <- k+1
      indices_gauss <- which(item == "gaussian")
      copy[indices_gauss] <- x$outputs$estimators$list_matrices[[k]][[1]]
    }

    if("multinomial" %in% item) {
      k <- k+1
      indices_mult <- which(item == "multinomial")
      copy[indices_mult] <- x$outputs$estimators$list_matrices[[k]][[1]]
    }

    parameters <- NULL
    transparameters <- c(x$outputs$estimators$vectors[[1]][[1]], unlist(copy))

  }

  # Collect all the information about the optimization:
  opt <- list(parameters = parameters,
              transparameters = transparameters,
              posterior = x$posterior,
              iterations = x$iterations,
              convergence = x$convergence,
              ng = x$ng,
              data = data,
              data_list = data_list,
              manifold_setup = manifold_setup,
              transform_setup = transform_setup,
              estimator_setup = estimator_setup,
              control_setup = control_setup,
              outputs = x$outputs,
              args = args)

  ## Collect information from the fitted model ##
  # Logarithm likelihood:
  loglik <- -x$f
  # Logarithm likelihood of each response pattern:
  logliks <- x$outputs$estimators$vectors[[1]][[3]]
  # Sum of logarithm likelihoods by response pattern:
  loglik_case <- n * logliks
  # Posterior probabilities:
  posterior <- x$posterior
  colnames(posterior) <- paste("P(", "Class", 1:nclasses, "|Y)", sep = "")

  ## Summary table with information for each response pattern ##
  # Estimated counts for each response pattern:
  estimated <- round(exp(logliks) * sum(n), 0)
  # Posterior classification:
  state <- apply(posterior, MARGIN = 1, FUN = which.max)
  # Data table of response patterns:
  Freqs <- cbind(Pattern = Y + 1,
                 Observed = n,
                 Estimated = estimated,
                 Posterior = posterior,
                 State = state,
                 loglik = logliks,
                 total_loglik = loglik_case)
  Freqs <- as.data.frame(Freqs)
  # Sort the patterns by increasing order:
  Freqs <- Freqs[do.call(order, as.data.frame(Freqs)), ]
  rownames(Freqs) <- paste("pattern", 1:nrow(Freqs), sep = "")

  # Allocate the parameters (logarithm scale) in the full model specification:
  ncomplete <- length(unlist(model$log_model))
  vec <- vector(length = ncomplete)
  if(!is.null(parameters)) {
    vec[args$indices_full_param_vector] <- parameters[args$indices_param_vector2[!is.na(args$indices_param_vector2)]]
    vec[args$indices_full_fixed_vector] <- args$full_fixed_vector
    parameters <- fill_list_with_vector(model$log_model, vec)
    parameters <- convert_all_to_numeric(parameters)
  }

  # Allocate the transformed parameters (probability scale) in the full model
  # specification:
  transformed_parameters <- vector(length = ncomplete)
  transformed_parameters[args$indices_full_transfixed_vector] <- args$full_transfixed_vector
  transformed_parameters[args$indices_full_transparam_vector] <- x$transparameters[args$indices_transparam_vector2[!is.na(args$indices_transparam_vector2)]]
  transformed_parameters <- fill_list_with_vector(model$prob_model, transformed_parameters)
  transformed_parameters <- convert_all_to_numeric(transformed_parameters)

  # Fill some results:
  posterior <- posterior[map2full, ]
  rownames(posterior) <- rownames(data)
  state <- state[map2full]
  loglik <- loglik
  loglik_case <- logliks[map2full]
  summary_table <- Freqs

  # Check the existence of gaussian items:
  gauss <- "gaussian" %in% item
  # Check the existence of multinomial items:
  multin <- "multinomial" %in% item

  # Additional outputs only for multinomial models:
  if(multin & !gauss) {

    classes <- x$outputs$estimators$vectors[[2]][[1]]

    # For multinomial items:
    indices <- which(item == "multinomial")
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

    ClassConditional <- multinomial_conditionals
    RespConditional <- RespConditional[indices]
    names(RespConditional) <- colnames(data)#paste("Item", 1:nitems)
    #probCat <- probCat

  }

  modelInfo <- list(item = item,
                    model = model$log_model,
                    prob_model = model$prob_model,
                    nobs = nobs,
                    nparam = nparam,
                    ntransparam = ntransparam,
                    npatterns = npatterns,
                    df = npatterns - nparam)

  #result$modelInfo <- modelInfo

  #result$opt <- opt
  #result$elapsed <- x$elapsed

  #class(result) <- "lca"

  llca <- new("llca",
                 version            = as.character( packageVersion('latent') ),
                 call               = mc, # matched call
                 timing             = x$elapsed, # timing information
                 modelInfo          = modelInfo, # modelInfo
                 Optim              = opt, # opt
                 parameters         = parameters,
                 transformed_pars   = transformed_parameters,
                 posterior          = list(posterior),
                 state              = list(state),
                 loglik             = loglik, # loglik values and info
                 loglik_case        = loglik_case,
                 summary_table      = summary_table,
                 ClassConditional   = ClassConditional,
                 RespConditional    = RespConditional,
                 probCat            = probCat
  )

  return(llca)

}


# Andres Chull, fit indices
# input output lavaan
# missing data CFA EFA LCA
# 1:nclasses
# Estimated values

# parameters package R
