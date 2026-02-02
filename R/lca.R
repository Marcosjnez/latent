# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 12/12/2025
#'
#' @title
#' Latent Class Analysis.
#' @description
#'
#' Estimate latent class models with gaussian and multinomial item models,
#' with or without covariates to predict class membership.
#'
#' @usage
#'
#' lca(data, nclasses = 2L, item = rep("gaussian", ncol(data)),
#'     X = NULL, penalties = TRUE, model = NULL, mimic = "LG",
#'     start = NULL, do.fit = TRUE, verbose = TRUE, control = NULL)
#'
#' @param data data frame or matrix.
#' @param nclasses Number of latent classes.
#' @param item Character vector with the model for each item (i.e., "gaussian" or "multinomial"). Defaults to "gaussian" for all the items.
#' @param X Matrix of covariates.
#' @param penalties Boolean or list of penalty terms for the parameters.
#' @param model List of parameter labels. See 'details' for more information.
#' @param mimic String. Replicate the output of other softwares. Use "LG" to replicate the output of LatentGOLD.
#' @param start List of starting values for the parameters. See 'details' for more information.
#' @param do.fit TRUE to fit the model and FALSE to return only the model setup. Defaults to TRUE.
#' @param control List of control parameters for the optimization algorithm. See 'details' for more information.
#' @param verbose Print information of model estimation. Defaults to FALSE.
#'
#' @details \code{lca} estimates models with categorical and continuous data.
#'
#' @return List with the following objects:
#' \item{version}{Version number of 'latent' when the model was estimated.}
#' \item{call}{Code used to estimate the model.}
#' \item{ModelInfo}{Model information.}
#' \item{Optim}{Output of the optimizer.}
#' \item{user_model}{Structure with most relevant parameters.}
#' \item{parameters}{Structure with all model parameters.}
#' \item{transparameters}{Structure with all transformed model parameters.}
#' \item{posterior}{Posterior probability of class membership.}
#' \item{state}{Class with highest posterior probability.}
#' \item{loglik}{Logarithm likelihood of the model.}
#' \item{penalized_loglik}{Logarithm likelihood + logarithm priors of the model.}
#' \item{loglik_case}{Logarithm likelihood of each pattern.}
#' \item{summary_table}{Table of summaries of the fitted model. Useful when all the items are 'multinomial'.}
#' \item{ClassConditional}{Parameters that are conditional on the class memberships.}
#' \item{ResponseConditional}{Probability of class memberships conditional on item response. Only available when all the items have a 'multinomial' likelihood.}
#' \item{probCat}{Marginal probability of an item response. Only available when all the items have a 'multinomial' likelihood.}
#'
#' @examples
#'
#' \dontrun{
#'
#' fit <- lca(data = gss82, nclasses = 3L,
#'            item = rep("multinomial", ncol(gss82)),
#'            penalties = TRUE, do.fit = TRUE)
#' fit@timing
#' fit@loglik # -2754.643
#' fit@penalized_loglik # -2759.507
#' fit@Optim$opt$iterations
#'
#' # Plot model fit info:
#' fit
#'
#' # Get fit indices:
#' getfit(fit)
#'
#' # Get a summary:
#' summary(fit)
#'
#' # Inspect model objects:
#' latInspect(fit, what = "coefs", digits = 3)
#' latInspect(fit, what = "classes", digits = 3)
#' latInspect(fit, what = "profile", digits = 3)
#' latInspect(fit, what = "posterior", digits = 3)
#' latInspect(fit, what = "table", digits = 3)
#' latInspect(fit, what = "pattern", digits = 3)
#'
#' # Get standard errors:
#' SE <- se(fit, type = "standard", digits = 4)
#' SE$table
#'
#' # Get confidence intervals:
#' CI <- ci(fit, type = "standard", confidence = 0.95, digits = 2)
#' CI$table
#' }
#'
#' @references
#'
#' None yet.
#'
#' @export
lca <- function(data, nclasses = 2L, item = rep("gaussian", ncol(data)),
                X = NULL, penalties = TRUE, model = NULL, mimic = "LG",
                start = NULL, do.fit = TRUE, control = NULL, verbose = TRUE) {

  # Fix the measurement part of the model if a previous llca fit is available:
  # This is usually done for two-step estimators of covariate coefficients.
  if(class(model) == "llca") {

    model_list <- model@transformed_pars
    model_list$beta <- NULL

    # Fit the model with covariates fixing the measurement part:
    fit <- lca(data, nclasses = nclasses, item = item,
               X = X, penalties = penalties, model = model_list,
               start = start, do.fit = do.fit, control = control,
               verbose = verbose)

    return(fit)

  }

  # Check control parameters:
  control$penalties <- penalties
  control$start <- start
  control <- lca_control(control)

  #### Initial input checks ####

  # Number of subjects:
  nobs <- nrow(data)
  # Number of items:
  nitems <- ncol(data)

  # Check that data is either a data.frame or a matrix:
  if(!is.data.frame(data) & !is.matrix(data)) {
    stop("data must be a matrix or data.frame")
  }

  # Check that item is a character vector with a string for each column of data:
  if(!is.character(item) || length(item) != nitems) {

    stop("item must be a character vector with as many elements as columns in data")

  }

  # Process the covariates:
  if(is.null(X)) {

    X <- matrix(1, nrow = nobs, ncol = 1L)
    colnames(X) <- "(Intercept)"

  } else {

    # Check that X is either a data.frame or a matrix:
    if(!is.data.frame(X) & !is.matrix(X)) {
      stop("data must be a matrix or data.frame")
    }

    if(nrow(X) != nobs) {
      stop("Number of cases in the data and covariates do not match")
    }

    # Transform characters into factors:
    X_df <- as.data.frame(X)
    X_df[] <- lapply(X_df, function(x) if (is.character(x)) factor(x) else x)

    # Create the design matrix:
    X <- model.matrix(~ . + 1, X_df)
    # Center the variables:
    # X[, -1] <- apply(X[, -1], MARGIN = 2, FUN = \(x) x-mean(x))

    # Put an underscore between the variable names and their level names:
    for (v in names(X_df)[sapply(X_df, is.factor)]) {
      i <- which(startsWith(colnames(X), v))
      colnames(X)[i] <- paste0(v, "_", make.names(levels(X_df[[v]]))[seq_along(i)])
    }

  }

  #### Process the data ####

  # Check if any column in "data" is modeled with the multinomial model and
  # transform it to a factor:

  condition <- item == "multinomial"
  factor_indices <- which(condition)
  # Transform into factors:
  data[, factor_indices] <- lapply(data[, factor_indices, drop = FALSE],
                                   FUN = factor)
  # Save category names for data modeled with the multinomial model:
  factor_names <- lapply(data[, factor_indices, drop = FALSE], levels)
  data[, factor_indices] <- lapply(data[, factor_indices, drop = FALSE],
                                   FUN = function(col) {
                                     if (is.factor(col)) as.integer(col) - 1L else col
                                   })

  # Get the number of possible response patterns:
  if(all(condition)) { # If all the items are multinomial...

    # Count the number of categories:
    Ks <- apply(data[, factor_indices, drop = FALSE], MARGIN = 2,
                FUN = \(x) length(unique(x[!is.na(x)])))
    # If all the items are multinomial, the number of possible patterns is:
    npossible_patterns <- min(prod(Ks)-1, nobs)
    # npossible_patterns <- min(prod(Ks), npatterns)

  } else { # If any item is not multinomial...

    npossible_patterns <- nobs

  }

  # Combine the data and the covariates:
  dt <- cbind(X, data)
  pcov <- ncol(X) # Number of covariates

  # Convert data to a data.table object:
  dt <- data.table::as.data.table(dt)

  ## Collect some information from the data ##
  counts_dt <- dt[, .(index = .I[1], count = .N), by = names(dt)]
  # Data matrix with the unique response patterns:
  patterns <- as.matrix(counts_dt[, names(dt), with = FALSE])[, -(1:pcov) , drop = FALSE]
  # Covariates with unique response patterns:
  cov_patterns2 <- as.matrix(counts_dt[, names(dt), with = FALSE])[, 1:pcov , drop = FALSE]
  # Number of unique response patterns:
  npatterns <- nrow(counts_dt)
  # Counts of each response pattern:
  weights <- counts_dt$count
  # Indices to map the original data to the matrix of unique patterns:
  full2short <- counts_dt$index
  # Indices to map the matrix of unique patterns to the original data:
  short2full <- match(do.call(paste, dt),
                      do.call(paste, counts_dt[, -c("index", "count"), with = FALSE]))

  ## Collect some information from the covariates ##
  cov_dt <- data.table::as.data.table(X)
  counts_cov_dt <- cov_dt[, .(index = .I[1], count = .N), by = names(cov_dt)]
  # Covariates with unique patterns:
  cov_patterns <- as.matrix(counts_cov_dt[, names(cov_dt), with = FALSE])
  # Number of unique covariate patterns:
  ncov_patterns <- nrow(cov_patterns)
  # Indices to map the original data to the matrix of unique patterns:
  cov_full2short <- counts_cov_dt$index
  # Indices to map the matrix of unique patterns to the original data:
  cov_short2full <- match(do.call(paste, cov_dt),
                          do.call(paste, counts_cov_dt[, -c("index", "count"),
                                                       with = FALSE]))

  # Put in a list the objects generated form the data:
  data_list <- vector("list")
  data_list$data <- data
  data_list$X <- X
  data_list$item <- item
  data_list$nobs <- nobs
  data_list$patterns <- patterns
  data_list$cov_patterns2 <- cov_patterns2
  data_list$npatterns <- npatterns
  data_list$npossible_patterns <- npossible_patterns
  data_list$nitems <- nitems
  data_list$weights <- weights
  data_list$full2short <- full2short
  data_list$short2full <- short2full
  data_list$item_names <- colnames(data)
  data_list$factor_indices <- factor_indices
  data_list$factor_names <- factor_names
  data_list$cov_patterns <- cov_patterns
  data_list$ncov_patterns <- ncov_patterns
  data_list$cov_full2short <- cov_full2short
  data_list$cov_short2full <- cov_short2full

  #### Initialize objects to store all the models ####

  model0 <- model
  NCLASSES <- nclasses
  nmodels <- length(NCLASSES)

  llca_list <- list()

  # Loop to create and fit the models for different nclasses:
  for(NK in 1:nmodels) {

    # NK <- 1L
    if(verbose) {

      print(paste0("Model nclasses = ", NCLASSES[NK]) )

    }

    # Ensure that nclasses is an integer:
    nclasses <- as.integer(NCLASSES[NK])
    # Check that nclasses is positive:
    if(nclasses < 1L) {
      stop("nclasses must be a positive integer")
    }

    ## store original call
    mc  <- match.call()

    #### Create the model ####

    # Get the model specification:
    full_model <- get_full_lca_model(data_list = data_list, nclasses = nclasses,
                                     item = item, model = model,
                                     control = control)
    list2env(full_model, envir = environment())

    # Get the short model specification (in logarithm and probability scale) with
    # labels for each parameter:
    short_model <- get_short_lca_model(data_list = data_list, nclasses = nclasses,
                                       item = item, lca_all = lca_all,
                                       model = model)
    list2env(short_model, envir = environment())

    #### Create the structures ####

    # Generate the structures for optimization:
    structures <- get_lca_structures(data_list = data_list,
                                     full_model = full_model,
                                     control = control)
    list2env(structures, envir = environment())

    #### Fit the model ####

    result <- vector("list") # Initialize the object to be returned

    # Model information:
    modelInfo <- list(nclasses = nclasses,
                      item = item,
                      nobs = nobs,
                      npatterns = npatterns,
                      npossible_patterns = npossible_patterns,
                      nparam = nparam,
                      dof = npossible_patterns - nparam,
                      ntrans = ntrans,
                      model = log_model,
                      prob_model = prob_model,
                      parameters_labels = parameters_labels,
                      transparameters_labels = transparameters_labels,
                      lca_param = lca_param,
                      lca_trans = lca_trans,
                      lca_all = lca_all,
                      control_manifold = control_manifold,
                      control_transform = control_transform,
                      control_estimator = control_estimator,
                      control = control)

    # Fit the model or just get the model specification:
    if(!do.fit) {

      llca_list[[NK]] <- new("llca",
                             version            = as.character( packageVersion('latent') ),
                             call               = mc, # matched call
                             timing             = numeric(), # timing information
                             data_list          = data_list,
                             modelInfo          = modelInfo,
                             Optim              = list(),
                             user_model         = list(),
                             parameters         = list(),
                             transformed_pars   = list(),
                             posterior          = matrix(),
                             state              = vector(),
                             loglik             = numeric(), # loglik values
                             penalized_loglik   = numeric(),
                             loglik_case        = numeric(),
                             summary_table      = data.frame(),
                             ClassConditional   = list(),
                             RespConditional    = list(),
                             probCat            = list()
      )

      next # Go to the next model with a different value of "nclasses"

    }

    if(control$opt == "em-lbfgs") {

      # Run first EM and then LBFGS:

      # Backup for LBFGS:
      maxit_lbfgs <- control$maxit
      eps_lbfgs <- control$eps

      # Setup for EM:
      control$opt <- "em"
      control$maxit <- control$maxit_em
      control$eps <- control$em_eps
      control$cores <- min(control$rstarts, control$cores)

      # Perform EM:
      xEM <- optimizer(control_manifold = control_manifold,
                       control_transform = control_transform,
                       control_estimator = control_estimator,
                       control_optimizer = control)

      # Reset everything for LBFGS:
      # Put as many rstarts in LBFGS as number of selected iterations in EM:
      rstarts_lbfgs <- control$pick
      control$rstarts <- rstarts_lbfgs
      control$pick <- 0L # Set to 0L to get a full output next time
      control$opt <- "lbfgs"
      control$maxit <- maxit_lbfgs # Reset number of iteration for LBFGS
      control$eps <- eps_lbfgs # Reset tolerance value for LBFGS

      # Update the initial parameters values for LBFGS:
      control$parameters <- control$transparameters <-
        vector("list", length = rstarts_lbfgs)
      for(i in 1:rstarts_lbfgs) {
        control$parameters[[i]] <- xEM[[i]]$parameters
        control$transparameters[[i]] <- xEM[[i]]$transparameters
      }

    }

    control$cores <- min(control$rstarts, control$cores)
    # Perform the optimization (fit the model):
    x <- optimizer(control_manifold = control_manifold,
                   control_transform = control_transform,
                   control_estimator = control_estimator,
                   control_optimizer = control)

    # Collect all the information about the optimization:

    Optim <- x

    #### Process the outputs ####

    # Logarithm likelihood:
    loglik <- x$outputs$estimators$doubles[[1]]
    penalized_loglik <- -x$f
    # Logarithm likelihood of each response pattern:
    loglik_case <- x$outputs$estimators$vectors[[1]][[1]]
    # Sum of logarithm likelihoods by response pattern:
    loglik_pattern <- weights * loglik_case

    # Create the transformed parameters:
    transformed_pars <- fill_list_with_vector(modelInfo$lca_all,
                                              Optim$transparameters)
    transformed_pars <- allnumeric(transformed_pars)

    ## Summary table with information for each response pattern ##
    # Estimated "counts" for each response pattern:
    estimated <- exp(loglik_case) * nobs
    # Posterior:
    posterior <- exp(matrix(x$outputs$estimators$matrices[[1]][[2]],
                            nrow = npatterns, ncol = nclasses))
    colnames(posterior) <- paste("P(", "Class", 1:nclasses, "|Y)", sep = "")
    # Posterior classification:
    state <- apply(posterior, MARGIN = 1, FUN = which.max)
    # Data table of response patterns:
    summary_table <- cbind(Pattern = patterns + 1,
                           Observed = weights,
                           Estimated = estimated,
                           Posterior = posterior,
                           State = state,
                           loglik_case = loglik_case,
                           loglik_pattern = loglik_pattern)
    summary_table <- as.data.frame(summary_table)
    # Sort the patterns by increasing order:
    summary_table <- summary_table[do.call(order, as.data.frame(summary_table)), ]
    rownames(summary_table) <- paste("pattern", 1:nrow(summary_table), sep = "")

    # Check the existence of gaussian items:
    gauss <- "gaussian" %in% item
    # Check the existence of multinomial items:
    multin <- "multinomial" %in% item

    # Create the parameter table:
    selection <- match(unlist(prob_model), transparameters_labels)
    out <- x$transparameters[selection]
    user_model <- fill_list_with_vector(prob_model, out)
    user_model <- allnumeric(user_model)

    selection <- match(unlist(log_model), transparameters_labels)
    out <- x$transparameters[selection]
    raw_model <- fill_list_with_vector(log_model, out)
    raw_model <- allnumeric(raw_model)

    ClassConditional <- user_model$items
    RespConditional <- probCat <- list() # Only for full multinomial models

    # Additional outputs for full multinomial models:
    if(all(item == "multinomial")) {

      classes <- colMeans(transformed_pars$class)
      conditionals <- user_model$items

      probCat <- lapply(conditionals, FUN = \(mat) {
        # Calculate P(y|X)*P(X), the joint probability:
        jointp <- t(mat) * classes
        # Calculate P(y), the denominator of the posterior:
        probCat <- colSums(jointp)
        return(probCat)
      })

      RespConditional <- lapply(conditionals, FUN = \(mat) {
        # Calculate P(y|X)*P(X), the joint probability:
        jointp <- t(mat) * classes
        # Calculate P(y), the denominator of the posterior:
        probCat <- colSums(jointp)
        # Calculate P(X|y) = P(y|X)*P(X)/P(y), the posterior:
        posterior <- t(jointp) / probCat
        return(posterior)
      })

      names(RespConditional) <- colnames(data)

    }

    #### Return ####

    loglik_case <- loglik_case[short2full]
    posterior <- posterior[short2full, , drop = FALSE]
    state <- state[short2full]
    rownames(posterior) <- names(state) <-
      names(loglik_case) <- rownames(data)

    if(control$opt == "em-lbfgs") {
      elapsed <- xEM$elapsed + x$elapsed
    } else {
      elapsed <- x$elapsed
    }

    llca_list[[NK]] <- new("llca",
                           version            = as.character( packageVersion('latent') ),
                           call               = mc, # matched call
                           timing             = elapsed, # timing information
                           data_list          = data_list,
                           modelInfo          = modelInfo,
                           Optim              = Optim,
                           user_model         = user_model,
                           parameters         = raw_model,
                           transformed_pars   = transformed_pars,
                           posterior          = posterior,
                           state              = state,
                           loglik             = loglik,
                           penalized_loglik   = penalized_loglik,
                           loglik_case        = loglik_case,
                           summary_table      = summary_table,
                           ClassConditional   = ClassConditional,
                           RespConditional    = RespConditional,
                           probCat            = probCat
    )


  }

  #### Return ####

  class(llca_list) <- "llcalist"

  if(nmodels == 1) {
    return(llca_list[[1]])
  } else{
    return(llca_list)
  }

}

