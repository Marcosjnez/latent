# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 22/05/2026
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
#'     X = NULL, penalties = TRUE, model = NULL,
#'     start = NULL, do.fit = TRUE, verbose = TRUE, control = NULL)
#'
#' @param data data frame or matrix.
#' @param nclasses Number of latent classes.
#' @param item Character vector with the model for each item (i.e., "gaussian" or "multinomial"). Defaults to "gaussian" for all the items.
#' @param X Matrix of covariates.
#' @param penalties Boolean or list of penalty terms for the parameters.
#' @param model List of parameter labels. See 'details' for more information.
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
lca <- function(data,
                nclasses = 1L,
                gaussian = NULL,
                multinomial = NULL,
                X = NULL,
                penalties = TRUE,
                model = NULL,
                start = NULL,
                control = NULL,
                do.fit = TRUE,
                verbose = TRUE) {

  ## Store original call:
  mc  <- match.call()
  args <- as.list(match.call(expand.dots = TRUE))[-1]

  # Check that data is either a data.frame or matrix:
  if(!is.data.frame(data) & !is.matrix(data)) {
    stop("data must be a matrix or data.frame")
  }

  # Check that nclasses is a (vector of) positive integer(s):
  if(any(nclasses < 1L)) {
    stop("nclasses must be a (vector of) positive integer(s)")
  }

  # Class enumeration: run a model for different number of classes
  # Returns a list of llca models:
  if(length(nclasses) > 1L) {

    result <- vector("list", length = length(nclasses))
    for(i in 1:length(nclasses)) {

      if(verbose) print(paste0("Model nclasses=", nclasses[i]))

      result[[i]] <- lca(data, nclasses = nclasses[i],
                         gaussian = gaussian, multinomial = multinomial,
                         X = X, penalties = penalties, model = model,
                         start = start, control = control,
                         do.fit = do.fit, verbose = verbose)

    }

    class(result) <- "llcalist"

    return(result)

  }

  #### Create the dataList ####

  dataList_and_control <- create_lca_dataList(data = data,
                                              X = X,
                                              gaussian = gaussian,
                                              multinomial = multinomial,
                                              model = model,
                                              penalties = penalties,
                                              control = control,
                                              start = start)
  list2env(dataList_and_control, envir = environment())

  #### Create the model ####

  # Fix the measurement part of the model if a previous llca fit was provided
  # model. This is usually done for two-step estimation:
  if(class(model) == "llca") {

    model <- model@parameters
    model$beta <- NULL

  }

  # Get the model specification:
  full_model <- create_lca_model(dataList = dataList,
                                 nclasses = nclasses,
                                 item = item,
                                 model = model,
                                 control = control)
  list2env(full_model, envir = environment())

  #### Create the structures ####

  # Generate the structures for optimization:
  modelInfo <- create_lca_modelInfo(dataList = dataList,
                                    full_model = full_model,
                                    control = control)

  #### Fit the model ####

  # Fit the model or just get the model specification:
  if(!do.fit) {

    result <- new("llca",
                  version            = as.character(packageVersion('latent')),
                  call               = mc,
                  timing             = numeric(),
                  dataList           = dataList,
                  modelInfo          = modelInfo,
                  Optim              = list(),
                  parameters         = list(),
                  transformed_pars   = list(),
                  extra              = list()
    )

    return(result)

  }

  modelInfo$control_optimizer$cores <- min(modelInfo$control_optimizer$rstarts,
                                           modelInfo$control_optimizer$cores)
  # Perform the optimization (fit the model):
  Optim <- optimizer(control_manifold = modelInfo$control_manifold,
                     control_transform = modelInfo$control_transform,
                     control_estimator = modelInfo$control_estimator,
                     control_optimizer = modelInfo$control_optimizer)

  # Collect all the information about the optimization:

  names(Optim$parameters) <- modelInfo$parameters_labels
  names(Optim$transparameters) <- modelInfo$transparameters_labels

  #### Process the outputs ####

  # Create the transformed parameters:
  transformed_pars <- fill_in(modelInfo$trans,
                              Optim$transparameters)

  parameters <- transformed_pars[names(modelInfo$param)]

  #### Process the fit information ####

  loss <- sum(unlist(lapply(Optim$outputs$estimators$doubles,
                            FUN = \(x) x[[1]])))
  loglik <- sum(unlist(lapply(Optim$outputs$estimators$doubles,
                              FUN = \(x) x[[4]])))
  penalty <- sum(unlist(lapply(Optim$outputs$estimators$doubles,
                               FUN = \(x) x[[7]])))
  penalized_loss <- loss + penalty
  penalized_loglik <- loglik - penalty

  #### latent object ####

  elapsed <- Optim$elapsed

  result <- new("llca",
                version            = as.character(packageVersion('latent')),
                call               = mc,
                timing             = elapsed,
                dataList           = dataList,
                modelInfo          = modelInfo,
                Optim              = Optim,
                parameters         = parameters,
                transformed_pars   = transformed_pars,
                extra              = list()
  )


  #### Return ####

  return(result)

}

create_lca_dataList <- function(data,
                                X = NULL,
                                gaussian = NULL,
                                multinomial = NULL,
                                model = NULL,
                                control = list(),
                                penalties = FALSE,
                                start = NULL) {

  # Matrix of predictors for latent class probabilities:
  if(is.character(X)) {
    X <- data[, X, drop = FALSE]
  }

  # Indicators:

  # If at least two items are gaussian, check which of them are involved in
  # covariance structures in the model argument so they are modeled with a
  # multivariate normal distribution. Create also a target matrix where TRUE
  # means that a given covariance is estimated and FALSE means that it is fixed
  # to zero. The diagonal of target should be FALSE because variances are
  # not estimated directly but parameterized as exp(log variances):
  mvgaussian <- NULL
  if (length(gaussian) >= 2L) {

    error_covs <- NULL

    if (!is.null(model)) {
      error_covs <- extract_cov_pairs(model)
    }

    if (!is.null(error_covs)) {

      error_covs$lhs <- as.character(error_covs$lhs)
      error_covs$rhs <- as.character(error_covs$rhs)

      # Keep only covariance pairs where both variables are gaussian:
      keep <- error_covs$lhs %in% gaussian &
        error_covs$rhs %in% gaussian &
        error_covs$lhs != error_covs$rhs

      error_covs <- error_covs[keep, , drop = FALSE]

      if (nrow(error_covs) > 0L) {

        candidate_mvgaussian <- unique(c(error_covs$lhs, error_covs$rhs))

        if (length(candidate_mvgaussian) >= 2L) {

          mvgaussian <- candidate_mvgaussian
          nmvgaussian <- length(mvgaussian)

          target <- matrix(
            FALSE,
            nrow = nmvgaussian,
            ncol = nmvgaussian,
            dimnames = list(mvgaussian, mvgaussian)
          )

          idx <- cbind(
            match(error_covs$lhs, mvgaussian),
            match(error_covs$rhs, mvgaussian)
          )

          target[idx] <- TRUE
          target[idx[, 2:1, drop = FALSE]] <- TRUE
          control$joint_gaussian$target <- target

          # Remove the variables in error_covs from gaussian:
          gaussian <- gaussian[!(gaussian %in% mvgaussian)]
          # These removed variables will be only present in mvgaussian

        }
      }
    }
  }

  # Count the number of items for each likelihood model:
  ngaussian <- length(gaussian)
  nmvgaussian <- length(mvgaussian)
  nmultinomial <- length(multinomial)

  # Select the subset of data with the variables that are used in the
  # measurement model (keeping the original ordering):
  model_labels <- rep(c("gaussian", "mvgaussian", "multinomial"),
                      times = c(ngaussian, nmvgaussian, nmultinomial))
  model_vector <- c(gaussian, mvgaussian, multinomial)
  idx <- match(colnames(data), model_vector); idx <- idx[!is.na(idx)]
  measurement <- data[, model_vector[idx]]

  # match each item with a likelihood model:
  item <- model_labels[idx]

  #### More input checks ####

  # Number of subjects:
  nobs <- nrow(measurement)
  # Number of items:
  nitems <- ncol(measurement)

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
  pX <- ncol(X) # Number of columns of the design matrix

  #### Process the data ####

  # Check if any column in "data" is modeled with the multinomial model and
  # transform it to a factor:
  condition <- item == "multinomial"
  factor_indices <- which(condition)
  # Transform multinomial variables into factors:
  measurement[, factor_indices] <- lapply(measurement[, factor_indices,
                                                      drop = FALSE],
                                          FUN = factor)
  # measurement contains the measurement variables

  # Save levels of factor variables:
  factor_names <- lapply(measurement[, factor_indices, drop = FALSE], levels)
  measurement_recoded <- measurement

  # Transform the factor variables into integers:
  measurement_recoded[, factor_indices] <- lapply(measurement[, factor_indices,
                                                       drop = FALSE],
                                           FUN = function(col) {
                                             if (is.factor(col)) as.integer(col) - 1L else col
                                           })
  # measurement_recoded contains the measurement variables wuth the factor
  # variables transformed into integers

  # Get the number of possible response patterns:
  if(all(condition)) { # If all the items are multinomial...

    # Count the number of categories:
    Ks <- apply(measurement[, factor_indices, drop = FALSE], MARGIN = 2,
                FUN = \(x) length(unique(x[!is.na(x)])))
    # If all the items are multinomial, the number of possible patterns is:
    npossible_patterns <- min(prod(Ks)-1L, nobs)
    # From technical manual of LatentGOLD 6.1 (p.68)

  } else { # If any item is not multinomial...

    npossible_patterns <- nobs

  }

  # Append the data and the covariates. This is necessary to find the unique
  # data patterns:
  dt <- data.table::as.data.table(cbind(X, measurement_recoded))
  # dt <- data.table::as.data.table(cbind(X, measurement))
  # Collect some information from the measurment + covariates:
  counts_dt <- dt[, .(index = .I[1], count = .N), by = names(dt)]
  # Measurement data with unique response patterns:
  patterns <- as.matrix(counts_dt[, names(dt), with = FALSE])[, -(1:pX) ,
                                                              drop = FALSE]
  # Covariate data with unique patterns:
  cov_patterns <- as.matrix(counts_dt[, names(dt)[seq_len(pX)], with = FALSE])

  # Number of unique response patterns:
  npatterns <- nrow(counts_dt)
  # Counts of each response pattern:
  weights <- counts_dt$count
  # Indices to map the full data to the unique patterns:
  full2short <- counts_dt$index
  # Indices to map unique patterns to the full data:
  short2full <- match(do.call(paste, dt),
                      do.call(paste, counts_dt[, -c("index", "count"),
                                               with = FALSE]))
  dt <- data.table::as.data.table(cbind(X, measurement))
  counts_dt <- dt[, .(index = .I[1], count = .N), by = names(dt)]
  # patterns_original <- as.matrix(counts_dt[, names(dt), with = FALSE])[, -(1:pX) ,
  #                                                                      drop = FALSE]
  patterns_original <- as.data.frame(counts_dt[, names(dt), with = FALSE])[, -(1:pX) ,
                                                                           drop = FALSE]

  #### Return ####

  # Put in a list the objects generated form the data:
  dataList <- vector("list")
  dataList$data <- data
  dataList$measurement <- measurement
  dataList$X <- X
  dataList$measurement_recoded <- measurement_recoded
  dataList$nitems <- nitems
  dataList$item <- item
  dataList$item_names <- colnames(measurement)
  dataList$nobs <- nobs
  dataList$patterns <- patterns
  dataList$cov_patterns <- cov_patterns
  dataList$npatterns <- npatterns
  dataList$npossible_patterns <- npossible_patterns
  dataList$weights <- weights
  dataList$full2short <- full2short
  dataList$short2full <- short2full
  dataList$factor_indices <- factor_indices
  dataList$factor_names <- factor_names
  dataList$args <- args
  dataList$patterns_original <- patterns_original

  # Check control parameters:
  control$penalties <- penalties # Either a logical or a list of named penalties
  control$start <- start # Optional named list with initial values
  control <- lca_control(control) # Check and update the control inputs and create defaults

  return(list(dataList = dataList, control = control))

}

create_lca_model <- function(dataList, nclasses, item,
                             model = NULL, control) {

  # Generate the model syntax and initial parameter values

  list2env(dataList, envir = environment())

  class_names <- paste("Class", 1:nclasses, sep = "")
  pattern_names <- paste("pattern", 1:npatterns, sep = "")

  #### Model for the transformed parameters ####

  p <- ncol(X) # Number of predictors
  pred_names <- colnames(X) # Names of predictors

  list_struct <- vector("list")
  k <- 1L

  # Model for the betas:

  list_struct[[k]] <- list(name = "beta",
                           type = "matrix",
                           dim = c(p, nclasses),
                           rownames = pred_names,
                           colnames = class_names)
  k <- k+1L

  list_struct[[k]] <- list(name = "theta",
                           type = "matrix",
                           dim = c(npatterns, nclasses),
                           rownames = pattern_names,
                           colnames = class_names)
  k <- k+1L

  list_struct[[k]] <- list(name = "class",
                           type = "matrix",
                           dim = c(npatterns, nclasses),
                           rownames = pattern_names,
                           colnames = class_names)
  k <- k+1L

  # list_struct[[k]] <- list(name = "loglik",
  #                          type = "array",
  #                          dim = c(npatterns, nitems, nclasses),
  #                          dimnames = list(pattern_names,
  #                                          item_names,
  #                                          class_names))

  list_struct[[k]] <- list(name = "loglik",
                           type = "matrix",
                           dim = c(npatterns, nclasses),
                           dimnames = list(pattern_names,
                                           class_names))
  k <- k+1L

  # Model for gaussian items:
  if(any(item == "gaussian")) {

    gauss <- which(item == "gaussian")
    Jgauss <- length(gauss) # Number of gaussian items

    for(j in 1:Jgauss) {

      list_struct[[k]] <- list(name = item_names[gauss[j]],
                               type = "matrix",
                               dim = c(3, nclasses),
                               rownames = c("mean",
                                            "stdv",
                                            "log(stdv)"),
                               colnames = class_names)
      k <- k+1L

    }

  }

  # Model for multinomial items:
  if(any(item == "multinomial")) {

    multinom <- which(item == "multinomial")
    Jmulti <- length(multinom) # Number of multinomial items

    K <- unlist(lapply(factor_names, FUN = length))
    if(any(K > 30)) {
      stop("You cannot use a multinomial likelihood to model a variable with more than 30 unique categories")
    }

    for(j in 1:Jmulti) {

      log_probs <- paste("log", item_names[multinom[j]], sep = "_")
      list_struct[[k]] <- list(name = log_probs,
                               type = "matrix",
                               dim = c(K[j], nclasses),
                               rownames = factor_names[[j]],
                               colnames = class_names)
      k <- k+1L

      list_struct[[k]] <- list(name = item_names[multinom[j]],
                               type = "matrix",
                               dim = c(K[j], nclasses),
                               rownames = factor_names[[j]],
                               colnames = class_names)
      k <- k+1L

    }

    # Define the joint probability parameters:
    error_covs <- NULL

    if (!is.null(model)) {
      error_covs <- extract_cov_pairs(model)
    }

    if (!is.null(error_covs)) {

      npairs <- nrow(error_covs)
      logpair_names <- pair_names <- pair_levels <- vector("list", length = npairs)
      logjoint_names <- joint_names <- vector("list", length = npairs)
      loginter_names <- vector("list", length = npairs)
      joint_vars <- matrix(NA, nrow = nrow(patterns), ncol = npairs)

      for(h in 1:npairs) {

        pair_names[[h]] <- c(error_covs[h, 1], error_covs[h, 2])
        levels1 <- factor_names[[pair_names[[h]][1]]]
        levels2 <- factor_names[[pair_names[[h]][2]]]
        logpair_names[[h]] <- paste("log_", pair_names[[h]], sep = "")
        pair_levels[[h]] <- paste(rep(levels1, times = length(levels2)),
                               rep(levels2, each = length(levels1)), sep = ".")
        logjoint_names[[h]] <- paste("log_", pair_names[[h]], collapse = ".",
                                  sep = "")
        joint_names[[h]] <- paste(pair_names[[h]], collapse = ".")

        list_struct[[k]] <- list(name = logjoint_names[[h]],
                                 type = "matrix",
                                 dim = c(length(pair_levels[[h]]), nclasses),
                                 rownames = pair_levels[[h]],
                                 colnames = paste("Class", 1:nclasses, sep = ""))
        k <- k+1L

        list_struct[[k]] <- list(name = joint_names[[h]],
                                 type = "matrix",
                                 dim = c(length(pair_levels[[h]]), nclasses),
                                 rownames = pair_levels[[h]],
                                 colnames = paste("Class", 1:nclasses, sep = ""))
        k <- k+1L

        joint_vars[, h] <- apply(patterns_original[, pair_names[[h]]], MARGIN = 1,
                                 FUN = paste, collapse = ".")


        loginter_names[[h]] <- paste("log_", pair_names[[h]], collapse = "x",
                                     sep = "")
        list_struct[[k]] <- list(name = loginter_names[[h]],
                                 type = "matrix",
                                 dim = c(length(levels1), length(levels2)),
                                 rownames = levels1,
                                 colnames = levels2)
        k <- k+1L

      }

      colnames(joint_vars) <- names(pair_levels) <- unlist(joint_names)
      names(logpair_names) <- unlist(logjoint_names)
      names(pair_names) <- unlist(joint_names)

      control$joint_multinomial <- list(
        npairs = npairs,
        logpair_names = logpair_names,
        pair_names = pair_names,
        pair_levels = pair_levels,
        logjoint_names = logjoint_names,
        joint_names = joint_names,
        joint_vars = joint_vars,
        loginter_names = loginter_names
      )

      remove <- unique(unlist(pair_names))
      multinomial_items <- item_names[multinom]
      keep <- multinomial_items[!(multinomial_items %in% remove)]
      new_probs <- unlist(c(keep, joint_names))

      y <- data.frame(patterns_original, joint_vars)[, new_probs]

      all_pair_levels <- c(factor_names[keep], pair_levels)
      new_K <- unlist(lapply(all_pair_levels, FUN = length))
      for(j in 1:length(new_probs)) {
        y[[new_probs[j]]] <- factor(y[[new_probs[j]]], ordered = TRUE,
                                    levels = all_pair_levels[[j]])
      }

      control$joint_multinomial$patterns <- y
      control$joint_multinomial$all_pair_levels <- all_pair_levels
      control$joint_multinomial$J <- ncol(y)
      control$joint_multinomial$K <- new_K

      y <- lapply(y,
                  FUN = function(col) {
                    if (is.factor(col)) as.integer(col) - 1L else col
                  })
      y <- as.matrix(as.data.frame(y))
      control$joint_multinomial$patterns_recoded <- y

    }

  }

  # Model for multivariate normal items:
  if(any(item == "mvgaussian")) {

    mvgauss <- which(item == "mvgaussian")
    Jmvgauss <- length(mvgauss) # Number of mvgauss items

    # mean_names <- paste("means|Class", 1:nclasses, sep = "")
    sigma_names <- paste("sigma|Class", 1:nclasses, sep = "")

    list_struct[[k]] <- list(name = "means",
                             type = "matrix",
                             dim = c(Jmvgauss, nclasses),
                             rownames = item_names[mvgauss],
                             colnames = paste("Class", 1:nclasses, sep = ""))
    k <- k+1L

    list_struct[[k]] <- list(name = "logsigma",
                             type = "matrix",
                             dim = c(Jmvgauss, nclasses),
                             rownames = item_names[mvgauss],
                             colnames = paste("Class", 1:nclasses, sep = ""))
    k <- k+1L

    # Sigma:
    for(j in 1:nclasses) {

      list_struct[[k]] <- list(name = sigma_names[j],
                               type = "matrix",
                               dim = c(Jmvgauss, Jmvgauss),
                               rownames = item_names[mvgauss],
                               colnames = item_names[mvgauss],
                               symmetric = TRUE)
      k <- k+1L

    }

  }

  # Create the full transparameter structure:
  trans <- create_parameters(list_struct)

  #### Model for the parameters ####

  param <- list()

  # Model for the betas:
  param$beta <- trans$beta
  param$beta[, 1] <- "0"

  # Model for gaussian items:
  if(any(item == "gaussian")) {

    cl <- length(param)
    idx_gauss <- (cl + 1L):(cl + Jgauss)
    param[idx_gauss] <- trans[item_names[gauss]]
    # Fix the standard deviations just to avoid that they are taken
    # as parameters and not as fixed parameters:
    param[idx_gauss] <- lapply(trans[item_names[gauss]],
                               FUN = \(x) {
                                 x[2, ] <- "1"
                                 return(x)
                               })
    names(param)[idx_gauss] <- item_names[gauss]

  }

  # Model for multinomial items:
  if(any(item == "multinomial")) {

    cl <- length(param)
    idx_multinom <- (cl + 1L):(cl + Jmulti)
    item_names_multinom <- paste("log", item_names[multinom], sep = "_")
    param[idx_multinom] <- trans[item_names_multinom]
    names(param)[idx_multinom] <- item_names_multinom

    # Fix the first item category to zero for identification after the
    # softmax transformation:
    param[idx_multinom] <- lapply(param[idx_multinom],
                                  FUN = \(x) {
                                    x[1, ] <- "0"; return(x)
                                  })

    if(!is.null(control$joint_multinomial)) {
      loginter_names <- unlist(control$joint_multinomial$loginter_names)
      param[loginter_names] <- trans[loginter_names]
      npairs <- control$joint_multinomial$npairs
      for(j in 1:npairs) {
        param[[loginter_names[[j]]]][1, ] <- "0"
        param[[loginter_names[[j]]]][, 1] <- "0"
        # param[[loginter_names[[j]]]][1, 1] <- "0"
      }
    }

  }

  # Model for multivariate gaussian items:
  if(any(item == "mvgaussian")) {

    param$means <- trans$means
    param$logsigma <- trans$logsigma
    # param[sigma_name] <- trans[sigma_names]
    ordering <- rownames(trans[[sigma_names[1]]])
    target <- control$joint_gaussian$target
    target <- target[ordering, ordering]
    for(j in 1:nclasses) {
      param[[sigma_names[j]]] <- diag(Jmvgauss)
      param[[sigma_names[j]]][target] <- trans[[sigma_names[j]]][target]
      # diag(param[[sigma_names[j]]]) <- diag(trans[[sigma_names[j]]])
      dimnames(param[[sigma_names[j]]]) <- dimnames(trans[[sigma_names[j]]])
    }

  }

  #### Fixed parameters and equality constraints ####

  # Replace the transformed parameter by custom values, if available:
  if(!is.null(model)) {

    # Replace the parameters by custom values:

    nm <- intersect(names(model), names(param))
    nm <- nm[!vapply(model[nm], is.null, logical(1))]
    param[nm] <- model[nm]

  }

  # Set equality constraints in trans:

  # Paste to trans only the alphabetical and alphanumerical elements in param,
  # excluding the numerical elements:
  for (nm in names(param)) {

    x <- param[[nm]]

    # TRUE if contains at least one alphabetic character
    keep <- grepl("[A-Za-z]", x)

    trans[[nm]][keep] <- x[keep]

  }

  #### Create the initial values for the parameters ####

  # Initial values for the log coefficients (betas):
  # betas are sampled from a normal distribution:
  init_param <- vector("list", length = control$rstarts)

  if(any(item == "gaussian")) {

    init_mean <- init_sd <- init_logsd <- list()
    for(j in 1:Jgauss) {

      init_mean[[j]] <- rep(mean(measurement_recoded[, gauss[j]], na.rm = TRUE),
                            times = nclasses)
      init_sd[[j]] <- rep(sd(measurement_recoded[, gauss[j]], na.rm = TRUE), times = nclasses)
      init_logsd[[j]] <- log(init_sd[[j]])

    }

  }

  pi_hat_list <- list()
  if(any(item == "multinomial")) {

    eta_hat_list <- list()
    for(j in 1:Jmulti) {

      int_vector <- measurement_recoded[, multinom[j]]
      int_vector <- int_vector[!is.na(int_vector)] # Remove missing values
      nsize <- length(int_vector)
      props <- count(int_vector, nsize, K[j]) / nsize
      pi_hat_list[[j]] <- props %*% t(rep(1, nclasses))
      log_props <- log(props)
      eta_hat_list[[j]] <- (log_props-log_props[1]) %*% t(rep(1, nclasses))

    }

    pi_hat_vector <- unlist(pi_hat_list)
    vars <- (1-pi_hat_vector)/(nobs*pi_hat_vector)
    sds <- sqrt(vars)
    Ks <- length(unlist(eta_hat_list))

  }

  if(any(item == "mvgaussian")) {

    init_mean <- matrix(rep(colMeans(measurement_recoded[, mvgauss], na.rm = TRUE),
                            times = nclasses), nrow = Jmvgauss, ncol = nclasses)
    init_logsigma <- matrix(rep(apply(measurement_recoded[, mvgauss], MARGIN = 2, FUN = var,
                                      na.rm = TRUE),
                                times = nclasses), nrow = Jmvgauss, ncol = nclasses)
    init_sigma <- vector("list", length = nclasses)
    for(j in 1:nclasses) {

      init_sigma[[j]] <- cov(measurement_recoded[, mvgauss], use = "pairwise.complete.obs")

    }

  }

  for(i in 1:control$rstarts) {

    init_param[[i]] <- vector("list", length = length(param))
    names(init_param[[i]]) <- names(param)

    # Initial values for betas:
    init_beta <- matrix(rnorm(p*nclasses), nrow = p, ncol = nclasses)
    init_param[[i]][["beta"]] <- init_beta
    init_param[[i]]$beta[, 1] <- 0
    dimnames(init_param[[i]]$beta) <- dimnames(param$beta)

    # Initial values for gaussian items:
    if(any(item == "gaussian")) {

      for(j in 1:Jgauss) {

        rmean <- init_mean[[j]] + rnorm(nclasses,
                                        mean = 0,
                                        sd = init_sd[[j]]/sqrt(nobs))

        init_param[[i]][[item_names[gauss[j]]]] <- rbind(rmean,
                                                         init_logsd[[j]],
                                                         init_sd[[j]])
        dimnames(init_param[[i]][[item_names[gauss[j]]]]) <-
          dimnames(param[[item_names[gauss[j]]]])

      }

    }

    # Initial values for multinomial items:
    if(any(item == "multinomial")) {

      item_names_multinom <- paste("log", item_names[multinom], sep = "_")

      for(j in 1:Jmulti) {

        sds <- (1-pi_hat_list[[j]])/(nobs*pi_hat_list[[j]])

        init_param[[i]][[item_names_multinom[j]]] <- eta_hat_list[[j]] +
          rnorm(length(eta_hat_list[[j]]), 0, sds)
        init_param[[i]][[item_names_multinom[j]]][1, ] <- 0

        dimnames(init_param[[i]][[item_names[item_names_multinom[j]]]]) <-
          dimnames(param[[item_names[item_names_multinom[j]]]])

      }

      if(!is.null(control$joint_multinomial)) {

        npairs <- control$joint_multinomial$npairs
        loginter_names <- control$joint_multinomial$loginter_names
        for(j in 1:npairs) {
          dims <- dim(trans[[loginter_names[[j]]]])
          init_param[[i]][[loginter_names[[j]]]] <- matrix(0, nrow = dims[1],
                                                           ncol = dims[2])
        }

      }

    }

    # Initial values for multivariate gaussian items:
    if(any(item == "mvgaussian")) {

      # TO-DO: ALLOW RANDOM STARTING VALUES

      init_param[[i]][["means"]] <- init_mean
      init_param[[i]][["logsigma"]] <- init_logsigma
      dimnames(init_param[[i]][["means"]]) <-
        dimnames(init_param[[i]][["logsigma"]]) <- dimnames(param[["means"]])

      for(j in 1:nclasses) {

        init_param[[i]][[sigma_names[j]]] <- init_sigma[[j]]
        dimnames(init_param[[i]][[sigma_names[j]]]) <- dimnames(param[[sigma_names[j]]])

      }

    }

  }

  #### Custom initial values ####

  # Replace initial starting values by custom starting values:

  if(!is.null(control$start)) {

    nm <- names(control$start)
    nm <- nm[!vapply(control$start, is.null, logical(1))]

    for (i in seq_len(control$rstarts)) {
      common_nm <- intersect(nm, names(init_param[[i]]))
      for (j in common_nm) {
        init_param[[i]][[j]] <- insert_object(init_param[[i]][[j]],
                                              control$start[[j]])
      }
    }

  }

  #### Return ####

  result <- list(param = param,
                 trans = trans,
                 init_param = init_param,
                 pi_hat_list = pi_hat_list,
                 control = control)

  return(result)

}

create_lca_modelInfo <- function(dataList, full_model, control) {

  # Generate control_manifold, control_transform, and control_estimator

  list2env(dataList, envir = environment())
  list2env(full_model, envir = environment())

  nclasses <- ncol(trans$class)

  #### Manifolds ####

  manifolds <- list(
    list(manifold = "euclidean", parameters = names(param))
  )

  control_manifold <- create_manifolds(manifolds = manifolds,
                                       structures = param)

  #### Transformations ####

  transforms <- list()
  k <- 1L

  # betas to thetas:
  transforms[[k]] <- list(transform = "column_space",
                          parameters_in = "beta",
                          parameters_out = "theta",
                          extra = list(X = cov_patterns))
  k <- k+1L

  # thetas to classes with softmax:
  for(s in 1:npatterns) {

    transforms[[k]] <- list(transform = "softmax",
                            parameters_in = list(trans$theta[s, ]),
                            parameters_out = list(trans$class[s, ]))
    k <- k + 1L

  }

  # Conditional item likelihoods (gaussian):
  if(any(item == "gaussian")) {

    gauss <- which(item == "gaussian")
    Jgauss <- length(gauss) # Number of gaussian items

    means <- c(do.call(rbind, lapply(trans[item_names[gauss]],
                                     FUN = \(x) x[1, ])))
    stdv <- c(do.call(rbind, lapply(trans[item_names[gauss]],
                                    FUN = \(x) x[2, ])))
    log_stdv <- c(do.call(rbind, lapply(trans[item_names[gauss]],
                                        FUN = \(x) x[3, ])))

    transforms[[k]] <- list(transform = "exponential",
                            parameters_in = list(log_stdv),
                            parameters_out = list(stdv))
    k <- k+1L

    y <- as.matrix(patterns[, gauss])
    # transforms[[k]] <- list(transform = "normal",
    #                         parameters_in = list(means, stdv),
    #                         parameters_out = list(trans$loglik[, gauss, ]),
    #                         extra = list(y = y, S = npatterns, J = Jgauss,
    #                                      I = nclasses))

    transforms[[k]] <- list(transform = "normal2",
                            parameters_in = list(means, stdv, trans$loglik),
                            parameters_out = list(trans$loglik),
                            extra = list(y = y, S = npatterns, J = Jgauss,
                                         I = nclasses))
    k <- k+1L

  }

  # Conditional item likelihoods (multinomial):
  if(any(item == "multinomial")) {

    # All of this is redundant in get_full_lca_model:
    multinom <- which(item == "multinomial")
    Jmulti <- length(multinom) # Number of multinomial items
    probs <- item_names[multinom]
    log_probs <- paste("log", item_names[multinom], sep = "_")

    # Softmax transformations:
    for(j in 1:Jmulti) {

      for(i in 1:nclasses) {

        transforms[[k]] <- list(transform = "softmax",
                                parameters_in = list(trans[[log_probs[j]]][, i]),
                                parameters_out = list(trans[[probs[j]]][, i]))
        k <- k+1L

      }

    }

    # Model for the joint probabilities:
    if(!is.null(control$joint_multinomial)) {

      logpair_names <- control$joint_multinomial$logpair_names
      pair_names <- control$joint_multinomial$pair_names
      logjoint_names <- control$joint_multinomial$logjoint_names
      joint_names <- control$joint_multinomial$joint_names
      loginter_names <- control$joint_multinomial$loginter_names

      for(j in 1:control$joint_multinomial$npairs) {

        inter_in <- c(trans[[loginter_names[[j]]]])

        for(i in 1:nclasses) {

          pairs1 <- trans[[logpair_names[[j]][1]]][, i]
          pairs2 <- trans[[logpair_names[[j]][2]]][, i]
          pairs_in1 <- rep(pairs1, times = length(pairs2))
          pairs_in2 <- rep(pairs2, each = length(pairs1))

          pars_out <- trans[[logjoint_names[[j]]]][, i]
          transforms[[k]] <- list(transform = "sum_vectors",
                                  parameters_in = list(pairs_in1,
                                                       pairs_in2,
                                                       inter_in),
                                  parameters_out = list(pars_out))

          k <- k+1L

          pars_in_softmax <- trans[[logjoint_names[[j]]]][, i]
          pars_out_softmax <- trans[[joint_names[[j]]]][, i]
          transforms[[k]] <- list(transform = "softmax",
                                  parameters_in = list(pars_in_softmax),
                                  parameters_out = list(pars_out_softmax))

          k <- k+1L

        }

      }

    }

    if(!is.null(control$joint_multinomial)) {

      y <- control$joint_multinomial$patterns_recoded
      new_K <- control$joint_multinomial$K
      newJ <- control$joint_multinomial$J

      transforms[[k]] <- list(transform = "multinomial2",
                              parameters_in = list(trans[colnames(y)],
                                                   trans$loglik),
                              parameters_out = list(trans$loglik),
                              extra = list(y = y, S = npatterns,
                                           J = newJ, I = nclasses,
                                           K = new_K))

      k <- k+1L

    } else {

      y <- as.matrix(patterns[, multinom])
      K <- unlist(lapply(dataList$factor_names, FUN = length))

      # transforms[[k]] <- list(transform = "multinomial",
      #                         parameters_in = list(trans[item_names[multinom]]),
      #                         parameters_out = list(trans$loglik[, multinom, ]),
      #                         extra = list(y = y, S = npatterns, J = Jmulti,
      #                                      I = nclasses, K = K))

      transforms[[k]] <- list(transform = "multinomial2",
                              parameters_in = list(trans[probs],
                                                   trans$loglik),
                              parameters_out = list(trans$loglik),
                              extra = list(y = y, S = npatterns, J = Jmulti,
                                           I = nclasses, K = K))

      k <- k+1L

    }

  }

  # Conditional item likelihoods (multivariate gaussian):
  if(any(item == "mvgaussian")) {

    mvgauss <- which(item == "mvgaussian")
    Jmvgauss <- length(mvgauss) # Number of multivariate gaussian items

    sigma_names <- paste("sigma|Class", 1:nclasses, sep = "")
    vars <- unname(do.call(c, lapply(trans[sigma_names], FUN = \(x) diag(x))))

    transforms[[k]] <- list(transform = "exponential",
                            parameters_in = list(trans$logsigma),
                            parameters_out = list(vars))
    k <- k+1L

    y <- as.matrix(patterns[, mvgauss])
    transforms[[k]] <- list(transform = "mvnormal2",
                            parameters_in = list(unlist(trans$means),
                                                 unlist(trans[sigma_names]),
                                                 trans$loglik),
                            parameters_out = list(trans$loglik),
                            extra = list(y = y, S = npatterns, J = Jmvgauss,
                                         I = nclasses))
    k <- k+1L

  }

  control_transform <- create_transforms(transforms = transforms,
                                         structures = trans)

  #### Estimators ####

  estimators <- list()
  G <- 1L

  # estimators[[G]] <- list(estimator = "lca",
  #                         parameters = c("class", "loglik"),
  #                         extra = list(S = npatterns,
  #                                      J = nitems,
  #                                      I = nclasses,
  #                                      weights = weights,
  #                                      double_names = "lca"))

  estimators[[G]] <- list(estimator = "lca2",
                          parameters = c("class", "loglik"),
                          extra = list(S = npatterns,
                                       I = nclasses,
                                       weights = weights,
                                       double_names = "lca"))
  G <- G + 1L

  # Choose whether using Bayes constants:
  if(control$reg) {

    # Bayes Constant for class probabilities:
    alpha <- control$penalties$class$alpha
    if(alpha != 0) {

      # Get the indices corresponding to the unique covariate patterns:
      dt_uniq_X <- data.table::as.data.table(cov_patterns)
      counts_X <- dt_uniq_X[, .(index = .I[1], count = .N), by = names(dt_uniq_X)]
      uniques <- counts_X$index
      U <- length(uniques)

      for(i in uniques) {

        estimators[[G]] <- list(estimator = "bayesconst1",
                                parameters = list(trans$class[i, ]),
                                extra = list(K = nclasses,
                                             alpha = alpha,
                                             U = U,
                                             N = nobs,
                                             double_names = paste("P(Class|pattern",
                                                                  i, ")", sep = "")))
        G <- G+1L

      }
    }

    # Bayes Constant for standard deviations:
    alpha <- control$penalties$sd$alpha
    if(any(item == "gaussian") & alpha != 0) {

      Y <- measurement_recoded[, gauss, drop = FALSE]
      # sigma_class <- split(trans$sigma, rep(1:nclasses, each = Jgauss))
      varshat <- apply(Y, MARGIN = 2, FUN = var, na.rm = TRUE)*(nobs-1)/nobs

      for(i in 1:nclasses) {

        stdv_by_class <- unlist(lapply(trans[item_names[gauss]],
                                       FUN = \(x) x[2, i]))
        estimators[[G]] <- list(estimator = "bayesconst3",
                                parameters = list(stdv_by_class),
                                extra = list(K = nclasses,
                                             varshat = varshat,
                                             alpha = alpha,
                                             N = nobs,
                                             double_names = paste("stdv|Class",
                                                                  i, sep = "")))
        G <- G+1L

      }

    }

    # Bayes Constant for error covariance matrices:
    alpha <- control$penalties$Sigma$alpha
    if(any(item == "mvgaussian") & alpha != 0) {

      Y <- measurement_recoded[, mvgauss, drop = FALSE]
      D <- diag(apply(Y, MARGIN = 2, FUN = var, na.rm = TRUE)*(nobs-1)/nobs)

      for(i in 1:nclasses) {

        estimators[[G]] <- list(estimator = "bayesconst4",
                                parameters = list(trans[[sigma_names[i]]]),
                                extra = list(K = nclasses,
                                             alpha = alpha,
                                             D = D,
                                             J = Jmvgauss,
                                             double_names = paste("Sigma|Class",
                                                                  i, sep = "")))
        G <- G+1L

      }

    }

    # Bayes Constant for multinomial probabilities:
    alpha <- control$penalties$prob$alpha
    if(any(item == "multinomial") & alpha != 0) {

      for(j in 1:Jmulti) {
        for(i in 1:nclasses) {

          pihat <- pi_hat_list[[j]][, i]
          probs_by_class <- unlist(lapply(trans[item_names[multinom[j]]],
                                          FUN = \(x) x[, i]))

          estimators[[G]] <- list(estimator = "bayesconst2",
                                  parameters = list(probs_by_class),
                                  extra = list(K = nclasses,
                                               pihat = pihat,
                                               alpha = alpha,
                                               N = nobs,
                                               double_names = paste("P(",
                                                            item_names[multinom[j]],
                                                            "|Class", i, ")",
                                                            sep = "")))
          G <- G+1L

        }
      }

    }

    # Gaussian regularization for coefficients:
    alpha <- control$penalties$beta$alpha
    if(alpha != 0) {

      p <- nrow(trans$beta)-1L
      q <- ncol(trans$beta)
      means <- matrix(0, nrow = p, ncol = q)
      sds <- apply(cov_patterns[, -1], MARGIN = 2, sd, na.rm = TRUE)
      sds <- matrix(sds, nrow = p, ncol = q) / alpha

      estimators[[G]] <- list(estimator = "gaussian_loglik",
                              parameters = list(trans$beta[-1, ]),
                              extra = list(means = means,
                                           sds = sds,
                                           alpha = alpha,
                                           N = nobs,
                                           double_names = "Gaussian_coeffs"))
      G <- G+1L


    }

    # Ridge regularization for coefficients:
    lambda <- control$penalties$beta$lambda
    power <- control$penalties$beta$power
    if(!is.null(lambda) & !is.null(power)) {

      if(lambda != 0 && power != 0) {

        estimators[[G]] <- list(estimator = "ridge",
                                parameters = list(trans$beta[-1, ]),
                                extra = list(lambda = lambda,
                                             power = power,
                                             N = nobs,
                                             double_names = paste("Ridge(", power,
                                                          ")_coeffs",
                                                          sep = "")))
        G <- G+1L

      }
    }

  }

  control_estimator <- create_estimators(estimators = estimators,
                                         structures = trans)

  #### Pass the initial values to vectors ####

  idx_transformed <- unlist(lapply(control_transform,
                                   FUN = \(x) unlist(x$indices_out)+1L))
  inits <- create_init(trans, param, init_param,
                       idx_transformed = idx_transformed, control)

  parameters <- inits$parameters
  parameters_labels <- names(parameters[[1]])
  nparam <- length(parameters_labels)

  transparameters <- inits$transparameters
  transparameters_labels <- names(transparameters[[1]])
  ntrans <- length(transparameters_labels)

  trans2param <- match(parameters_labels, transparameters_labels)

  #### Set up the optimizer ####

  # Create defaults for the control of the optimizer:
  control_optimizer <- control
  control_optimizer$parameters <- parameters
  control_optimizer$transparameters <- transparameters
  control_optimizer$init_param <- init_param
  control_optimizer$transparam2param <- trans2param-1L

  #### Return ####

  modelInfo <- list(param = param,
                    trans = trans,
                    nparam = nparam,
                    ntrans = ntrans,
                    parameters_labels = parameters_labels,
                    transparameters_labels = transparameters_labels,
                    dof = npossible_patterns - nparam,
                    control_manifold = control_manifold,
                    control_transform = control_transform,
                    control_estimator = control_estimator,
                    control_optimizer = control_optimizer)

  return(modelInfo)

}

lca_control <- function(control) {

  # Auxiliary function for lca.R

  # Control input

  if(isFALSE(control$penalties)) {

    control$reg <- FALSE

  } else if(isTRUE(control$penalties)) {

    control$reg <- TRUE

    control$penalties <- list(
      # beta  = list(alpha = 0, lambda = 0, power = 0),
      beta  = list(alpha = 0),
      class = list(alpha = 1),
      prob  = list(alpha = 1),
      sd    = list(alpha = 1),
      Sigma = list(alpha = 1)
    )

  } else if(is.list(control$penalties)) {

    control$reg <- TRUE

  } else {

    stop("penalties should be TRUE, FALSE, or a list")

  }

  if(is.null(control$step_maxit)) {
    control$step_maxit <- 30L
  } else if(control$step_maxit < 1L) {
    stop("step_maxit must be an integer greater than 0")
  }

  if(is.null(control$c1)) {
    control$c1 <- 0.5
  } else if(control$c1 < 0L) {
    stop("c1 must be a positive number, preferable lower than c2")
  }

  if(is.null(control$c2)) {
    control$c2 <- 0.5
  } else if(control$c2 < 0L) {
    stop("c2 must be a positive number, preferable larger than c1")
  }

  if(is.null(control$step_eps)) {
    control$step_eps <- 1e-09
  } else if(control$step_eps < 0) {
    stop("step_eps must be a positive number, preferable close to 0")
  }

  if(is.null(control$df_eps)) {
    control$df_eps <- 1e-09
  } else if(control$df_eps < 0) {
    stop("df_eps must be a positive number, preferable close to 0")
  }

  if(is.null(control$eps)) {
    control$eps <- 1e-05
  } else if(control$eps < 0) {
    stop("eps must be a positive number, preferable close to 0")
  }

  if(is.null(control$M)) {
    control$M <- 100L
  } else if(control$M < 0L) {
    stop("M must be a positive integer")
  }

  if(is.null(control$ss_fac)) {
    control$ss_fac <- 2
  } else if(control$ss_fac <= 1) {
    stop("ss_fac must be a positive integer larger than 1")
  }

  if(is.null(control$maxit)) {
    control$maxit <- 1000L
  } else if(control$maxit < 0L) {
    stop("maxit must be a positive integer")
  }

  if(is.null(control$cores)) {
    control$cores <- parallel::detectCores()-1L
  } else if(control$cores < 0L) {
    stop("cores must be a positive integer")
  }

  if(is.null(control$tcg_maxit)) {
    control$tcg_maxit <- 10L
  } else if(control$tcg_maxit < 0L) {
    stop("tcg_maxit must be a positive integer")
  }

  if(is.null(control$opt)) {
    control$opt <- "lbfgs"
  }

  if(control$opt == "newton") {
    if(is.null(control$tcg_maxit)) {
      control$tcg_maxit <- 100L
    }
  }

  if(is.null(control$rstarts)) {
    control$rstarts <- 16L
  } else if(control$rstarts < 0L) {
    stop("rstarts must be a positive integer")
  }

  return(control)

}

extract_cov_pairs <- function(model) {

  flatten_text <- function(x) {
    if (is.null(x)) {
      return(character(0))
    }

    if (is.list(x)) {
      return(unlist(lapply(x, flatten_text), use.names = FALSE))
    }

    out <- tryCatch(as.character(x), error = function(e) character(0))
    out[!is.na(out)]
  }

  txt <- flatten_text(model)

  if (length(txt) == 0L || !any(grepl("~~", txt, fixed = TRUE))) {
    return(NULL)
  }

  lines <- unlist(strsplit(txt, "\n", fixed = TRUE), use.names = FALSE)
  lines <- trimws(lines[grepl("~~", lines, fixed = TRUE)])

  if (length(lines) == 0L) {
    return(NULL)
  }

  pairs_list <- lapply(lines, function(line) {

    line <- sub("#.*$", "", line)
    line <- trimws(line)

    vars <- trimws(unlist(strsplit(line, "\\s*~~\\s*"), use.names = FALSE))
    vars <- vars[nzchar(vars)]

    if (length(vars) < 2L) {
      return(NULL)
    }

    out <- t(utils::combn(vars, 2))
    colnames(out) <- c("lhs", "rhs")
    out

  })

  pairs_list <- Filter(Negate(is.null), pairs_list)

  if (length(pairs_list) == 0L) {
    return(NULL)
  }

  pairs <- do.call(rbind, pairs_list)

  as.data.frame(pairs, stringsAsFactors = FALSE)

}
