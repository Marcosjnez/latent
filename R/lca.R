# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 10/06/2026
#'
#'
#' Estimate latent class models for continuous and categorical indicators, with
#' optional covariates for latent class membership probabilities.
#'
#' @description
#' \code{lca()} fits latent class models in which the measurement model may
#' contain Gaussian indicators, multinomial indicators, or both. Class membership
#' probabilities can be unconditional or modeled as a multinomial-logit function
#' of observed covariates. The function also supports class enumeration by
#' passing several values to \code{nclasses}.
#'
#' @usage
#' lca(data, nclasses = 1L, gaussian = NULL, multinomial = NULL,
#' X = NULL, Y = NULL, penalties = TRUE, model = NULL, start = NULL,
#' control = NULL, do.fit = TRUE, verbose = TRUE)
#'
#' @param data A \code{data.frame} or \code{matrix} containing the observed data.
#' Columns listed in \code{gaussian} and/or \code{multinomial} are used as
#' measurement indicators. Columns used in \code{X} may also be taken from
#' \code{data}.
#' @param nclasses A positive integer giving the number of latent classes. If a
#' vector of positive integers is supplied, one model is fitted for each value
#' and the result is returned as an object of class \code{"llcalist"}.
#' @param gaussian Optional character vector with the names of the variables in
#' \code{data} to be modeled as Gaussian indicators. These variables are
#' modeled with class-specific means and variances. If residual covariance
#' terms involving Gaussian indicators are specified in \code{model}, the
#' corresponding variables are modeled jointly with a class-specific
#' multivariate Gaussian distribution.
#' @param multinomial Optional character vector with the names of the variables
#' in \code{data} to be modeled as multinomial indicators. Variables are
#' internally converted to factors, and the first factor level is used as the
#' reference category for the multinomial logits. Variables with more than 30
#' observed categories are not allowed.
#' @param X Optional covariates for predicting latent class membership. It can be
#' \code{NULL}, a character vector with column names in \code{data}, or a
#' \code{data.frame}/\code{matrix} with the same number of rows as
#' \code{data}. An intercept is always included. Character covariates are
#' converted to factors and factors are expanded using \code{model.matrix()}.
#' @param penalties Logical value or named list controlling regularization.
#' If \code{FALSE}, no penalties are used. If \code{TRUE}, default penalties are
#' used. If a named list is supplied, missing penalty blocks are filled with
#' their defaults. Valid penalty blocks are \code{beta}, \code{class},
#' \code{prob}, \code{var}, and \code{Sigma}. The default is
#' \code{list(beta = list(alpha = 0), class = list(alpha = 1),
#' prob = list(alpha = 1), var = list(alpha = 1),
#' Sigma = list(alpha = 1))}. The \code{beta} block may also contain
#' \code{lambda} and \code{power} for ridge-type regularization.
#' @param model Optional model specification. This can be a named list used to
#' fix parameters, impose equality constraints, or provide custom parameter
#' labels. Names should match the internal parameter blocks, for example
#' \code{beta}, Gaussian item names, \code{log_<item>} for multinomial logits,
#' \code{means}, \code{logsigma}, or \code{sigma|Class<i>} for multivariate
#' Gaussian blocks. Character strings containing residual-dependency syntax
#' such as \code{"y1 ~~ y2"} are also used to identify residual covariances
#' among Gaussian indicators or residual associations among multinomial
#' indicators. If an object of class \code{"llca"} is supplied, its measurement
#' parameters are reused while the class-membership regression coefficients
#' are re-estimated.
#' @param start Optional named list of starting values. Names should correspond
#' to parameter blocks in the model. Supplied values replace the corresponding
#' default initial values, allowing partial specification of starting values.
#' @param control Optional list of optimizer and estimation controls. Common
#' entries include \code{rstarts} for the number of random starts,
#' \code{cores} for parallel computation, \code{maxit} for the maximum number
#' of optimizer iterations, \code{opt} for the optimizer type, and convergence
#' tolerances such as \code{eps}, \code{df_eps}, and \code{step_eps}. Missing
#' entries are replaced by internal defaults.
#' @param do.fit Logical. If \code{TRUE}, the model is estimated. If \code{FALSE},
#' the function returns an \code{"llca"} object containing the processed data,
#' model structure, and optimization setup, but without running the optimizer.
#' @param verbose Logical. If \code{TRUE}, progress information is printed,
#' especially when fitting several values of \code{nclasses}.
#'
#' @details
#' The measurement model is defined by the variables supplied through
#' \code{gaussian} and \code{multinomial}. Variables not listed there are ignored
#' by the measurement model unless they are used as covariates in \code{X}.
#'
#' For Gaussian indicators, the model estimates class-specific means and
#' variances. For multinomial indicators, the model estimates class-specific
#' category probabilities through a softmax parameterization with the first
#' category fixed as the reference.
#'
#' When covariates are supplied through \code{X}, class probabilities are modeled
#' through multinomial-logit coefficients stored in the \code{beta} parameter
#' block. The first class is the reference class, so its coefficients are fixed
#' to zero for identification.
#'
#' The function compresses the data into unique response/covariate patterns
#' before estimation and uses their frequencies as weights. This can make
#' estimation more efficient, especially for categorical data with repeated
#' response patterns.
#'
#' Residual dependencies can be requested by including covariance-style
#' statements such as \code{"y1 ~~ y2"} in \code{model}. For Gaussian indicators,
#' such terms define class-specific residual covariance parameters. For
#' multinomial indicators, they define joint categorical response structures
#' with residual interaction parameters.
#'
#' @return
#' If \code{length(nclasses) == 1}, an S4 object of class \code{"llca"} with
#' the following slots:
#' \describe{
#' \item{\code{version}}{Version of the \pkg{latent} package used to fit the model.}
#' \item{\code{call}}{The matched function call.}
#' \item{\code{timing}}{Elapsed optimization time. Empty when \code{do.fit = FALSE}.}
#' \item{\code{dataList}}{Processed data objects, including the measurement data,
#' covariate design matrix, unique response patterns, pattern weights, and
#' mappings between original rows and unique patterns.}
#' \item{\code{modelInfo}}{Internal model structure used by the optimizer,
#' including parameter labels, transformed-parameter labels, degrees of
#' freedom, manifolds, transformations, estimators, and optimizer controls.}
#' \item{\code{Optim}}{Raw output from the optimizer. Empty when
#' \code{do.fit = FALSE}.}
#' \item{\code{parameters}}{Estimated model parameters on the estimation scale,
#' organized by parameter block. Empty when \code{do.fit = FALSE}.}
#' \item{\code{transformed_pars}}{Estimated parameters after applying model
#' transformations, including class probabilities, item probabilities,
#' variances, standard deviations, and log-likelihood components. Empty when
#' \code{do.fit = FALSE}.}
#' \item{\code{extra}}{Additional information reserved for downstream methods.}
#' }
#'
#' If \code{nclasses} contains several values, a list of \code{"llca"} objects is
#' returned with class \code{"llcalist"}.
#'
#' @examples
#' \dontrun{
#' # Three-class model for categorical indicators
#' fit <- lca(
#' data = gss82,
#' nclasses = 3L,
#' multinomial = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"),
#' penalties = TRUE
#' )
#'
#' fit
#' summary(fit)
#' getfit(fit)
#'
#' latInspect(fit, what = "coefs", digits = 3)
#' latInspect(fit, what = "classes", digits = 3)
#' latInspect(fit, what = "profile", digits = 3)
#' latInspect(fit, what = "posterior", digits = 3)
#'
#' # Class enumeration
#' fits <- lca(
#' data = gss82,
#' nclasses = 1:4,
#' multinomial = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"),
#' penalties = TRUE
#' )
#'
#' # Latent class regression with covariates
#' fit_x <- lca(
#' data = empathy,
#' nclasses = 4L,
#' gaussian = c("ec1", "ec2", "ec3", "ec4", "ec5", "ec6"),
#' X = c("sex", "pt1", "pt2", "pt3", "pt4")
#' )
#'
#' }
#'
#' @references
#' None yet.
#'
#' @export
lca <- function(data,
                nclasses = 1L,
                gaussian = NULL,
                multinomial = NULL,
                X = NULL,
                Y = NULL,
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
                                              Y = Y,
                                              gaussian = gaussian,
                                              multinomial = multinomial,
                                              model = model,
                                              penalties = penalties,
                                              control = control,
                                              start = start)
  list2env(dataList_and_control, envir = environment())
  dataList$args <- args

  #### Create the model ####

  # Fix the measurement part of the model if a previous llca fit was provided
  # model. This is usually done for two-step estimation:
  if(class(model) == "llca") {

    control$model <- model
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
  names(Optim$g) <- modelInfo$parameters_labels
  names(Optim$rg) <- modelInfo$parameters_labels
  names(Optim$dir) <- modelInfo$parameters_labels

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
                                Y = NULL,
                                gaussian = NULL,
                                multinomial = NULL,
                                model = NULL,
                                control = list(),
                                penalties = FALSE,
                                start = NULL) {

  #### Process the covariates ####

  original_X <- X
  # Matrix of predictors for latent class probabilities:
  # X can be a character vector with the name of the predictors to be found in
  # data or a data.frame or matrix with named predictor variables
  X <- make_design_matrix(X = X, data = data)
  pX <- ncol(X) # Number of columns of the design matrix

  # Remove participants with missing data in any covariate:
  remove <- which(apply(X, MARGIN = 1, FUN = anyNA))

  if(length(remove) > 0) {

    missing <- colSums(is.na(X))
    missing <- missing[missing > 0]

    warning(
      "Missing data detected in covariates:\n",
      paste0(
        names(missing), ": ",
        missing,
        " missing observation",
        ifelse(missing == 1L, "", "s"),
        collapse = "\n"
      ),
      "\nIn total, ",
      length(remove),
      " observation",
      if (length(remove) == 1L) " was" else "s were",
      " removed."
    )

    X <- X[-remove, ]
    data <- data[-remove, ]
  }

  # Number of subjects:
  nobs <- nrow(data)

  #### Indicator types ####

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
  measurement_recoded <- measurement

  # match each item with a likelihood model:
  item <- model_labels[idx]
  item_names <- colnames(measurement)
  # Number of items:
  nitems <- ncol(measurement)

  # Check that item is a character vector with a string for each column of data:
  if(!is.character(item) || length(item) != nitems) {
    stop("item must be a character vector with as many elements as columns in data")
  }

  # Indicators:

  any_gaussian <- any(item == "gaussian")
  if(any_gaussian) {
    gaussian_names <- item_names[which(item == "gaussian")]
    Jgauss <- length(gaussian_names) # Number of gaussian items
    control$gaussian <- list(gaussian_names = gaussian_names,
                             Jgauss = Jgauss)
  }

  any_mvgaussian <- any(item == "mvgaussian")
  if(any_mvgaussian) {
    mvgaussian_names <- item_names[which(item == "mvgaussian")]
    Jmvgauss <- length(mvgaussian_names) # Number of mvgaussian items
    control$mvgaussian <- list(mvgaussian_names = mvgaussian_names,
                               Jmvgauss = Jmvgauss,
                               target = target)
  }

  any_multinomial <- any(item == "multinomial")
  if(any_multinomial) {
    multinomial_names <- item_names[which(item == "multinomial")]
    logmultinomial_names <- paste("log", multinomial_names, sep = "_")
    Jmultinom <- length(multinomial_names) # Number of multinomial items

    # Transform multinomial variables into factors:
    measurement[, multinomial_names] <- lapply(measurement[, multinomial_names,
                                                           drop = FALSE],
                                               FUN = factor)
    # Transform the factor variables into integers:
    measurement_recoded[, multinomial_names] <- lapply(measurement[, multinomial_names,
                                                                drop = FALSE],
                                                    FUN = function(col) {
                                                      if (is.factor(col)) as.integer(col) - 1L else col
                                                    })

    # Save levels of factor variables:
    multinomial_factor_levels <- lapply(measurement[, multinomial_names, drop = FALSE], levels)
    multinomial_factor_lengths <- unlist(lapply(multinomial_factor_levels, FUN = length))
    multinomial_reference <- lapply(multinomial_factor_levels, FUN = \(x) x[1])

    control$multinomial <- list(multinomial_names = multinomial_names,
                                logmultinomial_names = logmultinomial_names,
                                Jmultinom = Jmultinom,
                                multinomial_factor_levels = multinomial_factor_levels,
                                multinomial_factor_lengths = multinomial_factor_lengths,
                                multinomial_reference = multinomial_reference)

  }

  #### Possible response patterns and degrees of freedom ####

  # Get the number of possible response patterns:
  if(all(item == "multinomial")) { # If all the items are multinomial...

    # Count the number of categories:
    Ks <- apply(measurement[, multinomial_names, drop = FALSE], MARGIN = 2,
                FUN = \(x) length(unique(x[!is.na(x)])))
    # If all the items are multinomial, the number of possible patterns is:
    npossible_patterns <- min(prod(Ks)-1L, nobs)
    # From technical manual of LatentGOLD 6.1 (p.68)

  } else { # If any item is not multinomial...

    npossible_patterns <- nobs

  }

  #### Process the data ####

  # Outcomes:
  if(is.null(Y)) {
    control$outcomes <- FALSE
    Y_patterns <- NULL
  } else {
    control$outcomes <- TRUE
    if (is.character(Y)) {
      Y <- data[, Y, drop = FALSE]
    } else {
      # Check that X is either a data.frame or a matrix:
      if (!is.data.frame(Y) && !is.matrix(Y)) {
        stop("Y must be a character vector, matrix or data.frame")
      }
      if (nrow(Y) != nobs) {
        stop("Number of cases in the data and covariates does not match")
      }
    }
  }

  # Append the covariates, indicators, and distal outcomes.
  # This is necessary to find the unique data patterns. Later, we split again:
  if(control$outcomes) {
    dt <- data.table::as.data.table(cbind(X, measurement_recoded, Y))
  } else {
    dt <- data.table::as.data.table(cbind(X, measurement_recoded))
  }
  # Collect some information from the measurement + covariates:
  counts_dt <- dt[, .(index = .I[1], count = .N), by = names(dt)]
  # Number of unique response patterns:
  npatterns <- nrow(counts_dt)
  pattern_names <- paste("pattern", 1:npatterns, sep = "")
  # Measurement data with unique response patterns:
  patterns <- as.matrix(counts_dt[, colnames(measurement_recoded), with = FALSE])
  # Covariate data with unique patterns:
  cov_patterns <- as.matrix(counts_dt[, colnames(X), with = FALSE])
  # Distal outcomes with unique patterns:
  if(control$outcomes) {
    Y_patterns <- as.matrix(counts_dt[, colnames(Y), with = FALSE])
    rownames(Y_patterns) <- pattern_names
  }
  # Pattern names:
  rownames(patterns) <- rownames(cov_patterns) <- pattern_names
  # Counts of each response pattern:
  weights <- counts_dt$count
  # Indices to map the full data to the unique patterns:
  full2short <- counts_dt$index
  # Indices to map unique patterns to the full data:
  short2full <- match(do.call(paste, dt),
                      do.call(paste, counts_dt[, -c("index", "count"),
                                               with = FALSE]))

  # Compute the nonrecoded patterns of items (patterns_original):
  if(control$outcomes) {
    dt <- data.table::as.data.table(cbind(X, measurement, Y))
  } else {
    dt <- data.table::as.data.table(cbind(X, measurement))
  }
  counts_dt <- dt[, .(index = .I[1], count = .N), by = names(dt)]
  # patterns_original <- as.matrix(counts_dt[, colnames(measurement_recoded),, with = FALSE])
  patterns_original <- as.data.frame(counts_dt[, colnames(X), with = FALSE])
  rownames(patterns_original) <- pattern_names

  #### Store objects in dataList ####

  dataList <- vector("list")
  dataList$data <- data
  dataList$measurement <- measurement
  dataList$X <- X
  dataList$original_X <- original_X
  dataList$Y <- Y
  dataList$Y_patterns <- Y_patterns
  dataList$measurement_recoded <- measurement_recoded
  dataList$nitems <- nitems
  dataList$item <- item
  dataList$item_names <- item_names
  dataList$nobs <- nobs
  dataList$patterns <- patterns
  dataList$cov_patterns <- cov_patterns
  dataList$npatterns <- npatterns
  dataList$npossible_patterns <- npossible_patterns
  dataList$pattern_names <- pattern_names
  dataList$pX <- pX
  dataList$weights <- weights
  dataList$full2short <- full2short
  dataList$short2full <- short2full
  dataList$patterns_original <- patterns_original
  dataList$any_gaussian <- any_gaussian
  dataList$any_multinomial <- any_multinomial
  dataList$any_mvgaussian <- any_mvgaussian

  # Update and check control parameters:
  control$penalties <- penalties # Either a logical or a list of named penalties
  control$start <- start # Optional named list with initial values
  control <- lca_control(control) # Check and update the control inputs and create defaults

  #### Check for residual dependencies in multinomial items ####

  if(length(multinomial) >= 2L) {

    # Define the joint probability parameters:
    error_covs <- NULL

    if (!is.null(model)) {
      error_covs <- extract_cov_pairs(model)
    }

    if (!is.null(error_covs)) {

      error_covs$lhs <- as.character(error_covs$lhs)
      error_covs$rhs <- as.character(error_covs$rhs)

      # Keep only pairs where both variables are multinomial:
      keep <- error_covs$lhs %in% multinomial &
        error_covs$rhs %in% multinomial &
        error_covs$lhs != error_covs$rhs

      error_covs <- error_covs[keep, , drop = FALSE]

      if(nrow(error_covs) > 0L) {

        npairs <- nrow(error_covs)
        pairs <- vector("list", length = npairs)
        loginter_names <- vector("list", length = npairs)

        for(h in 1:npairs) {

          pairs[[h]] <- c(error_covs[h, 1], error_covs[h, 2])
          loginter_names[[h]] <- paste("log_", pairs[[h]], collapse = "x",
                                       sep = "")
        }

        names(pairs) <- unlist(lapply(pairs, FUN = paste, collapse = "."))

        removed_multinomial <- unique(unlist(pairs))
        keep <- multinomial_names[!(multinomial_names %in% removed_multinomial)]

        indep_pairs <- group_connected_pairs(error_covs)$pairs
        indep_pairs_list <- group_connected_pairs(error_covs)$elements
        names(indep_pairs) <- names(indep_pairs_list) <-
          lapply(indep_pairs_list, FUN = paste, collapse = ".")
        joints <- lapply(indep_pairs_list, FUN = \(x) {
          apply(patterns_original[, x], MARGIN = 1, FUN = paste, collapse = ".")
        })
        joints <- as.data.frame(joints)
        new_probs <- unlist(c(keep, colnames(joints)))
        new_J <- length(new_probs)
        patterns_mvmultinomial <- data.frame(patterns_original, joints)[, new_probs]
        joint_levels <- lapply(indep_pairs_list, FUN = \(x) {
          apply(expand.grid(multinomial_factor_levels[x]),
                MARGIN = 1, FUN = paste, collapse = ".")
        })
        new_levels <- c(multinomial_factor_levels[keep], joint_levels)
        new_K <- unlist(lapply(new_levels, FUN = length))
        for(j in 1:length(new_probs)) {
          patterns_mvmultinomial[[new_probs[j]]] <-
            factor(patterns_mvmultinomial[[new_probs[j]]], ordered = TRUE,
                   levels = new_levels[[j]])
        }
        njoints <- length(joint_levels)
        multinomial_factor_levels_numeric <- sapply(multinomial_factor_lengths,
                                                    FUN = seq_len)
        joint_orderings <- lapply(indep_pairs_list, FUN = \(x) {
          expand.grid(multinomial_factor_levels_numeric[x])
        })

        patterns_mvmultinomial_recoded <- lapply(patterns_mvmultinomial,
                    FUN = function(col) {
                      if (is.factor(col)) as.integer(col) - 1L else col
                    })
        patterns_mvmultinomial_recoded <-
          as.matrix(as.data.frame(patterns_mvmultinomial_recoded))

        nremoved_multinomial <- length(removed_multinomial)
        removed_logmultinomial <- paste("log_", removed_multinomial,
                                        "|reference", sep = "")
        removed_condmultinomial <- paste(removed_multinomial,
                                         "|reference", sep = "")
        removed_multinomial_factor_levels <- multinomial_factor_levels[removed_multinomial]
        removed_multinomial_factor_lengths <- multinomial_factor_lengths[removed_multinomial]

        control$mvmultinomial <- list(
          npairs = npairs,
          pairs = pairs,
          joints = joints,
          loginter_names = loginter_names,
          removed_multinomial = removed_multinomial,
          new_J = new_J,
          new_K = new_K,
          patterns_mvmultinomial = patterns_mvmultinomial,
          patterns_mvmultinomial_recoded = patterns_mvmultinomial_recoded,
          nremoved_multinomial = nremoved_multinomial,
          removed_logmultinomial = removed_logmultinomial,
          removed_condmultinomial = removed_condmultinomial,
          removed_multinomial_factor_levels = removed_multinomial_factor_levels,
          removed_multinomial_factor_lengths = removed_multinomial_factor_lengths,
          indep_pairs = indep_pairs,
          indep_pairs_list = indep_pairs_list,
          new_levels = new_levels,
          njoints = njoints,
          joint_orderings = joint_orderings
        )

      }

    }

  }

  #### Return ####

  return(list(dataList = dataList, control = control))

}

create_lca_model <- function(dataList, nclasses, item,
                             model = NULL, control) {

  # Generate the model syntax and initial parameter values

  list2env(dataList, envir = environment())

  class_names <- paste("Class", 1:nclasses, sep = "")
  control$class_names <- class_names

  #### Model for the transformed parameters ####

  pred_names <- colnames(X) # Names of predictors

  list_struct <- vector("list")
  k <- 1L

  # Model for the betas:

  list_struct[[k]] <- list(name = "beta",
                           type = "matrix",
                           dim = c(pX, nclasses),
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
  if(any_gaussian) {

    list2env(control$gaussian, envir = environment())

    for(name in gaussian_names) {

      list_struct[[k]] <- list(name = name,
                               type = "matrix",
                               dim = c(4L, nclasses),
                               rownames = c("mean",
                                            "var",
                                            "log(var)",
                                            "stdv"),
                               colnames = class_names)
      k <- k+1L

    }

  }

  # Model for multinomial items:
  if(any_multinomial) {

    list2env(control$multinomial, envir = environment())

    if(any(multinomial_factor_lengths > 30)) {
      stop("You cannot use a multinomial likelihood to model a variable with more than 30 unique categories")
    }

    has_mvmultinomial <- !is.null(control$mvmultinomial)
    removed_multinomial <- character(0L)

    if(has_mvmultinomial) {
      list2env(control$mvmultinomial, envir = environment())
    }

    for(j in seq_len(Jmultinom)) {

      # For ordinary multinomial items, the item probability matrix is obtained
      # directly from the item logits with a softmax transformation. For items
      # involved in residual dependencies, the original item probability matrix
      # is kept in the transformed object, but it is later computed as the
      # marginal probability obtained by summing over the corresponding joint
      # probability matrix.
      if(!(multinomial_names[j] %in% removed_multinomial)) {

        list_struct[[k]] <- list(name = logmultinomial_names[j],
                                 type = "matrix",
                                 dim = c(multinomial_factor_lengths[j], nclasses),
                                 rownames = multinomial_factor_levels[[j]],
                                 colnames = class_names)
        k <- k+1L

      }

      list_struct[[k]] <- list(name = multinomial_names[j],
                               type = "matrix",
                               dim = c(multinomial_factor_lengths[j], nclasses),
                               rownames = multinomial_factor_levels[[j]],
                               colnames = class_names)
      k <- k+1L

    }

    if(has_mvmultinomial) {

      # Conditional probability blocks for the multinomial variables involved
      # in residual dependencies. These are probabilities of the levels of one
      # variable conditional on the other variable being at its reference level.
      # They are not used as marginal item probabilities in the likelihood.
      for(j in seq_len(nremoved_multinomial)) {

        list_struct[[k]] <- list(name = removed_logmultinomial[j],
                                 type = "matrix",
                                 dim = c(removed_multinomial_factor_lengths[j],
                                         nclasses),
                                 rownames = removed_multinomial_factor_levels[[j]],
                                 colnames = class_names)
        k <- k+1L

        list_struct[[k]] <- list(name = removed_condmultinomial[j],
                                 type = "matrix",
                                 dim = c(removed_multinomial_factor_lengths[j],
                                         nclasses),
                                 rownames = removed_multinomial_factor_levels[[j]],
                                 colnames = class_names)
        k <- k+1L

      }

      for(h in seq_len(npairs)) {

        levels1 <- multinomial_factor_levels[[pairs[[h]][1]]]
        levels2 <- multinomial_factor_levels[[pairs[[h]][2]]]

        list_struct[[k]] <- list(name = loginter_names[[h]],
                                 type = "matrix",
                                 dim = c(length(levels1), length(levels2)),
                                 rownames = levels1,
                                 colnames = levels2)
        k <- k+1L

      }

      for(h in 1:njoints) {

        joint_name <- paste(indep_pairs_list[[h]], collapse = ".")
        logjoint_name <- paste("log_", joint_name, sep = "")

        list_struct[[k]] <- list(name = logjoint_name,
                                 type = "matrix",
                                 dim = c(new_K[joint_name], nclasses),
                                 rownames = new_levels[[joint_name]],
                                 colnames = class_names)
        k <- k+1L

        list_struct[[k]] <- list(name = joint_name,
                                 type = "matrix",
                                 dim = c(new_K[joint_name], nclasses),
                                 rownames = new_levels[[joint_name]],
                                 colnames = class_names)
        k <- k+1L

      }

    }

  }

  # Model for multivariate normal items:
  if(any_mvgaussian) {

    list2env(control$mvgaussian, envir = environment())

    # mean_names <- paste("means|Class", 1:nclasses, sep = "")
    sigma_names <- paste("sigma|Class", 1:nclasses, sep = "")

    list_struct[[k]] <- list(name = "means",
                             type = "matrix",
                             dim = c(Jmvgauss, nclasses),
                             rownames = mvgaussian_names,
                             colnames = class_names)
    k <- k+1L

    list_struct[[k]] <- list(name = "logsigma",
                             type = "matrix",
                             dim = c(Jmvgauss, nclasses),
                             rownames = mvgaussian_names,
                             colnames = class_names)
    k <- k+1L

    # Sigma:
    for(j in 1:nclasses) {

      list_struct[[k]] <- list(name = sigma_names[j],
                               type = "matrix",
                               dim = c(Jmvgauss, Jmvgauss),
                               rownames = mvgaussian_names,
                               colnames = mvgaussian_names,
                               symmetric = TRUE)
      k <- k+1L

    }

  }

  # Distal outcomes:
  if(control$outcomes) {

    list_struct[[k]] <- list(name = "distal_beta",
                             type = "matrix",
                             dim = c(nclasses, ncol(Y_patterns)),
                             rownames = class_names,
                             colnames = colnames(Y_patterns))
    k <- k+1L

  }

  # Create the full transparameter structure:
  trans <- create_parameters(list_struct)

  #### Model for the parameters ####

  param <- list()

  # Model for the betas:
  param$beta <- trans$beta
  param$beta[, 1] <- "0"

  # Model for gaussian items:
  if(any_gaussian) {

    cl <- length(param)
    idx_gauss <- (cl + 1L):(cl + Jgauss)
    param[idx_gauss] <- trans[gaussian_names]
    # Fix the standard deviations just to avoid that they are taken
    # as parameters and not as fixed parameters:
    param[idx_gauss] <- lapply(trans[gaussian_names],
                               FUN = \(x) {
                                 x[2, ] <- "1";
                                 x[4, ] <- "1"
                                 return(x)
                               })
    names(param)[idx_gauss] <- gaussian_names

  }

  # Model for multinomial items:
  if(any_multinomial) {

    multinomial_param_names <- multinomial_names
    logmultinomial_param_names <- logmultinomial_names

    if(!is.null(control$mvmultinomial)) {

      list2env(control$mvmultinomial, envir = environment())

      multinomial_param_names <- multinomial_names[
        !(multinomial_names %in% removed_multinomial)
      ]
      logmultinomial_param_names <- paste("log", multinomial_param_names,
                                          sep = "_")

    }

    if(length(multinomial_param_names) > 0L) {

      cl <- length(param)
      idx_multinom <- (cl + 1L):(cl + length(logmultinomial_param_names))
      param[idx_multinom] <- trans[logmultinomial_param_names]
      names(param)[idx_multinom] <- logmultinomial_param_names

      # Fix the first item category to zero for identification after the
      # softmax transformation:
      param[idx_multinom] <- lapply(param[idx_multinom],
                                    FUN = \(x) {
                                      x[1, ] <- "0"; return(x)
                                    })

    }

    if(!is.null(control$mvmultinomial)) {

      # Conditional logits for variables involved in residual dependencies.
      # These generate P(u_j | reference level of the paired variable), not the
      # marginal P(u_j). Marginal probabilities are computed later from the
      # joint probabilities using sum_vectors transformations.
      cl <- length(param)
      idx_cond <- (cl + 1L):(cl + nremoved_multinomial)
      param[idx_cond] <- trans[removed_logmultinomial]
      names(param)[idx_cond] <- removed_logmultinomial

      param[idx_cond] <- lapply(param[idx_cond],
                                FUN = \(x) {
                                  x[1, ] <- "0"; return(x)
                                })

      param[unlist(loginter_names)] <- trans[unlist(loginter_names)]
      for(j in 1:npairs) {
        param[[loginter_names[[j]]]][1, ] <- "0"
        param[[loginter_names[[j]]]][, 1] <- "0"
      }
    }

  }

  # Model for multivariate gaussian items:
  if(any_mvgaussian) {

    param$means <- trans$means
    param$logsigma <- trans$logsigma
    # param[sigma_name] <- trans[sigma_names]
    ordering <- rownames(trans[[sigma_names[1]]])
    target <- target[ordering, ordering]
    for(j in 1:nclasses) {
      param[[sigma_names[j]]] <- diag(Jmvgauss)
      param[[sigma_names[j]]][target] <- trans[[sigma_names[j]]][target]
      # diag(param[[sigma_names[j]]]) <- diag(trans[[sigma_names[j]]])
      dimnames(param[[sigma_names[j]]]) <- dimnames(trans[[sigma_names[j]]])
    }

  }

  # Distal outcomes:
  param$distal_beta <- trans$distal_beta

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

  if(any_gaussian) {

    init_mean <- init_var <- init_logvar <- list()
    j <- 1L
    for(name in gaussian_names) {

      init_mean[[j]] <- rep(mean(measurement_recoded[, name], na.rm = TRUE),
                            times = nclasses)
      init_var[[j]] <- rep(var(measurement_recoded[, name], na.rm = TRUE),
                          times = nclasses)
      init_logvar[[j]] <- log(init_var[[j]])
      j <- j+1L

    }

  }

  pi_hat_list <- list()
  if(any_multinomial) {

    eta_hat_list <- list()
    j <- 1L
    for(name in multinomial_names) {

      int_vector <- measurement_recoded[, name]
      int_vector <- int_vector[!is.na(int_vector)] # Remove missing values
      nsize <- length(int_vector)
      props <- count(int_vector, nsize, multinomial_factor_lengths[j]) / nsize
      pi_hat_list[[j]] <- props %*% t(rep(1, nclasses))
      log_props <- log(props)
      eta_hat_list[[j]] <- (log_props-log_props[1]) %*% t(rep(1, nclasses))
      j <- j+1L

    }

    pi_hat_vector <- unlist(pi_hat_list)
    vars <- (1-pi_hat_vector)/(nobs*pi_hat_vector)
    sds <- sqrt(vars)
    Ks <- length(unlist(eta_hat_list))

  }

  if(any_mvgaussian) {

    init_mean <- matrix(rep(colMeans(measurement_recoded[, mvgaussian_names],
                                     na.rm = TRUE),
                            times = nclasses), nrow = Jmvgauss, ncol = nclasses)
    init_logsigma <- matrix(rep(apply(measurement_recoded[, mvgaussian_names],
                                      MARGIN = 2,
                                      FUN = var,
                                      na.rm = TRUE),
                                times = nclasses), nrow = Jmvgauss,
                            ncol = nclasses)
    init_sigma <- vector("list", length = nclasses)
    for(j in 1:nclasses) {

      init_sigma[[j]] <- cov(measurement_recoded[, mvgaussian_names],
                             use = "pairwise.complete.obs")

    }

  }

  for(i in 1:control$rstarts) {

    init_param[[i]] <- vector("list", length = length(param))
    names(init_param[[i]]) <- names(param)

    # Initial values for betas:
    init_beta <- matrix(rnorm(pX*nclasses), nrow = pX, ncol = nclasses)
    init_param[[i]][["beta"]] <- init_beta
    init_param[[i]]$beta[, 1] <- 0
    dimnames(init_param[[i]]$beta) <- dimnames(param$beta)

    # Initial values for gaussian items:
    if(any_gaussian) {

      j <- 1L
      for(name in gaussian_names) {

        rmean <- init_mean[[j]] + rnorm(nclasses,
                                        mean = 0,
                                        sd = init_var[[j]]/sqrt(nobs))

        init_param[[i]][[name]] <- rbind(rmean,
                                         init_var[[j]],
                                         init_logvar[[j]],
                                         sqrt(init_var[[j]]))
        dimnames(init_param[[i]][[name]]) <-
          dimnames(param[[name]])
        j <- j+1L

      }

    }

    # Initial values for multinomial items:
    if(any_multinomial) {

      multinomial_init_names <- multinomial_names
      if(!is.null(control$mvmultinomial)) {
        multinomial_init_names <- multinomial_names[
          !(multinomial_names %in% removed_multinomial)
        ]
      }

      for(name in multinomial_init_names) {

        j <- match(name, multinomial_names)
        log_name <- logmultinomial_names[j]
        sds <- (1-pi_hat_list[[j]])/(nobs*pi_hat_list[[j]])

        init_param[[i]][[log_name]] <- eta_hat_list[[j]] +
          rnorm(length(eta_hat_list[[j]]), 0, sds)
        init_param[[i]][[log_name]][1, ] <- 0
        dimnames(init_param[[i]][[log_name]]) <- dimnames(param[[log_name]])

      }

      if(!is.null(control$mvmultinomial)) {

        for(name in removed_multinomial) {

          j <- match(name, multinomial_names)
          log_name <- paste("log_", name, "|reference", sep = "")
          sds <- (1-pi_hat_list[[j]])/(nobs*pi_hat_list[[j]])

          init_param[[i]][[log_name]] <- eta_hat_list[[j]] +
            rnorm(length(eta_hat_list[[j]]), 0, sds)
          init_param[[i]][[log_name]][1, ] <- 0
          dimnames(init_param[[i]][[log_name]]) <- dimnames(param[[log_name]])

        }

        for(j in 1:npairs) {
          dims <- dim(trans[[loginter_names[[j]]]])
          init_param[[i]][[loginter_names[[j]]]] <- matrix(0, nrow = dims[1],
                                                           ncol = dims[2])
          dimnames(init_param[[i]][[loginter_names[[j]]]]) <-
            dimnames(param[[loginter_names[[j]]]])
        }

      }

    }

    # Initial values for multivariate gaussian items:
    if(any_mvgaussian) {

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

    # Distal outcomes:
    if(control$outcomes) {
      init_param[[i]][["distal_beta"]] <- matrix(rnorm(nclasses*ncol(Y_patterns)),
                                                 nrow = nclasses, ncol = ncol(Y_patterns))
      dimnames(init_param[[i]]$distal_beta) <- dimnames(param$distal_beta)
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
  if(any_gaussian) {

    list2env(control$gaussian, envir = environment())

    means <- c(do.call(rbind, lapply(trans[gaussian_names],
                                     FUN = \(x) x[1, ])))
    vars <- c(do.call(rbind, lapply(trans[gaussian_names],
                                    FUN = \(x) x[2, ])))
    log_vars <- c(do.call(rbind, lapply(trans[gaussian_names],
                                        FUN = \(x) x[3, ])))
    stdv <- c(do.call(rbind, lapply(trans[gaussian_names],
                                    FUN = \(x) x[4, ])))

    transforms[[k]] <- list(transform = "exponential",
                            parameters_in = list(log_vars),
                            parameters_out = list(vars))
    k <- k+1L

    y <- as.matrix(patterns[, gaussian_names])
    # transforms[[k]] <- list(transform = "normal",
    #                         parameters_in = list(means, stdv),
    #                         parameters_out = list(trans$loglik[, gauss, ]),
    #                         extra = list(y = y, S = npatterns, J = Jgauss,
    #                                      I = nclasses))

    transforms[[k]] <- list(transform = "normal2",
                            parameters_in = list(means, vars, trans$loglik),
                            parameters_out = list(trans$loglik),
                            extra = list(y = y, S = npatterns, J = Jgauss,
                                         I = nclasses))
    k <- k+1L

    transforms[[k]] <- list(transform = "sqrt_vector",
                            parameters_in = list(vars),
                            parameters_out = list(stdv),
                            extra = list())
    k <- k+1L

  }

  # Conditional item likelihoods (multinomial):
  if(any_multinomial) {

    list2env(control$multinomial, envir = environment())

    has_mvmultinomial <- !is.null(control$mvmultinomial)
    removed_multinomial <- character(0L)

    if(has_mvmultinomial) {
      list2env(control$mvmultinomial, envir = environment())
    }

    multinomial_softmax_names <- multinomial_names[
      !(multinomial_names %in% removed_multinomial)
    ]

    # Softmax transformations for ordinary multinomial items. These directly
    # define marginal item probabilities.
    for(name in multinomial_softmax_names) {

      j <- match(name, multinomial_names)
      log_name <- logmultinomial_names[j]

      for(i in 1:nclasses) {

        transforms[[k]] <- list(transform = "softmax",
                                parameters_in = list(trans[[log_name]][, i]),
                                parameters_out = list(trans[[name]][, i]))
        k <- k+1L

      }

    }

    # Model for the joint probabilities:
    if(has_mvmultinomial) {

      # Softmax transformations for conditional probabilities. For a pair
      # u1 ~~ u2, these objects represent P(u1 | u2 = reference) and
      # P(u2 | u1 = reference). The original trans$u1 and trans$u2 objects are
      # reserved for the marginal probabilities computed from the joint table.
      for(j in seq_len(nremoved_multinomial)) {

        for(i in 1:nclasses) {

          transforms[[k]] <- list(transform = "softmax",
                                  parameters_in = list(trans[[removed_logmultinomial[j]]][, i]),
                                  parameters_out = list(trans[[removed_condmultinomial[j]]][, i]))
          k <- k+1L

        }

      }

      for(j in 1:njoints) {

        joint_name <- paste(indep_pairs_list[[j]], collapse = ".")
        logjoint_name <- paste("log_", joint_name, sep = "")

        mains_names <- sapply(indep_pairs_list[[j]],
                              FUN = \(x) paste("log_", x, "|reference",
                                               sep = ""))
        mains <- trans[mains_names]

        inters_lognames <- apply(indep_pairs[[j]], MARGIN = 1,
                              FUN = \(x) paste("log_", x, sep = "",
                                               collapse = "x"))
        inters <- trans[inters_lognames]

        orderings <- joint_orderings[[j]][, indep_pairs_list[[j]]]
        colnames(orderings) <- mains_names
        orderings2 <- orderings
        colnames(orderings2) <- indep_pairs_list[[j]]

        for(i in 1:nclasses) {

          mains_effects <- lapply(mains, FUN = \(x) x[, i])
          for(l in mains_names) {
            mains_effects[[l]] <- mains_effects[[l]][orderings[, l]]
          }

          inters_effects <- vector("list", length = length(inters_lognames))
          for(l in 1:nrow(indep_pairs[[j]])) {
            rows <- orderings2[, indep_pairs[[j]][l, 1]]
            cols <- orderings2[, indep_pairs[[j]][l, 2]]
            inters_effects[[l]] <- inters[[l]][cbind(rows, cols)]
          }


          pars_out <- trans[[logjoint_name]][, i]
          transforms[[k]] <- list(transform = "sum_vectors",
                                  parameters_in = c(mains_effects,
                                                    inters_effects),
                                  parameters_out = list(pars_out))

          k <- k+1L

          pars_in_softmax <- trans[[logjoint_name]][, i]
          pars_out_softmax <- trans[[joint_name]][, i]
          transforms[[k]] <- list(transform = "softmax",
                                  parameters_in = list(pars_in_softmax),
                                  parameters_out = list(pars_out_softmax))

          k <- k+1L

          # Marginal probabilities for the variables involved in the joint:
          for(lprob in indep_pairs_list[[j]]) {

            group_level <- split(pars_out_softmax, orderings2[, lprob])
            input <- do.call(Map, c(list(f = c), group_level))
            pars_out_softmax[orderings2[, lprob]]
            transforms[[k]] <- list(transform = "sum_vectors",
                                    parameters_in = input,
                                    parameters_out = list(trans[[lprob]][, i]))
            k <- k+1L

          }

        }

      }

      y <- patterns_mvmultinomial_recoded

      transforms[[k]] <- list(transform = "multinomial2",
                              parameters_in = list(trans[colnames(y)],
                                                   trans$loglik),
                              parameters_out = list(trans$loglik),
                              extra = list(y = y, S = npatterns,
                                           J = new_J, I = nclasses,
                                           K = new_K))

      k <- k+1L

    } else {

      y <- as.matrix(patterns[, multinomial_names])
      K <- unlist(lapply(multinomial_factor_levels, FUN = length))

      # transforms[[k]] <- list(transform = "multinomial",
      #                         parameters_in = list(trans[multinomial_names]),
      #                         parameters_out = list(trans$loglik[, multinom, ]),
      #                         extra = list(y = y, S = npatterns, J = Jmulti,
      #                                      I = nclasses, K = K))

      transforms[[k]] <- list(transform = "multinomial2",
                              parameters_in = list(trans[multinomial_names],
                                                   trans$loglik),
                              parameters_out = list(trans$loglik),
                              extra = list(y = y, S = npatterns, J = Jmultinom,
                                           I = nclasses, K = K))

      k <- k+1L

    }

  }

  # Conditional item likelihoods (multivariate gaussian):
  if(any_mvgaussian) {

    list2env(control$mvgaussian, envir = environment())

    sigma_names <- paste("sigma|Class", 1:nclasses, sep = "")
    vars <- unname(do.call(c, lapply(trans[sigma_names], FUN = \(x) diag(x))))

    transforms[[k]] <- list(transform = "exponential",
                            parameters_in = list(trans$logsigma),
                            parameters_out = list(vars))
    k <- k+1L

    y <- as.matrix(patterns[, mvgaussian_names])
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

  if(control$outcomes) {

    estimators[[G]] <- list(estimator = "lca_outcomes",
                            parameters = c("class", "loglik", "distal_beta"),
                            extra = list(S = npatterns,
                                         I = nclasses,
                                         weights = weights,
                                         Y = Y_patterns,
                                         double_names = "lca"))
    G <- G + 1L

  } else {

    estimators[[G]] <- list(estimator = "lca2",
                            parameters = c("class", "loglik"),
                            extra = list(S = npatterns,
                                         I = nclasses,
                                         weights = weights,
                                         double_names = "lca"))
    G <- G + 1L

  }

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
    alpha <- control$penalties$var$alpha
    if(any_gaussian & alpha != 0) {

      Y <- measurement_recoded[, gaussian_names, drop = FALSE]
      # sigma_class <- split(trans$sigma, rep(1:nclasses, each = Jgauss))
      varshat <- apply(Y, MARGIN = 2, FUN = var, na.rm = TRUE)*(nobs-1)/nobs

      for(i in 1:nclasses) {

        vars_by_class <- unlist(lapply(trans[gaussian_names],
                                       FUN = \(x) x[2, i]))
        estimators[[G]] <- list(estimator = "bayesconst3",
                                parameters = list(vars_by_class),
                                extra = list(K = nclasses,
                                             varshat = varshat,
                                             alpha = alpha,
                                             N = nobs,
                                             double_names = paste("vars|Class",
                                                                  i, sep = "")))
        G <- G+1L

      }

    }

    # Bayes Constant for error covariance matrices:
    alpha <- control$penalties$Sigma$alpha
    if(any_mvgaussian & alpha != 0) {

      Y <- measurement_recoded[, mvgaussian_names, drop = FALSE]
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
    if(any_multinomial & alpha != 0) {

      j <- 1L
      for(name in multinomial_names) {
        for(i in 1:nclasses) {

          pihat <- pi_hat_list[[j]][, i]
          probs_by_class <- unlist(lapply(trans[name],
                                          FUN = \(x) x[, i]))

          estimators[[G]] <- list(estimator = "bayesconst2",
                                  parameters = list(probs_by_class),
                                  extra = list(K = nclasses,
                                               pihat = pihat,
                                               alpha = alpha,
                                               N = nobs,
                                               double_names = paste("P(",
                                                                    name,
                                                                    "|Class", i, ")",
                                                                    sep = "")))
          G <- G+1L

        }
        j <- j+1L
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

  penalty_defaults <- list(
    # beta  = list(alpha = 0, lambda = 0, power = 0),
    beta  = list(alpha = 0),
    class = list(alpha = 1),
    prob  = list(alpha = 1),
    var   = list(alpha = 1),
    Sigma = list(alpha = 1)
  )

  if (isFALSE(control$penalties)) {

    control$reg <- FALSE

  } else if (isTRUE(control$penalties)) {

    control$reg <- TRUE
    control$penalties <- penalty_defaults

  } else if (is.list(control$penalties)) {

    control$reg <- TRUE

    # Check that all supplied penalty names are valid
    unknown_penalties <- setdiff(names(control$penalties), names(penalty_defaults))

    if (length(unknown_penalties) > 0L) {
      stop(
        "Unknown penalty name(s): ",
        paste(unknown_penalties, collapse = ", ")
      )
    }

    # Fill missing penalty objects, while keeping user-specified values
    control$penalties <- utils::modifyList(
      penalty_defaults,
      control$penalties
    )

  } else {

    stop("penalties should be TRUE, FALSE, or a named list")

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

  if(is.null(control$opt)) {
    control$opt <- "lbfgs"
  } else if(control$opt == "newton") {
    if(is.null(control$tcg_maxit)) {
      control$tcg_maxit <- 100L
    } else if(control$tcg_maxit < 0L) {
      stop("tcg_maxit must be a positive integer")
    }
  }

  if(is.null(control$rstarts)) {
    control$rstarts <- 16L
  } else if(control$rstarts < 0L) {
    stop("rstarts must be a positive integer")
  }

  return(control)

}

make_design_matrix <- function(X, data) {

  nobs <- nrow(data)
  make_X_matrix <- function(X_df) {

    # Transform characters into factors:
    char_cols <- vapply(X_df, is.character, logical(1L))
    X_df[char_cols] <- lapply(X_df[char_cols], factor)

    # if (anyNA(X_df)) {
    #   stop(
    #     "Missing values were found in the covariates. ",
    #     "Please use imputation or complete-case data."
    #   )
    # }

    # Create the design matrix:
    mf <- model.frame(~ . + 1, data = X_df, na.action = na.pass)
    X_mat  <- model.matrix(~ . + 1, data = mf)
    # X_mat <- model.matrix(~ . + 1, X_df)

    # Center the variables:
    # X_mat[, -1] <- scale(X_mat[, -1], center = TRUE, scale = FALSE)

    # Put an underscore between the variable names and their level names:
    fac_cols <- names(X_df)[vapply(X_df, is.factor, logical(1L))]

    for (v in fac_cols) {

      # Use assign attribute to identify the columns created by each term.
      # This is safer than startsWith(), because variable names may overlap.
      term_id <- match(v, attr(terms(~ . + 1, data = X_df), "term.labels"))
      i <- which(attr(X_mat, "assign") == term_id)

      if (length(i) > 0L) {
        old_names <- colnames(X_mat)[i]
        level_names <- sub(paste0("^", v), "", old_names)
        colnames(X_mat)[i] <- paste0(v, "_", level_names)
      }
    }

    X_mat
  }

  if (is.null(X)) {

    # Create just the intercept column:
    X <- matrix(
      1,
      nrow = nobs,
      ncol = 1L,
      dimnames = list(NULL, "(Intercept)")
    )

  } else {

    if (is.character(X)) {

      X <- data[, X, drop = FALSE]

    } else {

      # Check that X is either a data.frame or a matrix:
      if (!is.data.frame(X) && !is.matrix(X)) {
        stop("X must be a character vector, matrix or data.frame")
      }

      if (nrow(X) != nobs) {
        stop("Number of cases in the data and covariates does not match")
      }
    }

    X <- make_X_matrix(as.data.frame(X))
  }

  X
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

group_connected_pairs <- function(pairs) {

  if (!is.data.frame(pairs) || ncol(pairs) < 2L) {
    stop("`pairs` must be a data.frame with at least two columns.")
  }

  lhs <- as.character(pairs[[1L]])
  rhs <- as.character(pairs[[2L]])

  if (anyNA(lhs) || anyNA(rhs)) {
    stop("`pairs` cannot contain NA values.")
  }

  # Preserve order of first appearance
  elements <- unique(c(rbind(lhs, rhs)))

  # Create adjacency list
  adj <- setNames(vector("list", length(elements)), elements)

  for (i in seq_along(lhs)) {
    a <- lhs[i]
    b <- rhs[i]

    adj[[a]] <- unique(c(adj[[a]], b))
    adj[[b]] <- unique(c(adj[[b]], a))
  }

  visited <- setNames(rep(FALSE, length(elements)), elements)

  element_groups <- list()
  pair_groups <- list()

  for (el in elements) {

    if (visited[[el]]) next

    stack <- el
    group <- character()

    while (length(stack) > 0L) {

      current <- stack[[1L]]
      stack <- stack[-1L]

      if (visited[[current]]) next

      visited[[current]] <- TRUE
      group <- c(group, current)

      neighbours <- adj[[current]]
      neighbours <- neighbours[!visited[neighbours]]

      stack <- c(stack, neighbours)
    }

    # Pairs involved in this connected component
    pair_idx <- lhs %in% group | rhs %in% group

    element_groups[[length(element_groups) + 1L]] <- group
    pair_groups[[length(pair_groups) + 1L]] <- pairs[pair_idx, , drop = FALSE]
  }

  names(element_groups) <- paste0("group", seq_along(element_groups))
  names(pair_groups) <- paste0("group", seq_along(pair_groups))

  list(
    elements = element_groups,
    pairs = pair_groups
  )
}
