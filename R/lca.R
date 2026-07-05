# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 05/07/2026
#'
#' Latent Class Analysis
#'
#' Estimate latent class models with continuous and categorical indicators, with
#' optional covariates and distal outcomes.
#'
#' @description
#' \code{lca()} estimates latent class models in which the measurement model may
#' contain Gaussian indicators, multinomial indicators, or both. Class membership
#' probabilities can be unconditional or modeled as a softmax transformation
#' of observed covariates. Distal outcomes can also be included and modeled
#' conditionally on the latent classes. The function supports one-step
#' estimation, Bakk--Kuha two-step estimation, ML three-step estimation with
#' modal or proportional classification, and class enumeration by passing several
#' values to \code{nclasses}.
#'
#' @usage
#' lca(data, nclasses = 1L, gaussian = NULL, multinomial = NULL,
#'     covariates = NULL, outcomes = NULL, penalties = TRUE,
#'     model = NULL, weights = NULL, start = NULL, adjustment = "bk",
#'     classification = "modal", control = NULL, do.fit = TRUE,
#'     verbose = TRUE)
#'
#' @param data A \code{data.frame} containing the observed variables. Variables
#'   listed in \code{gaussian} and/or \code{multinomial} are used as indicators.
#'   Variables listed in \code{covariates} are used to predict class membership.
#'   Variables listed in \code{outcomes} are treated as distal outcomes.
#' @param nclasses A positive integer giving the number of latent classes. If a
#'   vector of positive integers is supplied, one model is fitted for each value
#'   and the result is returned as an object of class \code{"llcalist"}. This
#'   is usually termed as "class enumeration".
#' @param gaussian Optional character vector with the names of variables in
#'   \code{data} to be modeled as Gaussian indicators. These variables are
#'   modeled with class-specific means and variances. If residual covariance
#'   terms involving Gaussian indicators are specified in \code{model}, the
#'   corresponding variables are modeled jointly with a class-specific
#'   multivariate Gaussian distribution.
#' @param multinomial Optional character vector with the names of variables in
#'   \code{data} to be modeled as multinomial indicators. These variables are
#'   converted to factors internally. The first factor level is used as the
#'   reference category for the softmax transformation. Variables with more than 30
#'   observed categories are not allowed because estimated probabilities would be
#'   unreliable.
#' @param covariates Optional character vector with the names of variables in
#'   \code{data} used to predict latent class membership. Internally, an intercept
#'   is always included. Character covariates are converted to factors and factors
#'   are dummy-coded using \code{model.matrix()}. Observations with missing values
#'   in the covariates are removed before estimation.
#' @param outcomes Optional distal outcomes. It can be either a character vector
#'   with variable names in \code{data}, or a named list specifying the likelihood
#'   for each outcome. If a character vector is supplied, numeric variables are
#'   treated as Gaussian outcomes and non-numeric variables as multinomial
#'   outcomes. If a list is supplied, supported names are \code{gaussian} and
#'   \code{multinomial}. For example, \code{outcomes = c("y1", "y2", "z")}
#'   automatically assigns numeric outcomes to the Gaussian likelihood and
#'   non-numeric outcomes to the multinomial likelihood;
#'   \code{outcomes = list(gaussian = c("y1", "y2"))} adds Gaussian distal
#'   outcomes; and \code{outcomes = list(multinomial = "z")} adds a multinomial
#'   distal outcome. It is possible to predict different types of outcomes at a
#'   time using \code{outcomes = list(gaussian = c("y1", "y2"), multinomial = "z")}
#' @param penalties Logical value or named list controlling regularization. If
#'   \code{FALSE}, no penalties are used. If \code{TRUE}, default penalties are
#'   used. If a named list is supplied, missing penalty blocks are filled with
#'   their defaults. Valid penalty blocks are \code{beta}, \code{class},
#'   \code{prob}, \code{var}, and \code{Sigma}. The default is
#'   \code{list(beta = list(alpha = 0), class = list(alpha = 1),
#'   prob = list(alpha = 1), var = list(alpha = 1),
#'   Sigma = list(alpha = 1))}. The \code{beta} block may also contain
#'   \code{lambda} and \code{power} for ridge-type regularization.
#' @param model Optional model specification. This can be a named list used to
#'   fix parameters, impose equality constraints, or provide custom parameter
#'   labels. Names should match internal parameter blocks, for example
#'   \code{beta}, Gaussian item names, \code{log_<item>} for multinomial logs,
#'   \code{means}, \code{logsigma}, or \code{sigma|Class<i>} for multivariate
#'   Gaussian blocks. Character strings containing residual-dependency syntax
#'   such as \code{"y1 ~~ y2"} or \code{"u1 ~~ u2 ~~ u3"} are also used to
#'   identify residual covariances or residual associations. If an object of
#'   class \code{"llca"} is supplied, its measurement parameters are reused
#'   while the class-membership regression coefficients are re-estimated.
#' @param weights Optional numeric vector of observation weights, with one value
#'   per row of \code{data}.
#' @param start Optional named list of starting values. Names should correspond
#'   to parameter blocks in the model. Supplied values replace the corresponding
#'   default initial values, allowing partial specification of starting values.
#' @param adjustment Character string selecting the estimation strategy when
#'   \code{covariates} are supplied. Use \code{"none"} for one-step estimation,
#'   \code{"bk"} for the Bakk--Kuha two-step method, or \code{"ml"} for the
#'   ML three-step correction based on classification error. Defaults to
#'   \code{"bk"}.
#' @param classification Character string used when \code{adjustment = "ml"}.
#'   Use \code{"modal"} for modal assignment or \code{"prop"} for proportional
#'   assignment when estimating the classification-error matrix.
#' @param control Optional list of optimizer and estimation controls. Common
#'   entries include \code{rstarts} for the number of random starts,
#'   \code{cores} for parallel computation, \code{maxit} for the maximum number
#'   of optimizer iterations, \code{opt} for the optimizer type, and convergence
#'   tolerances such as \code{eps}, \code{df_eps}, and \code{step_eps}. Missing
#'   entries are replaced by internal defaults.
#' @param do.fit Logical. If \code{TRUE}, the model is estimated. If
#'   \code{FALSE}, the function returns an \code{"llca"} object containing the
#'   processed data, model structure, and optimization setup, but without running
#'   the optimizer.
#' @param verbose Logical. If \code{TRUE}, progress information is printed,
#'   especially when fitting several values of \code{nclasses}.
#'
#' @details
#' The measurement model is defined by the variables supplied through
#' \code{gaussian} and \code{multinomial}. Variables not listed there are ignored
#' by the measurement model unless they are used as covariates or distal
#' outcomes.
#'
#' For Gaussian indicators, the model estimates class-specific means and
#' variances. Variances are parameterized through log-variances and transformed
#' to the positive scale during optimization.
#'
#' For multinomial indicators, the model estimates class-specific category
#' probabilities through a softmax parameterization. The first category is fixed
#' to zero on the log scale for identification.
#'
#' When \code{covariates} are supplied and \code{adjustment = "none"}, class
#' probabilities are modeled in one step through multinomial-log coefficients
#' stored in the \code{beta} parameter block. The first class is the reference
#' class, so its coefficients are fixed to zero.
#'
#' When \code{adjustment = "bk"}, the measurement model is first estimated
#' without covariates. The structural model is then estimated with the
#' measurement parameters fixed. When \code{adjustment = "ml"}, the measurement
#' model is first estimated, cases are assigned to classes using either modal or
#' proportional classification, and the structural model is estimated using a
#' classification-error correction.
#'
#' The function compresses the data into unique response/covariate patterns and
#' uses their frequencies as pattern weights. If observation weights are supplied
#' through \code{weights}, they are incorporated into the pattern weights.
#'
#' Residual dependencies can be requested with covariance-style syntax in
#' \code{model}. For Gaussian indicators, terms such as \code{"y1 ~~ y2"}
#' define class-specific residual covariance parameters. For multinomial
#' indicators, dependency syntax defines joint categorical response blocks. For
#' example, \code{"u1 ~~ u2"} replaces the locally independent likelihood
#' contribution of \code{u1} and \code{u2} by a joint multinomial contribution
#' for \code{u1.u2}. Similarly, a connected dependency such as
#' \code{"u1 ~~ u2 ~~ u3"} defines one joint variable \code{u1.u2.u3}. Pairwise
#' interaction parameters are used to build the joint log-probabilities, and the
#' marginal item probabilities are recovered from the joint probabilities by
#' summing over the relevant joint levels.
#'
#' @return
#' If \code{length(nclasses) == 1} and \code{adjustment = "none"}, an S4 object
#' of class \code{"llca"} with the following slots:
#' \describe{
#'   \item{\code{version}}{Version of the \pkg{latent} package used to fit the model.}
#'   \item{\code{call}}{The matched function call.}
#'   \item{\code{timing}}{Elapsed optimization time. Empty when \code{do.fit = FALSE}.}
#'   \item{\code{dataList}}{Processed data objects, including measurement data,
#'   covariate design matrix, unique response patterns, pattern weights, and
#'   mappings between original rows and unique patterns.}
#'   \item{\code{modelInfo}}{Internal model structure used by the optimizer,
#'   including parameter labels, transformed-parameter labels, degrees of
#'   freedom, manifolds, transformations, estimators, and optimizer controls.}
#'   \item{\code{Optim}}{Raw output from the optimizer. Empty when
#'   \code{do.fit = FALSE}.}
#'   \item{\code{parameters}}{Estimated model parameters on the estimation scale,
#'   organized by parameter block. Empty when \code{do.fit = FALSE}.}
#'   \item{\code{transformed_pars}}{Estimated parameters after applying model
#'   transformations, including class probabilities, item probabilities,
#'   variances, standard deviations, joint probabilities, and log-likelihood
#'   components. Empty when \code{do.fit = FALSE}.}
#'   \item{\code{extra}}{Additional information reserved for downstream methods.}
#' }
#'
#' If \code{adjustment = "bk"} or \code{adjustment = "ml"}, the function returns
#' a list with two elements:
#' \describe{
#'   \item{\code{measurement}}{The measurement-model fit from the first step.}
#'   \item{\code{structural}}{The structural-model fit with covariates and/or
#'   distal outcomes.}
#' }
#'
#' If \code{nclasses} contains several values, a list of fitted objects is
#' returned with class \code{"llcalist"}.
#'
#' @examples
#' \dontrun{
#' # Three-class model for categorical indicators
#' fit <- lca(
#'   data = gss82,
#'   nclasses = 3L,
#'   multinomial = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"),
#'   penalties = TRUE
#' )
#'
#' fit
#' summary(fit)
#' getfit(fit)
#'
#' latInspect(fit, what = "coefs", digits = 3)
#' latInspect(fit, what = "profile", digits = 3)
#' latInspect(fit, what = "posterior", digits = 3)
#'
#' # Class enumeration
#' fits <- lca(
#'   data = gss82,
#'   nclasses = 1:4,
#'   multinomial = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"),
#'   penalties = TRUE
#' )
#'
#' # Latent class regression with Bakk--Kuha adjustment
#' fit_bk <- lca(
#'   data = empathy,
#'   nclasses = 4L,
#'   gaussian = c("ec1", "ec2", "ec3", "ec4", "ec5", "ec6"),
#'   covariates = c("sex", "pt1", "pt2", "pt3", "pt4"),
#'   outcomes = list(gaussian = c("pt5", "pt6")),
#'   adjustment = "bk"
#' )
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
                covariates = NULL,
                outcomes = NULL,
                model = NULL,
                weights = NULL,
                adjustment = "bk",
                classification = "modal",
                penalties = TRUE,
                start = NULL,
                control = NULL,
                do.fit = TRUE,
                verbose = TRUE) {

  #### Check input arguments ####

  # Check weights if available:
  if(!is.null(weights)) {

    if(!is.numeric(weights) ||
       any(weights < 0) ||
       !all(is.finite(weights))) {
      stop("weights must be a positive numeric vector.")
    }

    if(length(c(weights)) != nrow(data)) {
      stop("weights must be the same length than the number of rows in data.")
    }

  }

  # Check that data is a data.frame:
  if(!is.data.frame(data)) {
    stop("data must be a data.frame")
  }

  # Check that nclasses is a (vector of) positive integer(s):
  if(any(nclasses < 1L) ||
     !all(nclasses == as.integer(nclasses))) {
    stop("nclasses must be a (vector of) positive integer(s)")
  }

  # Reject empty measurement models:
  if(length(c(gaussian, multinomial, unlist(outcomes))) == 0L) {
    stop("At least one indicator or outcome variable must be supplied.")
  }

  # Check that the named indicators, covariates, and outcomes exist in data:
  # Indicators:
  indicators <- c(multinomial, gaussian)
  missing_indicators <- setdiff(indicators, colnames(data))
  if(length(missing_indicators) > 0) {
    stop(
      "The following indicators are missing in `data`: ",
      paste(missing_indicators, collapse = ", ")
    )
  }


  # Covariates:
  if (!is.null(covariates) && !is.character(covariates)) {
    stop("covariates must be NULL or a character vector with variable names in data.")
  }
  missing_covariates <- setdiff(covariates, colnames(data))
  if(length(missing_covariates) > 0) {
    stop(
      "The following covariates are missing in `data`: ",
      paste(missing_covariates, collapse = ", ")
    )
  }

  # Outcomes:
  missing_outcomes <- setdiff(unlist(outcomes), colnames(data))
  if(length(missing_outcomes) > 0) {
    stop(
      "The following outcomes are missing in `data`: ",
      paste(missing_outcomes, collapse = ", ")
    )
  }

  there_are_covariates <- !is.null(covariates)
  there_are_outcomes <- !is.null(outcomes)
  structural <- there_are_covariates || there_are_outcomes

  # Check overlapping variable roles

  outcome_names <- if(is.null(outcomes)) {
    character(0L)
  } else if(is.character(outcomes)) {
    outcomes
  } else {
    unlist(outcomes, use.names = FALSE)
  }

  all_roles <- c(gaussian, multinomial, covariates, outcome_names)
  all_roles <- all_roles[!is.na(all_roles)]

  duplicated_roles <- unique(all_roles[duplicated(all_roles)])

  if(length(duplicated_roles) > 0L) {
    stop(
      "Duplicated variable(s) found as indicator, covariate, or outcome: ",
      paste(duplicated_roles, collapse = ", ")
    )
  }

  # Check that control is either NULL or a list:
  if(is.null(control)) {
    control <- list()
  } else if(!is.list(control)) {
    stop("control must be NULL or a list.")
  }

  #### Class enumeration: run a model for different number of classes ####

  # Return a list of llca models:

  if(length(nclasses) > 1L) {

    result <- vector("list", length = length(nclasses))
    for(i in 1:length(nclasses)) {

      if(verbose) print(paste0("Model nclasses=", nclasses[i]))

      result[[i]] <- lca(data, nclasses = nclasses[i],
                         gaussian = gaussian, multinomial = multinomial,
                         covariates = covariates, outcomes = outcomes,
                         model = model, weights = weights,
                         adjustment = adjustment, classification = classification,
                         penalties = penalties, start = start,
                         control = control, do.fit = do.fit, verbose = verbose)

      names(result)[i] <- paste("nclasses=", nclasses[i], sep = "")

    }

    class(result) <- "llcalist"

    return(result)

  }

  #### Adjustment methods to fit models with covariates and distal outcomes ####

  adjustment <- tolower(adjustment) # lowercase

  # Return a list of lcca models:

  if(structural && adjustment != "none") {

    if(adjustment == "bk") { # Bakk and Kuha adjustment

      result <- lca_bakk_kuha(data = data,
                              nclasses = nclasses,
                              gaussian = gaussian,
                              multinomial = multinomial,
                              covariates = covariates,
                              outcomes = outcomes,
                              penalties = penalties,
                              model = model,
                              weights = weights,
                              start = start,
                              control = control,
                              do.fit = do.fit,
                              verbose = verbose)

    } else if(adjustment == "ml") { # ML LatentGold adjustment

      result <- lca_ml(data = data,
                       nclasses = nclasses,
                       gaussian = gaussian,
                       multinomial = multinomial,
                       covariates = covariates,
                       outcomes = outcomes,
                       penalties = penalties,
                       model = model,
                       weights = weights,
                       start = start,
                       control = control,
                       do.fit = do.fit,
                       verbose = verbose,
                       classification = classification)

    } else {
      stop("Unknown adjustment method")
    }

    return(result)

  }

  # If no structural part or adjustment method is used, then proceed...

  #### Store original call ####

  mc <- match.call(expand.dots = TRUE)
  args <- lapply(as.list(mc)[-1], eval, envir = parent.frame())

  #### Create the dataList ####

  dataList <- create_lca_dataList(data = data,
                                  nclasses = nclasses,
                                  covariates = covariates,
                                  outcomes = outcomes,
                                  gaussian = gaussian,
                                  multinomial = multinomial,
                                  model = model,
                                  weights = weights)
  list2env(dataList, envir = environment())
  dataList$args <- args

  #### Check control parameters ####

  # Update and check control parameters:
  control$penalties <- penalties # Either a logical or a list of named penalties
  control$start <- start # Optional named list with initial values
  control <- lca_control(control) # Check and update the control inputs and create defaults

  #### Create the model ####

  # Fix the measurement part of the model if a previous llca fit was provided in
  # the model argument. This is usually done for two or three-step estimation:

  control$model <- model

  if(inherits(model, "llca")) {

    model <- model@parameters
    model$beta <- NULL # Free the covariate coefficients

  } else if(is.list(model)) {

    # Find what objects in model are llca objects:
    llca_obj <- which(vapply(model, inherits, logical(1), what = "llca"))

    # Collect all the @parameters from the llca objects. If some of them are
    # repeated (i.e., share the same name), then keep the last one:
    if(length(llca_obj) > 0L) {
      model <- model[-llca_obj]
      for (obj in control$model[llca_obj]) {
        if (inherits(obj, "llca")) {
          model[names(obj@parameters)] <- obj@parameters
        }
      }
      model$beta <- NULL # Free the covariate coefficients
    }

  } else if(!is.null(model)) {
    stop("model must be a llca object or a list containing parameter blocks
         or llca objects")
  }
  # Here, the covariate coefficients are always freed if any llca object is
  # identified in the model argument

  # Get the model specification:
  full_model <- create_lca_model(dataList = dataList,
                                 nclasses = nclasses,
                                 model = model,
                                 control = control)
  list2env(full_model, envir = environment())

  #### Create the manifold, transformation, and estimator structures ####

  # Generate the structures for optimization:
  modelInfo <- create_lca_modelInfo(dataList = dataList,
                                    full_model = full_model,
                                    control = control)

  #### Fit the model ####

  if(!do.fit) { # Just get the model specification (empty model)

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

#### Function to create the dataList and modelInfo ####

create_lca_dataList <- function(data = NULL,
                                nclasses = NULL,
                                covariates = NULL,
                                outcomes = NULL,
                                gaussian = NULL,
                                multinomial = NULL,
                                model = NULL,
                                weights = NULL) {

  #### Out covariates, indicators, and outcomes names ####

  indicators_names <- intersect(c(multinomial, gaussian), colnames(data))
  covariates_names <- intersect(covariates, colnames(data))
  outcomes_names <- intersect(unlist(outcomes), colnames(data))
  variables <- c(indicators_names, covariates_names, outcomes_names)

  # Check for all missing variables:
  all_missing <- vapply(
    data[, variables, drop = FALSE],
    FUN = function(x) all(is.na(x)),
    FUN.VALUE = logical(1L)
  )
  if (any(all_missing)) {
    stop(
      "The following variable(s) have only missing values: ",
      paste(names(all_missing)[all_missing], collapse = ", ")
    )
  }

  #### Add each outcome to a likelihood model ####

  if(is.list(outcomes)) {
    if(is.null(names(outcomes)) || any(names(outcomes) == "")) {
      stop("All elements in outcomes must be named.")
    }
    # Check for unknown likelihoods:
    unknown <- setdiff(names(outcomes), c("gaussian", "multinomial"))
    if(length(unknown) > 0L) {
      stop(
        "Unknown outcome likelihood(s): ",
        paste(unknown, collapse = ", ")
      )
    }
  }

  # If outcomes is a character vector, transform it into a named likelihood list
  if(is.character(outcomes)) {

    out <- data[, outcomes, drop = FALSE]
    numeric_vars <- vapply(out, is.numeric, logical(1))
    outcomes <- list(
      gaussian = outcomes[numeric_vars],
      multinomial = outcomes[!numeric_vars]
    )

  }

  # Now outcomes is either NULL or a named list of likelihoods
  if(is.null(outcomes)) {
    outcomes$any_outcomes <- FALSE
    outcomes$outcomes_names <- NULL
  } else {
    outcomes$any_outcomes <- TRUE
    if(!is.null(outcomes$multinomial)) {
      multinomial <- c(multinomial, outcomes$multinomial)
      outcomes$outcomes_names <- c(outcomes$outcomes_names, outcomes$multinomial)
    }
    if(!is.null(outcomes$gaussian)) {
      gaussian <- c(gaussian, outcomes$gaussian)
      outcomes$outcomes_names <- c(outcomes$outcomes_names, outcomes$gaussian)
    }

    outcomes$outcomes_names <- intersect(outcomes$outcomes_names, colnames(data))

  }

  #### Create the design matrix from the covariates ####

  # Create the matrix of predictors for latent class probabilities:
  # covariates must be a character vector with the name of the predictors in data:
  design <- make_design_matrix(covariates = covariates_names,
                               data = data)
  pdesign <- ncol(design) # Number of columns of the design matrix

  # Remove participants with missing data in any covariate:
  remove <- which(apply(design, MARGIN = 1, FUN = anyNA))
  if(length(remove) > 0) {

    missing <- colSums(is.na(design))
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
      if(length(remove) == 1L) " was" else "s were",
      " removed."
    )

    design <- design[-remove, , drop = FALSE]
    data <- data[-remove, , drop = FALSE]
    weights <- weights[-remove]

  }

  if(nrow(data) == 0L) {
    stop("No observations remain after removing rows with missing covariates.")
  }

  # Check for all missing variables again:
  all_missing <- vapply(
    data[, variables, drop = FALSE],
    FUN = function(x) all(is.na(x)),
    FUN.VALUE = logical(1L)
  )
  if(any(all_missing)) {
    stop(
      "The following variable(s) have only missing values after removing rows
      with missing covariates: ",
      paste(names(all_missing)[all_missing], collapse = ", ")
    )
  }

  # Number of subjects after removing rows with missingness in the covariates:
  nobs <- nrow(data)

  #### Find the unique data patterns and compute the weights ####

  any_continuous <- FALSE # If any variable is continuous, don't find patterns
  if(length(gaussian) > 0L) any_continuous <- TRUE
  patterns <- make_patterns(data = data[, variables, drop = FALSE],
                            weights = weights,
                            any_continuous = any_continuous)

  design_patterns <- design[patterns$full2short, , drop = FALSE]

  #### Process the indicators ####

  # If at least two items are gaussian, check which of them are involved in
  # covariance structures in the model argument so they are modeled with a
  # multivariate normal distribution. Create also a target matrix where TRUE
  # means that a given covariance is estimated and FALSE means that it is fixed
  # to zero. The diagonal of target should be FALSE because variances are
  # not estimated directly but parameterized as exp(log variances):
  gaussian_objects <- check_mvgaussian(gaussian, model)
  # check_mvgaussian outputs gaussian, mvgaussian, and target:
  list2env(gaussian_objects, envir = environment())

  multinomial <- objects_multinomial_lca(data, patterns, multinomial, nclasses)
  gaussian <- objects_gaussian_lca(data, patterns, gaussian, nclasses)
  mvgaussian <- objects_mvgaussian_lca(data, patterns, mvgaussian, nclasses,
                                       target)

  list2env(multinomial, envir = environment())
  list2env(gaussian, envir = environment())
  list2env(mvgaussian, envir = environment())

  # # Count the number of items for each likelihood model:
  # nmultinomial <- multinomial$nmultinomial
  # ngaussian <- gaussian$ngaussian
  # nmvgaussian <- mvgaussian$nmvgaussian

  # Select the subset of data with the variables that are used in the
  # measurement model (keeping the original ordering):
  model_labels <- rep(c("multinomial", "gaussian", "mvgaussian"),
                      times = c(nmultinomial, ngaussian, nmvgaussian))
  model_vector <- c(multinomial_names, gaussian_names, mvgaussian_names)
  idx <- match(colnames(data), model_vector); idx <- idx[!is.na(idx)]
  # match each item with its respective likelihood model:
  variable_type <- model_labels[idx]
  names(variable_type) <- model_vector[idx]

  #### Possible response patterns and degrees of freedom ####

  # Get the number of possible response patterns:
  if(all(variable_type == "multinomial")) { # If all the items are multinomial...

    # If all the items are multinomial, the number of possible patterns is:
    npossible_patterns <- min(prod(multinomial_factor_lengths)-1L, nobs)
    # From technical manual of LatentGOLD 6.1 (p.68)

  } else { # If any item is not multinomial...

    npossible_patterns <- nobs

  }

  #### Check for residual dependencies in multinomial items ####

  mvmultinomial <- check_mvmultinomial(multinomial, model, patterns)

  #### Store objects in dataList ####

  dataList <- vector("list")
  dataList$data <- data
  dataList$design <- design
  dataList$variable_type <- variable_type
  dataList$nobs <- nobs
  dataList$cov_patterns <- design_patterns
  dataList$npatterns <- patterns$npatterns
  dataList$npossible_patterns <- npossible_patterns
  dataList$pattern_names <- patterns$pattern_names
  dataList$pdesign <- pdesign
  dataList$pattern_weights <- patterns$pattern_weights
  dataList$full2short <- patterns$full2short
  dataList$short2full <- patterns$short2full
  dataList$gaussian <- gaussian
  dataList$mvgaussian <- mvgaussian
  dataList$multinomial <- multinomial
  dataList$mvmultinomial <- mvmultinomial
  dataList$indicators_names <- indicators_names
  dataList$covariates_names <- covariates_names
  dataList$outcomes_names <- outcomes_names
  dataList$variables <- variables
  dataList$outcomes <- outcomes
  dataList$variable_type <- variable_type

  #### Return ####

  return(dataList)

}

#### Function to create the model structures ####

create_lca_model <- function(dataList, nclasses, model = NULL, control) {

  # Generate the model syntax and initial parameter values

  list2env(dataList, envir = environment())

  class_names <- paste("Class", 1:nclasses, sep = "")
  control$class_names <- class_names

  #### Model for the transformed parameters ####

  pred_names <- colnames(design) # Names of predictors

  list_struct <- vector("list")
  k <- 1L

  # Model for the betas:

  list_struct[[k]] <- list(name = "beta",
                           type = "matrix",
                           dim = c(pdesign, nclasses),
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

  list_struct <- c(list_struct,
                   model_gaussian_lca(nclasses, dataList),
                   model_mvgaussian_lca(nclasses, dataList),
                   model_multinomial_lca(nclasses, dataList))

  # Create the full transparameter structure:
  trans <- create_parameters(list_struct)

  #### Model for the parameters ####

  param <- list()

  # Model for the betas:
  param$beta <- trans$beta
  param$beta[, 1] <- "0"

  # Model for multinomial items:
  if(dataList$multinomial$any_multinomial) {

    list2env(dataList$multinomial, envir = environment())
    multinomial_param_names <- multinomial_names
    logmultinomial_param_names <- logmultinomial_names

    if(!is.null(dataList$mvmultinomial)) {

      list2env(dataList$mvmultinomial, envir = environment())

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

    if(!is.null(dataList$mvmultinomial)) {

      # Conditional logs for variables involved in residual dependencies.
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

  # Model for gaussian items:
  if(dataList$gaussian$any_gaussian) {

    list2env(dataList$gaussian, envir = environment())
    cl <- length(param)
    idx_gauss <- (cl + 1L):(cl + ngaussian)
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

  # Model for multivariate gaussian items:
  if(dataList$mvgaussian$any_mvgaussian) {

    list2env(dataList$mvgaussian, envir = environment())
    param$means <- trans$means
    param$logsigma <- trans$logsigma
    # param[sigma_name] <- trans[sigma_names]
    sigma_names <- paste("sigma|Class", 1:nclasses, sep = "")
    ordering <- rownames(trans[[sigma_names[1]]])
    target <- target[ordering, ordering]
    for(j in 1:nclasses) {
      param[[sigma_names[j]]] <- diag(nmvgaussian)
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

  init_gauss <- start_gaussian_lca(param, dataList, control$rstarts)
  init_mvgauss <- start_mvgaussian_lca(param, dataList, control$rstarts)
  init_multinom <- start_multinomial_lca(param, dataList, control$rstarts)

  for(i in 1:control$rstarts) {

    init_param[[i]] <- list()

    # Initial values for betas:
    init_param[[i]][["beta"]] <- matrix(rnorm(pdesign*nclasses),
                                        nrow = pdesign, ncol = nclasses)
    init_param[[i]]$beta[, 1] <- 0
    dimnames(init_param[[i]]$beta) <- dimnames(param$beta)

    init_param[[i]] <- c(init_param[[i]],
                         init_multinom[[i]],
                         init_gauss[[i]],
                         init_mvgauss[[i]])

    # Preserve the same order than in param:
    init_param[[i]] <- init_param[[i]][names(param)]

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
                 control = control)

  return(result)

}

#### Function to create the modelInfo ####

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

  # Conditional likelihoods:
  transforms <- c(transforms,
                  transformations_multinomial_lca(trans, dataList),
                  transformations_gaussian_lca(trans, dataList),
                  transformations_mvgaussian_lca(trans, dataList))

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
                                       weights = pattern_weights,
                                       double_names = "lca"))
  G <- G + 1L

  # Choose whether using Bayes constants:
  if(control$reg) {

    # Bayes Constant for class probabilities:
    alpha <- control$penalties$class$alpha
    estimators_bayes1 <- bayes1(trans, nclasses, alpha, nobs, cov_patterns)

    # Bayes Constant for standard deviations:
    alpha <- control$penalties$var$alpha
    estimators_bayes3 <- bayes3(trans, gaussian, alpha, nclasses, nobs)

    # Bayes Constant for error covariance matrices:
    alpha <- control$penalties$Sigma$alpha
    estimators_bayes4 <- bayes4(trans, mvgaussian, alpha, nclasses, nobs)

    # Bayes Constant for multinomial probabilities:
    alpha <- control$penalties$prob$alpha
    estimators_bayes2 <- bayes2(trans, multinomial, nclasses, alpha, nobs)

    # Gaussian regularization for coefficients:
    alpha <- control$penalties$beta$alpha
    estimators_gaussloglik <- gaussloglik(trans, alpha, nobs, cov_patterns)

    # Ridge regularization for coefficients:
    lambda <- control$penalties$beta$lambda
    power <- control$penalties$beta$power
    estimators_ridge <- ridge(trans, lambda, power, nobs)

    # add the regularizations:
    estimators <- c(estimators,
                    estimators_bayes1, estimators_bayes2,
                    estimators_bayes3, estimators_bayes4,
                    estimators_gaussloglik, estimators_ridge)

  }

  control_estimator <- create_estimators(estimators = estimators,
                                         structures = trans)

  #### Pass the initial values to vectors ####

  inits <- create_init(trans, param, init_param,
                       control_transform = control_transform, control)

  parameters <- inits$parameters
  parameters_labels <- inits$parameters_labels
  nparam <- inits$nparam

  transparameters <- inits$transparameters
  transparameters_labels <- inits$transparameters_labels
  ntrans <- inits$ntrans

  trans2param <- inits$trans2param

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

#### Function to create the control list of optimization parameters ####

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
    control$cores <- max(1L, parallel::detectCores()-1L)
  } else if(control$cores < 1L) {
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
  } else if(control$rstarts < 1L ||
            !all(control$rstarts == as.integer(control$rstarts))) {
    stop("rstarts must be a positive integer")
  }

  return(control)

}

#### Auxiliary functions for create_lca_dataList ####

make_design_matrix <- function(covariates, data) {

  nobs <- nrow(data)

  if(is.null(covariates) || length(covariates) == 0L) {

    X <- matrix(
      1,
      nrow = nobs,
      ncol = 1L,
      dimnames = list(NULL, "(Intercept)")
    )

  } else {

    if(!is.character(covariates)) {
      stop("covariates must be NULL or a character vector with variable names in data.")
    }

    X_df <- data[, covariates, drop = FALSE]

    # Transform characters into factors:
    char_cols <- vapply(X_df, is.character, logical(1L))
    X_df[char_cols] <- lapply(X_df[char_cols], factor)

    # Create design matrix:
    X <- model.matrix(~ . + 1, X_df)

    # Rename factor dummy columns:
    fac_cols <- names(X_df)[vapply(X_df, is.factor, logical(1L))]

    for(v in fac_cols) {

      term_id <- match(v, attr(terms(~ . + 1, data = X_df), "term.labels"))
      i <- which(attr(X, "assign") == term_id)

      if (length(i) > 0L) {
        old_names <- colnames(X)[i]
        level_names <- sub(paste0("^", v), "", old_names)
        colnames(X)[i] <- paste0(v, "_", level_names)
      }
    }
  }

  return(X)

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

  # Remove duplicated/reversed pairs:
  pairs <- t(apply(pairs, 1L, sort))
  colnames(pairs) <- c("lhs", "rhs")
  pairs <- unique(pairs)

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
    pair_idx <- lhs %in% group & rhs %in% group

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

check_mvgaussian <- function(gaussian, model) {

  # If at least two items are gaussian, check which of them are involved in
  # covariance structures in the model argument so they are modeled with a
  # multivariate normal distribution. Create also a target matrix where TRUE
  # means that a given covariance is estimated and FALSE means that it is fixed
  # to zero. The diagonal of target should be FALSE because variances are
  # not estimated directly but parameterized as exp(log variances):

  ngaussian <- length(gaussian)
  mvgaussian <- target <- NULL

  if(ngaussian > 1L) {

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

        if (length(candidate_mvgaussian) > 1L) {

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

  result <- list(gaussian = gaussian,
                 mvgaussian = mvgaussian,
                 target = target)

  return(result)

}

check_mvmultinomial <- function(multinomial, model, patterns) {

  mvmultinomial <- NULL

  if(multinomial$any_multinomial) {

    list2env(multinomial, envir = environment())

    # Define the joint probability parameters:
    error_covs <- NULL

    if (!is.null(model)) {
      error_covs <- extract_cov_pairs(model)
    }

    if (!is.null(error_covs)) {

      error_covs$lhs <- as.character(error_covs$lhs)
      error_covs$rhs <- as.character(error_covs$rhs)

      # Keep only pairs where both variables are multinomial:
      keep <- error_covs$lhs %in% multinomial_names &
        error_covs$rhs %in% multinomial_names &
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

        groups <- group_connected_pairs(error_covs)
        indep_pairs <- groups$pairs
        indep_pairs_list <- groups$elements
        joint_names <- vapply(indep_pairs_list, paste, character(1), collapse = ".")
        names(indep_pairs) <- names(indep_pairs_list) <- joint_names
        patterns_original <- patterns$patterns # dataList$patterns_original
        joints <- lapply(indep_pairs_list, FUN = \(x) {
          apply(patterns_original[, x], MARGIN = 1, FUN = paste, collapse = ".")
        })
        joints <- as.data.frame(joints)
        new_probs <- unlist(c(keep, colnames(joints)))
        new_J <- length(new_probs)
        patterns_mvmultinomial <- data.frame(patterns_original,
                                             joints)[, new_probs, drop = FALSE]
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
        multinomial_factor_levels_numeric <-
          setNames(lapply(multinomial_factor_lengths, seq_len),
                   names(multinomial_factor_lengths))
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

        mvmultinomial <- list(
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

  return(mvmultinomial)

}

make_patterns <- function(data, weights = NULL, any_continuous) {

  nobs <- nrow(data)

  if (any_continuous) {

    npatterns <- nobs
    pattern_names <- paste("pattern", seq_len(npatterns), sep = "")

    pattern_weights <- if (is.null(weights)) {
      rep(1, nobs)
    } else {
      weights
    }

    result <- list(
      patterns = data,
      npatterns = npatterns,
      pattern_names = pattern_names,
      pattern_weights = pattern_weights,
      full2short = seq_len(npatterns),
      short2full = seq_len(npatterns)
    )

    return(result)
  }

  dt <- data.table::as.data.table(data)
  data_names <- colnames(data)

  if (is.null(weights)) {

    counts_dt <- dt[
      ,
      .(
        index = .I[1],
        pattern_weights = .N
      ),
      by = data_names
    ]

  } else {

    dt[, ".weights" := weights]

    counts_dt <- dt[
      ,
      .(
        index = .I[1],
        pattern_weights = sum(.weights)
      ),
      by = data_names
    ]
  }

  npatterns <- nrow(counts_dt)
  pattern_names <- paste("pattern", seq_len(npatterns), sep = "")

  patterns <- as.matrix(counts_dt[, data_names, with = FALSE])
  rownames(patterns) <- pattern_names

  pattern_weights <- counts_dt$pattern_weights

  full2short <- counts_dt$index

  short2full <- match(
    do.call(paste, dt[, data_names, with = FALSE]),
    do.call(paste, counts_dt[, data_names, with = FALSE])
  )

  result <- list(
    patterns = patterns,
    npatterns = npatterns,
    pattern_names = pattern_names,
    pattern_weights = pattern_weights,
    full2short = full2short,
    short2full = short2full
  )

  return(result)
}

#### Structural after measurement methods ####

lca_bakk_kuha <- function(data,
                          nclasses = NULL,
                          gaussian = NULL,
                          multinomial = NULL,
                          covariates = NULL,
                          outcomes = NULL,
                          penalties = TRUE,
                          model = NULL,
                          weights = NULL,
                          start = NULL,
                          control = NULL,
                          do.fit = TRUE,
                          verbose = TRUE) {

  # This is the Bakk and Kuha method for estimating covariates:
  # https://bpspsychub.onlinelibrary.wiley.com/doi/10.1111/bmsp.12227

  #### Step 1: Measurement model ####

  fit1 <- lca(data = data,
              nclasses = nclasses,
              gaussian = gaussian,
              multinomial = multinomial,
              covariates = NULL,
              outcomes = NULL,
              penalties = penalties,
              model = model,
              weights = weights,
              start = start,
              adjustment = "none",
              control = control,
              do.fit = do.fit,
              verbose = verbose)

  #### Step 2: Fitting the structural part and fixing the measurement part ####

  fit2 <- lca(data = data,
              nclasses = nclasses,
              gaussian = gaussian,
              multinomial = multinomial,
              covariates = covariates,
              outcomes = outcomes,
              penalties = penalties,
              model = c(model, fit1),
              weights = weights,
              start = start,
              adjustment = "none",
              control = control,
              do.fit = do.fit,
              verbose = verbose)

  #### Return ####

  result <- list(measurement = fit1, structural = fit2)
  class(result) <- "llcalist"

  return(result)

}

lca_ml <- function(data,
                   nclasses = NULL,
                   gaussian = NULL,
                   multinomial = NULL,
                   covariates = NULL,
                   outcomes = NULL,
                   penalties = TRUE,
                   model = NULL,
                   weights = NULL,
                   start = NULL,
                   control = NULL,
                   do.fit = TRUE,
                   verbose = TRUE,
                   classification = "modal") {

  # This code runs the three-steps method of LatentGold with options:
  # Analysis: Covariates
  # Classification: Modal / Proportional
  # Adjustment: ML

  #### Step 1: Measurement model ####

  fit1 <- lca(data = data,
              nclasses = nclasses,
              gaussian = gaussian,
              multinomial = multinomial,
              covariates = NULL,
              outcomes = NULL,
              penalties = penalties,
              model = model,
              weights = weights,
              start = start,
              adjustment = "none",
              control = control,
              do.fit = do.fit,
              verbose = verbose)

  #### Step 2: Modal assignment ####

  if(classification == "modal") {

    data$states <- latInspect(fit1, what = "state")
    class_error <- latInspect(fit1, what = "classification")$class_error_modal

  } else if(classification == "prop") {

    N <- nrow(data)
    data <- data[rep(seq_len(N), each = nclasses), ] # Expand the dataset
    data$states <- factor(rep(seq_len(nclasses), times = N),
                                   levels = seq_len(nclasses))
    posterior <- latInspect(fit1, what = "posterior")
    weights <- as.vector(t(posterior)) # Original weights are overwritten
    class_error <- latInspect(fit1, what = "classification")$class_error_prop

  } else {
    stop("Unknown classification method")
  }

  log_class_error <- log(class_error)
  log_class_error <- apply(log_class_error, MARGIN = 1L,
                           FUN = \(x) x - x[1L])
  # t(apply(log_class_error, 2, soft, a=1)) / class_error

  #### Step 3: Fitting the structural part using the states ####

  fit2 <- lca(data = data,
              nclasses = nclasses,
              multinomial = "states",
              covariates = covariates,
              outcomes = outcomes,
              model = list(model, log_states = log_class_error),
              weights = weights,
              penalties = penalties,
              start = start,
              adjustment = "none",
              control = control,
              do.fit = do.fit,
              verbose = verbose)

  #### Return ####

  result <- list(measurement = fit1, structural = fit2)
  class(result) <- "llcalist"

  return(result)

}

#### Objects for items with a given conditional likelihood ####

objects_multinomial_lca <- function(data, patterns, multinomial, nclasses) {

  nmultinomial <- length(multinomial) # Number of multinomial items
  any_multinomial <- ifelse(nmultinomial > 0L, TRUE, FALSE)
  multinomial_names <- multinomial
  multinomial <- list(any_multinomial = any_multinomial,
                      nmultinomial = nmultinomial,
                      multinomial_names = multinomial_names)

  nobs <- nrow(data)

  # Transform multinomial variables into factors:

  if(any_multinomial) {

    logmultinomial_names <- paste("log", multinomial_names, sep = "_")

    # Transform multinomial variables into factors and drop unused levels:
    measurement <- lapply(data[, multinomial_names, drop = FALSE],
                          FUN = function(x) droplevels(factor(x)))

    # Save levels of factor variables:
    multinomial_factor_levels <- lapply(measurement, levels)
    multinomial_factor_lengths <- unlist(lapply(multinomial_factor_levels, FUN = length))
    multinomial_reference <- lapply(multinomial_factor_levels, FUN = \(x) x[1])

    problem <- multinomial_factor_lengths < 2L
    if(any(problem)) {
      stop("Multinomial variable(s) ",
           paste0("'", multinomial_names[problem], "'", collapse = ", "),
           " have fewer than two observed categories.")
    }

    # Transform the factor variables into integers:
    measurement <- as.data.frame(lapply(measurement, FUN = function(col) {
      if (is.factor(col)) as.integer(col) - 1L else col
    }))

    multinomial <- list(multinomial_names = multinomial_names,
                        logmultinomial_names = logmultinomial_names,
                        nmultinomial = nmultinomial,
                        multinomial_factor_levels = multinomial_factor_levels,
                        multinomial_factor_lengths = multinomial_factor_lengths,
                        multinomial_reference = multinomial_reference)
    # Check which multinomial items are ordered:
    idx <- unlist(lapply(data[, multinomial_names, drop = FALSE], FUN = is.ordered))
    ordered_multinomial <- multinomial_names[idx]

    # For starting values:
    pi_hat_list <- list()
    eta_hat_list <- list()
    j <- 1L
    for(name in multinomial_names) {

      int_vector <- measurement[, name]
      int_vector <- int_vector[!is.na(int_vector)] # Remove missing values
      nsize <- length(int_vector)
      props <- count(int_vector, nsize, multinomial_factor_lengths[j]) / nsize
      pi_hat_list[[j]] <- props %*% t(rep(1, nclasses))
      log_props <- log(props)
      eta_hat_list[[j]] <- (log_props-log_props[1]) %*% t(rep(1, nclasses))
      j <- j+1L

    }

    multinomial$eta_hat_list <- eta_hat_list
    multinomial$pi_hat_list <- pi_hat_list
    multinomial$pi_hat_vector <- pi_hat_vector <- unlist(pi_hat_list)
    multinomial$vars <- (1-pi_hat_vector)/(nobs*pi_hat_vector)
    multinomial$sds <- sqrt(multinomial$vars)
    multinomial$Ks <- length(unlist(eta_hat_list))
    multinomial$multinomial_names <- multinomial_names
    multinomial$patterns_multinomial <- as.matrix(
      measurement[patterns$full2short, multinomial_names, drop = FALSE])

  }
  multinomial$any_multinomial <- any_multinomial

  return(multinomial)

}

objects_gaussian_lca <- function(data, patterns, gaussian, nclasses) {

  ngaussian <- length(gaussian) # Number of gaussian items
  any_gaussian <- ifelse(ngaussian > 0L, TRUE, FALSE)
  gaussian_names <- gaussian
  gaussian <- list(any_gaussian = any_gaussian,
                   ngaussian = ngaussian,
                   gaussian_names = gaussian_names)

  if(any_gaussian) {

    init_mean <- init_var <- init_logvar <- list()

    j <- 1L
    for(name in gaussian_names) {

      init_mean[[j]] <- rep(mean(data[, name], na.rm = TRUE), times = nclasses)
      v <- stats::var(data[, name], na.rm = TRUE)
      if (!is.finite(v) || v <= 0) {
        stop("Gaussian variable '", name,
             "' has zero, undefined, or non-finite variance.")
      }
      init_var[[j]] <- rep(v, times = nclasses)
      init_logvar[[j]] <- log(init_var[[j]])
      j <- j+1L

    }

    gaussian$init_mean <- init_mean
    gaussian$init_var <- init_var
    gaussian$init_logvar <- init_logvar
    gaussian$patterns_gaussian <- as.matrix(
      data[patterns$full2short, gaussian_names, drop = FALSE])

  }

  return(gaussian)

}

objects_mvgaussian_lca <- function(data, patterns, mvgaussian, nclasses, target) {

  nmvgaussian <- length(mvgaussian) # Number of mvgaussian items
  any_mvgaussian <- ifelse(nmvgaussian > 0L, TRUE, FALSE)
  mvgaussian_names <- mvgaussian
  mvgaussian <- list(any_mvgaussian = any_mvgaussian,
                     nmvgaussian = nmvgaussian,
                     mvgaussian_names = mvgaussian_names,
                     target = target)

  if(any_mvgaussian) {

    init_sigma <- list()
    sigma_names <- paste("sigma|Class", 1:nclasses, sep = "")
    init_mvmean <- matrix(rep(colMeans(data[, mvgaussian_names, drop = FALSE],
                                       na.rm = TRUE),
                              times = nclasses), nrow = nmvgaussian,
                          ncol = nclasses)

    v <- apply(data[, mvgaussian_names, drop = FALSE], MARGIN = 2,
               FUN = stats::var,
               na.rm = TRUE)

    if(any(!is.finite(v) | v <= 0)) {
      stop("Multivariate Gaussian variable(s) have zero, undefined, or non-finite variance: ",
           paste(names(v)[!is.finite(v) | v <= 0], collapse = ", ")
      )
    }

    init_logsigma <- matrix(rep(log(v), times = nclasses),
                            nrow = nmvgaussian, ncol = nclasses )

    for(j in 1:nclasses) {
      init_sigma[[j]] <- cov(data[, mvgaussian_names, drop = FALSE],
                             use = "pairwise.complete.obs")
    }

    mvgaussian$init_mvmean <- init_mvmean
    mvgaussian$init_logsigma <- init_logsigma
    mvgaussian$init_sigma <- init_sigma
    mvgaussian$patterns_mvgaussian <- as.matrix(
      data[patterns$full2short, mvgaussian_names, drop = FALSE])

  }

  return(mvgaussian)

}

#### Functions for conditional likelihood models structures ####

model_multinomial_lca <- function(nclasses, dataList) {

  list2env(dataList, envir = environment())

  class_names <- paste("Class", 1:nclasses, sep = "")
  list_struct <- list()
  k <- 1L

  # Model for multinomial items:
  list2env(multinomial, envir = environment())
  if(any_multinomial) {

    if(any(multinomial_factor_lengths > 30)) {
      stop("You cannot use a multinomial likelihood to model a variable with more than 30 unique categories")
    }

    has_mvmultinomial <- !is.null(dataList$mvmultinomial)
    removed_multinomial <- character(0L)

    if(has_mvmultinomial) {
      list2env(dataList$mvmultinomial, envir = environment())
    }

    for(j in seq_len(nmultinomial)) {

      # For ordinary multinomial items, the item probability matrix is obtained
      # directly from the item logs with a softmax transformation. For items
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

  return(list_struct)

}

model_gaussian_lca <- function(nclasses, dataList) {

  list2env(dataList, envir = environment())

  class_names <- paste("Class", 1:nclasses, sep = "")
  list_struct <- list()
  k <- 1L

  # Model for gaussian items:
  list2env(gaussian, envir = environment())
  if(any_gaussian) {

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

  return(list_struct)

}

model_mvgaussian_lca <- function(nclasses, dataList) {

  list2env(dataList, envir = environment())

  class_names <- paste("Class", 1:nclasses, sep = "")
  sigma_names <- paste("sigma|Class", 1:nclasses, sep = "")
  list_struct <- list()
  k <- 1L

  # Model for multivariate normal items:
  list2env(mvgaussian, envir = environment())
  if(any_mvgaussian) {

    list_struct[[k]] <- list(name = "means",
                             type = "matrix",
                             dim = c(nmvgaussian, nclasses),
                             rownames = mvgaussian_names,
                             colnames = class_names)
    k <- k+1L

    list_struct[[k]] <- list(name = "logsigma",
                             type = "matrix",
                             dim = c(nmvgaussian, nclasses),
                             rownames = mvgaussian_names,
                             colnames = class_names)
    k <- k+1L

    # Sigma:
    for(j in 1:nclasses) {

      list_struct[[k]] <- list(name = sigma_names[j],
                               type = "matrix",
                               dim = c(nmvgaussian, nmvgaussian),
                               rownames = mvgaussian_names,
                               colnames = mvgaussian_names,
                               symmetric = TRUE)
      k <- k+1L

    }

  }

  return(list_struct)

}

#### Functions for starting values of conditional likelihood models ####

start_multinomial_lca <- function(param, dataList, rstarts) {

  list2env(dataList, envir = environment())
  list2env(multinomial, envir = environment())

  init_param <- vector("list", length = rstarts)

  for(i in 1:rstarts) {

    init_param[[i]] <- list()

    # Initial values for multinomial items:
    if(any_multinomial) {

      multinomial_init_names <- multinomial_names
      if(!is.null(dataList$mvmultinomial)) {
        list2env(mvmultinomial, envir = environment())
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

      if(!is.null(dataList$mvmultinomial)) {

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
          dims <- dim(param[[loginter_names[[j]]]])
          init_param[[i]][[loginter_names[[j]]]] <- matrix(0, nrow = dims[1],
                                                           ncol = dims[2])
          dimnames(init_param[[i]][[loginter_names[[j]]]]) <-
            dimnames(param[[loginter_names[[j]]]])
        }

      }

    }

  }

  return(init_param)

}

start_gaussian_lca <- function(param, dataList, rstarts) {

  list2env(dataList, envir = environment())
  list2env(gaussian, envir = environment())

  nclasses <- ncol(param$beta)
  init_param <- vector("list", length = rstarts)

  for(i in 1:rstarts) {

    init_param[[i]] <- list()
    # Initial values for gaussian items:
    j <- 1L

    if(any_gaussian) {

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

  }

  return(init_param)

}

start_mvgaussian_lca <- function(param, dataList, rstarts) {

  list2env(dataList, envir = environment())
  list2env(mvgaussian, envir = environment())

  nclasses <- ncol(param$beta)
  sigma_names <- paste("sigma|Class", 1:nclasses, sep = "")
  init_param <- vector("list", length = rstarts)

  for(i in 1:rstarts) {

    init_param[[i]] <- list()
    j <- 1L

    # Initial values for multivariate gaussian items:
    if(any_mvgaussian) {

      # TO-DO: ALLOW RANDOM STARTING VALUES

      init_param[[i]][["means"]] <- init_mvmean
      init_param[[i]][["logsigma"]] <- init_logsigma
      dimnames(init_param[[i]][["means"]]) <-
        dimnames(init_param[[i]][["logsigma"]]) <- dimnames(param[["means"]])

      for(j in 1:nclasses) {

        init_param[[i]][[sigma_names[j]]] <- init_sigma[[j]]
        dimnames(init_param[[i]][[sigma_names[j]]]) <- dimnames(param[[sigma_names[j]]])

      }

    }

  }

  return(init_param)

}

#### Functions for conditional likelihood transformations ####

transformations_multinomial_lca <- function(trans, dataList) {

  list2env(dataList, envir = environment())

  nclasses <- ncol(trans$beta)
  transforms <- list()
  k <- 1L

  # Conditional item likelihoods (multinomial):
  list2env(multinomial, envir = environment())
  if(multinomial$any_multinomial) {

    list2env(multinomial, envir = environment())

    has_mvmultinomial <- !is.null(dataList$mvmultinomial)
    removed_multinomial <- character(0L)

    if(has_mvmultinomial) {
      list2env(dataList$mvmultinomial, envir = environment())
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
            transforms[[k]] <- list(transform = "sum_vectors",
                                    parameters_in = input,
                                    parameters_out = list(trans[[lprob]][, i]))
            k <- k+1L

          }

        }

      }

      transforms[[k]] <- list(transform = "multinomial2",
                              parameters_in = list(trans[colnames(patterns_mvmultinomial_recoded)],
                                                   trans$loglik),
                              parameters_out = list(trans$loglik),
                              extra = list(y = patterns_mvmultinomial_recoded,
                                           S = npatterns, J = new_J,
                                           I = nclasses, K = new_K))

      k <- k+1L

    } else {

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
                              extra = list(y = patterns_multinomial, S = npatterns,
                                           J = nmultinomial, I = nclasses, K = K))

      k <- k+1L

    }

  }

  return(transforms)

}

transformations_gaussian_lca <- function(trans, dataList) {

  list2env(dataList, envir = environment())

  nclasses <- ncol(trans$beta)
  transforms <- list()
  k <- 1L

  # Conditional item likelihoods (gaussian):
  list2env(gaussian, envir = environment())
  if(any_gaussian) {

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

    # transforms[[k]] <- list(transform = "normal",
    #                         parameters_in = list(means, stdv),
    #                         parameters_out = list(trans$loglik[, gauss, ]),
    #                         extra = list(y = y, S = npatterns, J = Jgauss,
    #                                      I = nclasses))

    transforms[[k]] <- list(transform = "normal2",
                            parameters_in = list(means, vars, trans$loglik),
                            parameters_out = list(trans$loglik),
                            extra = list(y = patterns_gaussian, S = npatterns,
                                         J = ngaussian, I = nclasses))
    k <- k+1L

    transforms[[k]] <- list(transform = "sqrt_vector",
                            parameters_in = list(vars),
                            parameters_out = list(stdv),
                            extra = list())
    k <- k+1L

  }

  return(transforms)

}

transformations_mvgaussian_lca <- function(trans, dataList) {

  list2env(dataList, envir = environment())

  nclasses <- ncol(trans$beta)
  sigma_names <- paste("sigma|Class", 1:nclasses, sep = "")
  transforms <- list()
  k <- 1L

  # Conditional item likelihoods (multivariate gaussian):
  list2env(mvgaussian, envir = environment())
  if(any_mvgaussian) {

    vars <- unname(do.call(c, lapply(trans[sigma_names], FUN = \(x) diag(x))))

    transforms[[k]] <- list(transform = "exponential",
                            parameters_in = list(trans$logsigma),
                            parameters_out = list(vars))
    k <- k+1L

    transforms[[k]] <- list(transform = "mvnormal2",
                            parameters_in = list(unlist(trans$means),
                                                 unlist(trans[sigma_names]),
                                                 trans$loglik),
                            parameters_out = list(trans$loglik),
                            extra = list(y = patterns_mvgaussian, S = npatterns,
                                         J = nmvgaussian, I = nclasses))
    # k <- k+1L

  }

  return(transforms)

}

#### Functions for regularization ####

bayes1 <- function(trans, nclasses, alpha, nobs, cov_patterns) {

  estimators <- list()

  if(alpha != 0) {

    # Get the indices corresponding to the unique covariate patterns:
    dt_uniq_X <- data.table::as.data.table(cov_patterns)
    counts_X <- dt_uniq_X[, .(index = .I[1], count = .N), by = names(dt_uniq_X)]
    uniques <- counts_X$index
    U <- length(uniques)
    G <- 1L

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

  return(estimators)

}

bayes2 <- function(trans, multinomial, nclasses, alpha, nobs) {

  estimators <- list()

  if(multinomial$any_multinomial & alpha != 0) {

    list2env(multinomial, envir = environment())

    G <- j <- 1L
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

  return(estimators)

}

bayes3 <- function(trans, gaussian, alpha, nclasses, nobs) {

  estimators <- list()

  if(gaussian$any_gaussian & alpha != 0) {

    list2env(gaussian, envir = environment())
    Y <- patterns_gaussian[, gaussian_names, drop = FALSE]
    varshat <- apply(Y, MARGIN = 2, FUN = var, na.rm = TRUE)*(nobs-1)/nobs
    G <- 1L

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

  return(estimators)

}

bayes4 <- function(trans, mvgaussian, alpha, nclasses, nobs) {

  estimators <- list()

  if(mvgaussian$any_mvgaussian & alpha != 0) {

    list2env(mvgaussian, envir = environment())
    Y <- patterns_mvgaussian[, mvgaussian_names, drop = FALSE]
    D <- diag(apply(Y, MARGIN = 2, FUN = var, na.rm = TRUE)*(nobs-1)/nobs)
    sigma_names <- paste("sigma|Class", 1:nclasses, sep = "")
    G <- 1L

    for(i in 1:nclasses) {

      estimators[[G]] <- list(estimator = "bayesconst4",
                              parameters = list(trans[[sigma_names[i]]]),
                              extra = list(K = nclasses,
                                           alpha = alpha,
                                           D = D,
                                           J = nmvgaussian,
                                           double_names = paste("Sigma|Class",
                                                                i, sep = "")))
      G <- G+1L

    }

  }

  return(estimators)

}

ridge <- function(trans, lambda, power, nobs) {

  estimators <- list()

  if(!is.null(lambda) & !is.null(power)) {

    if(lambda != 0 && power != 0) {

      estimators[[1]] <- list(estimator = "ridge",
                              parameters = list(trans$beta[-1, ]),
                              extra = list(lambda = lambda,
                                           power = power,
                                           N = nobs,
                                           double_names = paste("Ridge(", power,
                                                                ")_coeffs",
                                                                sep = "")))

    }
  }

  return(estimators)

}

gaussloglik <- function(trans, alpha, nobs, cov_patterns) {

  estimators <- list()

  if(alpha != 0) {

    p <- nrow(trans$beta)-1L
    q <- ncol(trans$beta)
    means <- matrix(0, nrow = p, ncol = q)
    sds <- apply(cov_patterns[, -1], MARGIN = 2, sd, na.rm = TRUE)
    sds <- matrix(sds, nrow = p, ncol = q) / alpha

    estimators[[1]] <- list(estimator = "gaussian_loglik",
                            parameters = list(trans$beta[-1, ]),
                            extra = list(means = means,
                                         sds = sds,
                                         alpha = alpha,
                                         N = nobs,
                                         double_names = "Gaussian_coeffs"))

  }

  return(estimators)

}
