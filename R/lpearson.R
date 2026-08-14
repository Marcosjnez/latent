# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 13/08/2026
#'
#' Pearson Covariance or Correlation Matrix
#'
#' Estimate a covariance or correlation matrix and its asymptotic covariance
#' matrix from continuous data.
#'
#' @usage
#' lpearson(data, model = NULL, std.ov = FALSE,
#'          acov = "standard", likelihood = "normal",
#'          missing = "pairwise.complete.obs", do.fit = TRUE,
#'          message = FALSE, control = NULL, ...)
#'
#' @param data A data frame or matrix containing numeric observed variables.
#' @param model Optional model specification reserved for internal model setup.
#' @param std.ov Logical. If \code{TRUE}, the observed variables are standardized
#'   before the sample matrix is returned, so the off-diagonal elements represent
#'   correlations rather than covariances.
#' @param acov Character string selecting the asymptotic covariance estimator.
#'   Available options are \code{"standard"} and \code{"robust"}.
#' @param likelihood Character string controlling the covariance denominator.
#'   Use \code{"normal"} for the maximum-likelihood denominator \eqn{N} and
#'   \code{"wishart"} for the usual sample-covariance denominator \eqn{N-1}.
#' @param missing Character string passed to \code{stats::cov()} to control the
#'   handling of missing values. \code{"fiml"} is treated as
#'   \code{"pairwise.complete.obs"} because FIML is handled at the CFA level.
#' @param do.fit Logical. If \code{FALSE}, return the prepared but unfitted
#'   \code{"latent"} object.
#' @param message Logical. Print progress messages during estimation.
#' @param control Optional list of internal controls. Custom starting values can
#'   be supplied through \code{control$start}.
#' @param ... Additional arguments reserved for future extensions.
#'
#' @details
#' \code{lpearson()} computes the empirical covariance matrix directly rather
#' than estimating it numerically. With \code{std.ov = TRUE}, the sample matrix
#' is standardized before being stored in the fitted object.
#'
#' Standard asymptotic covariance matrices are computed under multivariate
#' normality. With \code{acov = "robust"}, fourth-moment information from the
#' observed data is used instead.
#'
#' @return An S4 object of class \code{"latent"}. The object contains the
#'   processed data in \code{dataList}, the parameter and optimization structures
#'   in \code{modelInfo}, the direct estimation output and asymptotic covariance
#'   information in \code{Optim}, and the estimated parameter structures in
#'   \code{parameters} and \code{transformed_pars}.
#'
#' @examples
#' \dontrun{
#' fit <- lpearson(data = HolzingerSwineford1939[, paste0("x", 1:9)])
#' fit_cor <- lpearson(data = HolzingerSwineford1939[, paste0("x", 1:9)],
#'                     std.ov = TRUE)
#' }
#'
#' @export
lpearson <- function(data,
                     model = NULL,
                     std.ov = FALSE,
                     acov = "standard",
                     likelihood = "normal",
                     missing = "pairwise.complete.obs",
                     do.fit = TRUE,
                     message = FALSE,
                     control = NULL,
                     ...) {

  #### Check input arguments ####

  if(!is.data.frame(data) && !is.matrix(data)) {
    stop("data must be a data.frame or matrix")
  }

  if(is.matrix(data)) {
    data <- as.data.frame(data)
  }

  if(nrow(data) == 0L || ncol(data) == 0L) {
    stop("data must contain at least one observation and one variable")
  }

  if(is.null(colnames(data)) ||
     any(colnames(data) == "") ||
     anyDuplicated(colnames(data))) {
    stop("data must have unique, non-empty column names")
  }

  numeric_data <- vapply(data, is.numeric, logical(1L))
  if(!all(numeric_data)) {
    stop("All variables in data must be numeric")
  }

  if(!is.logical(std.ov) || length(std.ov) != 1L || is.na(std.ov)) {
    stop("std.ov must be TRUE or FALSE")
  }

  if(!is.logical(do.fit) || length(do.fit) != 1L || is.na(do.fit)) {
    stop("do.fit must be TRUE or FALSE")
  }

  if(!is.logical(message) || length(message) != 1L || is.na(message)) {
    stop("message must be TRUE or FALSE")
  }

  acov <- match.arg(tolower(acov), c("standard", "robust"))
  likelihood <- match.arg(tolower(likelihood), c("normal", "wishart"))

  missing <- tolower(missing)
  if(missing == "fiml") {
    missing <- "pairwise.complete.obs"
  }

  missing_methods <- c("everything", "all.obs", "complete.obs",
                       "na.or.complete", "pairwise.complete.obs")

  if(!missing %in% missing_methods) {
    stop("Unknown missing-data method")
  }

  if(!is.null(control) && !is.list(control)) {
    stop("control must be NULL or a list")
  }

  if(is.null(control)) {
    control <- list()
  }

  #### Store original call ####

  mc <- match.call(expand.dots = TRUE)
  args <- lapply(as.list(mc)[-1L], eval, envir = parent.frame())

  #### Check control parameters ####

  control$std.ov <- std.ov
  control$acov <- acov
  control$likelihood <- likelihood
  control$missing <- missing
  control <- lpearson_control(control)

  #### Create the dataList ####

  dataList <- create_lpearson_dataList(data = data,
                                       control = control)
  dataList$args <- args

  #### Create the model ####

  full_model <- create_lpearson_model(dataList = dataList,
                                      model = model,
                                      control = control)

  #### Create the manifold, transformation, and estimator structures ####

  modelInfo <- create_lpearson_modelInfo(dataList = dataList,
                                         full_model = full_model,
                                         control = control)

  #### Fit the model ####

  if(!do.fit) {

    result <- new("latent",
                  version          = as.character(packageVersion("latent")),
                  call             = mc,
                  timing           = numeric(),
                  dataList         = dataList,
                  modelInfo        = modelInfo,
                  Optim            = list(),
                  parameters       = list(),
                  transformed_pars = list(),
                  extra            = list())

    #### Result ####

    return(result)

  }

  if(message) {
    print_lpearson_message("Fitting the model")
  }

  Optim <- fit_lpearson(dataList = dataList,
                        modelInfo = modelInfo)

  #### Process the outputs ####

  transformed_pars <- fill_in(modelInfo$trans,
                              Optim$transparameters)

  parameters <- transformed_pars[names(modelInfo$param)]

  #### Standard errors ####

  Optim$SE <- compute_se_lpearson(dataList = dataList,
                                  modelInfo = modelInfo,
                                  Optim = Optim,
                                  control = control)

  #### latent object ####

  result <- new("latent",
                version          = as.character(packageVersion("latent")),
                call             = mc,
                timing           = Optim$elapsed,
                dataList         = dataList,
                modelInfo        = modelInfo,
                Optim            = Optim,
                parameters       = parameters,
                transformed_pars = transformed_pars,
                extra            = list())

  #### Result ####

  return(result)

}

#### Function to create the dataList ####

create_lpearson_dataList <- function(data, control) {

  #### Sample covariance or correlation matrix ####

  sample_stats <- sample_lpearson(data = data,
                                  control = control)

  #### Result ####

  result <- list(data = data,
                 nobs = nrow(data),
                 nitems = ncol(data),
                 npatterns = nrow(data),
                 item_names = colnames(data),
                 S = sample_stats$S,
                 acov = sample_stats$acov,
                 likelihood = sample_stats$likelihood)

  return(result)

}

#### Function to create the model ####

create_lpearson_model <- function(dataList, model = NULL, control) {

  # Generate the model syntax and initial parameter values.

  #### Parameter-block names ####

  data_param <- create_lpearson_data_param(dataList = dataList,
                                            control = control)

  #### Model for the transformed parameters ####

  trans <- model_lpearson(dataList = dataList,
                          data_param = data_param,
                          control = control)

  #### Model for the parameters ####

  param <- constraints_lpearson(dataList = dataList,
                                data_param = data_param,
                                trans = trans,
                                control = control)

  #### Create the initial values for the parameters ####

  init_param <- start_lpearson(dataList = dataList,
                               data_param = data_param,
                               param = param,
                               trans = trans,
                               control = control)

  #### Custom initial values ####

  init_param <- custom_init_param(control$start, init_param)

  #### Result ####

  result <- list(param = param,
                 trans = trans,
                 init_param = init_param,
                 data_param = data_param,
                 control = control)

  return(result)

}

#### Function to create the modelInfo ####

create_lpearson_modelInfo <- function(dataList, full_model, control) {

  # Generate the manifold, transformation, estimator, and optimizer structures.

  list2env(full_model, envir = environment())

  #### Manifolds ####

  manifolds <- manifolds_lpearson(dataList = dataList,
                                  data_param = data_param,
                                  param = param,
                                  control = control)

  control_manifold <- create_manifolds(manifolds = manifolds,
                                       structures = param)

  #### Transformations ####

  transforms <- transformations_lpearson(dataList = dataList,
                                          data_param = data_param,
                                          trans = trans,
                                          control = control)

  control_transform <- create_transforms(transforms = transforms,
                                         structures = trans)

  #### Estimators ####

  estimators <- estimators_lpearson(dataList = dataList,
                                    data_param = data_param,
                                    trans = trans,
                                    control = control)

  control_estimator <- create_estimators(estimators = estimators,
                                         structures = trans)

  #### Pass the initial values to vectors ####

  inits <- create_init(trans, param, init_param,
                       control_transform = control_transform,
                       control = control)
  list2env(inits, envir = environment())

  #### Set up the optimizer ####

  control_optimizer <- control
  control_optimizer$parameters <- parameters
  control_optimizer$transparameters <- transparameters
  control_optimizer$init_param <- init_param
  control_optimizer$transparam2param <- trans2param-1L

  #### Result ####

  modelInfo <- list(param = param,
                    trans = trans,
                    nparam = nparam,
                    ntrans = ntrans,
                    parameters_labels = parameters_labels,
                    transparameters_labels = transparameters_labels,
                    dof = dataList$npatterns-nparam,
                    control_manifold = control_manifold,
                    control_transform = control_transform,
                    control_estimator = control_estimator,
                    control_optimizer = control_optimizer)

  return(modelInfo)

}

#### Function to create the control list of optimization parameters ####

lpearson_control <- function(control) {

  # Auxiliary function for lpearson.R.

  if(is.null(control$subfix)) {
    control$subfix <- ""
  }

  if(is.null(control$opt)) {
    control$opt <- "lbfgs"
  }

  if(is.null(control$step_maxit)) {
    control$step_maxit <- 30L
  } else if(control$step_maxit < 1L) {
    stop("step_maxit must be an integer greater than 0")
  }

  if(is.null(control$c1)) {
    control$c1 <- 0.5
  } else if(control$c1 < 0) {
    stop("c1 must be a positive number")
  }

  if(is.null(control$c2)) {
    control$c2 <- 0.5
  } else if(control$c2 < 0) {
    stop("c2 must be a positive number")
  }

  if(is.null(control$step_eps)) {
    control$step_eps <- 1e-09
  } else if(control$step_eps < 0) {
    stop("step_eps must be a positive number")
  }

  if(is.null(control$df_eps)) {
    control$df_eps <- 1e-09
  } else if(control$df_eps < 0) {
    stop("df_eps must be a positive number")
  }

  if(is.null(control$M)) {
    control$M <- 100L
  } else if(control$M < 1L) {
    stop("M must be a positive integer")
  }

  if(is.null(control$eps)) {
    control$eps <- 1e-04
  } else if(control$eps < 0) {
    stop("eps must be a positive number")
  }

  if(is.null(control$ss_fac)) {
    control$ss_fac <- 2
  } else if(control$ss_fac <= 1) {
    stop("ss_fac must be larger than 1")
  }

  if(is.null(control$maxit)) {
    control$maxit <- 1000L
  } else if(control$maxit < 1L) {
    stop("maxit must be a positive integer")
  }

  if(is.null(control$tcg_maxit)) {
    control$tcg_maxit <- 10L
  } else if(control$tcg_maxit < 1L) {
    stop("tcg_maxit must be a positive integer")
  }

  if(is.null(control$ss)) {
    control$ss <- 0.001
  } else if(control$ss <= 0) {
    stop("ss must be a positive number")
  }

  # lpearson uses direct estimation and therefore a single start and core.
  control$rstarts <- 1L
  control$cores <- 1L

  #### Result ####

  return(control)

}

#### Auxiliary functions for create_lpearson_dataList ####

sample_lpearson <- function(data, control) {

  nobs <- nrow(data)

  if(nobs < 2L) {

    S <- t(as.matrix(data)) %*% as.matrix(data)
    acov <- "standard"
    likelihood <- "wishart"

  } else {

    S <- stats::cov(data, use = control$missing)
    acov <- control$acov
    likelihood <- control$likelihood

  }

  rownames(S) <- colnames(S) <- colnames(data)

  if(control$std.ov) {

    if(any(!is.finite(diag(S))) || any(diag(S) <= 0)) {
      stop("Observed variables cannot be standardized because at least one variance is non-positive or non-finite")
    }

    inv_sqrtdiagS <- diag(1/sqrt(diag(S)), nrow = ncol(S))
    S <- inv_sqrtdiagS %*% S %*% inv_sqrtdiagS
    rownames(S) <- colnames(S) <- colnames(data)

  }

  # Preserve the existing lpearson convention: under the normal likelihood the
  # sample matrix uses N rather than N-1 in the denominator.
  if(likelihood == "normal") {
    S <- S*(nobs-1L)/nobs
  }

  #### Result ####

  result <- list(S = S,
                 acov = acov,
                 likelihood = likelihood)

  return(result)

}

#### Auxiliary functions for create_lpearson_model ####

create_lpearson_data_param <- function(dataList, control) {

  S_matrix <- paste0("S", control$subfix)

  #### Result ####

  result <- list(S_matrix = S_matrix)

  return(result)

}

model_lpearson <- function(dataList, data_param, control) {

  list2env(data_param, envir = environment())

  list_struct <- list(
    list(name = S_matrix,
         type = "matrix",
         dim = c(dataList$nitems, dataList$nitems),
         rownames = dataList$item_names,
         colnames = dataList$item_names,
         symmetric = TRUE)
  )

  #### Result ####

  result <- create_parameters(list_struct)

  return(result)

}

constraints_lpearson <- function(dataList, data_param, trans, control) {

  list2env(data_param, envir = environment())

  param <- trans

  if(control$std.ov) {
    diag(param[[S_matrix]]) <- "1"
  }

  #### Result ####

  return(param)

}

start_lpearson <- function(dataList, data_param, param, trans, control) {

  list2env(data_param, envir = environment())

  init_param <- vector("list", length = control$rstarts)

  for(rs in seq_len(control$rstarts)) {

    init_param[[rs]] <- list()
    init_param[[rs]][[S_matrix]] <- dataList$S
    dimnames(init_param[[rs]][[S_matrix]]) <- dimnames(trans[[S_matrix]])

  }

  #### Result ####

  return(init_param)

}

#### Auxiliary functions for create_lpearson_modelInfo ####

manifolds_lpearson <- function(dataList, data_param, param, control) {

  list2env(data_param, envir = environment())

  manifolds <- list(
    list(manifold = "euclidean",
         parameters = S_matrix)
  )

  #### Result ####

  return(manifolds)

}

transformations_lpearson <- function(dataList, data_param, trans, control) {

  #### Result ####

  return(list())

}

estimators_lpearson <- function(dataList, data_param, trans, control) {

  #### Result ####

  return(list())

}

#### Auxiliary functions for fitting and post-estimation ####

fit_lpearson <- function(dataList, modelInfo) {

  #### Identify the sample-matrix parameter block ####

  S_name <- names(modelInfo$trans)[startsWith(names(modelInfo$trans), "S")]

  if(length(S_name) != 1L) {
    stop("Unable to identify the Pearson covariance or correlation matrix")
  }

  S_labels <- c(modelInfo$trans[[S_name]])
  S_values <- c(dataList$S)

  #### Transformed parameters ####

  idx <- match(modelInfo$transparameters_labels, S_labels)

  if(anyNA(idx)) {
    stop("Unable to match the transformed Pearson parameters")
  }

  transparameters <- S_values[idx]
  names(transparameters) <- modelInfo$transparameters_labels

  #### Free parameters ####

  idx <- match(modelInfo$parameters_labels,
               modelInfo$transparameters_labels)

  if(anyNA(idx)) {
    stop("Unable to match the free Pearson parameters")
  }

  parameters <- transparameters[idx]
  names(parameters) <- modelInfo$parameters_labels

  #### Optimizer-style output ####

  Optim <- list(f = 0,
                elapsed = 0,
                parameters = parameters,
                transparameters = transparameters)

  #### Result ####

  return(Optim)

}

compute_se_lpearson <- function(dataList, modelInfo, Optim, control) {

  #### Asymptotic covariance ####

  if(dataList$acov == "standard") {

    ACOV <- asymptotic_normal(dataList$S,
                              cov = !control$std.ov,
                              diag = FALSE)

  } else if(dataList$acov == "robust") {

    ACOV <- asymptotic_general(as.matrix(dataList$data),
                               cov = !control$std.ov,
                               diag = FALSE)

  } else {

    stop("Unknown acov method")

  }

  if(nrow(ACOV) != length(modelInfo$parameters_labels) ||
     ncol(ACOV) != length(modelInfo$parameters_labels)) {
    stop("The asymptotic covariance matrix does not match the number of free Pearson parameters")
  }

  rownames(ACOV) <- colnames(ACOV) <- modelInfo$parameters_labels

  se <- sqrt(diag(ACOV)/dataList$nobs)
  names(se) <- modelInfo$parameters_labels

  #### Result ####

  result <- list(ACOV = ACOV,
                 se = se)

  return(result)

}

print_lpearson_message <- function(msg) {

  w <- nchar(msg)+4L

  cat("\n", "+", strrep("-", w), "+\n",
      "|  ", msg, "  |\n",
      "+", strrep("-", w), "+\n\n", sep = "")

  #### Result ####

  return(invisible(NULL))

}
