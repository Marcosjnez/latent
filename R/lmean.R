# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 15/08/2026
#'
#' Sample Means
#'
#' Estimate the sample means of a set of observed variables and their
#' asymptotic covariance matrix.
#'
#' @usage
#' lmean(data, model = NULL, std.ov = FALSE,
#'       do.fit = TRUE, message = FALSE,
#'       control = NULL, ...)
#'
#' @param data A data frame or matrix containing numeric observed variables.
#'   Column names are required and are used as the observed-variable labels.
#' @param model Optional model specification reserved for internal model setup.
#' @param std.ov Logical. If \code{TRUE}, the observed variables are treated as
#'   standardized and all mean parameters are fixed to zero. Their asymptotic
#'   covariance matrix and standard errors are consequently also zero.
#' @param do.fit Logical. If \code{TRUE}, compute the sample means. If
#'   \code{FALSE}, return the prepared but unfitted \code{"latent"} object.
#' @param message Logical. Print progress messages during estimation.
#' @param control Optional list of internal control parameters. Custom starting
#'   values can be supplied through \code{control$start}.
#' @param ... Additional arguments reserved for future extensions.
#'
#' @details
#' \code{lmean()} computes the arithmetic mean of each observed variable using
#' all available non-missing observations for that variable. When
#' \code{std.ov = FALSE}, the asymptotic covariance matrix of the sample means
#' is represented by the diagonal matrix of observed-variable variances, and
#' the reported standard errors are obtained by dividing this matrix by the
#' sample size.
#'
#' When \code{std.ov = TRUE}, the mean parameters are fixed to zero rather than
#' estimated. This is the appropriate mean structure for standardized observed
#' variables. Because these values are fixed, their asymptotic covariance
#' matrix and standard errors are set to zero.
#'
#' The function uses the same parameter/model infrastructure as the other
#' estimators in \pkg{latent}, even though the sample means themselves are
#' obtained directly rather than through numerical optimization.
#'
#' @return An S4 object of class \code{"latent"}. The object contains the
#'   processed data in \code{dataList}, the parameter and model structures in
#'   \code{modelInfo}, the estimated means and standard-error information in
#'   \code{Optim}, and the parameter-shaped results in \code{parameters} and
#'   \code{transformed_pars}.
#'
#' @examples
#' \dontrun{
#' fit <- lmean(data = HolzingerSwineford1939[, paste0("x", 1:9)])
#'
#' fit_std <- lmean(data = HolzingerSwineford1939[, paste0("x", 1:9)],
#'                  std.ov = TRUE)
#' }
#'
#' @export
lmean <- function(data,
                  model = NULL,
                  std.ov = FALSE,
                  do.fit = TRUE,
                  message = FALSE,
                  control = NULL,
                  ...) {

  #### Check input arguments ####

  if(!is.data.frame(data) && !is.matrix(data)) {
    stop("data must be a data.frame or matrix")
  }

  if(nrow(data) < 1L || ncol(data) < 1L) {
    stop("data must contain at least one observation and one variable")
  }

  if(is.data.frame(data)) {
    numeric_columns <- vapply(data, is.numeric, logical(1L))
    if(!all(numeric_columns)) {
      stop("All variables in data must be numeric")
    }
  } else if(!is.numeric(data)) {
    stop("data must be numeric")
  }

  if(is.null(colnames(data)) ||
     any(colnames(data) == "") ||
     anyDuplicated(colnames(data))) {
    stop("data must have unique, non-empty column names")
  }

  if(length(std.ov) != 1L || !is.logical(std.ov) || is.na(std.ov)) {
    stop("std.ov must be TRUE or FALSE")
  }

  if(length(do.fit) != 1L || !is.logical(do.fit) || is.na(do.fit)) {
    stop("do.fit must be TRUE or FALSE")
  }

  if(length(message) != 1L || !is.logical(message) || is.na(message)) {
    stop("message must be TRUE or FALSE")
  }

  #### Store original call ####

  mc <- match.call()

  #### Check control parameters ####

  control$std.ov <- std.ov
  control <- lmean_control(control)

  #### Create the dataList ####

  dataList <- create_lmean_dataList(data = data,
                                    control = control)

  #### Create the model ####

  full_model <- create_lmean_model(dataList = dataList,
                                   model = model,
                                   control = control)

  #### Create the manifold, transformation, and estimator structures ####

  modelInfo <- create_lmean_modelInfo(dataList = dataList,
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
    print_lmean_message("Fitting the model")
  }

  Optim <- fit_lmean(dataList = dataList,
                     modelInfo = modelInfo,
                     control = control)

  #### Standard errors ####

  Optim$SE <- compute_se_lmean(dataList = dataList,
                               modelInfo = modelInfo,
                               Optim = Optim,
                               control = control)

  #### Process the outputs ####

  transformed_pars <- fill_in(modelInfo$trans,
                              Optim$transparameters)

  parameters <- transformed_pars[names(modelInfo$param)]

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

create_lmean_dataList <- function(data, control) {

  #### Data information ####

  dataList <- list(
    data = data,
    nobs = nrow(data),
    nitems = ncol(data),
    npatterns = nrow(data),
    item_names = colnames(data)
  )

  #### Result ####

  return(dataList)

}

#### Function to create the model ####

create_lmean_model <- function(dataList, model, control) {

  data_param <- create_lmean_data_param(dataList = dataList,
                                        control = control)

  #### Model for the transformed parameters ####

  trans <- model_lmean(dataList = dataList,
                       data_param = data_param,
                       control = control)

  #### Model for the parameters ####

  param <- constraints_lmean(trans = trans,
                             dataList = dataList,
                             data_param = data_param,
                             control = control)

  #### Create the initial values for the parameters ####

  init_param <- start_lmean(param = param,
                            trans = trans,
                            dataList = dataList,
                            data_param = data_param,
                            control = control)

  #### Custom initial values ####

  init_param <- custom_init_param(control$start, init_param)

  #### Result ####

  result <- list(param = param,
                 trans = trans,
                 init_param = init_param,
                 data_param = data_param)

  return(result)

}

#### Function to create the data parameter names ####

create_lmean_data_param <- function(dataList, control) {

  means_vector <- paste("means", control$subfix, sep = "")

  #### Result ####

  result <- list(means_vector = means_vector)

  return(result)

}

#### Function to create the transformed-parameter model ####

model_lmean <- function(dataList, data_param, control) {

  means_vector <- data_param$means_vector

  list_struct <- list(
    list(name = means_vector,
         type = "matrix",
         dim = c(dataList$nitems, 1L),
         rownames = dataList$item_names,
         colnames = "intrcpt")
  )

  trans <- create_parameters(list_struct)

  #### Result ####

  return(trans)

}

#### Function to create the parameter constraints ####

constraints_lmean <- function(trans, dataList, data_param, control) {

  param <- trans

  #### Result ####

  return(param)

}

#### Function to create starting values ####

start_lmean <- function(param, trans, dataList, data_param, control) {

  means_vector <- data_param$means_vector

  if(control$std.ov) {
    means <- rep(0, dataList$nitems)
  } else {
    means <- colMeans(dataList$data, na.rm = TRUE)
  }

  init_param <- vector("list", length = control$rstarts)

  for(rs in seq_len(control$rstarts)) {

    init_param[[rs]] <- list()

    init_param[[rs]][[means_vector]] <-
      matrix(means,
             nrow = dataList$nitems,
             ncol = 1L,
             dimnames = dimnames(trans[[means_vector]]))

  }

  #### Result ####

  return(init_param)

}

#### Function to create the modelInfo ####

create_lmean_modelInfo <- function(dataList, full_model, control) {

  list2env(full_model, envir = environment())
  list2env(data_param, envir = environment())

  #### Manifolds ####

  manifolds <- manifolds_lmean(param = param,
                               dataList = dataList,
                               data_param = data_param,
                               control = control)

  control_manifold <- create_manifolds(manifolds = manifolds,
                                       structures = param)

  #### Transformations ####

  transforms <- transformations_lmean(trans = trans,
                                      dataList = dataList,
                                      data_param = data_param,
                                      control = control)

  control_transform <- create_transforms(transforms = transforms,
                                         structures = trans)

  #### Estimators ####

  estimators <- estimators_lmean(trans = trans,
                                 dataList = dataList,
                                 data_param = data_param,
                                 control = control)

  control_estimator <- create_estimators(estimators = estimators,
                                         structures = trans)

  #### Pass the initial values to vectors ####

  inits <- create_init(trans, param, init_param,
                       control_transform = control_transform, control)

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
                    dof = dataList$npatterns - nparam,
                    control_manifold = control_manifold,
                    control_transform = control_transform,
                    control_estimator = control_estimator,
                    control_optimizer = control_optimizer)

  return(modelInfo)

}

#### Function to create the manifolds ####

manifolds_lmean <- function(param, dataList, data_param, control) {

  manifolds <- list(
    list(manifold = "euclidean",
         parameters = data_param$means_vector)
  )

  #### Result ####

  return(manifolds)

}

#### Function to create the transformations ####

transformations_lmean <- function(trans, dataList, data_param, control) {

  transforms <- list()

  #### Result ####

  return(transforms)

}

#### Function to create the estimators ####

estimators_lmean <- function(trans, dataList, data_param, control) {

  estimators <- list()

  #### Result ####

  return(estimators)

}

#### Function to fit the model ####

fit_lmean <- function(dataList, modelInfo, control) {

  means_vector <- paste("means", control$subfix, sep = "")

  if(control$std.ov) {
    means <- rep(0, dataList$nitems)
  } else {
    means <- colMeans(dataList$data, na.rm = TRUE)
  }

  means <- matrix(means,
                  nrow = dataList$nitems,
                  ncol = 1L,
                  dimnames = dimnames(modelInfo$trans[[means_vector]]))

  values <- c(means)
  names(values) <- c(modelInfo$trans[[means_vector]])

  Optim <- list()
  Optim$f <- 0

  Optim$parameters <- values[modelInfo$parameters_labels]
  Optim$transparameters <- values[modelInfo$transparameters_labels]

  names(Optim$parameters) <- modelInfo$parameters_labels
  names(Optim$transparameters) <- modelInfo$transparameters_labels

  Optim$elapsed <- 0

  #### Result ####

  return(Optim)

}

#### Function to compute standard errors ####

compute_se_lmean <- function(dataList, modelInfo, Optim, control) {

  if(dataList$nobs < 2L || control$std.ov) {

    ACOV <- matrix(0,
                   nrow = dataList$nitems,
                   ncol = dataList$nitems)

  } else {

    variances <- apply(dataList$data, MARGIN = 2L,
                       FUN = var, na.rm = TRUE)

    ACOV <- diag(variances,
                 nrow = dataList$nitems,
                 ncol = dataList$nitems)

  }

  rownames(ACOV) <- colnames(ACOV) <- modelInfo$parameters_labels

  SE <- list(
    ACOV = ACOV,
    se = sqrt(diag(ACOV) / dataList$nobs)
  )

  names(SE$se) <- modelInfo$parameters_labels

  #### Result ####

  return(SE)

}

#### Function to create the control list ####

lmean_control <- function(control) {

  #### Optimizer defaults ####

  if(is.null(control$step_maxit)) {
    control$step_maxit <- 30L
  }

  if(is.null(control$c1)) {
    control$c1 <- 0.5
  }

  if(is.null(control$c2)) {
    control$c2 <- 0.5
  }

  if(is.null(control$step_eps)) {
    control$step_eps <- 1e-09
  }

  if(is.null(control$df_eps)) {
    control$df_eps <- 1e-09
  }

  if(is.null(control$M)) {
    control$M <- 100L
  }

  if(is.null(control$eps)) {
    control$eps <- 1e-04
  }

  if(is.null(control$ss_fac)) {
    control$ss_fac <- 2
  }

  if(is.null(control$maxit)) {
    control$maxit <- 1000L
  }

  if(is.null(control$rstarts)) {
    control$rstarts <- 1L
  }

  if(is.null(control$cores)) {
    control$cores <- 1L
  }

  if(is.null(control$tcg_maxit)) {
    control$tcg_maxit <- 10L
  }

  if(is.null(control$ss)) {
    control$ss <- 0.001
  }

  if(is.null(control$start)) {
    control$start <- NULL
  }

  if(is.null(control$subfix)) {
    control$subfix <- ""
  }

  #### Fixed rstarts and cores ####

  control$rstarts <- 1L
  control$cores <- 1L

  #### Result ####

  return(control)

}

#### Function to print progress messages ####

print_lmean_message <- function(msg) {

  w <- nchar(msg)+4L

  cat("\n", "+", strrep("-", w), "+\n",
      "|  ", msg, "  |\n",
      "+", strrep("-", w), "+\n\n", sep = "")

  #### Result ####

  return(invisible(NULL))

}
