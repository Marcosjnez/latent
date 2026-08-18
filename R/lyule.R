# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 17/08/2026
#'
#' Yule Correlation Matrix
#'
#' Estimate a matrix of pairwise Yule association coefficients and their
#' asymptotic standard errors from categorical data.
#'
#' @usage
#' lyule(data, model = NULL, do.fit = TRUE, control = NULL, ...)
#'
#' @param data A data frame or matrix containing categorical observed variables.
#'   Factor, integer, logical, and numeric variables are supported by the
#'   underlying Yule association routine.
#' @param model Optional model specification reserved for internal model setup.
#' @param do.fit Logical. If \code{FALSE}, return the prepared but unfitted
#'   \code{"latent"} object.
#' @param control Optional list of internal controls. Custom starting values can
#'   be supplied through \code{control$start}.
#' @param ... Additional arguments reserved for future extensions.
#'
#' @details
#' \code{lyule()} estimates the pairwise association matrix directly using the
#' Yule copula-based association routine implemented in \pkg{latent}. The
#' diagonal of the resulting matrix is fixed to one and only the off-diagonal
#' associations are treated as free parameters.
#'
#' Standard errors are obtained from the pairwise asymptotic standard errors
#' returned by the Yule association routine. The corresponding sampling
#' variance-covariance matrix is stored in \code{Optim$SE$VCOV}.
#'
#' @return An S4 object of class \code{"latent"}. The object contains the
#'   processed data in \code{dataList}, the parameter and optimization structures
#'   in \code{modelInfo}, the direct estimation output and asymptotic covariance
#'   information in \code{Optim}, and the estimated parameter structures in
#'   \code{parameters} and \code{transformed_pars}.
#'
#' @examples
#' \dontrun{
#' fit <- lyule(data = data.frame(
#'   x1 = factor(sample(1:3, 200, replace = TRUE)),
#'   x2 = factor(sample(1:4, 200, replace = TRUE)),
#'   x3 = factor(sample(1:2, 200, replace = TRUE))
#' ))
#' }
#'
#' @export
lyule <- function(data,
                  model = NULL,
                  do.fit = TRUE,
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

  supported <- vapply(data, FUN = function(x) {
    is.factor(x) || is.numeric(x) || is.integer(x) || is.logical(x)
  }, FUN.VALUE = logical(1L))

  if(!all(supported)) {
    stop("All variables in data must be factor, numeric, integer, or logical")
  }

  if(!is.logical(do.fit) || length(do.fit) != 1L || is.na(do.fit)) {
    stop("do.fit must be TRUE or FALSE")
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

  control <- lyule_control(control)

  #### Create the dataList ####

  dataList <- create_lyule_dataList(data = data,
                                    control = control)
  dataList$args <- args

  #### Create the model ####

  full_model <- create_lyule_model(dataList = dataList,
                                   model = model,
                                   control = control)

  #### Create the manifold, transformation, and estimator structures ####

  modelInfo <- create_lyule_modelInfo(dataList = dataList,
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

  fit_output <- fit_lyule(dataList = dataList,
                          modelInfo = modelInfo,
                          control = control)

  Optim <- fit_output$Optim

  #### Process the outputs ####

  transformed_pars <- fill_in(modelInfo$trans,
                              Optim$transparameters)

  parameters <- transformed_pars[names(modelInfo$param)]

  #### Standard errors ####

  Optim$SE <- compute_se_lyule(dataList = dataList,
                               modelInfo = modelInfo,
                               yule_output = fit_output$yule_output)

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

create_lyule_dataList <- function(data, control) {

  #### Result ####

  result <- list(data = data,
                 nobs = nrow(data),
                 nitems = ncol(data),
                 npatterns = nrow(data),
                 item_names = colnames(data))

  return(result)

}

#### Function to create the model ####

create_lyule_model <- function(dataList, model = NULL, control) {

  # Generate the model syntax and initial parameter values.

  #### Parameter-block names ####

  data_param <- create_lyule_data_param(dataList = dataList,
                                        control = control)

  #### Model for the transformed parameters ####

  trans <- model_lyule(dataList = dataList,
                       data_param = data_param,
                       control = control)

  #### Model for the parameters ####

  param <- constraints_lyule(dataList = dataList,
                             data_param = data_param,
                             trans = trans,
                             control = control)

  #### Create the initial values for the parameters ####

  init_param <- start_lyule(dataList = dataList,
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

create_lyule_modelInfo <- function(dataList, full_model, control) {

  # Generate the manifold, transformation, estimator, and optimizer structures.

  list2env(full_model, envir = environment())

  #### Manifolds ####

  manifolds <- manifolds_lyule(dataList = dataList,
                               data_param = data_param,
                               param = param,
                               control = control)

  control_manifold <- create_manifolds(manifolds = manifolds,
                                       structures = param)

  #### Transformations ####

  transforms <- transformations_lyule(dataList = dataList,
                                       data_param = data_param,
                                       trans = trans,
                                       control = control)

  control_transform <- create_transforms(transforms = transforms,
                                         structures = trans)

  #### Estimators ####

  estimators <- estimators_lyule(dataList = dataList,
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

lyule_control <- function(control) {

  # Auxiliary function for lyule.R.

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

  if(is.null(control$ipfp_maxit)) {
    control$ipfp_maxit <- 2000L
  } else if(control$ipfp_maxit < 1L) {
    stop("ipfp_maxit must be a positive integer")
  }

  if(is.null(control$ipfp_tol)) {
    control$ipfp_tol <- 1e-12
  } else if(control$ipfp_tol <= 0) {
    stop("ipfp_tol must be a positive number")
  }

  if(is.null(control$pinv_tol)) {
    control$pinv_tol <- 1e-12
  } else if(control$pinv_tol <= 0) {
    stop("pinv_tol must be a positive number")
  }

  # lyule uses direct estimation and therefore a single start and core.
  control$rstarts <- 1L
  control$cores <- 1L

  #### Result ####

  return(control)

}

#### Auxiliary functions for create_lyule_model ####

create_lyule_data_param <- function(dataList, control) {

  S_matrix <- paste0("S", control$subfix)

  #### Result ####

  result <- list(S_matrix = S_matrix)

  return(result)

}

model_lyule <- function(dataList, data_param, control) {

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

constraints_lyule <- function(dataList, data_param, trans, control) {

  list2env(data_param, envir = environment())

  param <- trans

  # Yule association matrices have unit diagonal.
  diag(param[[S_matrix]]) <- "1"

  #### Result ####

  return(param)

}

start_lyule <- function(dataList, data_param, param, trans, control) {

  list2env(data_param, envir = environment())

  init_param <- vector("list", length = control$rstarts)

  for(rs in seq_len(control$rstarts)) {

    init_param[[rs]] <- list()
    init_param[[rs]][[S_matrix]] <- diag(1, nrow = dataList$nitems,
                                         ncol = dataList$nitems)
    dimnames(init_param[[rs]][[S_matrix]]) <- dimnames(trans[[S_matrix]])

  }

  #### Result ####

  return(init_param)

}

#### Auxiliary functions for create_lyule_modelInfo ####

manifolds_lyule <- function(dataList, data_param, param, control) {

  list2env(data_param, envir = environment())

  manifolds <- list(
    list(manifold = "euclidean",
         parameters = S_matrix)
  )

  #### Result ####

  return(manifolds)

}

transformations_lyule <- function(dataList, data_param, trans, control) {

  #### Result ####

  return(list())

}

estimators_lyule <- function(dataList, data_param, trans, control) {

  #### Result ####

  return(list())

}

#### Auxiliary functions for fitting and post-estimation ####

fit_lyule <- function(dataList, modelInfo, control) {

  #### Estimate the Yule association matrix ####

  yule_output <- yule_cor_full_rcpp(dataList$data,
                                    ipfp_maxit = control$ipfp_maxit,
                                    ipfp_tol = control$ipfp_tol,
                                    pinv_tol = control$pinv_tol)

  S <- yule_output$cor
  rownames(S) <- colnames(S) <- dataList$item_names

  #### Identify the association-matrix parameter block ####

  S_name <- names(modelInfo$trans)[startsWith(names(modelInfo$trans), "S")]

  if(length(S_name) != 1L) {
    stop("Unable to identify the Yule association matrix")
  }

  S_labels <- c(modelInfo$trans[[S_name]])
  S_values <- c(S)

  #### Transformed parameters ####

  idx <- match(modelInfo$transparameters_labels, S_labels)

  if(anyNA(idx)) {
    stop("Unable to match the transformed Yule parameters")
  }

  transparameters <- S_values[idx]
  names(transparameters) <- modelInfo$transparameters_labels

  #### Free parameters ####

  idx <- match(modelInfo$parameters_labels,
               modelInfo$transparameters_labels)

  if(anyNA(idx)) {
    stop("Unable to match the free Yule parameters")
  }

  parameters <- transparameters[idx]
  names(parameters) <- modelInfo$parameters_labels

  #### Optimizer-style output ####

  Optim <- list(f = 0,
                elapsed = 0,
                parameters = parameters,
                transparameters = transparameters)

  #### Result ####

  result <- list(Optim = Optim,
                 yule_output = yule_output)

  return(result)

}

compute_se_lyule <- function(dataList, modelInfo, yule_output) {

  #### Off-diagonal standard errors ####

  lower_idx <- lower.tri(yule_output$cor, diag = FALSE)
  pair_se <- yule_output$se[lower_idx]

  if(length(pair_se) != length(modelInfo$parameters_labels)) {
    stop("The number of Yule standard errors does not match the number of free parameters")
  }

  #### Variance-covariance matrix ####

  # yule_cor_full_rcpp() already returns sampling standard errors, so no
  # additional multiplication or division by the sample size is required.
  VCOV <- diag(pair_se^2,
               nrow = length(pair_se),
               ncol = length(pair_se))

  rownames(VCOV) <- colnames(VCOV) <- modelInfo$parameters_labels

  se <- sqrt(diag(VCOV))
  names(se) <- modelInfo$parameters_labels

  #### Result ####

  result <- list(VCOV = VCOV,
                 se = se)

  return(result)

}
