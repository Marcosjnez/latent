#' @title
#' Maximum likelihood estimation of positive definite polychoric correlation matrices.
#'
#' @usage
#'
#' lpoly(data = NULL,
#' penalties = TRUE,
#' do.fit = TRUE,
#' control = NULL)
#'
#' @param data data frame or matrix with the raw data.
#' @param penalties Force a positive-definite solution. Defaults to TRUE.
#' @param do.fit TRUE to fit the model and FALSE to return only the model setup. Defaults to TRUE.
#' @param control List of control parameters for the optimization algorithm. See 'details' for more information.
#'
#' @details \code{lpoly} estimates positive-definite polychoric correlation matrices.
#'
#' @return List with the following objects:
#' \item{version}{Version number of 'latent' when the model was estimated.}
#' \item{call}{Code used to estimate the model.}
#' \item{ModelInfo}{Model information.}
#' \item{Optim}{Output of the optimizer.}
#' \item{parameters}{Structure with all model parameters.}
#' \item{transparameters}{Structure with all transformed model parameters.}
#' \item{loglik}{Logarithm likelihood of the model.}
#' \item{penalized_loglik}{Logarithm likelihood + logarithm priors of the model.}
#'
#' @examples
#'
#' \dontrun{
#' fit <- lpoly(data = values)
#'}
#'
#' @export
lpoly <- function(data,
                  model = NULL,
                  method = "crossprod",
                  penalties = TRUE,
                  do.fit = TRUE,
                  control = NULL,
                  ...) {

  dots <- list(...)

  ## store original call
  mc  <- match.call()

  # Check the arguments to control_optimizer and create defaults:
  control$penalties <- penalties
  control <- lpoly_control(control)

  #### Create the data_list ####

  data_list <- create_lpoly_datalist(data, control)
  list2env(data_list, envir = environment())

  #### Parameters of the model ####

  poly_param <- list()
  poly_param$taus <- vector("list", length = nitems)
  K <- vector(length = nitems)
  for(i in 1:nitems) {
    K[i] <- length(taus[[i]])-2L
    poly_param$taus[[i]] <- paste(".tau", i, ".", 1:K[i], sep = "")
  }

  if(method == "crossprod") {
    poly_param$X <- matrix(paste(".X", 1:(nitems*nitems), sep = ""),
                           nrow = nitems, ncol = nitems)
  } else {
    poly_param$R <- matrix(paste("r", 1:(nitems*nitems), sep = ""),
                           nrow = nitems, ncol = nitems)
    diag(poly_param$R) <- "1"
    poly_param$R[upper.tri(poly_param$R)] <- t(poly_param$R)[upper.tri(poly_param$R)]
  }

  # Create the model for the transformed parameters:

  poly_trans <- poly_param
  poly_trans$R <- matrix(paste("r", 1:(nitems*nitems), sep = ""),
                         nrow = nitems, ncol = nitems)
  poly_trans$R[upper.tri(poly_trans$R)] <- t(poly_trans$R)[upper.tri(poly_trans$R)]

  #### Fix parameters ####

  if(!is.null(model$taus)) {
    poly_param$taus <- model$taus
  }

  #### Arrange labels ####

  # Arrange parameter labels:
  vector_param <- unname(unique(unlist(poly_param)))

  # Select the unique, nonnumeric labels:
  nonfixed_pars <- which(is.na(suppressWarnings(as.numeric(vector_param))))
  parameters_labels <- vector_param[nonfixed_pars]
  nparam <- length(parameters_labels)

  # Arrange transparameter labels:
  vector_trans <- unname(unlist(poly_trans))
  transparameters_labels <- unique(vector_trans)
  ntrans <- length(transparameters_labels)

  #### Relate the transformed parameters to the parameters ####

  param2trans <- match(transparameters_labels, parameters_labels)
  param2trans <- param2trans[!is.na(param2trans)]
  # Relate the parameters to the transformed parameters:
  trans2param <- match(parameters_labels, transparameters_labels)

  #### Create the initial values for the parameters ####

  # Thresholds without the infite values:
  threslds <- lapply(taus, FUN = \(x) x[!is.infinite(x)])
  init_param <- list()
  init_param$taus <- threslds # Initial values for the thresholds
  if(method == "crossprod") {
    # Initial values for the square root of correlations:
    init_param$X <- real_sqrtmat(data_list$R)
    init_trans <- init_param
    # Initial values for the polychoric correlations:
    init_trans$R <- crossprod(init_trans$X)
  } else {
    init_param$R <- data_list$R
    init_trans <- init_param
  }

  #### Create the vectors of parameters and transformed parameters ####

  parameters <- transparameters <- vector("list", length = control$rstarts)
  # Indices of the unique transparameters in init_trans:
  trans_inds <- match(transparameters_labels, vector_trans)
  init_inds <- match(parameters_labels, vector_trans)

  transparameters <- unlist(init_trans)[trans_inds]
  names(transparameters) <- transparameters_labels
  parameters <- unlist(init_trans)[init_inds]
  names(parameters) <- parameters_labels

  #### Set up the optimizer ####

  # Create defaults for the control of the optimizer:
  control$parameters <- list(parameters)
  control$transparameters <- list(transparameters)
  control$param2transparam <- param2trans-1L
  control$transparam2param <- trans2param-1L

  #### Manifolds ####

  if(method == "crossprod") {

    manifolds <- list(
      list(manifold = "euclidean", parameters = "taus"),
      list(manifold = "oblq", parameters = "X",
           extra = list(p = nitems, q = nitems))
    )

  } else {

    manifolds <- list(
      list(manifold = "euclidean", parameters = "taus"),
      list(manifold = "euclidean", parameters = "R",
           extra = list(p = nitems, q = nitems))
    )

  }

  control_manifold <- create_manifolds(manifolds = manifolds,
                                       structures = poly_param)

  #### Transformations ####

  transform <- list()

  lower_indices <- which(lower.tri(poly_trans$R, diag = TRUE))
  dots$p <- nitems

  if(method == "crossprod") {

    transforms <- list(
      list(transform = "crossprod",
           parameters_in = "X",
           parameters_out = "R",
           extra = dots)
    )

  } else {

    transforms <- list()

  }

  control_transform <- create_transforms(transforms = transforms,
                                         structures = poly_trans)

  #### Estimators ####

  estimators <- list()

  estimators[[1]] <- list(estimator = "polycor",
                          parameters = c(list(poly_trans$R),
                                         poly_trans$taus),
                          extra = list(n = data_list$n,
                                       p = nitems,
                                       N = data_list$nobs))

  if(control$reg) {

    estimators[[2]] <- list(estimator = "logdetmat",
                            parameters = "R",
                            extra = list(lower_indices = lower_indices-1L,
                                         p = nitems,
                                         logdetw = control$penalties$logdet$w))

  }

  control_estimator <- create_estimators(estimators = estimators,
                                         structures = poly_trans)

  #### Collect all the model information ####

  # Model information:
  nparam <- 0.5*nitems*(nitems-1) + sum(K)
  modelInfo <- list(nobs = nobs,
                    nparam = nparam,
                    npatterns = npatterns,
                    dof = npatterns - nparam,
                    ntrans = ntrans,
                    parameters_labels = parameters_labels,
                    transparameters_labels = transparameters_labels,
                    poly_param = poly_param,
                    poly_trans = poly_trans,
                    init_param = init_param,
                    init_trans = init_trans,
                    control_manifold = control_manifold,
                    control_transform = control_transform,
                    control_estimator = control_estimator,
                    control = control)

  if(!do.fit) {

    result <- new("lpoly",
                  version            = as.character(packageVersion('latent')),
                  call               = mc,
                  timing             = numeric(),
                  modelInfo          = modelInfo,
                  Optim              = list(),
                  dataList           = data_list,
                  parameters         = list(),
                  transformed_pars   = list(),
                  loglik             = numeric(),
                  penalized_loglik   = numeric(),
                  loss               = numeric(),
                  penalized_loss     = numeric()
    )

    return(result)

  }

  #### Fit the model ####

  Optim <- optimizer(control_manifold, control_transform,
                     control_estimator, control)
  names(Optim$parameters) <- modelInfo$parameters_labels
  names(Optim$transparameters) <- modelInfo$transparameters_labels

  # Collect all the information about the optimization:
  elapsed <- Optim$elapsed

  #### Estimated model structures ####

  # Create the structures of untransformed parameters:
  parameters <- fill_in(modelInfo$poly_param, Optim$parameters)
  # FIXED PARAMETERS?

  # Create the structures of transformed parameters:
  transformed_pars <- fill_in(modelInfo$poly_trans, Optim$transparameters)

  #### Process the fit information ####

  all_estimators <- unlist(lapply(modelInfo$control_estimator, FUN = \(x) x$estimator))
  index_estimator <- which(all_estimators == "polycor")
  loss <- Optim$outputs$estimators$doubles[[index_estimator]][1]
  penalized_loss <- Optim$f
  loglik <- Optim$outputs$estimators$doubles[[index_estimator]][2]
  penalized_loglik <- sum(unlist(lapply(Optim$outputs$estimators$doubles,
                                        FUN = \(x) x[[2]])))

  #### Return ####

  result <- new("lpoly",
                version            = as.character( packageVersion('latent') ),
                call               = mc, # matched call
                timing             = elapsed, # timing information
                modelInfo          = modelInfo, # modelInfo
                Optim              = Optim, # Optim
                dataList           = data_list, # All data information
                parameters         = parameters,
                transformed_pars   = transformed_pars,
                loglik             = loglik, # loglik values
                penalized_loglik   = penalized_loglik,
                loss               = loss,
                penalized_loss     = penalized_loss
  )

  return(result)

}

asymptotic_poly <- function(fit, model = NULL) {

  control_manifold <- fit@modelInfo$control_manifold
  control_transform <- fit@modelInfo$control_transform
  control_estimator <- fit@modelInfo$control_estimator
  control_optimizer <- fit@modelInfo$control

  control_optimizer$parameters[[1]] <- fit@Optim$parameters
  control_optimizer$transparameters[[1]] <- fit@Optim$transparameters

  x <- get_hess(control_manifold, control_transform,
                control_estimator, control_optimizer,
                cores = parallel::detectCores())
  ACOV <- solve(x$h)

  return(ACOV)

}

create_lpoly_datalist <- function(data, control) {

  # Estimate the polychoric correlations without positive-definite constraints:
  polychorics <- polyfast(as.matrix(data))
  R <- polychorics$correlation # Polychoric correlation matrix
  taus <- polychorics$thresholds # Thresholds
  n <- polychorics$contingency_tables # Contingency tables
  nobs <- nrow(data)
  nitems <- ncol(data)
  item_label <- colnames(data)
  npatterns <- length(unlist(n))

  data_list <- vector("list")
  data_list$data <- data
  data_list$nobs <- nobs
  data_list$nitems <- nitems
  data_list$npatterns <- npatterns
  data_list$positive <- control$reg
  data_list$item_label <- item_label
  data_list$n <- n
  data_list$R <- R
  data_list$taus <- taus

  return(data_list)

}
