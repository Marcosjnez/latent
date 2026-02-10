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
                  penalties = TRUE,
                  do.fit = TRUE,
                  control = NULL) {

  ## store original call
  mc  <- match.call()

  # Check the arguments to control_optimizer and create defaults:
  control$penalties <- penalties
  control <- lpoly_control(control)

  #### Create the data_list ####

  # Estimate the polychoric correlations without positive-definite constraints:
  polychorics <- polyfast(as.matrix(data))
  R <- polychorics$correlation # Polychoric correlation matrix
  taus <- polychorics$thresholds # Thresholds
  n <- polychorics$contingency_tables # Contingency tables
  p <- ncol(data) # Number of items
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

  #### Parameters of the model ####

  poly_param <- list()
  poly_param$taus <- vector("list", length = p)
  K <- vector(length = p)
  for(i in 1:p) {
    K[i] <- length(taus[[i]])-2L
    poly_param$taus[[i]] <- paste(".tau", i, ".", 1:K[i], sep = "")
  }
  poly_param$X <- matrix(paste(".X", 1:(p*p), sep = ""), nrow = p, ncol = p)

  # Create the model for the transformed parameters:

  poly_trans <- poly_param
  poly_trans$R <- matrix(paste("r", 1:(p*p), sep = ""), nrow = p, ncol = p)
  # diag(poly_trans$R) <- "1"
  poly_trans$R[upper.tri(poly_trans$R)] <- t(poly_trans$R)[upper.tri(poly_trans$R)]

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

  #### Create the vectors of parameters and transformed parameters ####

  parameters <- transparameters <- vector("list", length = control$rstarts)
  # Indices of the unique transparameters in init_trans:
  trans_inds <- match(transparameters_labels, vector_trans)
  init_inds <- match(parameters_labels, vector_trans)

  #### Create the initial values for the parameters ####

  # Thresholds without the infite values:
  threslds <- lapply(taus, FUN = \(x) x[!is.infinite(x)])
  init_param <- list()
  init_param$taus <- threslds # Initial values for the thresholds
  # Initial values for the square root of correlations:
  init_param$X <- real_sqrtmat(R)
  init_trans <- init_param
  # Initial values for the polychoric correlations:
  init_trans$R <- crossprod(init_trans$X)

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

  #### Structures ####

  # Manifolds:
  control_manifold <- list()
  indices <- match(unlist(poly_param$taus), parameters_labels)
  labels <- parameters_labels[indices]
  control_manifold[[1]] <- list(manifold = "euclidean",
                                parameters = labels,
                                indices = list(indices-1L))

  indices <- match(c(poly_param$X), parameters_labels)
  labels <- parameters_labels[indices]
  control_manifold[[2]] <- list(manifold = "oblq",
                                parameters = labels,
                                indices = list(indices-1L),
                                q = p)

  # Transformations:
  control_transform <- list()
  lower_indices <- which(lower.tri(poly_trans$R, diag = TRUE))

  labels_in <- c(poly_trans$X)
  labels_out <- poly_trans$R[lower_indices]
  indices_in <- match(labels_in, transparameters_labels)
  indices_out <- match(labels_out, transparameters_labels)
  control_transform[[1]] <- list(transform = "crossprod",
                                 labels_in = labels_in,
                                 indices_in = list(indices_in-1L),
                                 labels_out = labels_out,
                                 indices_out = list(indices_out-1L),
                                 p = p)

  # Estimators:
  indices_taus <- vector("list", length = p)
  for(i in 1:p) {
    indices_taus[[i]] <- match(poly_trans$taus[[i]], transparameters_labels)-1L
  }
  indices_R <- match(poly_trans$R[lower_indices], transparameters_labels)-1L

  labels <- c(unlist(poly_trans$taus), poly_trans$R[lower_indices])
  indices <- list(match(labels, transparameters_labels)-1L)
  control_estimator <- list()
  control_estimator[[1]] <- list(estimator = "polycor",
                                 labels = labels,
                                 indices = indices,
                                 indices_taus = indices_taus,
                                 indices_R = indices_R,
                                 n = n,
                                 p = p,
                                 N = data_list$nobs)

  if(control$reg) {

    labels <- poly_trans$R[lower_indices]
    indices <- match(labels, transparameters_labels)
    control_estimator[[2]] <- list(estimator = "logdetmat",
                                   labels = labels,
                                   indices = list(indices-1L),
                                   lower_indices = lower_indices-1L,
                                   p = p,
                                   logdetw = control$penalties$logdet$w)

  }

  #### Collect all the model information ####

  # Model information:
  nparam <- 0.5*p*(p-1) + sum(K)
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
                  version            = as.character( packageVersion('latent') ),
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

  # Collect all the information about the optimization:
  elapsed <- Optim$elapsed

  # x$f
  # x$iterations
  # x$ng
  # polys <- matrix(fit$outputs$estimators$matrices[[1]][[1]], 16, 16)
  # det(polys)

  #### Estimated model structures ####

  # Create the structures of untransformed parameters:
  indices_pars <- match(modelInfo$parameters_labels,
                        unlist(modelInfo$poly_param))
  vv <- rep(0, times = length(unlist(modelInfo$poly_param)))
  vv[indices_pars] <- Optim$parameters
  parameters <- fill_list_with_vector(modelInfo$poly_param, vv)
  parameters <- allnumeric(parameters)
  # FIXED PARAMETERS?

  # Create the structures of transformed parameters:
  indices_trans <- match(modelInfo$transparameters_labels,
                         unlist(modelInfo$poly_trans))
  vv <- rep(0, times = length(unlist(modelInfo$poly_trans)))
  vv[indices_trans] <- Optim$transparameters
  transformed_pars <- fill_list_with_vector(modelInfo$poly_trans, vv)
  transformed_pars <- allnumeric(transformed_pars)

  #### Process the fit information ####

  all_estimators <- unlist(lapply(modelInfo$control_estimator, FUN = \(x) x$estimator))
  index_estimator <- which(all_estimators == "polycor")
  loss <- Optim$outputs$estimators$doubles[[index_estimator]][1]
  penalized_loss <- Optim$f
  loglik <- Optim$outputs$estimators$doubles[[index_estimator]][2]
  penalized_loglik <- sum(unlist(lapply(Optim$outputs$estimators$doubles, FUN = \(x) x[[2]])))

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
