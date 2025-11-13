poly_positive <- function(data = NULL, R = NULL, taus = NULL, n = NULL,
                          penalties = TRUE,
                          control = list(opt = "lbfgs", maxit = 100,
                                         eps = 1e-04)) {

  penalties = TRUE
  control = list(opt = "lbfgs", maxit = 100,
                 eps = 1e-04)
  # data: Full data
  # R: nonpositive-definite correlation matrix from polyfast
  # taus: thresholds from polyfast
  # n: contingency tables from polyfast

  control$rstarts <- 1L
  control$cores <- 1L
  control$penalties <- TRUE
  control$ss <- 0.001 # TAKE CARE OF TAUS WHEN THEY ARE ESTIMATED
  # Step sizes should be small so taus are not far from sensible bounds
  control <- lcfa_control(control)

  if(is.null(data)) {
    if(is.null(R) || is.null(taus) || is.null(n)) {
      stop("If data = NULL, then you must provide R, taus, and n instead.")
    }
  } else {
    polychorics <- polyfast(data)
    R <- polychorics$correlation
    taus <- polychorics$thresholds
    n <- polychorics$contingency_tables # This is the important thing for fitting
  }

  p <- nrow(R)

  #### Create the initial values for the parameters ####

  threslds <- lapply(taus, FUN = \(x) x[!is.infinite(x)])
  init_param <- list()
  init_param$taus <- threslds # Initial values for the thresholds
  init_param$X <- real_sqrtmat(R) # Initial values for the correlations
  init_trans <- init_param
  init_trans$R <- crossprod(init_trans$X)

  #### Parameters of the model ####

  poly_param <- list()
  poly_param$taus <- vector("list", length = p)
  for(i in 1:p) {
    K <- length(taus[[i]])-2L
    poly_param$taus[[i]] <- paste(".tau", i, ".", 1:K, sep = "")
  }
  poly_param$X <- matrix(paste(".X", 1:(p*p), sep = ""), nrow = p, ncol = p)
  # diag(poly_param$R) <- "1"
  # poly_param$R[upper.tri(poly_param$R)] <- t(poly_param$R)[upper.tri(poly_param$R)]

  # Create the model for the transformed parameters:

  poly_trans <- poly_param
  poly_trans$R <- matrix(paste("r", 1:(p*p), sep = ""), nrow = p, ncol = p)
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
  indices <- list()
  indices[[1]] <- match(labels, transparameters_labels)-1L
  control_estimator <- list()
  control_estimator[[1]] <- list(estimator = "polycor",
                                 labels = labels,
                                 indices = indices,
                                 indices_taus = indices_taus,
                                 indices_R = indices_R,
                                 n = n,
                                 p = p)

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

  fit <- optimizer(control_manifold, control_transform,
                   control_estimator, control)
  fit$f
  fit$iterations
  fit$ng
  # polys <- matrix(fit$outputs$estimators$matrices[[1]][[1]], 16, 16)
  # det(polys)

  result <- list(fit = fit,
                 control_manifold = control_manifold,
                 control_transform = control_transform,
                 control_estimator = control_estimator,
                 control = control)

  return(result)

}
