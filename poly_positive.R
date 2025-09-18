poly_positive <- function(data = NULL, R = NULL, taus = NULL, n = NULL,
                          control = list(opt = "lbfgs", maxit = 1000,
                                         eps = 1e-04)) {

  # data: Full data
  # R: non-positive correlation matrix from polyfast
  # taus: thresholds from polyfast
  # n: contingency tables from polyfast

  control$rstarts <- 1L
  control$cores <- 1L

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
  threslds <- unlist(taus)
  threslds <- threslds[!is.infinite(threslds)]
  # X <- roblq(p, p)
  X <- real_sqrtmat(R) # Initial values for the correlations
  init_param <- c(threslds, X)

  # Create the model:

  poly_param <- list()
  poly_param$taus <- vector("list", length = p)
  for(i in 1:p) {
    K <- length(taus[[i]])-2L
    poly_param$taus[[i]] <- paste(".tau", i, ".", 1:K, sep = "")
  }
  poly_param$R <- matrix(paste(".poly", 1:(p*p), sep = ""), nrow = p, ncol = p)
  # diag(poly_param$R) <- "1"
  # poly_param$R[upper.tri(poly_param$R)] <- t(poly_param$R)[upper.tri(poly_param$R)]
  vector_param <- unname(unlist(poly_param))

  indicator <- is.na(suppressWarnings(as.numeric(vector_param)))
  fixed_indices <- which(!indicator)
  fixed_values <- as.numeric(vector_param[fixed_indices])
  parameters_indices <- which(indicator)
  parameters_labels <- unique(vector_param[parameters_indices])
  nparam <- length(parameters_labels)

  # Created the model for the transformed parameters:

  poly_trans <- list()
  poly_trans$taus <- vector("list", length = p)
  for(i in 1:p) {
    K <- length(taus[[i]])-2L
    poly_trans$taus[[i]] <- paste("t.tau", i, ".", 1:K, sep = "")
  }
  poly_trans$R <- matrix(paste("t.poly", 1:(p*p), sep = ""), nrow = p, ncol = p)
  vector_trans <- unname(unlist(poly_trans))
  transparameters_labels <- unique(c(parameters_labels, vector_trans))
  parameters <- init_param[match(parameters_labels, vector_param)]
  transparameters <- init_param
  transparameters[fixed_indices] <- fixed_values

  # Relate the transformed parameters to the parameters:
  param2trans <- match(transparameters_labels, parameters_labels)
  param2trans <- param2trans[!is.na(param2trans)]
  # Relate the parameters to the transformed parameters:
  trans2param <- match(parameters_labels, transparameters_labels)
  transparameters <- c(parameters, transparameters)

  # Manifolds:
  control_manifold <- list()
  indices <- which(parameters_labels %in% unlist(poly_param$taus))
  labels <- parameters_labels[indices]
  control_manifold[[1]] <- list(manifold = "euclidean",
                                parameters = labels,
                                indices = list(indices-1L))

  indices <- which(parameters_labels %in% c(poly_param$R))
  labels <- parameters_labels[indices]
  control_manifold[[2]] <- list(manifold = "oblq",
                                parameters = labels,
                                indices = list(indices-1L),
                                q = p)

  # Transformations:
  control_transform <- list()

  all_param <- unlist(poly_param$taus)
  all_trans <- unlist(poly_trans$taus)
  positions <- which(all_param %in% parameters_labels)
  labels_in <- all_param[positions]
  labels_out <- all_trans[positions]
  indices_in <- match(labels_in, transparameters_labels)
  indices_out <- match(labels_out, transparameters_labels)
  control_transform[[1]] <- list(transform = "identity",
                                 labels_in = labels_in,
                                 indices_in = list(indices_in-1L),
                                 labels_out = labels_out,
                                 indices_out = list(indices_out-1L))

  # all_param <- c(poly_param$R)
  # all_trans <- c(poly_trans$R)
  # positions <- which(all_param %in% parameters_labels)
  # labels_in <- all_param[positions]
  # labels_out <- all_trans[positions]
  labels_in <- c(poly_param$R)
  labels_out <- c(poly_trans$R)
  indices_in <- match(labels_in, transparameters_labels)
  indices_out <- match(labels_out, transparameters_labels)
  control_transform[[2]] <- list(transform = "crossprod",
                                 labels_in = labels_in,
                                 indices_in = list(indices_in-1L),
                                 labels_out = labels_out,
                                 indices_out = list(indices_out-1L),
                                 p = nrow(poly_param$R),
                                 q = ncol(poly_param$R))

  # Estimators:
  indices <- list()
  indices[[1]] <- match(vector_trans, transparameters_labels)-1L
  k <- 2
  for(i in 1:p) {
    indices[[k]] <- match(poly_trans$taus[[i]], vector_trans)-1L
    k <- k+1
  }
  indices[[k]] <- match(c(poly_trans$R), vector_trans)-1L
  labels <- transparameters_labels[indices[[1]]+1L]
  control_estimator <- list()
  control_estimator[[1]] <- list(estimator = "polycor",
                                 labels = labels,
                                 indices = indices,
                                 n = n,
                                 p = p)

  # Control:
  control <- cfast_control(control)
  control$parameters <- list(parameters)
  control$transparameters <- list(transparameters)
  control$param2transparam <- param2trans-1L
  control$transparam2param <- trans2param-1L
  control$ss <- 0.001 # TAKE CARE OF TAUS WHEN THEY ARE ESTIMATED

  fit <- optimizer(control_manifold, control_transform,
                   control_estimator, control)

  result <- list(fit = fit,
                 control_manifold = control_manifold,
                 control_transform = control_transform,
                 control_estimator = control_estimator,
                 control = control)

  return(result)

}
