# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 01/11/2025

get_full_efa_model <- function(data_list, model, control = NULL) {

  # Generate the model syntax and initial parameter values

  list2env(data_list, envir = environment())

  # Initialize the objects to store the initial parameters:
  cfa_param <- cfa_trans <- vector("list", length = ngroups)
  fixed <- fixed_values <- nonfixed <- vector("list", length = ngroups)

  #### Parameters of the model for each group ####

  # Initialize the target matrices for positive-definite constraints:
  target_psi <- target_theta <- targets <- vector("list", length = ngroups)
  rest <- 0L

  for(i in 1:ngroups) {

    # Get the positions of parameters and fixed values:

    nonfixed[[i]] <- lapply(model[[i]], FUN = \(x) {
      which(is.na(suppressWarnings(as.numeric(x))))
    })

    fixed[[i]] <- lapply(model[[i]], FUN = \(x) {
      which(!is.na(suppressWarnings(as.numeric(x))))
    })

    fixed_values[[i]] <- lapply(model[[i]], FUN = \(x) {
      numerals <- suppressWarnings(as.numeric(x))
      inds <- which(!is.na(numerals))
      return(numerals[inds])
    })

    # Transformed parameters:

    # Unrotated lambda:
    ulambda_labels <- paste("g", i, ".ulambda[", rep(1:nitems[[i]], times = nfactors[[i]]),
                           ",", rep(1:nfactors[[i]], each = nitems[[i]]), "]", sep = "")
    ulambda_labels[nonfixed[[i]]$lambda] <- model[[i]]$lambda[nonfixed[[i]]$lambda]
    cfa_trans[[i]]$ulambda <- matrix(ulambda_labels, nrow = nitems[[i]], ncol = nfactors[[i]])

    # X:
    X_labels <- paste("g", i, ".X[", rep(1:nfactors[[i]], times = nfactors[[i]]),
                      ",", rep(1:nfactors[[i]], each = nfactors[[i]]), "]", sep = "")
    cfa_trans[[i]]$X <- matrix(X_labels, nrow = nfactors[[i]], ncol = nfactors[[i]])

    # Create additional parameters if there are positive-definite constraints:
    if(positive) {

      # Theta:
      pj_theta_labels <- paste("g", i, ".pj_theta[", rep(1:nitems[[i]], times = nitems[[i]]),
                               ",", rep(1:nitems[[i]], each = nitems[[i]]), "]", sep = "")
      cfa_trans[[i]]$pj_theta <- matrix(pj_theta_labels, nrow = nitems[[i]], ncol = nitems[[i]])

    }

    # Theta:
    theta_labels <- paste("g", i, ".theta[", rep(1:nitems[[i]], times = nitems[[i]]),
                          ",", rep(1:nitems[[i]], each = nitems[[i]]), "]", sep = "")
    theta_labels[nonfixed[[i]]$theta] <- model[[i]]$theta[nonfixed[[i]]$theta]
    cfa_trans[[i]]$theta <- matrix(theta_labels, nrow = nitems[[i]], ncol = nitems[[i]])
    # Force symmetry:
    cfa_trans[[i]]$theta[upper.tri(cfa_trans[[i]]$theta)] <- t(cfa_trans[[i]]$theta)[upper.tri(cfa_trans[[i]]$theta)]

    if(!orthogonal) {
      # Xinv:
      Xinv_labels <- paste("g", i, ".Xinv[", rep(1:nfactors[[i]], times = nfactors[[i]]),
                           ",", rep(1:nfactors[[i]], each = nfactors[[i]]), "]", sep = "")
      cfa_trans[[i]]$Xinv <- matrix(Xinv_labels, nrow = nfactors[[i]], ncol = nfactors[[i]])
    }

    # Rotated lambda:
    lambda_labels <- paste("g", i, ".lambda[", rep(1:nitems[[i]], times = nfactors[[i]]),
                           ",", rep(1:nfactors[[i]], each = nitems[[i]]), "]", sep = "")
    cfa_trans[[i]]$lambda <- matrix(lambda_labels, nrow = nitems[[i]], ncol = nfactors[[i]])

    # Psi:
    psi_labels <- paste("g", i, ".psi[", rep(1:nfactors[[i]], times = nfactors[[i]]),
                        ",", rep(1:nfactors[[i]], each = nfactors[[i]]), "]", sep = "")
    psi_labels[nonfixed[[i]]$psi] <- model[[i]]$psi[nonfixed[[i]]$psi]
    cfa_trans[[i]]$psi <- matrix(psi_labels, nrow = nfactors[[i]], ncol = nfactors[[i]])
    # Force symmetry:
    cfa_trans[[i]]$psi[upper.tri(cfa_trans[[i]]$psi)] <- t(cfa_trans[[i]]$psi)[upper.tri(cfa_trans[[i]]$psi)]

    # Model matrix:
    model_labels <- paste("g", i, ".model[", rep(1:nitems[[i]], times = nitems[[i]]),
                          ",", rep(1:nitems[[i]], each = nitems[[i]]), "]", sep = "")
    cfa_trans[[i]]$model <- matrix(model_labels, nrow = nitems[[i]], ncol = nitems[[i]])
    # Force symmetry:
    cfa_trans[[i]]$model[upper.tri(cfa_trans[[i]]$model)] <- t(cfa_trans[[i]]$model)[upper.tri(cfa_trans[[i]]$model)]

    # Untransformed parameters:

    cfa_param[[i]]$ulambda <- cfa_trans[[i]]$ulambda
    # Insert fixed values in the model:
    cfa_param[[i]]$ulambda[fixed[[i]]$lambda] <- model[[i]]$lambda[fixed[[i]]$lambda]

    cfa_param[[i]]$X <- cfa_trans[[i]]$X

    if(positive) {

      # Theta:
      cfa_param[[i]]$pj_theta <- cfa_trans[[i]]$pj_theta

    } else {

      # Theta:
      cfa_param[[i]]$theta <- cfa_trans[[i]]$theta
      # Insert fixed values in the model:
      cfa_param[[i]]$theta[fixed[[i]]$theta] <- model[[i]]$theta[fixed[[i]]$theta]

    }

    # Create the target matrices for positive-definite constraints:
    if(positive) {

      target_theta[[i]] <- matrix(0, nrow = nitems[[i]], ncol = nitems[[i]])
      target_theta[[i]][nonfixed[[i]]$theta] <- 1

      p <- nitems[[i]]
      lower_theta <- lower.tri(diag(p), diag = TRUE)
      targets[[i]] <- unlist(c(target_theta[[i]][lower_theta]))
      rest <- rest + 0.5*p*(p-1) + sum(targets[[i]] == 0)

    }

  }

  #### Arrange labels ####

  # Arrange parameter labels:
  vector_param <- unname(unique(unlist(cfa_param)))

  # Select the unique, nonnumeric labels:
  nonfixed_pars <- which(is.na(suppressWarnings(as.numeric(vector_param))))
  parameters_labels <- vector_param[nonfixed_pars]
  nparam <- length(parameters_labels)

  # Arrange transparameter labels:
  vector_trans <- unname(unlist(cfa_trans))
  transparameters_labels <- unique(vector_trans)
  ntrans <- length(transparameters_labels)

  #### Relate the transformed parameters to the parameters ####

  param2trans <- match(transparameters_labels, parameters_labels)
  param2trans <- param2trans[!is.na(param2trans)]
  # Relate the parameters to the transformed parameters:
  trans2param <- match(parameters_labels, transparameters_labels)

  #### Create the initial values for the parameters ####

  # Collect the unique nontransformed parameters and the unique transformed parameters:

  init_trans <- vector("list", length = control$rstarts)

  for(rs in 1:control$rstarts) {

    init_trans[[rs]] <- vector("list", length = ngroups)

    for(i in 1:ngroups) {

      init_trans[[rs]][[i]]$ulambda <- rorth(nitems[[i]], nfactors[[i]])
      init_trans[[rs]][[i]]$ulambda[fixed[[i]]$ulambda] <- fixed_values[[i]]$ulambda

      # init_trans[[rs]][[i]]$X <- rorth(nfactors[[i]], nfactors[[i]])
      init_trans[[rs]][[i]]$X <- roblq(nfactors[[i]], nfactors[[i]])

      if(positive) {

        init_trans[[rs]][[i]]$pj_theta <- rpoblq(nitems[[i]], nitems[[i]],
                                                 target = target_theta[[i]])
        init_trans[[rs]][[i]]$theta <- crossprod(init_trans[[rs]][[i]]$pj_theta)

      } else {

        U <- diag(1/diag(solve(correl[[i]]$R)))
        init_trans[[rs]][[i]]$theta <- U
        init_trans[[rs]][[i]]$theta[fixed[[i]]$theta] <- fixed_values[[i]]$theta

      }

      if(orthogonal) {
        init_trans[[rs]][[i]]$lambda <- init_trans[[rs]][[i]]$ulambda %*%
                                        init_trans[[rs]][[i]]$X
      } else {
        init_trans[[rs]][[i]]$Xinv <- solve(init_trans[[rs]][[i]]$X)
        init_trans[[rs]][[i]]$lambda <- init_trans[[rs]][[i]]$ulambda %*%
                                        t(init_trans[[rs]][[i]]$Xinv)
      }

      init_trans[[rs]][[i]]$psi <- crossprod(init_trans[[rs]][[i]]$X)

      Lambda <- init_trans[[rs]][[i]]$lambda
      Phi <- init_trans[[rs]][[i]]$psi
      Theta <- init_trans[[rs]][[i]]$theta
      init_trans[[rs]][[i]]$model <- Lambda %*% Phi %*% t(Lambda) + Theta

    }

  }

  #### Create the vectors of parameters and transformed parameters ####

  parameters <- transparameters <- vector("list", length = control$rstarts)
  # Indices of the unique transparameters in init_trans:
  trans_inds <- match(transparameters_labels, vector_trans)
  init_inds <- match(parameters_labels, vector_trans)

  for(rs in 1:control$rstarts) {

    transparameters[[rs]] <- unlist(init_trans[[rs]])[trans_inds]
    names(transparameters[[rs]]) <- transparameters_labels
    parameters[[rs]] <- unlist(init_trans[[rs]])[init_inds]
    names(parameters[[rs]]) <- parameters_labels

  }

  #### Set up the optimizer ####

  # Create defaults for the control of the optimizer:
  control$parameters <- parameters
  control$transparameters <- transparameters
  control$param2transparam <- param2trans-1L
  control$transparam2param <- trans2param-1L

  #### Return ####

  result <- list(parameters_labels = parameters_labels,
                 nparam = nparam,
                 transparameters_labels = transparameters_labels,
                 ntrans = ntrans,
                 cfa_param = cfa_param,
                 cfa_trans = cfa_trans,
                 target_psi = target_psi,
                 target_theta = target_theta,
                 fixed = fixed,
                 nonfixed = nonfixed,
                 init_trans = init_trans,
                 rest = rest,
                 control = control)

  return(result)

}

get_efa_structures <- function(data_list, full_model, control) {

  # Generate control_manifold, control_transform, and control_estimator

  list2env(data_list, envir = environment())
  list2env(full_model, envir = environment())

  #### Manifolds ####

  control_manifold <- list()

  # Euclidean for unrotated lambdas:
  ulambdas <- unlist(lapply(cfa_param, FUN = \(x) x$ulambda))
  indices_ulambda <- match(unique(ulambdas), parameters_labels)
  indices_ulambda <- indices_ulambda[!is.na(indices_ulambda)]
  labels <- parameters_labels[indices_ulambda]
  indices <- list(indices_ulambda-1L)
  control_manifold[[1]] <- list(manifold = "euclidean",
                                parameters = labels,
                                indices = indices)
  k <- 2L

  for(i in 1:ngroups) {

    if(orthogonal) {
      manifoldX <- "orth"
    } else {
      manifoldX <- "oblq"
    }

    # Orthogonal/Oblique X in each group:
    Xs <- c(cfa_param[[i]]$X)
    indices_Xs <- match(unique(Xs), parameters_labels)
    indices_Xs <- indices_Xs[!is.na(indices_Xs)]
    labels <- parameters_labels[indices_Xs]
    indices <- list(indices_Xs-1L)
    control_manifold[[k]] <- list(manifold = manifoldX,
                                  parameters = labels,
                                  indices = indices,
                                  q = nfactors[[i]])
    k <- k+1L

  }

  # Positive constraints for theta in each group:
  if(positive) {

    for(i in 1:ngroups) {

      indices_pj_theta <- match(unique(c(cfa_param[[i]]$pj_theta)),
                                parameters_labels)
      labels <- parameters_labels[indices_pj_theta]
      indices <- list(indices_pj_theta-1L)
      control_manifold[[k]] <- list(manifold = "poblq",
                                    parameters = labels,
                                    indices = indices,
                                    target = target_theta[[i]])
      k <- k+1L

    }

  } else {

    # Unconstrained thetas:
    thetas <- unlist(lapply(cfa_param, FUN = \(x) x$theta))
    indices_thetas <- match(unique(thetas), parameters_labels)
    indices_thetas <- indices_thetas[!is.na(indices_thetas)]
    labels <- parameters_labels[indices_thetas]
    indices <- list(indices_thetas-1L)
    control_manifold[[k]] <- list(manifold = "euclidean",
                                  parameters = labels,
                                  indices = indices)
    k <- k+1L

  }

  #### Transformations ####

  control_transform <- list()
  k <- 1L

  for(i in 1:ngroups) {

    if(orthogonal) {

      # Rotated lambda:
      labels_ulambda <- c(cfa_trans[[i]]$ulambda)
      labels_X <- c(cfa_trans[[i]]$X)
      labels_in <- list(labels_ulambda, labels_X)
      labels_out <- c(cfa_trans[[i]]$lambda)
      indices_in_ulambda <- match(labels_ulambda, transparameters_labels)-1L
      indices_in_X <- match(labels_X, transparameters_labels)-1L
      indices_in <- list(indices_in_ulambda, indices_in_X)
      indices_out <- list(match(labels_out, transparameters_labels)-1L)
      control_transform[[k]] <- list(transform = "XY",
                                     labels_in = labels_in,
                                     indices_in = indices_in,
                                     labels_out = labels_out,
                                     indices_out = indices_out,
                                     p = nitems[[i]],
                                     q = nfactors[[i]])
      k <- k+1L

    } else {

      # Inverse of X:
      labels_in <- c(cfa_trans[[i]]$X)
      labels_out <- c(cfa_trans[[i]]$Xinv)
      indices_in <- list(match(labels_in, transparameters_labels)-1L)
      indices_out <- list(match(labels_out, transparameters_labels)-1L)
      control_transform[[k]] <- list(transform = "matrix_inverse",
                                     labels_in = labels_in,
                                     indices_in = indices_in,
                                     labels_out = labels_out,
                                     indices_out = indices_out,
                                     p = nfactors[[i]])
      k <- k+1L

      # Rotated lambda:
      labels_ulambda <- c(cfa_trans[[i]]$ulambda)
      labels_Xinv <- c(cfa_trans[[i]]$Xinv)
      labels_in <- list(labels_ulambda, labels_Xinv)
      labels_out <- c(cfa_trans[[i]]$lambda)
      indices_in_ulambda <- match(labels_ulambda, transparameters_labels)-1L
      indices_in_Xinv <- match(labels_Xinv, transparameters_labels)-1L
      indices_in <- list(indices_in_ulambda, indices_in_Xinv)
      indices_out <- list(match(labels_out, transparameters_labels)-1L)
      control_transform[[k]] <- list(transform = "XYt",
                                     labels_in = labels_in,
                                     indices_in = indices_in,
                                     labels_out = labels_out,
                                     indices_out = indices_out,
                                     p = nitems[[i]],
                                     q = nfactors[[i]])
      k <- k+1L

    }

    # Latent covariances:
    lower_psi <- lower.tri(cfa_trans[[i]]$psi, diag = TRUE)
    labels_in <- c(cfa_trans[[i]]$X)
    indices_in <- list(match(labels_in, transparameters_labels)-1L)
    labels_out <- c(cfa_trans[[i]]$psi[lower_psi])
    indices_out <- list(match(labels_out, transparameters_labels)-1L)
    control_transform[[k]] <- list(transform = "crossprod",
                                   labels_in = labels_in,
                                   indices_in = indices_in,
                                   labels_out = labels_out,
                                   indices_out = indices_out,
                                   p = nfactors[[i]],
                                   q = nfactors[[i]])
    k <- k+1L

    # Error covariances:
    if(positive) {

      lower_theta <- lower.tri(cfa_trans[[i]]$theta, diag = TRUE)
      labels_in <- c(cfa_trans[[i]]$pj_theta)
      labels_out <- c(cfa_trans[[i]]$theta[lower_theta])
      indices_in <- list(match(labels_in, transparameters_labels)-1L)
      indices_out <- list(match(labels_out, transparameters_labels)-1L)
      control_transform[[k]] <- list(transform = "crossprod",
                                     labels_in = labels_in,
                                     indices_in = indices_in,
                                     labels_out = labels_out,
                                     indices_out = indices_out,
                                     p = nitems[[i]],
                                     q = nitems[[i]])
      k <- k+1L

    }

    # Model matrix correlation:
    lower_theta <- lower.tri(cfa_trans[[i]]$theta, diag = TRUE)
    indices_lambda <- match(cfa_trans[[i]]$lambda,
                            transparameters_labels)
    indices_psi <- match(cfa_trans[[i]]$psi[lower_psi],
                         transparameters_labels)
    indices_theta <- match(cfa_trans[[i]]$theta[lower_theta],
                           transparameters_labels)
    indices_all <- c(indices_lambda, indices_psi, indices_theta)
    indices_in <- list(indices = indices_all-1L,
                       indices_lambda = indices_lambda-1L,
                       indices_psi = indices_psi-1L,
                       indices_theta = indices_theta-1L)
    labels_in <- transparameters_labels[indices_all]

    lower_diag <- lower.tri(cfa_trans[[i]]$model, diag = TRUE)
    indices_model <- match(cfa_trans[[i]]$model[lower_diag],
                           transparameters_labels)
    labels_out <- transparameters_labels[indices_model]
    indices_out <- list(indices = indices_model-1L)

    control_transform[[k]] <- list(transform = "factor_cor",
                                   labels_in = labels_in,
                                   indices_in = indices_in,
                                   labels_out = labels_out,
                                   indices_out = indices_out,
                                   p = nitems[[i]],
                                   q = nfactors[[i]])
    k <- k+1L

  }

  #### Estimators ####

  control_estimator <- list()
  k <- 1L

  for(i in 1:ngroups) {

    # Rotation criteria:
    labels <- c(cfa_trans[[i]]$lambda)
    indices <- list(match(labels, transparameters_labels) - 1L)
    control_estimator[[k]] <- list(estimator = "oblimin",
                                   labels = labels,
                                   indices = indices,
                                   gamma = 0.00,
                                   p = nitems[[i]],
                                   q = nfactors[[i]])
    k <- k+1L

    lower_diag <- lower.tri(cfa_trans[[i]]$model, diag = TRUE)
    indices_model <- match(cfa_trans[[i]]$model[lower_diag],
                         transparameters_labels)
    labels <- transparameters_labels[indices_model]
    indices <- list(indices = indices_model-1L)

    estimator <- tolower(estimator)
    if(estimator == "uls" || estimator == "dwls") {
      cfa_estimator <- "cfa_dwls"
    } else if(estimator == "ml") {
      cfa_estimator <- "cfa_ml"
    } else if(estimator == "ml2") {
      cfa_estimator <- "cfa_ml2"
    } else if(estimator == "mlR") {
      cfa_estimator <- "cfa_ml_R"
    }

    control_estimator[[k]] <- list(estimator = cfa_estimator,
                                   labels = labels,
                                   indices = indices,
                                   R = correl[[i]]$R,
                                   W = correl[[i]]$W,
                                   w = nobs[[i]] / sum(unlist(nobs)),
                                   q = nfactors[[i]],
                                   n = nobs[[i]])
    k <- k+1L

  }

  if(positive & control$reg) {

    for(i in 1:ngroups) {

      # For the psi matrix:

      lower_indices <- which(lower.tri(cfa_trans[[i]]$psi, diag = TRUE))
      labels <- cfa_trans[[i]]$psi[lower_indices]
      indices <- match(labels, transparameters_labels)
      control_estimator[[k]] <- list(estimator = "logdetmat",
                                     labels = labels,
                                     indices = list(indices-1L),
                                     lower_indices = lower_indices-1L,
                                     p = nrow(cfa_trans[[i]]$psi),
                                     logdetw = control$penalties$logdet$w)
      k <- k+1L

      # For the theta matrix:

      lower_indices <- which(lower.tri(cfa_trans[[i]]$theta, diag = TRUE))
      labels <- cfa_trans[[i]]$theta[lower_indices]
      indices <- match(labels, transparameters_labels)
      control_estimator[[k]] <- list(estimator = "logdetmat",
                                     labels = labels,
                                     indices = list(indices-1L),
                                     lower_indices = lower_indices-1L,
                                     p = nrow(cfa_trans[[i]]$theta),
                                     logdetw = control$penalties$logdet$w)
      k <- k+1L

    }

  }

  #### Return ####

  result <- list(control_manifold = control_manifold,
                 control_transform = control_transform,
                 control_estimator = control_estimator)

  return(result)

}
