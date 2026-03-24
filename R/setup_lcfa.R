# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 16/03/2026

get_full_cfa_model <- function(data_list, model, control = NULL) {

  # Generate the model syntax and initial parameter values

  list2env(data_list, envir = environment())

  # Initialize the objects to store the initial parameters:
  param <- trans <- vector("list")
  fixed <- fixed_values <- nonfixed <- vector("list")

  #### Parameters of the model for each group ####

  # Initialize the target matrices for positive-definite constraints:
  target_psi <- target_theta <- targets <- vector("list", length = ngroups)
  rest <- 0L

  S_group <- paste("S.g", 1:ngroups, sep = "")
  lambda_group <- paste("lambda.g", 1:ngroups, sep = "")
  psi_group <- paste("psi.g", 1:ngroups, sep = "")
  theta_group <- paste("theta.g", 1:ngroups, sep = "")
  xpsi_group <- paste("xpsi.g", 1:ngroups, sep = "")
  xtheta_group <- paste("xtheta.g", 1:ngroups, sep = "")
  model_group <- paste("model.g", 1:ngroups, sep = "")

  # Transformed parameters:
  list_struct <- vector("list")
  k <- 1L
  for(i in 1:ngroups) {

    # Lambda:
    list_struct[[k]] <- list(name = lambda_group[i],
                             type = "matrix",
                             dim = c(nitems[[i]], nfactors[[i]]),
                             rownames = item_label[[i]],
                             colnames = factor_label[[i]])
    k <- k+1L

    # Create additional parameters if there are positive-definite constraints:
    if(positive) {

      # Theta:
      list_struct[[k]] <- list(name = xtheta_group[i],
                               type = "matrix",
                               dim = c(nitems[[i]], nitems[[i]]),
                               rownames = item_label[[i]],
                               colnames = item_label[[i]])
      k <- k+1L

      # Psi:
      list_struct[[k]] <- list(name = xpsi_group[i],
                               type = "matrix",
                               dim = c(nfactors[[i]], nfactors[[i]]),
                               rownames = factor_label[[i]],
                               colnames = factor_label[[i]])
      k <- k+1L

    }

    # Theta:
    list_struct[[k]] <- list(name = theta_group[i],
                             type = "matrix",
                             dim = c(nitems[[i]], nitems[[i]]),
                             rownames = item_label[[i]],
                             colnames = item_label[[i]],
                             symmetric = TRUE)
    k <- k+1L

    # Psi:
    list_struct[[k]] <- list(name = psi_group[i],
                             type = "matrix",
                             dim = c(nfactors[[i]], nfactors[[i]]),
                             rownames = factor_label[[i]],
                             colnames = factor_label[[i]],
                             symmetric = TRUE)
    k <- k+1L

    # Model matrix:
    list_struct[[k]] <- list(name = model_group[i],
                             type = "matrix",
                             dim = c(nitems[[i]], nitems[[i]]),
                             rownames = item_label[[i]],
                             colnames = item_label[[i]],
                             symmetric = TRUE)
    k <- k+1L

    # S:
    list_struct[[k]] <- list(name = S_group[i],
                             type = "matrix",
                             dim = c(nitems[[i]], nitems[[i]]),
                             rownames = item_label[[i]],
                             colnames = item_label[[i]],
                             symmetric = TRUE)
    k <- k+1L

  }

  trans <- create_parameters(list_struct)

  # Replace latent labels by lavaan labels:
  for(i in 1:ngroups) {

    # Get the positions of parameters and fixed values:

    group_i <- c(lambda_group[i], psi_group[i], theta_group[i])

    nonfixed[group_i] <- lapply(model[[i]], FUN = \(x) {
      which(is.na(suppressWarnings(as.numeric(x))))
    })

    fixed[group_i] <- lapply(model[[i]], FUN = \(x) {
      which(!is.na(suppressWarnings(as.numeric(x))))
    })

    fixed_values[group_i] <- lapply(model[[i]], FUN = \(x) {
      numerals <- suppressWarnings(as.numeric(x))
      inds <- which(!is.na(numerals))
      return(numerals[inds])
    })

    trans[[lambda_group[i]]][nonfixed[[lambda_group[i]]]] <- model[[i]]$lambda[nonfixed[[lambda_group[i]]]]
    trans[[theta_group[i]]][nonfixed[[theta_group[i]]]] <- model[[i]]$theta[nonfixed[[theta_group[i]]]]
    trans[[psi_group[i]]][nonfixed[[psi_group[i]]]] <- model[[i]]$psi[nonfixed[[psi_group[i]]]]

  }

  # Untransformed parameters:

  for(i in 1:ngroups) {

    param[[lambda_group[i]]] <- trans[[lambda_group[i]]]
    # Insert fixed values in the model:
    param[[lambda_group[i]]][fixed[[lambda_group[i]]]] <- model[[i]]$lambda[fixed[[lambda_group[i]]]]

    if(positive) {

      # Theta:
      param[[xtheta_group[i]]] <- trans[[xtheta_group[i]]]
      # Psi:
      param[[xpsi_group[i]]] <- trans[[xpsi_group[i]]]

    } else {

      # Theta:
      param[[theta_group[i]]] <- trans[[theta_group[i]]]
      # Insert fixed values in the model:
      param[[theta_group[i]]][fixed[[theta_group[i]]]] <- model[[i]]$theta[fixed[[theta_group[i]]]]
      if(control$deltaparam) {
        diag(param[[theta_group[i]]]) <- "1"
      }

      # Psi:
      param[[psi_group[i]]] <- trans[[psi_group[i]]]
      # Insert fixed values in the model:
      param[[psi_group[i]]][fixed[[psi_group[i]]]] <- model[[i]]$psi[fixed[[psi_group[i]]]]

    }

    # S:
    if(control$free_S) {
      param[[S_group[i]]] <- trans[[S_group[i]]]
    } else {
      param[[S_group[i]]] <- correl[[i]]$R
    }

    if(!control$free_S_diag) {
      diag(param[[S_group[i]]]) <- "1"
    }

    # Create the target matrices for positive-definite constraints:
    if(positive) {

      target_theta[[i]] <- matrix(0, nrow = nitems[[i]], ncol = nitems[[i]])
      target_theta[[i]][nonfixed[[theta_group[i]]]] <- 1

      target_psi[[i]] <- matrix(0, nrow = nfactors[[i]], ncol = nfactors[[i]])
      target_psi[[i]][nonfixed[[psi_group[i]]]] <- 1

      q <- nfactors[[i]]
      p <- nitems[[i]]
      lower_theta <- lower.tri(diag(p), diag = TRUE)
      lower_psi <- lower.tri(diag(q), diag = TRUE)
      targets[[i]] <- unlist(c(target_theta[[i]][lower_theta],
                               target_psi[[i]][lower_psi]))
      rest <- rest + 0.5*q*(q-1) + 0.5*p*(p-1) + sum(targets[[i]] == 0)

    }
  }

  #### Arrange labels ####

  # Arrange parameter labels:
  vector_param <- unname(unique(unlist(param)))

  # Select the unique, nonnumeric labels:
  nonfixed_pars <- which(is.na(suppressWarnings(as.numeric(vector_param))))
  parameters_labels <- vector_param[nonfixed_pars]
  nparam <- length(parameters_labels)

  # Arrange transparameter labels:
  vector_trans <- unname(unlist(trans))
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

    init_trans[[rs]] <- vector("list")

    for(i in 1:ngroups) {

      init_trans[[rs]][[lambda_group[i]]] <- rorth(nitems[[i]], nfactors[[i]])
      init_trans[[rs]][[lambda_group[i]]][fixed[[lambda_group[i]]]] <- fixed_values[[lambda_group[i]]]

      if(positive) {

        init_trans[[rs]][[xtheta_group[i]]] <- rpoblq(nitems[[i]], nitems[[i]], constraints = target_theta[[i]])
        init_trans[[rs]][[xpsi_group[i]]] <- rpoblq(nfactors[[i]], nfactors[[i]], constraints = target_psi[[i]])

        init_trans[[rs]][[theta_group[i]]] <- crossprod(init_trans[[rs]][[xtheta_group[i]]])
        init_trans[[rs]][[psi_group[i]]] <- crossprod(init_trans[[rs]][[xpsi_group[i]]])

      } else {

        U <- diag(1/diag(solve(correl[[i]]$R)))
        init_trans[[rs]][[theta_group[i]]] <- U
        init_trans[[rs]][[theta_group[i]]][fixed[[theta_group[i]]]] <- fixed_values[[theta_group[i]]]

        P <- diag(nfactors[[i]])
        init_trans[[rs]][[psi_group[i]]] <- P
        init_trans[[rs]][[psi_group[i]]][fixed[[psi_group[i]]]] <- fixed_values[[psi_group[i]]]

      }

      Lambda <- init_trans[[rs]][[lambda_group[i]]]
      Theta <- init_trans[[rs]][[theta_group[i]]]
      Psi <- init_trans[[rs]][[psi_group[i]]]
      init_trans[[rs]][[model_group[i]]] <- Lambda %*% Psi %*% t(Lambda) + Theta
      init_trans[[rs]][[S_group[i]]] <- correl[[i]]$R

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
  control$init <- init_trans
  control$param2transparam <- param2trans-1L
  control$transparam2param <- trans2param-1L

  #### Return ####

  result <- list(parameters_labels = parameters_labels,
                 nparam = nparam,
                 transparameters_labels = transparameters_labels,
                 ntrans = ntrans,
                 param = param,
                 trans = trans,
                 target_psi = target_psi,
                 target_theta = target_theta,
                 fixed = fixed,
                 nonfixed = nonfixed,
                 rest = rest,
                 control = control)

  return(result)

}

get_cfa_structures <- function(data_list, full_model, control) {

  # Generate control_manifold, control_transform, and control_estimator

  list2env(data_list, envir = environment())
  list2env(full_model, envir = environment())

  S_group <- paste("S.g", 1:ngroups, sep = "")
  lambda_group <- paste("lambda.g", 1:ngroups, sep = "")
  psi_group <- paste("psi.g", 1:ngroups, sep = "")
  theta_group <- paste("theta.g", 1:ngroups, sep = "")
  xpsi_group <- paste("xpsi.g", 1:ngroups, sep = "")
  xtheta_group <- paste("xtheta.g", 1:ngroups, sep = "")
  model_group <- paste("model.g", 1:ngroups, sep = "")

  #### Manifolds ####

  manifolds <- list()
  k <- 1L

  for(i in 1:ngroups) {

    manifolds[[k]] <- list(manifold = "euclidean",
                           parameters = c(lambda_group[i],
                                          S_group[i]))
    k <- k+1L

    if(positive) {

      manifolds[[k]] <- list(manifold = "poblq",
                             parameters = xpsi_group[i],
                             extra = list(
                               p = nfactors[[i]],
                               q = nfactors[[i]],
                               constraints = target_psi[[i]]))
      k <- k+1L

      manifolds[[k]] <- list(manifold = "poblq",
                             parameters = xtheta_group[i],
                             extra = list(
                               p = nitems[[i]],
                               q = nitems[[i]],
                               constraints = target_theta[[i]]))
      k <- k+1L

    } else {

      manifolds[[k]] <- list(manifold = "euclidean",
                             parameters = psi_group[i])
      k <- k+1L

      manifolds[[k]] <- list(manifold = "euclidean",
                             parameters = theta_group[i])
      k <- k+1L

    }

  }

  control_manifold <- create_manifolds(manifolds = manifolds,
                                       structures = param)

  #### Transformations ####

  transforms <- list()
  dots <- list()
  k <- 1L

  for(i in 1:ngroups) {

    if(positive) {

      lower_psi <- lower.tri(trans[[psi_group[i]]], diag = TRUE)
      lower_theta <- lower.tri(trans[[theta_group[i]]], diag = TRUE)

      dots$p <- nrow(trans[[psi_group[i]]])
      transforms[[k]] <- list(transform = "crossprod",
                             parameters_in = xpsi_group[i],
                             parameters_out = psi_group[i],
                             extra = dots)
      k <- k+1L

      dots$p <- nrow(trans[[theta_group[i]]])
      transforms[[k]] <- list(transform = "crossprod",
                              parameters_in = xtheta_group[i],
                              parameters_out = theta_group[i],
                              extra = dots)
      k <- k+1L

    }

    if(control$deltaparam) {

      dots$p <- nrow(trans[[theta_group[i]]])
      dots$q <- nrow(trans[[psi_group[i]]])
      transforms[[k]] <- list(transform = "deltaparam",
                              parameters_in = c(lambda_group[i],
                                                psi_group[i]),
                              parameters_out = list(diag(trans[[theta_group[i]]])),
                              extra = dots)
      k <- k+1L

    }

    # Model matrix correlation:
    lower_psi <- lower.tri(trans[[psi_group[i]]], diag = TRUE)
    lower_theta <- lower.tri(trans[[theta_group[i]]], diag = TRUE)
    lower_diag <- lower.tri(trans[[model_group[i]]], diag = TRUE)

    dots$p <- nrow(correl[[i]]$R)
    dots$q <- nrow(trans[[psi_group[i]]])
    transforms[[k]] <- list(transform = "factor_cor",
                            parameters_in = c(lambda_group[i],
                                              psi_group[i],
                                              theta_group[i]),
                            parameters_out = model_group[i],
                            extra = dots)
    k <- k+1L

  }

  control_transform <- get_transforms(transforms = transforms,
                                      structures = trans)

  #### Estimators ####

  estimators <- list()
  k <- 1L

  for(i in 1:ngroups) {

    estimator <- tolower(estimator)
    cfa_estimator <- switch(estimator,
                            uls = "cfa_dwls",
                            dwls  = "cfa_dwls",
                            ulsr = "cfa_dwls_error",
                            dwlsr  = "cfa_dwls_error",
                            ml = "cfa_fml",
                            fml  = "cfa_fml",
                            mlr = "cfa_fml_error",
                            fmlr  = "cfa_fml_error",
                            stop("Unknown estimator: ", estimator)
    )

    lower_diag <- lower.tri(trans[[model_group[i]]], diag = TRUE)
    estimators[[k]] <- list(estimator = cfa_estimator,
                            parameters = c(model_group[i], S_group[i]),
                            extra = list(R = correl[[i]]$R,
                                         W = correl[[i]]$W,
                                         w = nobs[[i]] / sum(unlist(nobs)),
                                         q = nrow(trans[[psi_group[i]]]),
                                         p = nitems[[i]],
                                         n = nobs[[i]]))
    k <- k+1L

  }

  if(positive & control$reg) {

    for(i in 1:ngroups) {

      # For the psi matrix:

      lower_indices <- which(lower.tri(trans[[psi_group[i]]], diag = TRUE))
      estimators[[k]] <- list(estimator = "logdetmat",
                              parameters = psi_group[i],
                              extra = list(lower_indices = lower_indices-1L,
                                           p = nrow(trans[[psi_group[i]]]),
                                           logdetw = control$penalties$logdet$w))
      k <- k+1L

      # For the theta matrix:

      lower_indices <- which(lower.tri(trans[[theta_group[i]]], diag = TRUE))
      estimators[[k]] <- list(estimator = "logdetmat",
                              parameters = theta_group[i],
                              extra = list(lower_indices = lower_indices-1L,
                                           p = nrow(trans[[theta_group[i]]]),
                                           logdetw = control$penalties$logdet$w))
      k <- k+1L

    }

  }

  control_estimator <- get_estimators(estimators = estimators,
                                      structures = trans)

  #### Return ####

  result <- list(control_manifold = control_manifold,
                 control_transform = control_transform,
                 control_estimator = control_estimator)

  return(result)

}
