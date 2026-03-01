# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 01/11/2025

get_full_cfa_model <- function(data_list, model, control = NULL) {

  # Generate the model syntax and initial parameter values

  list2env(data_list, envir = environment())

  # Initialize the objects to store the initial parameters:
  cfa_param <- cfa_trans <- vector("list")
  fixed <- fixed_values <- nonfixed <- vector("list")

  #### Parameters of the model for each group ####

  # Initialize the target matrices for positive-definite constraints:
  target_psi <- target_theta <- targets <- vector("list", length = ngroups)
  rest <- 0L

  lambda_group <- paste("lambda.group", 1:ngroups, sep = "")
  psi_group <- paste("psi.group", 1:ngroups, sep = "")
  theta_group <- paste("theta.group", 1:ngroups, sep = "")
  pj_psi_group <- paste("pj_psi.group", 1:ngroups, sep = "")
  pj_theta_group <- paste("pj_theta.group", 1:ngroups, sep = "")
  model_group <- paste("model.group", 1:ngroups, sep = "")

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

    # Transformed parameters:

    # Lambda:
    lambda_labels <- paste("g", i, ".lambda[", rep(1:nitems[[i]], times = nfactors[[i]]),
                           ",", rep(1:nfactors[[i]], each = nitems[[i]]), "]", sep = "")
    lambda_labels[nonfixed[[lambda_group[i]]]] <- model[[i]]$lambda[nonfixed[[lambda_group[i]]]]
    cfa_trans[[lambda_group[i]]] <- matrix(lambda_labels, nrow = nitems[[i]], ncol = nfactors[[i]])

    # Create additional parameters if there are positive-definite constraints:
    if(positive) {

      # Psi:
      pj_psi_labels <- paste("g", i, ".pj_psi[", rep(1:nfactors[[i]], times = nfactors[[i]]),
                             ",", rep(1:nfactors[[i]], each = nfactors[[i]]), "]", sep = "")
      cfa_trans[[pj_psi_group[i]]] <- matrix(pj_psi_labels, nrow = nfactors[[i]], ncol = nfactors[[i]])

      # Theta:
      pj_theta_labels <- paste("g", i, ".pj_theta[", rep(1:nitems[[i]], times = nitems[[i]]),
                               ",", rep(1:nitems[[i]], each = nitems[[i]]), "]", sep = "")
      cfa_trans[[pj_theta_group[i]]] <- matrix(pj_theta_labels, nrow = nitems[[i]], ncol = nitems[[i]])

    }

    # Psi:
    psi_labels <- paste("g", i, ".psi[", rep(1:nfactors[[i]], times = nfactors[[i]]),
                        ",", rep(1:nfactors[[i]], each = nfactors[[i]]), "]", sep = "")
    psi_labels[nonfixed[[psi_group[i]]]] <- model[[i]]$psi[nonfixed[[psi_group[i]]]]
    cfa_trans[[psi_group[i]]] <- matrix(psi_labels, nrow = nfactors[[i]], ncol = nfactors[[i]])
    # Force symmetry:
    cfa_trans[[psi_group[i]]][upper.tri(cfa_trans[[psi_group[i]]])] <- t(cfa_trans[[psi_group[i]]])[upper.tri(cfa_trans[[psi_group[i]]])]

    # Theta:
    theta_labels <- paste("g", i, ".theta[", rep(1:nitems[[i]], times = nitems[[i]]),
                          ",", rep(1:nitems[[i]], each = nitems[[i]]), "]", sep = "")
    theta_labels[nonfixed[[theta_group[i]]]] <- model[[i]]$theta[nonfixed[[theta_group[i]]]]
    cfa_trans[[theta_group[i]]] <- matrix(theta_labels, nrow = nitems[[i]], ncol = nitems[[i]])
    # Force symmetry:
    cfa_trans[[theta_group[i]]][upper.tri(cfa_trans[[theta_group[i]]])] <- t(cfa_trans[[theta_group[i]]])[upper.tri(cfa_trans[[theta_group[i]]])]

    # Model matrix:
    model_labels <- paste("g", i, ".model[", rep(1:nitems[[i]], times = nitems[[i]]),
                          ",", rep(1:nitems[[i]], each = nitems[[i]]), "]", sep = "")
    cfa_trans[[model_group[i]]] <- matrix(model_labels, nrow = nitems[[i]], ncol = nitems[[i]])
    # Force symmetry:
    cfa_trans[[model_group[i]]][upper.tri(cfa_trans[[model_group[i]]])] <- t(cfa_trans[[model_group[i]]])[upper.tri(cfa_trans[[model_group[i]]])]

    # Untransformed parameters:

    cfa_param[[lambda_group[i]]] <- cfa_trans[[lambda_group[i]]]
    # Insert fixed values in the model:
    cfa_param[[lambda_group[i]]][fixed[[lambda_group[i]]]] <- model[[i]]$lambda[fixed[[lambda_group[i]]]]

    if(positive) {

      # Psi:
      cfa_param[[pj_psi_group[i]]] <- cfa_trans[[pj_psi_group[i]]]
      # Theta:
      cfa_param[[pj_theta_group[i]]] <- cfa_trans[[pj_theta_group[i]]]

    } else {

      # Psi:
      cfa_param[[psi_group[i]]] <- cfa_trans[[psi_group[i]]]
      # Insert fixed values in the model:
      cfa_param[[psi_group[i]]][fixed[[psi_group[i]]]] <- model[[i]]$psi[fixed[[psi_group[i]]]]

      # Theta:
      cfa_param[[theta_group[i]]] <- cfa_trans[[theta_group[i]]]
      # Insert fixed values in the model:
      cfa_param[[theta_group[i]]][fixed[[theta_group[i]]]] <- model[[i]]$theta[fixed[[theta_group[i]]]]

    }

    # Create the target matrices for positive-definite constraints:
    if(positive) {

      target_psi[[i]] <- matrix(0, nrow = nfactors[[i]], ncol = nfactors[[i]])
      target_psi[[i]][nonfixed[[psi_group[i]]]] <- 1
      target_theta[[i]] <- matrix(0, nrow = nitems[[i]], ncol = nitems[[i]])
      target_theta[[i]][nonfixed[[theta_group[i]]]] <- 1

      q <- nfactors[[i]]
      p <- nitems[[i]]
      lower_psi <- lower.tri(diag(q), diag = TRUE)
      lower_theta <- lower.tri(diag(p), diag = TRUE)
      targets[[i]] <- unlist(c(target_psi[[i]][lower_psi],
                               target_theta[[i]][lower_theta]))
      rest <- rest + 0.5*q*(q-1) + 0.5*p*(p-1) + sum(targets[[i]] == 0)

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

    init_trans[[rs]] <- vector("list")

    for(i in 1:ngroups) {

      init_trans[[rs]][[lambda_group[i]]] <- rorth(nitems[[i]], nfactors[[i]])
      init_trans[[rs]][[lambda_group[i]]][fixed[[lambda_group[i]]]] <- fixed_values[[lambda_group[i]]]

      if(positive) {

        init_trans[[rs]][[pj_psi_group[i]]] <- rpoblq(nfactors[[i]], nfactors[[i]], constraints = target_psi[[i]])
        init_trans[[rs]][[pj_theta_group[i]]] <- rpoblq(nitems[[i]], nitems[[i]], constraints = target_theta[[i]])

        init_trans[[rs]][[psi_group[i]]] <- crossprod(init_trans[[rs]][[pj_psi_group[i]]])

        init_trans[[rs]][[theta_group[i]]] <- crossprod(init_trans[[rs]][[pj_theta_group[i]]])

      } else {

        P <- diag(nfactors[[i]])
        init_trans[[rs]][[psi_group[i]]] <- P
        init_trans[[rs]][[psi_group[i]]][fixed[[psi_group[i]]]] <- fixed_values[[psi_group[i]]]

        U <- diag(1/diag(solve(correl[[i]]$R)))
        init_trans[[rs]][[theta_group[i]]] <- U
        init_trans[[rs]][[theta_group[i]]][fixed[[theta_group[i]]]] <- fixed_values[[theta_group[i]]]

      }

      Lambda <- init_trans[[rs]][[lambda_group[i]]]
      Phi <- init_trans[[rs]][[psi_group[i]]]
      Theta <- init_trans[[rs]][[theta_group[i]]]
      init_trans[[rs]][[model_group[i]]] <- Lambda %*% Phi %*% t(Lambda) + Theta

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

get_cfa_structures <- function(data_list, full_model, control) {

  # Generate control_manifold, control_transform, and control_estimator

  list2env(data_list, envir = environment())
  list2env(full_model, envir = environment())

  lambda_group <- paste("lambda.group", 1:ngroups, sep = "")
  psi_group <- paste("psi.group", 1:ngroups, sep = "")
  theta_group <- paste("theta.group", 1:ngroups, sep = "")
  pj_psi_group <- paste("pj_psi.group", 1:ngroups, sep = "")
  pj_theta_group <- paste("pj_theta.group", 1:ngroups, sep = "")
  model_group <- paste("model.group", 1:ngroups, sep = "")

  #### Manifolds ####

  # control_manifold <- list()
  mani_and_labs <- list()
  k <- 1L

  for(i in 1:ngroups) {

    mani_and_labs[[k]] <- list("euclidean", cfa_param[[lambda_group[i]]])
    k <- k+1L

    if(positive) {

      mani_and_labs[[k]] <- list("poblq",
                                 cfa_param[[pj_psi_group[i]]],
                                 p = nfactors[[i]],
                                 q = nfactors[[i]],
                                 constraints = target_psi[[i]])
      k <- k+1L

      mani_and_labs[[k]] <- list("poblq",
                                 cfa_param[[pj_theta_group[i]]],
                                 p = nitems[[i]],
                                 q = nitems[[i]],
                                 constraints = target_theta[[i]])
      k <- k+1L

    } else {

      mani_and_labs[[k]] <- list("euclidean",
                                 cfa_param[[psi_group[i]]])
      k <- k+1L

      mani_and_labs[[k]] <- list("euclidean",
                                 cfa_param[[theta_group[i]]])
      k <- k+1L

    }

  }

  control_manifold <- create_manifolds(manifolds_and_labels = mani_and_labs,
                                       param_structures = cfa_param)

  #### Transformations ####

  trans_and_labs <- list()
  dots <- list()
  k <- 1L

  for(i in 1:ngroups) {

    if(positive) {

      lower_psi <- lower.tri(cfa_trans[[psi_group[i]]], diag = TRUE)
      lower_theta <- lower.tri(cfa_trans[[theta_group[i]]], diag = TRUE)

      dots$p <- nrow(cfa_trans[[psi_group[i]]])
      trans_and_labs[[k]] <- extra_transforms(transform = "crossprod",
                                              labels_in = list(cfa_trans[[pj_psi_group[i]]]),
                                              labels_out = list(cfa_trans[[psi_group[i]]][lower_psi]),
                                              dots)
      k <- k+1L

      dots$p <- nrow(cfa_trans[[theta_group[i]]])
      trans_and_labs[[k]] <- extra_transforms(transform = "crossprod",
                                              labels_in = list(cfa_trans[[pj_theta_group[i]]]),
                                              labels_out = list(cfa_trans[[theta_group[i]]][lower_theta]),
                                              dots)
      k <- k+1L

    }

    # Model matrix correlation:
    lower_psi <- lower.tri(cfa_trans[[psi_group[i]]], diag = TRUE)
    lower_theta <- lower.tri(cfa_trans[[theta_group[i]]], diag = TRUE)
    lower_diag <- lower.tri(cfa_trans[[model_group[i]]], diag = TRUE)

    dots$p <- nrow(correl[[i]]$R)
    dots$q <- nrow(cfa_trans[[psi_group[i]]])
    trans_and_labs[[k]] <- extra_transforms(transform = "factor_cor",
                                            labels_in = list(cfa_trans[[lambda_group[i]]],
                                                             cfa_trans[[psi_group[i]]][lower_psi],
                                                             cfa_trans[[theta_group[i]]][lower_theta]),
                                            labels_out = list(cfa_trans[[model_group[i]]][lower_diag]),
                                            dots)
    k <- k+1L

  }

  control_transform <- create_transforms(transforms_and_labels = trans_and_labs,
                                         param_structures = cfa_trans)

  #### Estimators ####

  control_estimator <- list()
  k <- 1L

  for(i in 1:ngroups) {

    lower_diag <- lower.tri(cfa_trans[[model_group[i]]], diag = TRUE)
    indices <- list(match(cfa_trans[[model_group[i]]][lower_diag],
                     transparameters_labels)-1L)

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
                                   # labels = labels,
                                   indices = indices,
                                   R = correl[[i]]$R,
                                   W = correl[[i]]$W,
                                   w = nobs[[i]] / sum(unlist(nobs)),
                                   q = nrow(cfa_trans[[psi_group[i]]]),
                                   n = nobs[[i]])
    k <- k+1L

  }

  if(positive & control$reg) {

    for(i in 1:ngroups) {

      # For the psi matrix:

      lower_indices <- which(lower.tri(cfa_trans[[psi_group[i]]], diag = TRUE))
      labels <- cfa_trans[[psi_group[i]]][lower_indices]
      indices <- list(match(labels, transparameters_labels)-1L)
      control_estimator[[k]] <- list(estimator = "logdetmat",
                                     # labels = labels,
                                     indices = indices,
                                     lower_indices = lower_indices-1L,
                                     p = nrow(cfa_trans[[psi_group[i]]]),
                                     logdetw = control$penalties$logdet$w)
      k <- k+1L

      # For the theta matrix:

      lower_indices <- which(lower.tri(cfa_trans[[theta_group[i]]], diag = TRUE))
      labels <- cfa_trans[[theta_group[i]]][lower_indices]
      indices <- list(match(labels, transparameters_labels)-1L)
      control_estimator[[k]] <- list(estimator = "logdetmat",
                                     # labels = labels,
                                     indices = indices,
                                     lower_indices = lower_indices-1L,
                                     p = nrow(cfa_trans[[theta_group[i]]]),
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
