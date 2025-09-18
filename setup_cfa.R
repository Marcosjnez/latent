# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 11/09/2025

get_full_cfa_model <- function(data_list, model = NULL, control = NULL) {

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

    # Lambda:
    lambda_labels <- paste("g", i, ".lambda[", rep(1:nitems, times = nfactors),
                           ",", rep(1:nfactors, each = nitems), "]", sep = "")
    lambda_labels[nonfixed[[i]]$lambda] <- model[[i]]$lambda[nonfixed[[i]]$lambda]
    cfa_trans[[i]]$lambda <- matrix(lambda_labels, nrow = nitems, ncol = nfactors)

    # Create additional parameters if there are positive-definite constraints:
    if(positive) {

      # Psi:
      pj_psi_labels <- paste("g", i, ".pj_psi[", rep(1:nfactors, times = nfactors),
                             ",", rep(1:nfactors, each = nfactors), "]", sep = "")
      cfa_trans[[i]]$pj_psi <- matrix(pj_psi_labels, nrow = nfactors, ncol = nfactors)

      # Theta:
      pj_theta_labels <- paste("g", i, ".pj_theta[", rep(1:nitems, times = nitems),
                               ",", rep(1:nitems, each = nitems), "]", sep = "")
      cfa_trans[[i]]$pj_theta <- matrix(pj_theta_labels, nrow = nitems, ncol = nitems)

    }

    # Psi:
    psi_labels <- paste("g", i, ".psi[", rep(1:nfactors, times = nfactors),
                        ",", rep(1:nfactors, each = nfactors), "]", sep = "")
    psi_labels[nonfixed[[i]]$psi] <- model[[i]]$psi[nonfixed[[i]]$psi]
    cfa_trans[[i]]$psi <- matrix(psi_labels, nrow = nfactors, ncol = nfactors)
    # Force symmetry:
    # cfa_trans[[i]]$psi[upper.tri(cfa_trans[[i]]$psi)] <- cfa_trans[[i]]$psi[lower.tri(cfa_trans[[i]]$psi)]

    # Theta:
    theta_labels <- paste("g", i, ".theta[", rep(1:nitems, times = nitems),
                          ",", rep(1:nitems, each = nitems), "]", sep = "")
    theta_labels[nonfixed[[i]]$theta] <- model[[i]]$theta[nonfixed[[i]]$theta]
    cfa_trans[[i]]$theta <- matrix(theta_labels, nrow = nitems, ncol = nitems)
    # Force symmetry:
    # cfa_trans[[i]]$theta[upper.tri(cfa_trans[[i]]$theta)] <- cfa_trans[[i]]$theta[lower.tri(cfa_trans[[i]]$theta)]

    # Untransformed parameters:

    cfa_param[[i]]$lambda <- cfa_trans[[i]]$lambda
    # Insert fixed values in the model:
    cfa_param[[i]]$lambda[fixed[[i]]$lambda] <- model[[i]]$lambda[fixed[[i]]$lambda]

    if(positive) {

      # Psi:
      cfa_param[[i]]$pj_psi <- cfa_trans[[i]]$pj_psi
      # Theta:
      cfa_param[[i]]$pj_theta <- cfa_trans[[i]]$pj_theta

    } else {

      # Psi:
      cfa_param[[i]]$psi <- cfa_trans[[i]]$psi
      # Force symmetry:
      # cfa_param[[i]]$psi[upper.tri(cfa_param[[i]]$psi)] <- cfa_param[[i]]$psi[lower.tri(cfa_param[[i]]$psi)]
      # Insert fixed values in the model:
      cfa_param[[i]]$psi[fixed[[i]]$psi] <- model[[i]]$psi[fixed[[i]]$psi]

      # Theta:
      cfa_param[[i]]$theta <- cfa_trans[[i]]$theta
      # Force symmetry:
      # cfa_param[[i]]$theta[upper.tri(cfa_param[[i]]$theta)] <- cfa_param[[i]]$theta[lower.tri(cfa_param[[i]]$theta)]
      # Insert fixed values in the model:
      cfa_param[[i]]$theta[fixed[[i]]$theta] <- model[[i]]$theta[fixed[[i]]$theta]

    }

    # Create the target matrices for positive-definite constraints:
    if(positive) {

      target_psi[[i]] <- matrix(0, nrow = nfactors, ncol = nfactors)
      target_psi[[i]][nonfixed[[i]]$psi] <- 1
      target_theta[[i]] <- matrix(0, nrow = nitems, ncol = nitems)
      target_theta[[i]][nonfixed[[i]]$theta] <- 1

      q <- nfactors
      p <- nitems
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

  init_trans <- vector("list", length = ngroups)

  for(i in 1:ngroups) {

    init_trans[[i]]$lambda <- rorth(nitems, nfactors)
    init_trans[[i]]$lambda[fixed[[i]]$lambda] <- fixed_values[[i]]$lambda

    if(positive) {

      init_trans[[i]]$pj_psi <- rpoblq(nfactors, nfactors, target = target_psi[[i]])
      init_trans[[i]]$pj_theta <- rpoblq(nitems, nitems, target = target_theta[[i]])

      init_trans[[i]]$psi <- crossprod(init_trans[[i]]$pj_psi)

      init_trans[[i]]$theta <- crossprod(init_trans[[i]]$pj_theta)

    } else {

      P <- diag(nfactors)
      init_trans[[i]]$psi <- P
      init_trans[[i]]$psi[fixed[[i]]$psi] <- fixed_values[[i]]$psi

      U <- diag(1/diag(solve(correl[[i]]$R)))
      init_trans[[i]]$theta <- U
      init_trans[[i]]$theta[fixed[[i]]$theta] <- fixed_values[[i]]$theta

    }

  }

  # Create the vectors of parameters and transformed parameters:
  parameters <- transparameters <- vector("list", length = control$rstarts)
  # Indices of the unique transparameters in init_trans:
  trans_inds <- match(transparameters_labels, vector_trans)
  init_inds <- match(parameters_labels, vector_trans)
  for(i in 1:control$rstarts) {

    transparameters[[i]] <- unlist(init_trans)[trans_inds]
    names(transparameters[[i]]) <- transparameters_labels
    parameters[[i]] <- unlist(init_trans)[init_inds]
    names(parameters[[i]]) <- parameters_labels

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

  #### Manifolds ####

  control_manifold <- list()

  if(!positive) {

    all <- unique(unlist(cfa_param))
    indices_all <- match(all, parameters_labels)
    indices_all <- indices_all[!is.na(indices_all)]
    labels <- parameters_labels[indices_all]
    indices <- list(indices_all-1L)
    control_manifold[[1]] <- list(manifold = "euclidean",
                                  parameters = labels,
                                  indices = indices)

  } else {

    lambdas <- unlist(lapply(cfa_param, FUN = \(x) x$lambda))
    indices_lambda <- match(unique(lambdas), parameters_labels)
    indices_lambda <- indices_lambda[!is.na(indices_lambda)]
    labels <- parameters_labels[indices_lambda]
    indices <- list(indices_lambda-1L)
    control_manifold[[1]] <- list(manifold = "euclidean",
                                  parameters = labels,
                                  indices = indices)
    k <- 2L

    for(i in 1:ngroups) {

      indices_pj_psi <- match(unique(c(cfa_param[[i]]$pj_psi)),
                              parameters_labels)
      labels <- parameters_labels[indices_pj_psi]
      indices <- list(indices_pj_psi-1L)
      control_manifold[[k]] <- list(manifold = "poblq",
                                    parameters = labels,
                                    indices = indices,
                                    target = target_psi[[i]])
      k <- k+1L

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

  }

  #### Transformations ####

  control_transform <- list()
  k <- 1L

  for(i in 1:ngroups) {

    if(positive) {

      lower_psi <- lower.tri(cfa_trans[[i]]$psi, diag = TRUE)
      lower_theta <- lower.tri(cfa_trans[[i]]$theta, diag = TRUE)

      # The result of the "crossprod" transformation is always symmetric,
      # so we take only the diagonal and lower diagonal elements

      labels_in <- c(cfa_trans[[i]]$pj_psi)
      labels_out <- c(cfa_trans[[i]]$psi[lower_psi])
      indices_in <- list(match(labels_in, transparameters_labels)-1L)
      indices_out <- list(match(labels_out, transparameters_labels)-1L)
      control_transform[[k]] <- list(transform = "crossprod",
                                     labels_in = labels_in,
                                     indices_in = indices_in,
                                     labels_out = labels_out,
                                     indices_out = indices_out,
                                     p = nrow(cfa_trans[[i]]$psi),
                                     q = ncol(cfa_trans[[i]]$psi))
      k <- k+1L

      labels_in <- c(cfa_trans[[i]]$pj_theta)
      labels_out <- c(cfa_trans[[i]]$theta[lower_theta])
      indices_in <- list(match(labels_in, transparameters_labels)-1L)
      indices_out <- list(match(labels_out, transparameters_labels)-1L)
      control_transform[[k]] <- list(transform = "crossprod",
                                     labels_in = labels_in,
                                     indices_in = indices_in,
                                     labels_out = labels_out,
                                     indices_out = indices_out,
                                     p = nrow(cfa_trans[[i]]$theta),
                                     q = ncol(cfa_trans[[i]]$theta))
      k <- k+1L

    }

  }

  #### Estimators ####

  control_estimator <- list()
  k <- 1L

  for(i in 1:ngroups) {

    lower_psi <- lower.tri(cfa_trans[[i]]$psi, diag = TRUE)
    lower_theta <- lower.tri(cfa_trans[[i]]$theta, diag = TRUE)

    indices_lambda <- match(cfa_trans[[i]]$lambda,
                            transparameters_labels)
    indices_psi <- match(cfa_trans[[i]]$psi[lower_psi],
                         transparameters_labels)
    indices_theta <- match(cfa_trans[[i]]$theta[lower_theta],
                           transparameters_labels)
    indices_all <- c(indices_lambda, indices_psi, indices_theta)

    indices <- list(indices = indices_all-1L,
                    indices_lambda = indices_lambda-1L,
                    indices_psi = indices_psi-1L,
                    indices_theta = indices_theta-1L)
    labels <- transparameters_labels[indices_all]

    estimator <- tolower(estimator)
    if(estimator == "uls" || estimator == "dwls") {
      cfa_estimator <- "cfa_dwls"
    } else if(estimator == "ml") {
      cfa_estimator <- "cfa_ml"
    }

    control_estimator[[k]] <- list(estimator = cfa_estimator,
                                   labels = labels,
                                   indices = indices,
                                   R = correl[[i]]$R,
                                   W = correl[[i]]$W,
                                   nfactors = nrow(cfa_trans[[i]]$psi))
    k <- k+1L

  }

  #### Return ####

  result <- list(control_manifold = control_manifold,
                 control_transform = control_transform,
                 control_estimator = control_estimator)

  return(result)

}
