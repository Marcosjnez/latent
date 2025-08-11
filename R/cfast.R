# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 19/07/2025
#'
#' @title
#' Fit a Confirmatory Factor Analysis (CFA) model with lavaan syntax.
#'
#' @usage
#'
#' cfast(data, model = NULL, cor = "pearson", estimator = "ml",
#' group = NULL, sample.cov = NULL, nobs = NULL,
#' missing = "pairwise.complete.obs", W = NULL, std.lv = FALSE,
#' positive = FALSE, do.fit = TRUE, control = NULL)
#'
#' @examples
#'
#' \dontrun{
#' # The famous Holzinger and Swineford (1939) example
#' HS.model <- ' visual  =~ x1 + x2 + x3
#'               textual =~ x4 + x5 + x6
#'               speed   =~ x7 + x8 + x9 '
#'
#' fit <- cfast(model = HS.model, data = HolzingerSwineford1939)
#' summary(fit, fit.measures = TRUE)
#'}
#'
#' @export
cfast <- function(data, model = NULL, cor = "pearson", estimator = "uls",
                  group = NULL, sample.cov = NULL, nobs = NULL,
                  missing = "pairwise.complete.obs", W = NULL, std.lv = FALSE,
                  positive = FALSE, do.fit = TRUE, control = NULL) {

  # Extract the lavaan model:
  extract_fit <- lavaan::cfa(model = model, data = data,
                             sample.cov = sample.cov, std.lv = std.lv,
                             do.fit = FALSE, group = group)
  item_names <- unique(extract_fit@ParTable$rhs[extract_fit@ParTable$op == "=~"])
  # Model for the parameters:
  matrices_param <- getmodel_cfa(extract_fit)

  if(!is.list(matrices_param[[1]])) {
    matrices_param <- list(matrices_param)
  }
  ngroups <- length(matrices_param)

  # Collect the parameter labels:
  vector_matrices_param <- unname(unlist(matrices_param))
  all_paramlabels <- is.na(suppressWarnings(as.numeric(vector_matrices_param)))
  parameter_labels <- unique(vector_matrices_param[all_paramlabels])
  # Prepare the vector of fixed parameters and their indices:
  fixed_indices <- which(!all_paramlabels)
  fixed_values <- as.numeric(vector_matrices_param[fixed_indices])

  if(positive) { # Redefine the model to find a positive-semidefinite solution

    targetpsi <- vector("list", length = ngroups)
    targettheta <- vector("list", length = ngroups)

    for(i in 1:ngroups) {

      p <- nrow(matrices_param[[i]]$lambda)
      q <- ncol(matrices_param[[i]]$lambda)

      param_psi <- suppressWarnings(as.numeric(matrices_param[[i]]$psi))
      param_theta <- suppressWarnings(as.numeric(matrices_param[[i]]$theta))

      # Generate the indicator matrix of cells that are freely estimated. These
      # matrices are passed to the pobl manifold:
      targetpsi[[i]] <- matrix(as.numeric(is.na(param_psi)), nrow = q, ncol = q)
      targettheta[[i]] <- matrix(as.numeric(is.na(param_theta)), nrow = p, ncol = p)

      # The model parameters must be modified if psi and theta are constrained
      # to be positive-semidefinite:
      psi_labels <- matrix(paste("g", i, ".pj.psi", 1:(q*q), sep = ""), nrow = q, ncol = q)
      matrices_param[[i]]$psi <- psi_labels
      theta_labels <- matrix(paste("g", i, ".pj.theta", 1:(p*p), sep = ""), nrow = p, ncol = p)
      matrices_param[[i]]$theta <- theta_labels

    }

    # Update the parameter labels:
    vector_matrices_param <- unname(unlist(matrices_param))
    all_paramlabels <- is.na(suppressWarnings(as.numeric(vector_matrices_param)))
    parameter_labels <- unique(vector_matrices_param[all_paramlabels])

  }

  # Get the correlation and weight matrices:
  correl <- vector("list", length = ngroups)
  if(ngroups > 1) {
    data_split <- split(data, data[[group]])
    for(i in 1:ngroups) {
      correl[[i]] <- correlation(data_split[[i]], item_names, cor,
                                 estimator, missing)
    }
  } else {
    correl[[1]] <- correlation(data, item_names, cor,
                               estimator, missing)
  }

  # Build the model for the transformed parameters:
  matrices_trans <- matrices_param
  init <- vector("list", length = ngroups)

  for(i in 1:ngroups) {

    p <- nrow(matrices_trans[[i]]$lambda)
    q <- ncol(matrices_trans[[i]]$lambda)
    # lambda:
    lambda_labels <- matrix(paste("g", i, ".lambda", 1:(p*q), sep = ""), nrow = p, ncol = q)
    matrices_trans[[i]]$lambda <- lambda_labels

    # psi:
    psi_labels <- matrix(paste("g", i, ".psi", 1:(q*q), sep = ""), nrow = q, ncol = q)
    matrices_trans[[i]]$psi <- psi_labels

    # theta:
    theta_labels <- matrix(paste("g", i, ".theta", 1:(p*p), sep = ""), nrow = p, ncol = p)
    matrices_trans[[i]]$theta <- theta_labels

    if(positive) { # Force a positive-semidefinite solution?

      init[[i]]$lambda <- rorth(p, q)
      init[[i]]$psi <- rpoblq(q, q, targetpsi[[i]])
      init[[i]]$theta <- rpoblq(p, p, targettheta[[i]])

    } else {

      u <- 1/diag(solve(correl[[i]]$R))
      init[[i]]$lambda <- rorth(p, q)
      init[[i]]$psi <- diag(q)
      init[[i]]$theta <- diag(u)

    }

  }

  # Collect the transformed parameters labels (including the parameter labels):
  vector_matrices_trans <- unname(unlist(matrices_trans))
  transparameter_labels <- unique(c(parameter_labels, vector_matrices_trans))

  # Initialize the vector of parameters
  nparameters <- length(parameter_labels)
  ntransparameters <- length(parameter_labels)

  # Relate the transformed parameters to the parameters:
  param2trans <- match(transparameter_labels, parameter_labels)
  param2trans <- param2trans[!is.na(param2trans)]

  # Relate the parameters to the transformed parameters:
  trans2param <- match(parameter_labels, transparameter_labels)

  # Put the initial parameter values in their place:
  transparameters <- unlist(init)
  indices <- match(parameter_labels, vector_matrices_param)
  parameters <- transparameters[indices]
  names(parameters) <- parameter_labels
  # Fix the transformed parameters according to the user:
  transparameters[fixed_indices] <- fixed_values
  transparameters <- c(parameters, transparameters)

  # Check the arguments to control_optimizer and create defaults:
  control <- cfast_control(control)

  # Create the list of manifolds:
  control_manifold <- list()
  k <- 1L
  if(positive) {

    for(i in 1:ngroups) {

      indices_lambda <- match(unique(c(matrices_param[[i]]$lambda)),
                              parameter_labels)
      indices_lambda <- indices_lambda[!is.na(indices_lambda)]
      labels <- parameter_labels[indices_lambda]
      control_manifold[[k]] <- list(manifold = "euclidean",
                                    parameters = labels,
                                    indices = list(indices_lambda-1L))
      k <- k+1L

      indices_psi <- match(unique(c(matrices_param[[i]]$psi)),
                           parameter_labels)
      indices_psi <- indices_psi[!is.na(indices_psi)]
      labels <- parameter_labels[indices_psi]
      control_manifold[[k]] <- list(manifold = "poblq",
                                    parameters = labels,
                                    indices = list(indices_psi-1L),
                                    target = targetpsi[[i]])
      k <- k+1L

      indices_theta <- match(unique(c(matrices_param[[i]]$theta)),
                             parameter_labels)
      indices_theta <- indices_theta[!is.na(indices_theta)]
      labels <- parameter_labels[indices_theta]
      control_manifold[[k]] <- list(manifold = "poblq",
                                    parameters = labels,
                                    indices = list(indices_theta-1L),
                                    target = targettheta[[i]])
      k <- k+1L

    }

  } else {

    indices <- 1:nparameters
    labels <- parameter_labels[indices]
    control_manifold[[k]] <- list(manifold = "euclidean",
                                  parameters = labels,
                                  indices = list(indices-1L))
    k <- k+1L

  }

  # Create the list of transformations:
  control_transform <- list()
  k <- 1L
  for(i in 1:ngroups) {
    if(positive) {

      positions <- which(matrices_param[[i]]$lambda %in% parameter_labels)
      labels_in <- matrices_param[[i]]$lambda[positions]
      labels_out <- matrices_trans[[i]]$lambda[positions]
      indices_in <- match(labels_in, transparameter_labels)
      indices_out <- match(labels_out, transparameter_labels)
      control_transform[[k]] <- list(transform = "identity",
                                     labels_in = labels_in,
                                     indices_in = list(indices_in-1L),
                                     labels_out = labels_out,
                                     indices_out = list(indices_out-1L))
      k <- k+1L

      labels_in <- c(matrices_param[[i]]$psi)
      labels_out <- c(matrices_trans[[i]]$psi)
      indices_in <- match(labels_in, transparameter_labels)
      indices_out <- match(labels_out, transparameter_labels)
      control_transform[[k]] <- list(transform = "crossprod",
                                     labels_in = labels_in,
                                     indices_in = list(indices_in-1L),
                                     labels_out = labels_out,
                                     indices_out = list(indices_out-1L),
                                     p = nrow(matrices_trans[[i]]$psi),
                                     q = ncol(matrices_trans[[i]]$psi))
      k <- k+1L

      labels_in <- c(matrices_param[[i]]$theta)
      labels_out <- c(matrices_trans[[i]]$theta)
      indices_in <- match(labels_in, transparameter_labels)
      indices_out <- match(labels_out, transparameter_labels)
      control_transform[[k]] <- list(transform = "crossprod",
                                     labels_in = labels_in,
                                     indices_in = list(indices_in-1L),
                                     labels_out = labels_out,
                                     indices_out = list(indices_out-1L),
                                     p = nrow(matrices_trans[[i]]$theta),
                                     q = ncol(matrices_trans[[i]]$theta))
      k <- k+1L

    } else {

      all_param <- unname(unlist(matrices_param[[i]]))
      all_trans <- unname(unlist(matrices_trans[[i]]))

      positions <- which(all_param %in% parameter_labels)
      labels_in <- all_param[positions]
      labels_out <- all_trans[positions]
      indices_in <- match(labels_in, transparameter_labels)
      indices_out <- match(labels_out, transparameter_labels)
      control_transform[[k]] <- list(transform = "identity",
                                     labels_in = labels_in,
                                     indices_in = list(indices_in-1L),
                                     labels_out = labels_out,
                                     indices_out = list(indices_out-1L))
      k <- k+1L

    }
  }

  control_estimator <- list()
  k <- 1L
  for(i in 1:ngroups) {

    all_param <- unname(unlist(matrices_param[[i]]))
    all_trans <- unname(unlist(matrices_trans[[i]]))

    indices_all <- match(all_trans, transparameter_labels)

    indices_lambda <- match(matrices_trans[[i]]$lambda,
                            transparameter_labels[indices_all])

    indices_psi <- match(matrices_trans[[i]]$psi,
                         transparameter_labels[indices_all])

    indices_theta <- match(matrices_trans[[i]]$theta,
                           transparameter_labels[indices_all])

    indices <- list(indices = indices_all-1L,
                    indices_lambda = indices_lambda-1L,
                    indices_psi = indices_psi-1L,
                    indices_theta = indices_theta-1L)
    labels <- transparameter_labels[indices_all]

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
                                   nfactors = nrow(matrices_trans[[i]]$psi))
    k <- k+1L

  }

  control$parameters <- list(parameters)
  control$transparameters <- list(transparameters)
  control$param2transparam <- param2trans-1L
  control$transparam2param <- trans2param-1L

  if(do.fit) {

    x <- optimizer(control_manifold = control_manifold,
                   control_transform = control_transform,
                   control_estimator = control_estimator,
                   control_optimizer = control)

  } else {

    x <- vector("list")
    x$control_manifold <- control_manifold
    x$control_transform <- control_transform
    x$control_estimator <- control_estimator
    x$control <- control
    x$parameters <- parameters
    x$transparameters <- transparameters
    x$model <- matrices_param
    x$model_full <- matrices_trans
    x$parameter_labels <- parameter_labels
    x$transparameter_labels <- transparameter_labels

    return(x)

  }

  x$control_manifold <- control_manifold
  x$control_transform <- control_transform
  x$control_estimator <- control_estimator
  x$control <- control

  outputs <- vector("list", length = ngroups)

  for(i in 1:ngroups) {

    p <- nrow(x$control_estimator[[i]]$R)
    q <- control_estimator[[i]]$nfactors

    # Arrange lambda parameter estimates:
    outputs[[i]]$lambda <- matrix(x$outputs$estimators$matrices[[i]][[1]], p, q)

    # Arrange psi parameter estimates:
    outputs[[i]]$psi <- matrix(x$outputs$estimators$matrices[[i]][[2]], q, q)

    # Arrange theta parameter estimates:
    outputs[[i]]$theta <- matrix(x$outputs$estimators$matrices[[i]][[3]], p, p)
    # uniquenesses_hat[[i]] <- diag(psi_hat[[i]])

    # Model matrix:
    outputs[[i]]$model <- matrix(x$outputs$estimators$matrices[[i]][[4]], p, p)

    # Residual matrix:
    outputs[[i]]$residuals <- matrix(x$outputs$estimators$matrices[[i]][[5]], p, p)

    # Weight matrix:
    if(x$control_estimator[[i]]$estimator == "uls" ||
       x$control_estimator[[i]]$estimator == "dwls") {
      outputs[[i]]$W <- matrix(x$outputs$estimators$matrices[[i]][[6]], p, p)
    }

    # Uniquenesses:
    outputs[[i]]$uniquenesses <- c(x$outputs$estimators$vectors[[i]][[1]])

  }

  x$model <- matrices_param
  x$model_full <- matrices_trans
  x$outputs <- outputs

  return(x)

}

