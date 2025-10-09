# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 09/10/2025
#'
#' @title
#' Standard Errors
#' @description
#'
#' Compute standard errors.
#'
#' @usage
#'
#' se(fit)
#'
#' @param fit model fitted with lca.
#' @param confidence Coverage of the confidence interval.
#'
#' @details Compute standard errors.
#'
#' @return List with the following objects:
#' \item{vcov}{Variance-covariance matrix between the parameters.}
#' \item{se}{Standard errors.}
#' \item{SE}{Standard errors in the model list.}
#'
#' @references
#'
#' None yet.
#'
#' @method se llca
#' @export
se.llca <- function(fit, type = "standard", model = "user", digits = 2) {

  # Select the parameters to display according to model type:

  if(model == "user") {

    model <- fit@modelInfo$prob_model
    est <- fit@user_model

  } else if(model == "model") {

    model <- fit@modelInfo$model
    est <- fit@parameters

  } else {
    stop("Unknown model")
  }

  if(type == "standard") {

    SE <- standard_se(fit = fit)
    SE$B <- matrix(, nrow = 0, ncol = 0) # Empty matrix

  } else if(type == "robust") {

    # Use H to compute vcov = solve(H) %*% B %*% solve(H):
    SE <- robust_se(fit = fit)

  } else {
    stop("Unknown type")
  }

  # Select the parameter labels for the table:
  mylabels <- unlist(model)
  selection <- match(mylabels, fit@modelInfo$transparameters_labels)

  # Standard errors:
  se <- SE$se[selection]
  names(se) <- mylabels

  # Tables:
  table_se <- fill_list_with_vector(model, se)
  table_se <- allnumeric(table_se)
  table <- combine_est_se(est, table_se, digits = digits)

  # Return:
  result <- list()
  result$table <- table
  result$table_se <- table_se
  result$se <- c(SE$se)
  result$vcov <- SE$vcov
  result$B <- SE$B

  return(result)

}

standard_se <- function(fit) {

  # Compute the variance-covariance matrix of the parameters:

  control_manifold <- fit@modelInfo$control_manifold
  control_transform <- fit@modelInfo$control_transform
  control_estimator <- fit@modelInfo$control_estimator
  control_optimizer <- fit@modelInfo$control
  control_optimizer$parameters[[1]] <- fit@Optim$opt$parameters
  control_optimizer$transparameters[[1]] <- fit@Optim$opt$transparameters

  # Calculate the Hessian matrix using numerical approximations:
  G <- function(parameters) {

    control_optimizer$parameters[[1]] <- parameters
    g <- get_grad(control_manifold = control_manifold,
                  control_transform = control_transform,
                  control_estimator = control_estimator,
                  control_optimizer = control_optimizer)$g

    return(g)

  }

  H <- numDeriv::jacobian(func = G, x = fit@Optim$opt$parameters)
  H <- 0.5*(H + t(H)) # Force symmetry

  # Get the variance-covariance matrix between the parameters:
  result <- get_vcov(control_manifold = control_manifold,
                     control_transform = control_transform,
                     control_estimator = control_estimator,
                     control_optimizer = control_optimizer,
                     H = H)

  # Name the parameters:
  result$se <- as.vector(result$se)
  names(result$se) <- colnames(result$vcov) <- rownames(result$vcov) <-
    fit@modelInfo$transparameters_labels

  colnames(H) <- rownames(H) <- fit@modelInfo$parameters_labels
  result$h <- H

  return(result)

}

robust_se <- function(fit) {

  # Rearrange control_estimator to get the gradient per response pattern:

  control_manifold <- fit@modelInfo$control_manifold
  control_transform <- fit@modelInfo$control_transform
  control_estimator <- fit@modelInfo$control_estimator
  control_optimizer <- fit@modelInfo$control
  control_optimizer$parameters[[1]] <- fit@Optim$opt$parameters
  control_optimizer$transparameters[[1]] <- fit@Optim$opt$transparameters

  # Calculate the Hessian matrix using numerical approximations:
  G <- function(parameters) {

    control_optimizer$parameters[[1]] <- parameters
    g <- get_grad(control_manifold = control_manifold,
                  control_transform = control_transform,
                  control_estimator = control_estimator,
                  control_optimizer = control_optimizer)$g

    return(g)

  }

  H <- numDeriv::jacobian(func = G, x = fit@Optim$opt$parameters)
  H <- 0.5*(H + t(H)) # Force symmetry

  transparameters_labels <- fit@modelInfo$transparameters_labels
  lca_trans <- fit@modelInfo$lca_trans
  full_loglik <- lca_trans$loglik
  weights <- fit@Optim$data_list$weights
  npatterns <- fit@Optim$data_list$npatterns
  nitems <- fit@Optim$data_list$nitems
  nparam <- fit@modelInfo$nparam
  nobs <- fit@Optim$data_list$nobs

  control_estimator <- control_estimator[-1]
  nclasses <- ncol(fit@modelInfo$lca_trans$class)
  K <- length(control_estimator)

  for(s in 1:npatterns) {

    labels <- c(lca_trans$class[s, ], full_loglik[s,,])
    all_indices <- match(labels, transparameters_labels)
    indices_classes <- match(lca_trans$class[s, ], transparameters_labels)
    indices_items <- match(full_loglik[s,,], transparameters_labels)
    indices_theta <- match(lca_trans$theta[s, ], transparameters_labels)
    indices <- list(all_indices-1L, indices_classes-1L, indices_items-1L,
                    indices_theta-1L)
    SJ <- 1L*nitems
    hess_indices <- lapply(0:(nclasses - 1), function(i) {
      nclasses + seq(1 + i * SJ, (i + 1) * SJ)-1L })

    control_estimator[[K+s]] <- list(estimator = "lca",
                                     labels = labels,
                                     indices = indices,
                                     S = 1L,
                                     J = nitems,
                                     I = nclasses,
                                     weights = weights[s],
                                     hess_indices = hess_indices)

  }

  # nest <- K + npatterns
  # weights <- c(rep(1, times = K), weights)
  nest <- npatterns
  f <- vector(length = nest)
  g <- matrix(NA, nrow = nest, ncol = nparam)
  B <- matrix(0, nrow = nparam, ncol = nparam)
  for(s in 1:nest) {

    computations <- grad_comp(control_manifold = control_manifold,
                              control_transform = control_transform,
                              control_estimator = control_estimator[-(1:K)][s],
                              # control_estimator = control_estimator[s],
                              control_optimizer = control_optimizer,
                              compute = "g")
    f[s] <- computations$f
    g[s, ] <- computations$g
    B <- B + weights[s]*(g[s, ]/weights[s]) %*% t(g[s, ]/weights[s])

  }
  # sum(f)
  # fit@penalized_loglik
  # colSums(g)
  # c(fit@Optim$opt$rg)
  B <- B*nobs/(nobs-1)
  # eigen(B)$values

  # Update the hessian:
  newH <- H %*% solve(B) %*% H

  # Recompute the variance-covariance matrix with this new hessian:
  control_manifold <- fit@modelInfo$control_manifold
  control_transform <- fit@modelInfo$control_transform
  control_estimator <- fit@modelInfo$control_estimator
  control_optimizer <- fit@modelInfo$control
  control_optimizer$parameters[[1]] <- fit@Optim$opt$parameters
  control_optimizer$transparameters[[1]] <- fit@Optim$opt$transparameters

  # Get the variance-covariance matrix between the parameters:
  result <- get_vcov(control_manifold = control_manifold,
                     control_transform = control_transform,
                     control_estimator = control_estimator,
                     control_optimizer = control_optimizer,
                     H = newH)

  # Add B to the results:
  result$B <- B

  # Name the parameters:
  names(result$se) <- colnames(result$vcov) <- rownames(result$vcov) <-
    fit@modelInfo$transparameters_labels

  colnames(newH) <- rownames(newH) <-
    colnames(result$B) <- rownames(result$B) <-
    fit@modelInfo$transparameters_labels <- fit@modelInfo$parameters_labels

  result$H <- H
  result$HBinvH <- newH

  return(result)

}

ci <- function(fit, type = "standard", model = "user",
               confidence = 0.95, digits = 2) {

  if(type == "standard") {

    SE <- standard_se(fit = fit)
    SE$B <- matrix(, nrow = 0, ncol = 0) # Empty matrix

  } else if(type == "robust") {

    # Use H to compute vcov = solve(H) %*% B %*% solve(H):
    SE <- robust_se(fit = fit)

  } else {
    stop("Unknown type")
  }

  se <- c(SE$se) # standard errors
  names(se) <- fit@modelInfo$transparameters_labels

  # Number of total parameters and transformed parameters:
  ntrans <- length(fit@Optim$opt$transparameters)
  # Initialize a vector of degrees of freedom:
  ps <- rep(1, times = ntrans)
  # slot for extracting the degrees of freedom from the transformations:
  slot <- 2
  # Get the indices of each transformed parameter:
  indices <- unlist(lapply(fit@modelInfo$control_transform,
                           FUN = \(x) x$indices_out[[1]]+1L))
  # Update the degrees of freedom of each transformed parameter:
  ps[indices] <- unlist(lapply(fit@Optim$opt$outputs$transformations$vectors,
                               FUN = \(x) x[[slot]]))

  # Compute confidence intervals for each transformed parameter:
  x <- c(fit@Optim$opt$transparameters)
  lower <- x - sqrt(qchisq(confidence, df = ps)) * se
  upper <- x + sqrt(qchisq(confidence, df = ps)) * se
  names(lower) <- names(upper) <- fit@modelInfo$transparameters_labels

  # Get confidence limits for the user model or raw model parameters:

  # Select the parameters according to model type:

  if(model == "user") {

    model <- fit@modelInfo$prob_model
    est <- fit@user_model

  } else if(model == "model") {

    model <- fit@modelInfo$model
    est <- fit@parameters

  }

  # Select the parameter labels for the table:
  mylabels <- unlist(model)
  selection <- match(mylabels, fit@modelInfo$transparameters_labels)

  # Tables:
  lower_ci <- fill_list_with_vector(model, lower[selection])
  lower_ci <- allnumeric(lower_ci)
  upper_ci <- fill_list_with_vector(model, upper[selection])
  upper_ci <- allnumeric(upper_ci)

  table <- combine_est_ci(lower_ci, est, upper_ci, digits = digits)

  # Return:
  result <- list()
  result$table <- table
  result$lower_table <- lower_ci
  result$upper_table <- upper_ci
  result$lower <- lower
  result$upper <- upper
  result$se <- se

  return(result)

}
