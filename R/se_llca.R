# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 13/10/2025
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
se.llca <- function(fit, type = "standard", model = "model", digits = 3) {

  # Select the parameters to display according to model type:

  if(model == "user") {

    model <- fit@modelInfo$prob_model
    est <- fit@user_model
    fit@modelInfo$control$minimal_se <- FALSE
    fit@modelInfo$control$se_names <- fit@modelInfo$transparameters_labels

  } else if(model == "model") {

    model <- fit@modelInfo$model
    est <- fit@parameters
    fit@modelInfo$control$minimal_se <- TRUE
    fit@modelInfo$control$se_names <- fit@modelInfo$parameters_labels

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
  selection <- match(mylabels, fit@modelInfo$control$se_names)
  se <- SE$se[selection]
  se[is.na(se)] <- 0
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
  result$H <- SE$H
  result$HBinvH <- SE$HBinvH

  return(result)

}

#' @method se llcalist
#' @export
se.llcalist <- function(model, digits = 3) {

  nmodels <- length(model)
  out <- vector("list", length = nmodels)
  for(i in 1:nmodels) {

    out[[i]] <- se.llca(model[[i]], type = type, model = model,
                        digits = digits)
    names(out)[i] <- paste("nclasses = ", model[[i]]@modelInfo$nclasses,
                           sep = "")

  }

  class(out) <- "se.llcalist"

  return(out)

}

standard_se <- function(fit) {

  # Compute the variance-covariance matrix of the parameters:

  control_manifold <- fit@modelInfo$control_manifold
  control_transform <- fit@modelInfo$control_transform
  control_estimator <- fit@modelInfo$control_estimator
  control_optimizer <- fit@modelInfo$control
  control_optimizer$parameters[[1]] <- fit@Optim$parameters
  control_optimizer$transparameters[[1]] <- fit@Optim$transparameters

  H <- get_hess(control_manifold = control_manifold,
                control_transform = control_transform,
                control_estimator = control_estimator,
                control_optimizer = control_optimizer)$h
  # For numerical stability in optimization, the loglik is divided by N
  # So multiply the Hessian by N:
  H <- H*fit@data_list$nobs[[1]]

  # Get the variance-covariance matrix between the parameters:
  result <- get_vcov(control_manifold = control_manifold,
                     control_transform = control_transform,
                     control_estimator = control_estimator,
                     control_optimizer = control_optimizer,
                     H = H)

  # Name the parameters:
  colnames(H) <- rownames(H) <- fit@modelInfo$parameters_labels
  result$H <- H

  result$se <- as.vector(result$se)
  names(result$se) <- colnames(result$vcov) <- rownames(result$vcov) <-
    fit@modelInfo$control$se_names

  # Return:
  return(result)

}

robust_se <- function(fit) {

  # Rearrange control_estimator to get the gradient per response pattern:

  control_manifold <- fit@modelInfo$control_manifold
  control_transform <- fit@modelInfo$control_transform
  control_estimator <- fit@modelInfo$control_estimator
  control_optimizer <- fit@modelInfo$control
  control_optimizer$parameters[[1]] <- fit@Optim$parameters
  control_optimizer$transparameters[[1]] <- fit@Optim$transparameters

  H <- get_hess(control_manifold = control_manifold,
                control_transform = control_transform,
                control_estimator = control_estimator,
                control_optimizer = control_optimizer)$h
  # For numerical stability in optimization, the loglik is divided by N
  # So multiply the Hessian by N:
  H <- H*fit@data_list$nobs

  #### Collect the gradient by response pattern ####

  transparameters_labels <- fit@modelInfo$transparameters_labels
  lca_all <- fit@modelInfo$lca_all
  full_loglik <- lca_all$loglik
  weights <- fit@data_list$weights
  npatterns <- fit@data_list$npatterns
  nitems <- fit@data_list$nitems
  nparam <- fit@modelInfo$nparam
  nobs <- fit@data_list$nobs

  control_estimator <- control_estimator[-1]
  nclasses <- ncol(fit@modelInfo$lca_all$class)
  K <- length(control_estimator)

  for(s in 1:npatterns) {

    indices <- list(match(lca_all$class[s, ], transparameters_labels)-1L,
                    match(full_loglik[s,,], transparameters_labels)-1L)

    control_estimator[[K+s]] <- list(estimator = "lca",
                                     # labels = labels,
                                     indices = indices,
                                     S = 1L,
                                     J = nitems,
                                     I = nclasses,
                                     weights = weights[s])

  }

  # nest <- K + npatterns
  # weights <- c(rep(1, times = K), weights)
  nest <- npatterns
  f <- vector(length = nest)
  g <- matrix(NA, nrow = nest, ncol = nparam)
  B <- matrix(0, nrow = nparam, ncol = nparam)
  for(s in 1:nest) {

    computations <- get_grad(control_manifold = control_manifold,
                             control_transform = control_transform,
                             control_estimator = control_estimator[-(1:K)][s],
                             control_optimizer = control_optimizer)
    f[s] <- computations$f
    g[s, ] <- computations$g
    B <- B + weights[s]*(g[s, ]/weights[s]) %*% t(g[s, ]/weights[s])

  }

  # Update the hessian:
  B <- B*nobs/(nobs-1)
  newH <- H %*% solve(B) %*% H

  # Recompute the variance-covariance matrix with this new hessian:
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
    fit@modelInfo$control$se_names

  colnames(newH) <- rownames(newH) <-
    colnames(result$B) <- rownames(result$B) <-
    fit@modelInfo$parameters_labels

  result$H <- H
  result$HBinvH <- newH

  return(result)

}

ci <- function(fit, type = "standard", model = "model",
               confidence = 0.95, digits = 3) {

  # Select the parameters to display according to model type:

  if(model == "user") {

    model <- fit@modelInfo$prob_model
    est <- fit@user_model
    fit@modelInfo$control$minimal_se <- FALSE
    fit@modelInfo$control$se_names <- fit@modelInfo$transparameters_labels

    x <- c(fit@Optim$transparameters)
    # Number of total parameters and transformed parameters:
    ntrans <- length(fit@Optim$transparameters)
    # Initialize a vector of degrees of freedom:
    ps <- rep(1, times = ntrans)
    # slot for extracting the degrees of freedom from the transformations:
    slot <- 2
    # Get the indices of each transformed parameter:
    indices <- unlist(lapply(fit@modelInfo$control_transform,
                             FUN = \(x) x$indices_out[[1]]+1L))
    # Update the degrees of freedom of each transformed parameter:
    ps[indices] <- unlist(lapply(fit@Optim$outputs$transformations$vectors,
                                 FUN = \(x) x[[slot]]))

  } else if(model == "model") {

    model <- fit@modelInfo$model
    est <- fit@parameters
    fit@modelInfo$control$minimal_se <- TRUE
    fit@modelInfo$control$se_names <- fit@modelInfo$parameters_labels
    x <- c(fit@Optim$parameters)
    ps <- rep(1L, fit@modelInfo$nparam)

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
  selection <- match(mylabels, fit@modelInfo$control$se_names)
  se <- SE$se[selection]
  se[is.na(se)] <- 0
  names(se) <- mylabels

  x <- x[selection]
  x[is.na(x)] <- 0
  ps <- ps[selection]
  ps[is.na(ps)] <- 0

  # Compute confidence intervals for each transformed parameter:
  lower <- x - sqrt(qchisq(confidence, df = ps)) * se
  upper <- x + sqrt(qchisq(confidence, df = ps)) * se
  names(lower) <- names(upper) <- mylabels

  # Get confidence limits for the user model or raw model parameters:

  # Tables:
  lower_ci <- fill_list_with_vector(model, lower)
  lower_ci <- allnumeric(lower_ci)
  upper_ci <- fill_list_with_vector(model, upper)
  upper_ci <- allnumeric(upper_ci)

  table <- combine_est_ci(lower_ci, est, upper_ci, digits = digits)

  # Return:
  result <- list()
  result$table <- table
  result$lower_table <- lower_ci
  result$upper_table <- upper_ci
  result$lower <- lower
  result$upper <- upper
  result$se <- c(SE$se)
  result$vcov <- SE$vcov
  result$B <- SE$B

  return(result)

}
