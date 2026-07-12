# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 19/06/2026
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
se.llca <- function(fit, type = "standard", digits = 3) {

  # Select the parameters to display according to model type:

  fit@modelInfo$control_optimizer$minimal_se <- TRUE

  if(class(fit@modelInfo$control_optimizer$model) == "llca") {

    SE <- se_twostep(fit2 = fit, type = type)

  } else {

    if(type == "standard") {

      SE <- standard_se(fit = fit)
      SE$B <- matrix(, nrow = 0, ncol = 0) # Empty matrix

    } else if(type == "robust") {

      # Use H to compute vcov = solve(H) %*% B %*% solve(H):
      SE <- robust_se_LG(fit = fit)

    } else {
      stop("Unknown type")
    }

  }

  est <- fill_in(fit@modelInfo$trans[names(fit@modelInfo$param)],
                 c(fit@Optim$parameters, fit@Optim$transparameters),
                 miss = NA)
  table_se <- fill_in(fit@modelInfo$trans[names(fit@modelInfo$param)],
                      SE$se, miss = NA)

  table <- combine_est_se(est, table_se, digits = digits)

  # Return:
  result <- list()
  result$table <- table
  result$table_se <- table_se
  result$se <- c(SE$se)
  result$vcov <- SE$vcov
  result$B <- SE$B
  result$H <- SE$H
  result$newH <- SE$newH

  return(result)

}

#' @method se llcalist
#' @export
se.llcalist <- function(model, type = "standard", digits = 3) {

  nmodels <- length(model)
  out <- vector("list", length = nmodels)
  for(i in 1:nmodels) {

    out[[i]] <- se.llca(model[[i]], type = type, digits = digits)
    names(out)[i] <- paste("nclasses=",
                           ncol(model[[i]]@modelInfo$param$beta),
                           sep = "")

  }

  class(out) <- "se.llcalist"

  return(out)

}

robust_se_LG <- function(fit) {

  # Robust standard errors in LatentGold's style
  # Rearrange control_estimator to get the gradient per response pattern:

  control_manifold <- fit@modelInfo$control_manifold
  control_transform <- fit@modelInfo$control_transform
  control_estimator <- fit@modelInfo$control_estimator
  control_optimizer <- fit@modelInfo$control_optimizer
  control_optimizer$parameters[[1]] <- fit@Optim$parameters
  control_optimizer$transparameters[[1]] <- fit@Optim$transparameters

  H <- get_hess(control_manifold = control_manifold,
                control_transform = control_transform,
                control_estimator = control_estimator,
                control_optimizer = control_optimizer,
                cores = parallel::detectCores())$h

  #### Collect the gradient by response pattern ####

  transparameters_labels <- fit@modelInfo$transparameters_labels
  trans <- fit@modelInfo$trans
  full_loglik <- trans$loglik
  pattern_weights <- fit@dataList$pattern_weights
  npatterns <- fit@dataList$npatterns
  nitems <- fit@dataList$nitems
  nparam <- fit@modelInfo$nparam
  nobs <- fit@dataList$nobs

  control_estimator <- control_estimator[-1]
  nclasses <- ncol(fit@modelInfo$trans$class)
  K <- length(control_estimator)

  for(s in 1:npatterns) {

    indices <- list(match(trans$class[s, ], transparameters_labels)-1L,
                    match(full_loglik[s,,], transparameters_labels)-1L)

    control_estimator[[K+s]] <- list(estimator = "lca",
                                     # labels = labels,
                                     indices = indices,
                                     S = 1L,
                                     J = nitems,
                                     I = nclasses,
                                     pattern_weights = pattern_weights[s])

  }

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
    B <- B + pattern_weights[s]*(g[s, ]/pattern_weights[s]) %*% t(g[s, ]/pattern_weights[s])

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
    fit@modelInfo$control_optimizer$se_names

  colnames(newH) <- rownames(newH) <-
    colnames(result$B) <- rownames(result$B) <-
    fit@modelInfo$parameters_labels

  result$H <- H
  result$newH <- newH

  return(result)

}

ci <- function(fit, type = "standard", confidence = 0.95, digits = 3) {

  # Select the parameters to display according to model type:

  fit@modelInfo$control_optimizer$minimal_se <- TRUE

  if(class(fit@modelInfo$original_model) == "llca") {

    SE <- se_twostep(fit2 = fit, type = type)

  } else {

    if(type == "standard") {

      SE <- standard_se(fit = fit)
      SE$B <- matrix(, nrow = 0, ncol = 0) # Empty matrix

    } else if(type == "robust") {

      # Use H to compute vcov = solve(H) %*% B %*% solve(H):
      SE <- robust_se_LG(fit = fit)

    } else {
      stop("Unknown type")
    }

  }

  x <- c(fit@Optim$parameters)
  ps <- rep(1L, fit@modelInfo$nparam)

  # Compute confidence intervals for each transformed parameter:
  lower <- x - sqrt(qchisq(confidence, df = ps)) * SE$se
  upper <- x + sqrt(qchisq(confidence, df = ps)) * SE$se
  names(lower) <- names(upper) <- fit@modelInfo$parameters_labels

  # Get confidence limits for the user model or raw model parameters:
  est <- fill_in(fit@modelInfo$trans[names(fit@modelInfo$param)],
                 c(fit@Optim$parameters, fit@Optim$transparameters),
                 miss = NA)

  # Tables:
  lower_ci <- fill_in(fit@modelInfo$param, lower)
  upper_ci <- fill_in(fit@modelInfo$param, upper)

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
  result$H <- SE$H

  return(result)

}

se_twostep <- function(fit2, type = "standard") {

  fit1 <- fit2@modelInfo$control_optimizer$model

  # Variance covariance of step 1:
  VCOV1 <- se(fit1, type = type)
  # Variance covariance of step 2:
  fit2@modelInfo$control_optimizer$model <- NULL # Avoid infinite recursion
  VCOV2 <- se(fit2, type = type)

  args <- fit2@dataList$args
  args$model <- NULL
  args$do.fit <- FALSE
  args$adjustment <- "none"
  fit <- do.call(lca, args)

  # Get the hessian matrix:
  control_manifold <- fit@modelInfo$control_manifold
  control_transform <- fit@modelInfo$control_transform
  control_estimator <- fit@modelInfo$control_estimator
  control_optimizer <- fit@modelInfo$control_optimizer

  parameters <- fit2@Optim$transparameters[fit@modelInfo$parameters_labels]
  transparameters <- fit2@Optim$transparameters[fit@modelInfo$transparameters_labels]
  control_optimizer$parameters[[1]] <- parameters
  control_optimizer$transparameters[[1]] <- transparameters
  x <- get_hess(control_manifold, control_transform,
                control_estimator, control_optimizer)
  colnames(x$h) <- rownames(x$h) <- fit@modelInfo$parameters_labels

  # Get the second-order devarivates between fixed and estimated parameters:
  model_pars <- fit@modelInfo$parameters_labels %in% fit2@modelInfo$parameters_labels
  nuisance_pars <- !model_pars
  df2_dparamdR <- x$h[nuisance_pars, model_pars]

  # Pick the right coefficients from the VCOV of the first step:
  ACOV <- VCOV1$vcov[fit@modelInfo$parameters_labels[nuisance_pars],
                     fit@modelInfo$parameters_labels[nuisance_pars]]
  # Ham of sandwich estimator:
  B <- t(df2_dparamdR) %*% ACOV %*% df2_dparamdR
  VCOV2$B <- B

  # Get the hessian matrix of second-step model:
  H_inv <- solve(VCOV2$H)
  VCOV2$vcov <- H_inv %*% B %*% H_inv

  # Update the standard errors of the model parameters in the VCOV2 object:
  VCOV2$se <- sqrt(diag(VCOV2$vcov))
  names(VCOV2$se) <- fit2@modelInfo$parameters_labels

  VCOV <- list()
  name_nuisance_pars <- fit@modelInfo$parameters_labels[nuisance_pars]
  VCOV$vcov <- block_diag(list(VCOV1$vcov[name_nuisance_pars, name_nuisance_pars],
                               VCOV2$vcov))
  VCOV$se <- sqrt(diag(VCOV$vcov))

  return(VCOV)

}
