# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 28/11/2025
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
#' @method se lcfa
#' @export
se.lcfa <- function(fit, type = "standard", digits = 5) {

  # Compute lavaan standard errors:
  if(fit@call$mimic == "lavaan") {

    result <- se_information_lavaan(fit = fit, type = type,
                                    model = model, digits = digits)
    return(result)

  }

  # Select the parameters to display according to model type:

  model <- fit@modelInfo$cfa_param
  est <- fit@parameters

  # Check if ml or uls are used:

  estimators <- unlist(lapply(fit@modelInfo$control_estimator,
                              FUN = \(x) x$estimator))
  if(all(estimators == "cfa_ml")) {
    all_ml <- TRUE
  } else {
    all_ml <- FALSE
  }

  if(type == "standard") {

    if(all_ml) {
      SE <- standard_se(fit = fit)
    } else {
      SE <- robust_general(fit = fit)
    }

  } else if(type == "robust") {

    if(all_ml) {
      # Use H to compute vcov = solve(H) %*% B %*% solve(H):
      SE <- robust_se(fit = fit)
    } else {
      SE <- robust_general(fit = fit)
    }

  } else {
    stop("Unknown type")
  }

  # Select the parameter labels for the table:
  # mylabels <- unlist(model)
  mylabels <- fit@modelInfo$parameters_labels
  selection <- match(mylabels, fit@modelInfo$transparameters_labels)
  VCOV <- SE$vcov[selection, ][, selection]

  # Standard errors:
  vector_se <- SE$se[selection]
  names(vector_se) <- mylabels

  # Select the parameter labels for the table:
  mylabels <- unlist(model)
  selection <- match(mylabels, fit@modelInfo$transparameters_labels)

  # Standard errors:
  se <- SE$se[selection]
  names(se) <- mylabels

  # Tables:
  table_se <- fill_list_with_vector(model, round(se, digits = digits))
  table_se <- allnumeric(table_se)
  # table <- combine_est_se(est, table_se, digits = digits)

  # Return:
  result <- list()
  # result$table <- table
  result$table_se <- table_se
  result$se <- c(vector_se)
  result$vcov <- VCOV
  result$B <- SE$B

  return(result)

}

dlambda_drhat <- function(lambda, psi) {

  # derivative of Lambda wrt Rhat

  p <- nrow(lambda)
  q <- ncol(lambda)
  g1 <- (lambda %*% psi) %x% diag(p)
  g21 <- diag(p) %x% (lambda %*% psi)
  g2 <- g21 %*% dxt(p, q)
  g <- g1 + g2

  return(g)

}
dpsi_drhat <- function(lambda, q) {

  gpsi <- lambda %x% lambda
  gpsi <- 2*gpsi
  diag(gpsi) <- 2*diag(gpsi)

  return(gpsi)

}
dtheta_drhat <- function(p) {

  gtheta <- diag(p) %x% diag(p)
  gtheta <- 0.5*gtheta
  diag(gtheta) <- 2*diag(gtheta)

  return(gtheta)

}

se_information_lavaan <- function(fit, type = "standard", model = "model",
                                  digits = 3) {

  ngroups <- fit@Optim$data_list$ngroups
  se <- VCOV <- table_se <- vector("list", length = ngroups)
  names(se) <- names(VCOV) <- names(table_se) <-
    fit@Optim$data_list$group_label
  fit@Optim$data_list$item_label

  for(i in 1:ngroups) {

    # Compute the model matrix:
    N <- fit@modelInfo$nobs[[i]]
    p <- fit@Optim$data_list$nitems[[i]]
    q <- fit@Optim$data_list$nfactors[[i]]
    lambda <- fit@transformed_pars[[i]]$lambda
    psi <- fit@transformed_pars[[i]]$psi
    theta <- fit@transformed_pars[[i]]$theta
    rhat <- fit@transformed_pars[[i]]$model
    rhat_inv <- solve(rhat)

    # Get the jacobian:
    select <- match(fit@modelInfo$transparameters_labels,
                    unlist(fit@modelInfo$cfa_param[[i]]))
    select <- select[!is.na(select)]

    dl_drhat <- dlambda_drhat(lambda, psi)
    dp_drhat <- dpsi_drhat(lambda, q)
    dt_drhat <- dtheta_drhat(p)
    J <- cbind(dl_drhat, dp_drhat, dt_drhat)
    J <- J[, select]
    # numDeriv::jacobian()

    # From p.13 of
    # https://bpspsychub.onlinelibrary.wiley.com/doi/pdf/10.1111/bmsp.70003
    exp_inf <- 0.5*rhat_inv %x% rhat_inv
    Iexp <- t(J) %*% exp_inf %*% J
    VCOV[[i]] <- solve(Iexp)
    se[[i]] <- sqrt(diag(VCOV[[i]]))
    se[[i]] <- round(se[[i]], digits = digits)

    # Build the tables:
    select <- match(unlist(fit@modelInfo$cfa_param[[i]]),
                    fit@modelInfo$transparameters_labels)
    select <- select[!is.na(select)]
    select <- unique(select)

    param_labels <- fit@modelInfo$transparameters_labels[select]

    select_se_lambda <- match(fit@modelInfo$cfa_param[[i]]$lambda,
                              param_labels)
    select_se_lambda <- select_se_lambda[!is.na(select_se_lambda)]
    select_lambda <- match(param_labels,
                           fit@modelInfo$cfa_param[[i]]$lambda)
    select_lambda <- select_lambda[!is.na(select_lambda)]
    lambda_se <- matrix(0, nrow = p, ncol = q)
    lambda_se[select_lambda] <- se[[i]][select_se_lambda]

    lower_psi <- lower.tri(diag(q), diag = TRUE)
    select_se_psi <- match(fit@modelInfo$cfa_param[[i]]$psi[lower_psi],
                           param_labels)
    select_se_psi <- select_se_psi[!is.na(select_se_psi)]
    select_psi <- match(param_labels,
                        fit@modelInfo$cfa_param[[i]]$psi)
    select_psi <- select_psi[!is.na(select_psi)]
    select_psi <- unique(select_psi)
    psi_se <- matrix(0, nrow = q, ncol = q)
    psi_se[select_psi] <- se[[i]][select_se_psi]

    lower_theta <- lower.tri(diag(p), diag = TRUE)
    select_se_theta <- match(fit@modelInfo$cfa_param[[i]]$theta[lower_theta],
                             param_labels)
    select_se_theta <- select_se_theta[!is.na(select_se_theta)]
    select_theta <- match(param_labels,
                          fit@modelInfo$cfa_param[[i]]$theta)
    select_theta <- select_theta[!is.na(select_theta)]
    select_theta <- unique(select_theta)
    theta_se <- matrix(0, nrow = p, ncol = p)
    theta_se[select_theta] <- se[[i]][select_se_theta]

    table_se[[i]] <- list(lambda = lambda_se,
                          psi = psi_se,
                          theta = theta_se)

    # Select the parameter labels for the table:
    mylabels <- unlist(model)
    selection <- match(mylabels, fit@modelInfo$transparameters_labels)

  }

  # Return:
  result <- list()
  result$table_se <- table_se
  result$se <- se
  result$vcov <- VCOV

  return(result)

}

robust_general <- function(fit) {

  SE <- standard_se(fit = fit)
  N <- sum(unlist(fit@data_list$nobs))

  # Get the asymptotic correlation matrix:
  ACOV <- asymptotic_normal(fit@data_list$correl[[1]]$R)
  # ACOV <- asymptotic_general(fit@data_list$correl[[1]]$R)

  # Get the jacobian:
  lambda <- fit@transformed_pars[[1]]$lambda
  psi <- fit@transformed_pars[[1]]$psi
  psi[upper.tri(psi)] <- t(psi)[upper.tri(psi)]
  theta <- fit@transformed_pars[[1]]$theta
  theta[upper.tri(theta)] <- t(theta)[upper.tri(theta)]

  select <- match(fit@modelInfo$transparameters_labels,
                  unlist(fit@modelInfo$cfa_param[[1]]))
  select <- select[!is.na(select)]
  dl_drhat <- dlambda_drhat(lambda, psi)
  dp_drhat <- dpsi_drhat(lambda, nrow(psi))
  dt_drhat <- dtheta_drhat(nrow(theta))
  J <- cbind(dl_drhat, dp_drhat, dt_drhat)
  J <- J[, select]

  # XXX <- gLPS_uls(diag(9), lambda, psi)
  # dim(XXX)

  # control_manifold <- fit@modelInfo$control_manifold
  # control_transform <- fit@modelInfo$control_transform
  # control_estimator <- fit@modelInfo$control_estimator
  # control_optimizer <- fit@modelInfo$control
  # control_optimizer$parameters[[1]] <- fit@Optim$parameters
  # control_optimizer$transparameters[[1]] <- fit@Optim$transparameters
  # JAC <- get_jacob(control_manifold = control_manifold,
  #                  control_transform = control_transform,
  #                  control_estimator = control_estimator,
  #                  control_optimizer = control_optimizer)
  # JACOB <- matrix(JAC$matrices[[1]][[1]], nrow = 45, ncol = 78)
  # mylabels <- fit@modelInfo$parameters_labels
  # selection <- match(mylabels, fit@modelInfo$transparameters_labels)
  # JACOB <- JACOB[, selection]
  # lower_inds <- which(lower.tri(theta, diag = TRUE))
  # B <- t(JACOB) %*% (ACOV/N)[lower_inds, ][, lower_inds] %*% JACOB
  # H_inv <- solve(SE$h)
  # VAR <- H_inv %*% B %*% H_inv
  # SE$se[selection] <- 2*sqrt(diag(VAR))

  # Sandwich se estimator:
  B <- t(J) %*% (ACOV/N) %*% J
  H_inv <- solve(SE$h)
  VAR <- H_inv %*% B %*% H_inv

  mylabels <- fit@modelInfo$parameters_labels
  selection <- match(mylabels, fit@modelInfo$transparameters_labels)
  SE$se[selection] <- sqrt(diag(VAR))

  # SE$B <- matrix(, nrow = 0, ncol = 0) # Empty matrix
  SE$B <- B

  return(SE)

}

# f <- function(x, p, q) {
#
#   pq <- p*q
#   qq <- q*q
#   x1 <- x[1:pq]
#   x2 <- x[1:(q*(q+1)/2) + pq]
#   x3 <- x[1:(p*(p+1)/2) + pq + (q*(q+1)/2)]
#   lambda <- matrix(x1, p, q)
#   psi <- matrix(0, q, q)
#   psi[lower.tri(psi, diag = TRUE)] <- x2
#   psi[upper.tri(psi)] <- t(psi)[upper.tri(psi)]
#   theta <- matrix(0, p, p)
#   theta[lower.tri(theta, diag = TRUE)] <- x3
#   theta[upper.tri(theta)] <- t(theta)[upper.tri(theta)]
#
#   Rhat <- lambda %*% psi %*% t(lambda) + theta
#   return(c(Rhat))
#
# }
#
# x1 <- 1:pq
# x2 <- 1:(q*(q+1)/2) + pq
# x3 <- 1:(p*(p+1)/2) + pq + (q*(q+1)/2)
# lambda <- latInspect(fit, what = "loadings", digits = 6)[[1]]
# psi <- latInspect(fit, what = "psi", digits = 6)[[1]]
# theta <- latInspect(fit, what = "theta", digits = 6)[[1]]
# p <- nrow(lambda)
# q <- ncol(lambda)
# dl_drhat <- dlambda_drhat(lambda, psi)
# dp_drhat <- dpsi_drhat(lambda, p)
# dt_drhat <- dtheta_drhat(p)
#
# x <- c(lambda,
#        psi[lower.tri(psi, diag = TRUE)],
#        theta[lower.tri(theta, diag = TRUE)])
# f(x, p, q)
# J2 <- numDeriv::jacobian(func = f, x = x, p = p, q = q)
# max(abs(J2[, x1] - dl_drhat))
# max(abs(J2[, x2] - dp_drhat[, lower.tri(psi, diag = TRUE)]))
# max(abs(J2[, x3] - dt_drhat[, lower.tri(theta, diag = TRUE)]))

