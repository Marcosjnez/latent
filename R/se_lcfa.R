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

  # Select the model parameters:
  model <- fit@modelInfo$cfa_param
  est <- fit@parameters

  # Collect all the estimators:
  estimators <- unlist(lapply(fit@modelInfo$control_estimator,
                              FUN = \(x) x$estimator))
  # Check if all estimators are ml:
  all_ml <- FALSE
  if(all(estimators == "cfa_ml")) all_ml <- TRUE

  if(type == "information") {
    if(all_ml) {
      SE <- standard_se(fit = fit)
    } else {
      stop("The information method is only available for maximum likelihood")
    }
  } else if(type == "standard" || type == "robust") {
    SE <- general_se(fit = fit, type = type)
  } else {
    stop("Unknown type of standard error estimation")
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

# Derivatives of model matrix wrt parameters:
drhat_dlambda <- function(lambda, psi) {

  # derivative of Lambda wrt Rhat

  p <- nrow(lambda)
  q <- ncol(lambda)
  g1 <- (lambda %*% psi) %x% diag(p)
  g21 <- diag(p) %x% (lambda %*% psi)
  g2 <- g21 %*% dxt(p, q)
  g <- g1 + g2

  return(g)

}
drhat_dpsi <- function(lambda, q) {

  gpsi <- lambda %x% lambda
  gpsi <- 2*gpsi
  diag(gpsi) <- 2*diag(gpsi)

  return(gpsi)

}
drhat_dtheta <- function(p) {

  gtheta <- diag(p) %x% diag(p)
  gtheta <- 0.5*gtheta
  diag(gtheta) <- 2*diag(gtheta)

  return(gtheta)

}

# Second-order derivatives between parameters and model matrix wrt loss:
df2_dtheta_dwls <- function(lambda, psi, theta, w) {

  q <- nrow(psi)
  p <- nrow(theta)
  drhat_dl <- drhat_dlambda(lambda, psi)
  drhat_dp <- drhat_dpsi(lambda, q)
  drhat_dt <- drhat_dtheta(p)
  J <- -w*cbind(drhat_dl, drhat_dp, drhat_dt)

  return(J)

}
df2_dtheta_ml <- function(lambda, psi, theta, w, n) {

  q <- nrow(psi)
  p <- nrow(theta)
  drhat_dl <- drhat_dlambda(lambda, psi)
  drhat_dp <- drhat_dpsi(lambda, q)
  drhat_dt <- drhat_dtheta(p)
  drhat_dtheta <- cbind(drhat_dl, drhat_dp, drhat_dt)

  rhat <- lambda %*% psi %*% t(lambda) + theta
  rhat_inv <- solve(rhat)
  df2_drhatdR <- -rhat_inv %x% rhat_inv
  df2_dthetadR <- 0.5*n*df2_drhatdR %*% drhat_dtheta

  return(df2_dthetadR)

}
df2_dtheta_ml2 <- function(lambda, psi, theta, w) {

  df2_dlambdadS <- function(lambda, psi, theta) {

    p <- nrow(lambda)
    q <- ncol(lambda)
    indexes_p <- c()
    for(i in 1:p) indexes_p[i] <- (i-1) * p + i

    Rhat <- lambda %*% psi %*% t(lambda) + theta
    Rhat_inv <- solve(Rhat)
    LP <- lambda %*% psi

    dRi_res_Ri_dS <- -2*Rhat_inv %x% Rhat_inv
    dRi_res_Ri_dS <- dRi_res_Ri_dS + dRi_res_Ri_dS %*% dxt(p, p)
    h <- t(LP) %x% diag(p) %*% dRi_res_Ri_dS
    h[, indexes_p] <- 0.5*h[, indexes_p]
    # h <- h[, which(lower.tri(theta, diag = TRUE))]

    return(t(h))

  }
  gf2_dpsiPdS <- function(lambda, psi, theta) {

    p <- nrow(lambda)
    q <- ncol(lambda)
    indexes_p <- c()
    for(i in 1:p) indexes_p[i] <- (i-1) * p + i
    indexes_q <- c()
    for(i in 1:q) indexes_q[i] <- (i-1) * q + i

    Rhat <- lambda %*% psi %*% t(lambda) + theta
    Rhat_inv <- solve(Rhat)

    dRi_res_Ri_dS <- -2*Rhat_inv %x% Rhat_inv
    dRi_res_Ri_dS <- dRi_res_Ri_dS + dRi_res_Ri_dS %*% dxt(p, p)
    gPS <- t(lambda) %x% t(lambda) %*% dRi_res_Ri_dS
    gPS[indexes_q, ] <- 0.5*gPS[indexes_q, ]
    gPS[, indexes_p] <- 0.5*gPS[, indexes_p]
    # gPS <- gPS[which(lower.tri(psi, diag = TRUE)),
    #            which(lower.tri(theta, diag = TRUE))]

    return(t(gPS))

  }
  df2_dthetadS <- function(lambda, psi, theta) {

    p <- nrow(lambda)
    q <- ncol(lambda)
    Rhat <- lambda %*% psi %*% t(lambda) + theta
    Rhat_inv <- solve(Rhat)

    dRi_res_Ri_dS <- -2*Rhat_inv %x% Rhat_inv
    dRi_res_Ri_dS <- dRi_res_Ri_dS + dRi_res_Ri_dS %*% dxt(p, p)

    indexes <- seq(1, p*p, by = p+1)
    gUS <- dRi_res_Ri_dS
    gUS[indexes, ] <- 0.5*gUS[indexes, ]
    gUS[, indexes] <- 0.5*gUS[, indexes]
    # gUS <- gUS[lower.tri(theta, diag = TRUE),
    #            lower.tri(theta, diag = TRUE)]

    return(t(gUS))

  }

  J <- 0.5*w*cbind(df2_dlambdadS(lambda, psi, theta),
               gf2_dpsiPdS(lambda, psi, theta),
               df2_dthetadS(lambda, psi, theta))

  return(J)

}

# Sandwhich estimator of standard errors:
general_se <- function(fit, type = "standard") {

  ngroups <- fit@data_list$ngroups
  if(ngroups > 1) stop("Standard errors are not available for multigroup models")

  # Collect all the estimators:
  estimators <- unlist(lapply(fit@modelInfo$control_estimator,
                              FUN = \(x) x$estimator))
  # Check if all the estimators are ml:
  all_ml <- all_ml2 <- all_dwls <- FALSE
  if(all(estimators == "cfa_ml")) all_ml <- TRUE
  if(all(estimators == "cfa_ml2")) all_ml2 <- TRUE
  if(all(estimators == "cfa_dwls")) all_dwls <- TRUE

  # Get the hessian matrix:
  SE <- standard_se(fit = fit)

  # Get the sample size and weight scalar:
  N <- fit@data_list$nobs[[1]]
  w <- fit@modelInfo$control_estimator[[1]]$w

  # Get the matrix of second-order derivatives between R and the parameters:
  lambda <- fit@transformed_pars[[1]]$lambda
  psi <- fit@transformed_pars[[1]]$psi
  psi[upper.tri(psi)] <- t(psi)[upper.tri(psi)]
  theta <- fit@transformed_pars[[1]]$theta
  theta[upper.tri(theta)] <- t(theta)[upper.tri(theta)]
  if(all_ml) { # SHOULD COMPUTE THIS AUTOMATICALY IN C++
    df2_dparamdR <- df2_dtheta_ml(lambda, psi, theta, w, N)
  } else if(all_ml2) {
    df2_dparamdR <- df2_dtheta_ml2(lambda, psi, theta, w)
  } else if(all_dwls) {
    df2_dparamdR <- df2_dtheta_dwls(lambda, psi, theta, w)
  } else {
    stop("All the estimators should be the same")
  }

  # Select the model parameters from df2_dparamdR:
  select <- match(fit@modelInfo$transparameters_labels,
                  unlist(fit@modelInfo$cfa_param[[1]]))
  select <- select[!is.na(select)]
  df2_dparamdR <- df2_dparamdR[, select]

  # Get the asymptotic correlation matrix:
  if(type == "standard") {
    ACOV <- asymptotic_normal(fit@data_list$correl[[1]]$R)
  } else if(type == "robust") {
    # CHECK FOR RAW DATA AVAILABILITY
    ACOV <- asymptotic_general(fit@data_list$data[[1]])
  } else {
    stop("Unknown type")
  }

  # Sandwich estimator of standard errors:
  B <- t(df2_dparamdR) %*% (ACOV/N) %*% df2_dparamdR
  H_inv <- solve(SE$h)
  VAR <- H_inv %*% B %*% H_inv

  # Update the standard errors of the model parameters in the SE object:
  mylabels <- fit@modelInfo$parameters_labels
  selection <- match(mylabels, fit@modelInfo$transparameters_labels)
  SE$se[selection] <- sqrt(diag(VAR))

  # Store the ham of the sandwich:
  SE$B <- B

  # UPDATE SE$VCOV AND ALL OTHER ELEMENTS IN SE$se

  return(SE)

}

# se_information_lavaan <- function(fit, type = "standard", model = "model",
#                                   digits = 3) {
#
#   ngroups <- fit@Optim$data_list$ngroups
#   se <- VCOV <- table_se <- vector("list", length = ngroups)
#   names(se) <- names(VCOV) <- names(table_se) <-
#     fit@Optim$data_list$group_label
#   fit@Optim$data_list$item_label
#
#   for(i in 1:ngroups) {
#
#     # Compute the model matrix:
#     N <- fit@modelInfo$nobs[[i]]
#     p <- fit@Optim$data_list$nitems[[i]]
#     q <- fit@Optim$data_list$nfactors[[i]]
#     lambda <- fit@transformed_pars[[i]]$lambda
#     psi <- fit@transformed_pars[[i]]$psi
#     theta <- fit@transformed_pars[[i]]$theta
#     rhat <- fit@transformed_pars[[i]]$model
#     rhat_inv <- solve(rhat)
#
#     # Get the jacobian:
#     select <- match(fit@modelInfo$transparameters_labels,
#                     unlist(fit@modelInfo$cfa_param[[i]]))
#     select <- select[!is.na(select)]
#
#     dl_drhat <- dlambda_drhat(lambda, psi)
#     dp_drhat <- dpsi_drhat(lambda, q)
#     dt_drhat <- dtheta_drhat(p)
#     J <- cbind(dl_drhat, dp_drhat, dt_drhat)
#     J <- J[, select]
#     # numDeriv::jacobian()
#
#     # From p.13 of
#     # https://bpspsychub.onlinelibrary.wiley.com/doi/pdf/10.1111/bmsp.70003
#     exp_inf <- 0.5*rhat_inv %x% rhat_inv
#     Iexp <- t(J) %*% exp_inf %*% J
#     VCOV[[i]] <- solve(Iexp)
#     se[[i]] <- sqrt(diag(VCOV[[i]]))
#     se[[i]] <- round(se[[i]], digits = digits)
#
#     # Build the tables:
#     select <- match(unlist(fit@modelInfo$cfa_param[[i]]),
#                     fit@modelInfo$transparameters_labels)
#     select <- select[!is.na(select)]
#     select <- unique(select)
#
#     param_labels <- fit@modelInfo$transparameters_labels[select]
#
#     select_se_lambda <- match(fit@modelInfo$cfa_param[[i]]$lambda,
#                               param_labels)
#     select_se_lambda <- select_se_lambda[!is.na(select_se_lambda)]
#     select_lambda <- match(param_labels,
#                            fit@modelInfo$cfa_param[[i]]$lambda)
#     select_lambda <- select_lambda[!is.na(select_lambda)]
#     lambda_se <- matrix(0, nrow = p, ncol = q)
#     lambda_se[select_lambda] <- se[[i]][select_se_lambda]
#
#     lower_psi <- lower.tri(diag(q), diag = TRUE)
#     select_se_psi <- match(fit@modelInfo$cfa_param[[i]]$psi[lower_psi],
#                            param_labels)
#     select_se_psi <- select_se_psi[!is.na(select_se_psi)]
#     select_psi <- match(param_labels,
#                         fit@modelInfo$cfa_param[[i]]$psi)
#     select_psi <- select_psi[!is.na(select_psi)]
#     select_psi <- unique(select_psi)
#     psi_se <- matrix(0, nrow = q, ncol = q)
#     psi_se[select_psi] <- se[[i]][select_se_psi]
#
#     lower_theta <- lower.tri(diag(p), diag = TRUE)
#     select_se_theta <- match(fit@modelInfo$cfa_param[[i]]$theta[lower_theta],
#                              param_labels)
#     select_se_theta <- select_se_theta[!is.na(select_se_theta)]
#     select_theta <- match(param_labels,
#                           fit@modelInfo$cfa_param[[i]]$theta)
#     select_theta <- select_theta[!is.na(select_theta)]
#     select_theta <- unique(select_theta)
#     theta_se <- matrix(0, nrow = p, ncol = p)
#     theta_se[select_theta] <- se[[i]][select_se_theta]
#
#     table_se[[i]] <- list(lambda = lambda_se,
#                           psi = psi_se,
#                           theta = theta_se)
#
#     # Select the parameter labels for the table:
#     mylabels <- unlist(model)
#     selection <- match(mylabels, fit@modelInfo$transparameters_labels)
#
#   }
#
#   # Return:
#   result <- list()
#   result$table_se <- table_se
#   result$se <- se
#   result$vcov <- VCOV
#
#   return(result)
#
# }
