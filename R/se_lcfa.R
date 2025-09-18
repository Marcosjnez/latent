# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 29/08/2025
#'
#' @title
#' Standard Errors
#' @description
#'
#' Compute standard errors.
#'
#' @usage
#'
#' se_cfa(fit)
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
se.lcfa <- function(fit, type = "standard", model = "user", digits = 2) {

  # x <- standard_se(fit)
  # x$se

  # control_manifold <- fit@Optim$control_manifold
  # control_transform <- fit@Optim$control_transform
  # control_estimator <- fit@Optim$control_estimator
  # control_optimizer <- fit@Optim$control
  # control_optimizer$parameters[[1]] <- fit@Optim$opt$parameters
  # control_optimizer$transparameters[[1]] <- fit@Optim$opt$transparameters
  # computations <- grad_comp(control_manifold = control_manifold,
  #                           control_transform = control_transform,
  #                           control_estimator = control_estimator,
  #                           control_optimizer = control_optimizer,
  #                           compute = "h")

  # Compute the model matrix:
  N <- fit@modelInfo$nobs
  p <- fit@Optim$data_list$nitems
  q <- fit@Optim$data_list$nfactors
  lambda <- fit@transformed_pars[[1]]$lambda
  psi <- fit@transformed_pars[[1]]$psi
  theta <- fit@transformed_pars[[1]]$theta
  rhat <- lambda %*% psi %*% t(lambda) + theta
  rhat_inv <- solve(rhat)

  # Get the jacobian:
  select <- match(fit@modelInfo$transparameters_labels,
                  unlist(fit@modelInfo$cfa_param))
  select <- select[!is.na(select)]

  dl_drhat <- dlambda_drhat(lambda, psi)
  dp_drhat <- dpsi_drhat(lambda, q)
  dt_drhat <- dtheta_drhat(p)
  J <- cbind(dl_drhat, dp_drhat, dt_drhat)
  J <- J[, select]

  exp_inf <- 0.5*rhat_inv %x% rhat_inv
  # exp_inf <- computations$h
  Iexp <- t(J) %*% exp_inf %*% J
  VCOV <- solve(Iexp) / N
  se <- sqrt(diag(VCOV))


  # Return:
  result <- list()
  # result$table <- table
  # result$table_se <- table_se
  result$se <- se
  result$vcov <- VCOV

  return(result)

}

dlambda_drhat <- function(lambda, phi) {

  # derivative of Lambda wrt Rhat

  p <- nrow(lambda)
  q <- ncol(lambda)
  g1 <- (lambda %*% phi) %x% diag(p)
  g21 <- diag(p) %x% (lambda %*% phi)
  g2 <- g21 %*% dxt(p, q)
  g <- g1 + g2

  return(g)

}
dpsi_drhat <- function(lambda, q) {

  g <- lambda %x% lambda

  return(g)

}
dtheta_drhat <- function(p) {

  gtheta <- diag(p) %x% diag(p)

  return(gtheta)

}

# f <- function(x, p, q) {
#
#   pq <- p*q
#   qq <- q*q
#   lambda <- matrix(x[1:pq], p, q)
#   x1 <- x[-(1:pq)]
#   psi <- matrix(x1[1:qq], q, q)
#   x2 <- x1[-(1:qq)]
#   theta <- matrix(x2, p, p)
#
#   Rhat <- lambda %*% psi %*% t(lambda) + theta
#   return(c(Rhat))
#
# }
# p <- nrow(lambda)
# q <- ncol(lambda)
# x <- c(lambda, psi, theta)
# f(x, p, q)
# J2 <- numDeriv::jacobian(func = f, x = x, p = p, q = q)
# max(abs(J2[, 1:27] - dl_drhat))
# max(abs(J2[, -(1:27)][, 1:9] - dp_drhat))
# max(abs(J2[, -(1:36)] - dt_drhat))

