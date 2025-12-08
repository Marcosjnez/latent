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
se.lcfa <- function(fit, type = "standard", model = "model", digits = 5) {

  # Compute lavaan standard errors:
  if(fit@call$mimic == "lavaan") {

    result <- se_information_lavaan(fit = fit, type = type,
                                    model = model, digits = digits)
    return(result)

  }

  # Select the parameters to display according to model type:

  if(model == "user") {

    model <- fit@modelInfo$cfa_trans
    est <- fit@transformed_pars

  } else if(model == "model") {

    model <- fit@modelInfo$cfa_param
    est <- fit@parameters

  } else {
    stop("Unknown model")
  }

  if(type == "standard") {

    SE <- standard_se(fit = fit)

    # # Get the asymptotic correlation matrix:
    # ACOV <- asymptotic_normal(fit@Optim$data_list$correl[[1]]$R)
    #
    # # Get the jacobian:
    # lambda <- fit@transformed_pars[[1]]$lambda
    # psi <- fit@transformed_pars[[1]]$psi
    # theta <- fit@transformed_pars[[1]]$theta
    # select <- match(fit@modelInfo$transparameters_labels,
    #                 unlist(fit@modelInfo$cfa_param[[1]]))
    # select <- select[!is.na(select)]
    # dl_drhat <- dlambda_drhat(lambda, psi)
    # dp_drhat <- dpsi_drhat(lambda, nrow(psi))
    # dt_drhat <- dtheta_drhat(nrow(theta))
    # J <- cbind(dl_drhat, dp_drhat, dt_drhat)
    # J <- J[, select]
    #
    # # Sandwich se estimator:
    # # lower_ind <- which(lower.tri(theta, diag = TRUE))
    # # ACOV <- ACOV[lower_ind, ]; ACOV <- ACOV[, lower_ind]
    # # J <- J[lower_ind, ]
    # B <- t(J) %*% ACOV %*% J
    # H_inv <- solve(SE$h)
    # VAR <- H_inv %*% B %*% H_inv
    #
    # mylabels <- fit@modelInfo$parameters_labels
    # selection <- match(mylabels, fit@modelInfo$transparameters_labels)
    # SE$se[selection] <- sqrt(diag(VAR))

    SE$B <- matrix(, nrow = 0, ncol = 0) # Empty matrix

  } else if(type == "robust") {

    # Use H to compute vcov = solve(H) %*% B %*% solve(H):
    SE <- robust_se(fit = fit)

  } else {
    stop("Unknown type")
  }

  N <- sum(unlist(fit@Optim$data_list$nobs))
  # Select the parameter labels for the table:
  # mylabels <- unlist(model)
  mylabels <- fit@modelInfo$parameters_labels
  selection <- match(mylabels, fit@modelInfo$transparameters_labels)
  VCOV <- SE$vcov[selection, ][, selection]

  # Standard errors:
  vector_se <- SE$se[selection]/sqrt(N)
  names(vector_se) <- mylabels

  # Select the parameter labels for the table:
  mylabels <- unlist(model)
  selection <- match(mylabels, fit@modelInfo$transparameters_labels)

  # Standard errors:
  se <- SE$se[selection]/sqrt(N)
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
  gpsi <- 0.5*gpsi
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
# lambda <- latInspect(fit, what = "loadings", digits = 3)[[1]]
# psi <- latInspect(fit, what = "psi", digits = 3)[[1]]
# theta <- latInspect(fit, what = "theta", digits = 3)[[1]]
# dl_drhat <- dlambda_drhat(lambda, psi)
# dp_drhat <- dpsi_drhat(lambda, p)
# dt_drhat <- dtheta_drhat(p)
#
# p <- nrow(lambda)
# q <- ncol(lambda)
# x <- c(lambda, psi, theta)
# f(x, p, q)
# J2 <- numDeriv::jacobian(func = f, x = x, p = p, q = q)
# max(abs(J2[, 1:27] - dl_drhat))
# max(abs(J2[, -(1:27)][, 1:9] - dp_drhat))
# max(abs(J2[, -(1:36)] - dt_drhat))

