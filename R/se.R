# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 26/08/2025
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
#' @export
se <- function(fit, type = "user", digits = 2) {

  control_manifold <- fit@Optim$control_manifold
  control_transform <- fit@Optim$control_transform
  control_estimator <- fit@Optim$control_estimator
  control_optimizer <- fit@Optim$control

  control_optimizer$parameters[[1]] <- fit@Optim$opt$parameters
  computations <- grad_comp(control_manifold = control_manifold,
                            control_transform = control_transform,
                            control_estimator = control_estimator,
                            control_optimizer = control_optimizer,
                            compute = "outcomes", eps = 1e-04)

  se <- computations$se # standard erros
  names(se) <- fit@Optim$opt$transparameters_labels

  SE <- standard_se(fit, computations, type = type, digits = digits)
  # user_se <- Visser_se(fit, computations, digits = digits)

  result <- list(table = SE$table,
                 table_se = SE$table_se,
                 se = se)

  return(result)

}

standard_se <- function(fit, computations, type = "user",
                        digits = 2) {

  # Select the parameters according to model type:

  if(type == "user") {

    model <- fit@modelInfo$prob_model
    est <- fit@transformed_pars

  } else if(type == "model") {

    model <- fit@modelInfo$model
    est <- fit@parameters

  }

  mylabels <- unlist(model)
  selection <- match(mylabels, fit@Optim$opt$transparameters_labels)

  # Standard errors:
  se <- computations$se[selection]
  # names(se) <- mylabels

  # Tables:
  table_se <- fill_list_with_vector(model, se)
  table_se <- allnumeric(table_se)
  table <- combine_est_se(est, table_se, digits = digits)

  # Return:
  result <- list()
  result$table <- table
  result$table_se <- table_se
  # result$se <- se

  return(result)

}

Visser_se <- function(fit, computations, digits = 2) {

  mylabels <- c(fit@Optim$opt$lca_trans$class,
                unlist(fit@Optim$opt$lca_trans$peta),
                fit@Optim$opt$lca_trans$mu,
                fit@Optim$opt$lca_trans$sigma)

  # Parameter types:
  nprob <- length(c(fit@Optim$opt$lca_trans$class,
                    unlist(fit@Optim$opt$lca_trans$peta)))
  nlinear <- length(fit@Optim$opt$lca_trans$mu)
  nsd <- length(fit@Optim$opt$lca_trans$sigma)
  types <- rep(c("prob", "linear", "sd"), times = c(nprob, nlinear, nsd))

  nparams <- length(mylabels)
  selection <- match(mylabels, fit@Optim$opt$transparameters_labels)
  # computations$se[selection]
  x <- fit@Optim$opt$transparameters[selection] # Parameters
  H <- computations$hess[selection, selection]
  # eigen(H)$values
  rownames(H) <- colnames(H) <- mylabels
  dconstraints <- computations$dconstr[selection, ]

  # Remove probabilities close to 0 or 1 to avoid numerical problems:
  condition_remove <- types == "prob" & (x > 0.999 | x < 0.001)
  remove <- which(condition_remove)
  keep <- which(!condition_remove)
  dconstraints <- dconstraints[keep, , drop = FALSE]
  H <- H[keep, keep]
  se <- vector(length = length(mylabels))
  C <- matrix(NA, nrow = nparams, ncol = nparams)

  # Remove empty columns in the constraints matrix:
  condition_remove <- colSums(dconstraints) < 1
  keep_constr <- which(!condition_remove)
  dconstraints <- dconstraints[, keep_constr, drop = FALSE]

  # The next formula is from the following paper:
  # Visser, I., Raijmakers, M.E.J. and Molenaar, P.C.M. (2000),
  # Confidence intervals for hidden Markov model parameters.
  # British Journal of Mathematical and Statistical Psychology, 53: 317-327.
  # https://doi.org/10.1348/000711000159240

  # Variance-covariance matrix:
  K <- t(dconstraints)
  D <- H + t(K) %*% K
  D_inv <- solve(D)
  C[keep, keep] <- D_inv - D_inv %*% t(K) %*% solve(K %*% D_inv %*% t(K)) %*% K %*% D_inv
  se[keep] <- sqrt(diag(C[keep, keep]))
  names(se) <- mylabels
  colnames(C) <- rownames(C) <- mylabels

  # Model:
  est <- fit@transformed_pars

  # Standard errors:
  reorder <- match(unlist(fit@modelInfo$prob_model), mylabels)
  se_model <- fill_list_with_vector(fit@modelInfo$prob_model, se[reorder])
  se_model <- allnumeric(se_model)
  est_se <- combine_est_se(est, se_model, digits = digits)
  C <- C[reorder, reorder]

  result <- list()
  result$est_se <- est_se
  result$se_model <- se_model
  result$se <- se
  result$vcov <- C

  return(result)

}

ci <- function(fit, type = "user", confidence = 0.95, digits = 2) {

  # Compute standard errors:
  control_manifold <- fit@Optim$control_manifold
  control_transform <- fit@Optim$control_transform
  control_estimator <- fit@Optim$control_estimator
  control_optimizer <- fit@Optim$control

  control_optimizer$parameters[[1]] <- fit@Optim$opt$parameters
  computations <- grad_comp(control_manifold = control_manifold,
                            control_transform = control_transform,
                            control_estimator = control_estimator,
                            control_optimizer = control_optimizer,
                            compute = "outcomes")

  se <- computations$se # standard errors
  names(se) <- fit@Optim$opt$transparameters_labels

  # Number of total parameters and transformed parameters:
  ntrans <- length(fit@Optim$opt$transparameters)
  # Initialize a vector of degrees of freedom:
  ps <- rep(1, times = ntrans)
  # slot for extracting the degrees of freedom from the transformations:
  slot <- 2
  # Get the indices of each transformed parameter:
  indices <- unlist(lapply(fit@Optim$control_transform,
                           FUN = \(x) x$indices_out[[1]]+1L))
  # Update the degrees of freedom of each transformed parameter:
  ps[indices] <- unlist(lapply(fit@Optim$opt$outputs$transformations$vectors,
                               FUN = \(x) x[[slot]]))

  # Compute confidence intervals for each transformed parameter:
  x <- fit@Optim$opt$transparameters
  lower <- x - sqrt(qchisq(confidence, df = ps)) * se
  upper <- x + sqrt(qchisq(confidence, df = ps)) * se
  names(lower) <- names(upper) <- fit@Optim$opt$transparameters_labels

  # Get confidence limits for the user model or raw model parameters:

  # Select the parameters according to model type:

  if(type == "user") {

    model <- fit@modelInfo$prob_model
    est <- fit@transformed_pars

  } else if(type == "model") {

    model <- fit@modelInfo$model
    est <- fit@parameters

  }

  mylabels <- unlist(model)
  selection <- match(mylabels, fit@Optim$opt$transparameters_labels)

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
