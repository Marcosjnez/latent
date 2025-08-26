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
se <- function(fit, confidence = 0.95, digits = 3) {

  control_manifold <- fit@Optim$control_manifold
  control_transform <- fit@Optim$control_transform
  control_estimator <- fit@Optim$control_estimator
  control_optimizer <- fit@Optim$control

  control_optimizer$parameters[[1]] <- fit@Optim$opt$parameters
  computations <- grad_comp(control_manifold = control_manifold,
                            control_transform = control_transform,
                            control_estimator = control_estimator,
                            control_optimizer = control_optimizer,
                            compute = "dconstr", eps = 1e-04)

  raw_model <- standard_se(fit, computations, type = "model",
                           confidence = confidence, digits = digits)
  user_model <- standard_se(fit, computations, type = "user",
                            confidence = confidence, digits = digits)
  # user_model <- Visser_se(fit, computations,
  #                         confidence = confidence, digits = digits)

  se <- raw_model$se
  nparam <- length(se)
  labels <- rownames(raw_model$se)
  selection <- match(labels, fit@Optim$opt$transparameters_labels)
  x <- fit@Optim$opt$transparameters[selection]
  types <- rep("linear", times = nparam)
  ci <- sapply(1:nparam, FUN = \(i) conf(x[i], se[i],
                                         confidence = confidence,
                                         type = types[i]))

  rownames(ci) <- c("lower", "upper")
  colnames(ci) <- labels
  est <- fit@parameters
  reorder <- match(unlist(fit@modelInfo$model), labels)
  ci_lower <- fill_list_with_vector(fit@modelInfo$model,
                                    ci[1, reorder])
  ci_lower <- allnumeric(ci_lower)
  ci_upper <- fill_list_with_vector(fit@modelInfo$model,
                                    ci[2, reorder])
  ci_upper <- allnumeric(ci_upper)
  est_ci <- combine_est_ci(ci_lower, est, ci_upper, digits = digits)

  result <- list(user_model = user_model,
                 raw_model = raw_model)

  # result$est_ci <- est_ci
  # result$ci_lower <- ci_lower
  # result$ci_upper <- ci_upper

  return(result)

}

standard_se <- function(fit, computations, type = "user",
                        confidence = 0.95, digits = 3) {

  if(type == "user") {

    all_labels <- c(fit@Optim$opt$lca_trans$class,
                    unlist(fit@Optim$opt$lca_trans$peta),
                    fit@Optim$opt$lca_trans$mu,
                    fit@Optim$opt$lca_trans$sigma)
    all_labels <- unname(all_labels)
    # Model:
    est <- fit@transformed_pars
    reorder <- match(unlist(fit@modelInfo$prob_model), all_labels)

  } else if(type == "model") {

    all_labels <- c(fit@Optim$opt$lca_trans$theta,
                    unlist(fit@Optim$opt$lca_trans$eta),
                    fit@Optim$opt$lca_trans$mu,
                    fit@Optim$opt$lca_trans$s)
    all_labels <- unname(all_labels)
    # Model:
    est <- fit@parameters
    reorder <- match(unlist(fit@modelInfo$model), all_labels)

  } else {

    stop("Unkown type")

  }

  nparam <- length(all_labels)

  trans_labels <- fit@Optim$opt$transparameters_labels
  selection <- match(all_labels, trans_labels)

  C <- computations$vcov[selection, selection]
  rownames(C) <- colnames(C) <- all_labels

  # Variance-covariance matrix:
  se <- computations$se[selection]
  names(se) <- all_labels

  # Standard errors:
  se_model <- fill_list_with_vector(fit@modelInfo$model, se[reorder])
  se_model <- allnumeric(se_model)
  est_se <- combine_est_se(est, se_model, digits = digits)

  # Return:
  result <- list()
  result$est_se <- est_se
  result$se_model <- se_model
  result$se <- se
  result$vcov <- C

  return(result)

}

Visser_se <- function(fit, computations,
                      confidence = 0.95, digits = 3) {

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

conf <- function(x, se, confidence = 0.95, type = "linear") {

  z <- stats::qnorm((1+confidence)/2)

  if(type == "prob") {

    logit_prob <- stats::qlogis(x)
    logit_se <- se / (x * (1 - x))  # Delta method

    logit_lower <- logit_prob - z * logit_se
    logit_upper <- logit_prob + z * logit_se

    # Back-transform to probability scale
    lower <- stats::plogis(logit_lower)
    upper <- stats::plogis(logit_upper)

  } else if(type == "linear") {

    lower <- x - z * se
    upper <- x + z * se

  } else if(type == "sd") {

    se <- se/x
    lower <- exp(log(x) - z * se)
    upper <- exp(log(x) + z * se)

  }


  return(c(lower, upper))

}
