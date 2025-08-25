# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 25/08/2025
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
se <- function(fit, confidence = 0.95, digits = 2) {

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

  mylabels <- c(fit@Optim$opt$lca_trans$class,
                unlist(fit@Optim$opt$lca_trans$peta),
                fit@Optim$opt$lca_trans$mu,
                fit@Optim$opt$lca_trans$sigma)
  nparams <- length(mylabels)
  selection <- match(mylabels, fit@Optim$opt$transparameters_labels)
  x <- fit@Optim$opt$transparameters[selection] # Parameters
  H <- computations$hess[selection, selection]
  rownames(H) <- colnames(H) <- mylabels
  dconstraints <- computations$dconstr[selection, ]

  # Parameter types:
  nprob <- length(c(fit@Optim$opt$lca_trans$class,
                    unlist(fit@Optim$opt$lca_trans$peta)))
  nlinear <- length(fit@Optim$opt$lca_trans$mu)
  nsd <- length(fit@Optim$opt$lca_trans$sigma)

  types <- rep(c("prob", "linear", "sd"), times = c(nprob, nlinear, nsd))

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

  K <- t(dconstraints)
  D <- H + t(K) %*% K
  D_inv <- solve(D)
  C[keep, keep] <- D_inv - D_inv %*% t(K) %*% solve(K %*% D_inv %*% t(K)) %*% K %*% D_inv
  se[keep] <- sqrt(diag(C[keep, keep]))
  names(se) <- mylabels
  colnames(C) <- rownames(C) <- mylabels

  reorder <- match(unlist(fit@modelInfo$prob_model), mylabels)
  se_user_model <- fill_list_with_vector(fit@modelInfo$prob_model, se[reorder])
  se_user_model <- allnumeric(se_user_model)

  ci <- sapply(1:nparams, FUN = \(i) conf(x[i], se[i], confidence = confidence,
                                         type = types[i]))
  rownames(ci) <- c("lower", "upper")
  colnames(ci) <- mylabels

  strings <- format(round(ci[, reorder], digits = digits), nsmall = digits)
  ci_ordered <- apply(strings, MARGIN = 2,
                      FUN = \(x) paste(x[1], x[2], sep = "-"))
  ci_user_model <- fill_list_with_vector(fit@modelInfo$prob_model, ci_ordered)

  user_model <- fit@transformed_pars
  est_se <- combine_est_se(user_model, se_user_model)

  result <- list()
  result$se <- se_user_model
  result$vcov <- C[reorder, reorder]
  result$ci <- ci_user_model
  result$est_se <- est_se

  return(result)

}

conf <- function(x, se, confidence, type) {

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

    lower <- exp(log(x) - z * se)
    upper <- exp(log(x) + z * se)

  }


  return(c(lower, upper))

}

