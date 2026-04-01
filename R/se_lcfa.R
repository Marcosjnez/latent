# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 01/04/2026
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

  SE <- general_se(fit = fit, type = type)

  # Create the tables of parameters with standard errors:
  table_se <- fill_in(fit@modelInfo$param, SE$se, miss = 0)

  # Return:
  result <- list()
  result$table_se <- table_se
  result$se <- c(SE$se)
  result$vcov <- SE$vcov
  result$H <- SE$H
  result$B <- SE$B

  return(result)

}

block_diag <- function(mats) {
  if (!is.list(mats) || length(mats) == 0L) {
    stop("`mats` must be a non-empty list of square matrices.")
  }

  dims <- vapply(mats, function(M) {
    if (!is.matrix(M)) stop("All elements of `mats` must be matrices.")
    nr <- nrow(M); nc <- ncol(M)
    if (nr != nc) stop("All matrices must be square.")
    nr
  }, integer(1))

  n_tot <- sum(dims)
  out <- matrix(0, n_tot, n_tot)

  idx <- 0L
  for (k in seq_along(mats)) {
    d <- dims[k]
    r <- (idx + 1L):(idx + d)
    out[r, r] <- mats[[k]]
    idx <- idx + d
  }

  out
}

general_se <- function(fit, type = "standard") {

  ngroups <- fit@data_list$ngroups

  # Get the matrix of second-order derivatives between R and the parameters:
  VAR <- vector("list", length = ngroups)
  for(i in 1:ngroups) {

    # Get the asymptotic correlation matrix:
    ACOVij <- vector("list", length = fit@data_list$correl[[i]]$npatterns)
    for(j in 1:fit@data_list$correl[[i]]$npatterns) {
      ACOVij[[j]] <- fit@data_list$correl[[i]][[j]]$ACOV /
        fit@data_list$correl[[i]][[j]]$nobs
    }

    ACOVi <- block_diag(ACOVij)
    VAR[[i]] <- ACOVi

  }

  ACOV <- block_diag(VAR)

  args <- fit@data_list$args
  args$control <- fit@modelInfo$control_optimizer
  args$control$free_S <- TRUE
  args$control$free_S_diag <- TRUE
  args$do.fit <- FALSE
  fit2 <- do.call(lcfa, args) # AVOID RECOMPUTING CORRELATIONS HERE

  parameters <- fit@Optim$transparameters[fit2@modelInfo$parameters_labels]
  transparameters <- fit@Optim$transparameters[fit2@modelInfo$transparameters_labels]

  control_manifold <- fit2@modelInfo$control_manifold
  control_transform <- fit2@modelInfo$control_transform
  control_estimator <- fit2@modelInfo$control_estimator
  control_optimizer <- fit2@modelInfo$control_optimizer
  control_optimizer$parameters[[1]] <- parameters
  control_optimizer$transparameters[[1]] <- transparameters
  x <- get_hess(control_manifold, control_transform,
                control_estimator, control_optimizer)
  colnames(x$h) <- rownames(x$h) <- fit2@modelInfo$parameters_labels

  model_pars <- fit2@modelInfo$parameters_labels %in%
    fit@modelInfo$parameters_labels
  nuisance_pars <- !model_pars
  df2_dparamdR <- x$h[nuisance_pars, model_pars]

  # Ham of sandwich estimator:
  B <- t(df2_dparamdR) %*% ACOV %*% df2_dparamdR

  # Get the hessian matrix of second-step model:
  SE <- standard_se(fit = fit)
  H_inv <- solve(SE$H)
  SE$vcov <- H_inv %*% B %*% H_inv

  # Update the standard errors of the model parameters in the SE object:
  SE$se <- sqrt(diag(SE$vcov))
  names(SE$se) <- fit@modelInfo$parameters_labels

  return(SE)

}
