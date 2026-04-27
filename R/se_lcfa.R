# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 27/04/2026
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
    stop("`mats` must be a non-empty list of square matrices or nested lists.")
  }

  flatten_mats <- function(x) {
    if (is.matrix(x)) {
      list(x)
    } else if (is.list(x)) {
      unlist(lapply(x, flatten_mats), recursive = FALSE)
    } else {
      stop("All elements must be matrices or lists containing matrices.")
    }
  }

  mats <- flatten_mats(mats)

  if (length(mats) == 0L) {
    stop("No matrices found in `mats`.")
  }

  dims <- vapply(mats, function(M) {
    nr <- nrow(M)
    nc <- ncol(M)
    if (nr != nc) stop("All matrices must be square.")
    nr
  }, integer(1))

  n_tot <- sum(dims)

  rn <- unlist(Map(function(M, d) {
    if (is.null(rownames(M))) rep(NA_character_, d) else rownames(M)
  }, mats, dims), use.names = FALSE)

  cn <- unlist(Map(function(M, d) {
    if (is.null(colnames(M))) rep(NA_character_, d) else colnames(M)
  }, mats, dims), use.names = FALSE)

  out <- matrix(
    0,
    nrow = n_tot,
    ncol = n_tot,
    dimnames = list(rn, cn)
  )

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

  ngroups <- fit@dataList$ngroups

  # Get the asymptotic variance-covariance matrix:
  ACOV_covij <- lapply(fit@dataList$fit_cov, FUN = \(grp) {

    if(length(grp@extra) == 0) {
      object <- list(grp)
    } else {
      object <- grp@extra
    }

    return(lapply(object, FUN = \(ij) ij@Optim$SE$ACOV / ij@dataList$nobs[[1]]))

  })
  ACOV_meansij <- lapply(fit@dataList$fit_means, FUN = \(grp) {

    if(length(grp@extra) == 0) {
      object <- list(grp)
    } else {
      object <- grp@extra
    }

    return(lapply(object, FUN = \(ij) ij@Optim$SE$ACOV / ij@dataList$nobs[[1]]))

  })

  if(fit@modelInfo$control_optimizer$meanstructure) {
    ACOV <- block_diag(c(ACOV_meansij, ACOV_covij))
  } else {
    ACOV <- block_diag(ACOV_covij)
  }

  args <- fit@dataList$args
  # args$control <- fit@modelInfo$control_optimizer
  args$control <- list()
  args$control$free_S <- TRUE
  args$control$free_M <- TRUE
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
                control_estimator, control_optimizer,
                cores = 1L)
  colnames(x$h) <- rownames(x$h) <- fit2@modelInfo$parameters_labels

  model_pars <- fit2@modelInfo$parameters_labels %in%
    fit@modelInfo$parameters_labels
  nuisance_pars <- !model_pars
  df2_dparamdR <- x$h[nuisance_pars, model_pars]

  # Ham of sandwich estimator:
  inter <- intersect(rownames(ACOV), rownames(df2_dparamdR))
  df2_dparamdR <- df2_dparamdR[inter, ]
  ACOV <- ACOV[inter, inter]
  # cbind(rownames(df2_dparamdR), rownames(ACOV))
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
