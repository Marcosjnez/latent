# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 15/08/2026
#'
#' LatentGold-Style Robust Variance-Covariance Matrix
#'
#' Construct an information matrix corresponding to a LatentGold-style
#' sandwich covariance estimator for a fitted latent class model.
#'
#' @param fit A fitted object of class \code{"llca"}.
#'
#' @details
#' Let \eqn{H} denote the Hessian and \eqn{B} the empirical covariance matrix
#' of score contributions. The robust covariance matrix is
#' \deqn{H^{-1} B H^{-1}.}
#'
#' Because \code{\link{vcov.latent}} expects an information matrix that is
#' subsequently inverted, this function returns the equivalent matrix
#' \deqn{H B^{-1} H.}
#'
#' @return
#' A symmetric numeric matrix representing the robust information matrix.
#'
#' @keywords internal
#' @export
robust.llca <- function(fit) {

  if(fit@dataList$nobs <= 1L) {
    stop("Robust standard errors require more than one observation.")
  }

  #### Compute the Hessian ####

  fit@modelInfo$control_optimizer$parameters[[1]] <-
    fit@Optim$parameters
  fit@modelInfo$control_optimizer$transparameters[[1]] <-
    fit@Optim$transparameters

  H <- hessian(fit)

  #### Collect the gradient by response pattern ####

  transparameters_labels <- fit@modelInfo$transparameters_labels
  pattern_weights <- fit@dataList$pattern_weights
  npatterns <- fit@dataList$npatterns
  nparam <- fit@modelInfo$nparam
  nobs <- fit@dataList$nobs
  nclasses <- fit@dataList$nclasses

  # Remove the lca estimator but keep everything else:
  control_estimator <- fit@modelInfo$control_estimator[-1L]
  K <- length(control_estimator)

  pattern_struct <- vector("list", length = npatterns)
  for(s in seq_len(npatterns)) {

    # Put back the lca estimator but for each response pattern:
    pattern_struct[[s]] <- list(estimator = "lca",
                                parameters = list(fit@modelInfo$trans$class[s, ],
                                                  fit@modelInfo$trans$loglik[s, ]),
                                extra = list(S = 1L,
                                             I = nclasses,
                                             weights = pattern_weights[s]))
  }

  pattern_estimators <- create_estimators(estimators = pattern_struct,
                                          structures = fit@modelInfo$trans)
  control_estimator <- c(control_estimator, pattern_estimators)

  #### Compute the B matrix ####

  control_manifold <- fit@modelInfo$control_manifold
  control_transform <- fit@modelInfo$control_transform
  control_optimizer <- fit@modelInfo$control_optimizer
  B <- matrix(0, nrow = nparam, ncol = nparam)
  for(s in seq_len(npatterns)) {
    idx <- c(seq_len(K), K+s)
    computations <- get_grad(control_manifold = control_manifold,
                             control_transform = control_transform,
                             control_estimator = control_estimator[idx],
                             control_optimizer = control_optimizer)
    gradient <- computations$g/pattern_weights[s]
    B <- B + pattern_weights[s] * gradient %*% t(gradient)
  }

  B <- B*nobs/(nobs-1L)

  #### Sandwhich estimator ####

  Hinv <- solve(H)
  VCOV <- Hinv %*% B %*% Hinv
  VCOV <- (VCOV + t(VCOV))/2

  #### Return ####

  result <- list(H = H, VCOV = VCOV, se = sqrt(diag(VCOV)))

  return(result)

}
