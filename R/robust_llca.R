# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 25/08/2026
#'
#' LatentGold-Style Robust Variance-Covariance Matrix
#'
#' Construct a sandwich covariance estimator for a fitted latent class model.
#'
#' @param fit A fitted object of class \code{"llca"}.
#'
#' @return A list containing the Hessian, empirical score covariance, robust
#'   variance-covariance matrix, standard errors, and method label.
#'
#' @method robust llca
#' @export
robust.llca <- function(fit) {

  if(fit@dataList$nobs <= 1L) {
    stop("Robust standard errors require more than one observation.")
  }

  #### Compute the Hessian ####

  fit@modelInfo$control_optimizer$parameters[[1L]] <-
    fit@Optim$parameters
  fit@modelInfo$control_optimizer$transparameters[[1L]] <-
    fit@Optim$transparameters

  labels <- fit@modelInfo$parameters_labels
  H <- hessian(fit)
  H <- validate_covariance_matrix(H, labels = labels,
                                  object_name = "LCA Hessian")
  Hinv <- invert_hessian_latent(
    fit = fit,
    H = H,
    labels = labels,
    object_name = "LCA Hessian"
  )

  #### Collect the gradient by response pattern ####

  pattern_weights <- fit@dataList$pattern_weights
  npatterns <- fit@dataList$npatterns
  nparam <- fit@modelInfo$nparam
  nobs <- fit@dataList$nobs
  nclasses <- fit@dataList$nclasses

  if(any(!is.finite(pattern_weights)) ||
     any(pattern_weights <= 0)) {
    stop("Response-pattern weights must be positive and finite.")
  }

  # Remove the aggregate LCA estimator but keep all other estimators.
  control_estimator <- fit@modelInfo$control_estimator[-1L]
  K <- length(control_estimator)

  pattern_struct <- vector("list", length = npatterns)

  for(s in seq_len(npatterns)) {

    pattern_struct[[s]] <- list(
      estimator = "lca",
      parameters = list(fit@modelInfo$trans$class[s, ],
                        fit@modelInfo$trans$loglik[s, ]),
      extra = list(S = 1L,
                   I = nclasses,
                   weights = pattern_weights[s])
    )

  }

  pattern_estimators <- create_estimators(
    estimators = pattern_struct,
    structures = fit@modelInfo$trans
  )
  control_estimator <- c(control_estimator, pattern_estimators)

  #### Empirical score covariance ####

  control_manifold <- fit@modelInfo$control_manifold
  control_transform <- fit@modelInfo$control_transform
  control_optimizer <- fit@modelInfo$control_optimizer

  B <- matrix(0,
              nrow = nparam,
              ncol = nparam,
              dimnames = list(labels, labels))

  for(s in seq_len(npatterns)) {

    idx <- c(seq_len(K), K+s)

    computations <- get_grad(
      control_manifold = control_manifold,
      control_transform = control_transform,
      control_estimator = control_estimator[idx],
      control_optimizer = control_optimizer
    )

    gradient <- computations$g/pattern_weights[s]
    B <- B + pattern_weights[s]*tcrossprod(gradient)

  }

  B <- B*nobs/(nobs-1L)
  B <- validate_covariance_matrix(B, labels = labels,
                                  object_name = "empirical score covariance")

  #### Sandwich estimator ####

  VCOV <- Hinv %*% B %*% Hinv
  VCOV <- validate_covariance_matrix(
    VCOV,
    labels = labels,
    object_name = "robust LCA variance-covariance matrix"
  )

  se <- standard_errors_from_vcov(
    VCOV,
    object_name = "robust LCA variance-covariance matrix"
  )

  #### Result ####

  result <- list(H = H,
                 B = B,
                 VCOV = VCOV,
                 se = se,
                 type = "robust")

  return(result)

}
