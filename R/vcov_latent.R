# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 21/08/2026
#'
#' Variance-Covariance Matrix for Latent Objects
#'
#' @param fit A fitted object inheriting from class \code{"latent"}.
#' @param v Variance-covariance matrix of the freely estimated, untransformed
#'   parameters.
#' @param parameters Optional parameter specification identifying the parameters
#'   or transformed parameters to return.
#'
#' @return A list containing the selected variance-covariance matrix and
#'   standard errors.
#'
#' @details
#' Transformation covariances are propagated incrementally in C++. The internal
#' Jacobian workspace used by that computation is not returned because it is not
#' the cumulative Jacobian exposed by \code{jacobian.latent()}.
#'
#' @method vcov latent
#' @export
vcov.latent <- function(fit, v, parameters = NULL) {

  #### Check inputs ####

  if(!inherits(fit, "latent")) {
    stop("fit must inherit from class 'latent'.")
  }

  labels <- fit@modelInfo$parameters_labels

  v <- validate_covariance_matrix(
    v,
    labels = labels,
    object_name = "untransformed variance-covariance matrix"
  )

  if(is.null(parameters)) {

    parameters <- labels

  } else {

    selected_parameters <- unique(unlist(parameters))

    if(!all(selected_parameters %in%
            fit@modelInfo$transparameters_labels)) {
      stop("Unknown parameters.")
    }

  }

  #### Fitted parameter values ####

  fit@modelInfo$control_optimizer$parameters[[1L]] <-
    fit@Optim$parameters
  fit@modelInfo$control_optimizer$transparameters[[1L]] <-
    fit@Optim$transparameters
  fit@modelInfo$control_optimizer$idx_transforms <-
    trans_depends(fit@modelInfo, parameters)

  #### Delta-method covariance ####

  VCOV <- get_vcov(
    control_manifold = fit@modelInfo$control_manifold,
    control_transform = fit@modelInfo$control_transform,
    control_estimator = fit@modelInfo$control_estimator,
    control_optimizer = fit@modelInfo$control_optimizer,
    vcov = v
  )

  trans_labels <- fit@modelInfo$transparameters_labels
  rownames(VCOV$vcov) <- colnames(VCOV$vcov) <- trans_labels

  selected_parameters <- unique(unlist(parameters))
  selected_idx <- match(selected_parameters, trans_labels)

  if(anyNA(selected_idx)) {
    stop("The selected parameters could not be matched to the transformed ",
         "parameter coordinates.")
  }

  VCOV$vcov <- as.matrix(
    VCOV$vcov[selected_idx, selected_idx, drop = FALSE]
  )
  VCOV$vcov <- validate_covariance_matrix(
    VCOV$vcov,
    labels = selected_parameters,
    object_name = "transformed variance-covariance matrix"
  )

  VCOV$se <- standard_errors_from_vcov(
    VCOV$vcov,
    object_name = "transformed variance-covariance matrix"
  )

  # get_vcov() uses x.jacob only as an internal workspace. It is not the
  # cumulative Jacobian returned by jacobian.latent().
  VCOV$jacob <- NULL

  #### Result ####

  return(VCOV)

}
