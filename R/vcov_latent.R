# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 17/08/2026
#'
#' Variance-Covariance Matrix for Latent Objects
#'
#' @param fit A fitted object inheriting from class \code{"latent"}.
#' @param v Variance-covariance matrix of the freely estimated, untransformed
#'   parameters.
#' @param parameters Optional parameter specification identifying the parameters
#'   or transformed parameters to return.
#'
#' @return A list containing the selected variance-covariance matrix, standard
#'   errors, and transformation Jacobian.
#'
#' @method vcov latent
#' @export
vcov.latent <- function(fit, v, parameters = NULL) {

  #### Check the untransformed covariance matrix ####

  if(!is.matrix(v)) {
    v <- as.matrix(v)
  }

  labels <- fit@modelInfo$parameters_labels
  nparam <- length(labels)

  if(nrow(v) != nparam ||
     ncol(v) != nparam ||
     !isSymmetric(v)) {
    stop("v should be a square, symmetric matrix of dimensions ",
         nparam, "x", nparam)
  }

  if(!is.null(rownames(v)) && !is.null(colnames(v)) &&
     all(labels %in% rownames(v)) &&
     all(labels %in% colnames(v))) {
    v <- v[labels, labels, drop = FALSE]
  }

  rownames(v) <- colnames(v) <- labels

  #### Parameters ####

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

  VCOV <- get_vcov(control_manifold = fit@modelInfo$control_manifold,
                   control_transform = fit@modelInfo$control_transform,
                   control_estimator = fit@modelInfo$control_estimator,
                   control_optimizer = fit@modelInfo$control_optimizer,
                   vcov = v)

  rownames(VCOV$vcov) <- colnames(VCOV$vcov) <- names(VCOV$se) <-
    fit@modelInfo$transparameters_labels
  rownames(VCOV$jacob) <- colnames(VCOV$jacob) <-
    fit@modelInfo$transparameters_labels

  selected_parameters <- unique(unlist(parameters))
  selected_idx <- match(selected_parameters,
                        fit@modelInfo$transparameters_labels)

  VCOV$vcov <- as.matrix(VCOV$vcov[selected_idx, selected_idx, drop = FALSE])
  VCOV$vcov <- (VCOV$vcov+t(VCOV$vcov))/2
  VCOV$se <- sqrt(diag(VCOV$vcov))
  names(VCOV$se) <- selected_parameters
  VCOV$jacob <- VCOV$jacob[selected_idx, selected_idx, drop = FALSE]

  #### Result ####

  return(VCOV)

}
