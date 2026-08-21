# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 21/08/2026
#'
#' Jacobian Matrix for Latent Models
#'
#' Compute the cumulative Jacobian from the freely estimated parameters to
#' selected transformed parameters.
#'
#' @param fit A fitted object inheriting from class \code{"latent"}.
#' @param parameters Optional parameter specification identifying the transformed
#'   parameters whose derivatives should be returned.
#'
#' @return A numeric matrix whose rows correspond to selected transformed
#'   parameters and whose columns correspond to the freely estimated parameters.
#'
#' @details
#' The cumulative Jacobian is computed only when this method is called. Local
#' transformation Jacobians are composed in dependency order, and the matrix is
#' stored relative to the freely estimated parameter vector rather than the
#' complete transformed-parameter vector.
#'
#' @method jacobian latent
#' @export
jacobian.latent <- function(fit, parameters = NULL) {

  #### Check inputs ####

  if(!inherits(fit, "latent")) {
    stop("fit must inherit from class 'latent'.")
  }

  if(length(fit@Optim$transparameters) == 0L) {
    stop("The latent object has not been fitted.")
  }

  if(is.null(parameters)) {

    parameters <- fit@modelInfo$trans[names(fit@modelInfo$param)]

  } else {

    selected_parameters <- unique(unlist(parameters))

    if(!all(selected_parameters %in%
            fit@modelInfo$transparameters_labels)) {
      stop("Unknown parameters.")
    }

  }

  selected_parameters <- unique(unlist(parameters))

  #### Fitted parameter values ####

  fit@modelInfo$control_optimizer$parameters[[1L]] <-
    fit@Optim$parameters
  fit@modelInfo$control_optimizer$transparameters[[1L]] <-
    fit@Optim$transparameters
  fit@modelInfo$control_optimizer$idx_transforms <-
    trans_depends(fit@modelInfo, parameters)

  #### Cumulative Jacobian ####

  result <- get_jacob(
    fit@modelInfo$control_manifold,
    fit@modelInfo$control_transform,
    fit@modelInfo$control_estimator,
    fit@modelInfo$control_optimizer
  )$jacob

  trans_labels <- fit@modelInfo$transparameters_labels
  free_labels <- fit@modelInfo$parameters_labels

  if(nrow(result) != length(trans_labels) ||
     ncol(result) != length(free_labels)) {
    stop("The cumulative Jacobian returned by C++ has incompatible dimensions.")
  }

  rownames(result) <- trans_labels
  colnames(result) <- free_labels

  selected_idx <- match(selected_parameters, trans_labels)

  if(anyNA(selected_idx)) {
    stop("The selected parameters could not be matched to the Jacobian rows.")
  }

  result <- as.matrix(
    result[selected_idx, , drop = FALSE]
  )
  rownames(result) <- selected_parameters
  colnames(result) <- free_labels

  #### Result ####

  return(result)

}
