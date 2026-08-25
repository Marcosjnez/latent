# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 25/08/2026
#'
#' Jacobian Matrix for Latent Models
#'
#' Compute the dependency Jacobian among transformed parameters.
#'
#' @param fit A fitted object inheriting from class \code{"latent"}.
#' @param parameters Optional parameter specification identifying the transformed
#'   parameters whose Jacobian submatrix should be returned. If \code{NULL}, the
#'   complete transformed-parameter Jacobian is returned.
#'
#' @return A sparse square matrix. Rows and columns correspond to the selected
#'   transformed parameters, or to all transformed parameters when
#'   \code{parameters = NULL}.
#'
#' @details
#' The complete Jacobian has one row and one column for every transformed
#' parameter. Its diagonal is the identity, while off-diagonal entries describe
#' direct and transitive dependencies induced by the sequence of parameter
#' transformations.
#'
#' This square dependency Jacobian is intended for inspecting relationships
#' among transformed parameters. It differs from the conventional delta-method
#' Jacobian, whose columns contain only freely estimated parameters.
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

  trans_labels <- fit@modelInfo$transparameters_labels

  if(is.null(parameters)) {

    selected_parameters <- trans_labels

  } else {

    selected_parameters <- unique(unlist(parameters))

    if(!all(selected_parameters %in% trans_labels)) {
      stop("Unknown parameters.")
    }

  }

  #### Fitted parameter values ####

  fit@modelInfo$control_optimizer$parameters[[1L]] <-
    fit@Optim$parameters
  fit@modelInfo$control_optimizer$transparameters[[1L]] <-
    fit@Optim$transparameters
  fit@modelInfo$control_optimizer$idx_transforms <-
    trans_depends(fit@modelInfo, selected_parameters)

  #### Full dependency Jacobian ####

  result <- get_jacob(
    fit@modelInfo$control_manifold,
    fit@modelInfo$control_transform,
    fit@modelInfo$control_estimator,
    fit@modelInfo$control_optimizer
  )$jacob

  if(nrow(result) != length(trans_labels) ||
     ncol(result) != length(trans_labels)) {
    stop("The dependency Jacobian returned by C++ has incompatible dimensions.")
  }

  rownames(result) <- colnames(result) <- trans_labels

  #### Selected submatrix ####

  selected_idx <- match(selected_parameters, trans_labels)

  if(anyNA(selected_idx)) {
    stop("The selected parameters could not be matched to the Jacobian.")
  }

  result <- result[selected_idx, selected_idx, drop = FALSE]
  rownames(result) <- colnames(result) <- selected_parameters

  #### Result ####

  return(result)

}
