# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 15/08/2026
#'
#' Constraint Derivatives for Latent Models
#'
#' Compute derivatives of transformation-induced constraints for selected
#' parameters of a fitted latent variable model.
#'
#' @description
#' Some parameter transformations imply constraints on their outputs. For
#' example, probabilities produced by a softmax transformation sum to one.
#' \code{constraints_derivs.latent()} extracts the derivatives associated with
#' such constraints for selected transformed parameters.
#'
#' @param fit A fitted object inheriting from class \code{"latent"}.
#' @param parameters Optional parameter specification identifying the transformed
#'   parameters for which constraint derivatives should be returned. Parameter
#'   labels must occur in \code{fit@modelInfo$transparameters_labels}. If
#'   \code{NULL}, the parameter blocks corresponding to the fitted model
#'   parameters are used.
#'
#' @return
#' A numeric matrix. Rows correspond to the selected transformed parameters and
#' columns correspond to transformation-induced constraints.
#'
#' @details
#' Only transformations on which the requested parameters depend are evaluated.
#' Transformations that do not impose an explicit constraint do not contribute
#' a constraint column.
#'
#' @seealso
#' \code{\link{constraints_derivs}}, \code{\link{jacobian.latent}},
#' \code{\link{vcov.latent}}
#'
#' @method constraints_derivs latent
#' @export
constraints_derivs.latent <- function(fit, parameters = NULL) {

  if(is.null(parameters)) {
    parameters <- fit@modelInfo$parameters_labels
  } else if(!any(unlist(parameters) %in% fit@modelInfo$transparameters_labels)) {
    stop("Unknown parameters.")
  }

  fit@modelInfo$control_optimizer$parameters[[1]] <- fit@Optim$parameters
  fit@modelInfo$control_optimizer$transparameters[[1]] <- fit@Optim$transparameters

  fit@modelInfo$control_optimizer$idx_transforms <- #c(0, 1)
    trans_depends(fit@modelInfo, parameters)

  dconstr <- get_dconstr(fit@modelInfo$control_manifold,
                         fit@modelInfo$control_transform,
                         fit@modelInfo$control_estimator,
                         fit@modelInfo$control_optimizer)$dconstr
  rownames(dconstr) <- fit@modelInfo$transparameters_labels
  idx_transforms <- fit@modelInfo$control_optimizer$idx_transforms

  # Not every transformation contributes with just one column of constraints, so
  # don't put names to columns (do not run this code):
  # colnames(dconstr) <- sapply(fit@modelInfo$control_transform[idx_transforms+1L],
  #                             FUN = \(x) x$transform)

  selected_parameters <- unique(unlist(parameters))
  selected_idx <- match(selected_parameters, fit@modelInfo$transparameters_labels)
  dconstr <- dconstr[selected_idx, , drop = FALSE]

  return(dconstr)

}
