# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 12/08/2026
#'
#' Jacobian Matrix for Latent Models
#'
#' Compute the Jacobian matrix associated with parameter transformations in a
#' fitted latent variable model.
#'
#' @description
#' \code{jacobian.latent()} evaluates the derivatives of transformed parameters
#' at the fitted parameter estimates. The Jacobian is used by \pkg{latent} to
#' propagate covariance matrices through parameter transformations using the
#' delta method.
#'
#' @param fit A fitted object inheriting from class \code{"latent"}.
#' @param parameters Optional parameter specification identifying the transformed
#'   parameters for which derivatives should be returned. Parameter labels must
#'   occur in \code{fit@modelInfo$transparameters_labels}. If \code{NULL}, the
#'   parameter blocks corresponding to the fitted model parameters are used.
#'
#' @return
#' A numeric matrix containing the Jacobian for the selected parameters. Row and
#' column names correspond to transformed-parameter labels.
#'
#' @details
#' Only transformations required to obtain the selected parameters are
#' evaluated. Dependencies between transformations are identified recursively,
#' so parameters obtained through several successive transformations are
#' supported.
#'
#' The same Jacobian machinery is used by \code{\link{vcov.latent}} to apply
#' the delta method to transformed parameters.
#'
#' @seealso
#' \code{\link{jacobian}}, \code{\link{hessian.latent}},
#' \code{\link{vcov.latent}}, \code{\link{constraints_derivs.latent}}
#'
#' @method jacobian latent
#' @export
jacobian.latent <- function(fit, parameters = NULL) {

  if(is.null(parameters)) {
    parameters <- fit@modelInfo$trans[names(fit@modelInfo$param)]
  } else if(!any(unlist(parameters) %in% fit@modelInfo$transparameters_labels)) {
    stop("Unknown parameters.")
  }

  fit@modelInfo$control_optimizer$parameters[[1]] <- fit@Optim$parameters
  fit@modelInfo$control_optimizer$transparameters[[1]] <- fit@Optim$transparameters
  fit@modelInfo$control_optimizer$idx_transforms <- trans_depends(fit, parameters)

  jacob <- get_jacob(fit@modelInfo$control_manifold,
                     fit@modelInfo$control_transform,
                     fit@modelInfo$control_estimator,
                     fit@modelInfo$control_optimizer)$jacob
  rownames(jacob) <- colnames(jacob) <- fit@modelInfo$transparameters_labels

  selected_parameters <- unique(unlist(parameters))
  selected_idx <- match(selected_parameters, fit@modelInfo$transparameters_labels)
  jacob <- jacob[selected_idx, selected_idx, drop = FALSE]

  return(jacob)

}
