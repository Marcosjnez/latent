# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 12/08/2026
#'
#' Hessian Matrix for Latent Models
#'
#' Compute the Hessian matrix of the objective function evaluated at the
#' parameter estimates of a fitted \code{"latent"} object.
#'
#' @description
#' The Hessian is computed with respect to the freely estimated parameters on
#' their optimization scale. Transformed parameters are therefore not included
#' as separate rows or columns.
#'
#' @param fit A fitted object inheriting from class \code{"latent"}.
#'
#' @return
#' A symmetric numeric matrix with one row and one column for each freely
#' estimated parameter. Row and column names correspond to
#' \code{fit@modelInfo$parameters_labels}.
#'
#' @details
#' The Hessian is evaluated at the parameter estimates stored in
#' \code{fit@Optim}. For models estimated by an alternative optimization
#' criterion, such as expectation-maximization, the estimator stored after
#' fitting determines the objective whose derivatives are evaluated.
#'
#' The Hessian is used internally by \code{\link{vcov.latent}} to obtain
#' covariance matrices and standard errors.
#'
#' @seealso
#' \code{\link{hessian}}, \code{\link{vcov.latent}},
#' \code{\link{jacobian.latent}}
#'
#' @method hessian latent
#' @export
hessian.latent <- function(fit) {

  fit@modelInfo$control_optimizer$parameters[[1]] <- fit@Optim$parameters
  fit@modelInfo$control_optimizer$transparameters[[1]] <- fit@Optim$transparameters

  H <- get_hess(fit@modelInfo$control_manifold,
                fit@modelInfo$control_transform,
                fit@modelInfo$control_estimator,
                fit@modelInfo$control_optimizer)$h

  rownames(H) <- colnames(H) <- fit@modelInfo$parameters_labels

  return(H)

}
