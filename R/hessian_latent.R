# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 05/09/2026
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
#' @param riemannian Logical. If \code{FALSE}, return the Euclidean Hessian. If
#'   \code{TRUE}, return the ambient constrained inverse operator
#'   \eqn{P=T H_R^{-1}T^\top}, where \eqn{H_R} is the tangent-coordinate
#'   Riemannian Hessian and \eqn{T} is the tangent-space basis.
#'
#' @return
#' A symmetric numeric matrix with one row and one column for each freely
#' estimated parameter. Row and column names correspond to
#' \code{fit@modelInfo$parameters_labels}. When \code{riemannian = FALSE}, the
#' matrix is the Euclidean Hessian. When \code{riemannian = TRUE}, it is the
#' ambient constrained inverse operator \eqn{P}.
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

  # fit@modelInfo$control_optimizer$parameters[[1]] <- fit@Optim$parameters
  # fit@modelInfo$control_optimizer$transparameters[[1]] <- fit@Optim$transparameters

  # H <- get_hess(fit@modelInfo$control_manifold,
  #               fit@modelInfo$control_transform,
  #               fit@modelInfo$control_estimator,
  #               fit@modelInfo$control_optimizer)$h
  H <- get_hess(fit)$h

  rownames(H) <- colnames(H) <- fit@modelInfo$parameters_labels

  return(H)

}
