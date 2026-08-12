#' Jacobian Matrix
#'
#' Compute Jacobian matrices for transformed model parameters.
#'
#' @description
#' \code{jacobian()} is a generic function for obtaining the derivatives
#' associated with parameter transformations in fitted latent variable models.
#'
#' @param x A fitted model object.
#' @param ... Additional arguments passed to methods.
#'
#' @return
#' A numeric matrix containing derivatives of transformed model parameters.
#'
#' @seealso
#' \code{\link{hessian}}, \code{\link{vcov}},
#' \code{\link{constraints_derivs}}
#'
#' @export
jacobian <- function(x, ...) {
  UseMethod("jacobian")
}
