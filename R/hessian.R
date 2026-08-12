#' Hessian Matrix
#'
#' Compute the Hessian matrix of a fitted latent variable model.
#'
#' @description
#' \code{hessian()} is a generic function for obtaining the matrix of
#' second-order derivatives of the objective function with respect to the
#' freely estimated model parameters.
#'
#' @param x A fitted model object.
#' @param ... Additional arguments passed to methods.
#'
#' @return
#' A numeric matrix containing the Hessian of the objective function with
#' respect to the freely estimated parameters.
#'
#' @seealso
#' \code{\link{vcov}}, \code{\link{jacobian}},
#' \code{\link{constraints_derivs}}
#'
#' @export
hessian <- function(x, ...) {
  UseMethod("hessian")
}
