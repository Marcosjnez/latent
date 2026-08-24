#' Derivatives of Parameter Constraints
#'
#' Compute first and second derivatives associated with constraints on model
#' parameters and transformed parameters.
#'
#' @description
#' \code{constraints_derivs()} is a generic function for obtaining first and
#' second derivatives of parameter constraints introduced by manifolds and
#' model transformations.
#'
#' @param x A fitted model object.
#' @param ... Additional arguments passed to methods.
#'
#' @return
#' A list containing the first constraint derivatives and the block-diagonal
#' matrix of second constraint derivatives.
#'
#' @seealso
#' \code{\link{jacobian}}, \code{\link{hessian}},
#' \code{\link{vcov}}
#'
#' @export
constraints_derivs <- function(x, ...) {
  UseMethod("constraints_derivs")
}
