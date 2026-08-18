# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 18/08/2026
#'
#' Robust Covariance for Multistep Models
#'
#' The covariance estimators of preceding sample statistics are selected when
#' the multistep model is fitted. This method therefore returns the same
#' multistep propagation as \code{information.multistep()}.
#'
#' @param fit A fitted object inheriting from class \code{"multistep"}.
#'
#' @return The result of \code{information.multistep(fit)}.
#'
#' @method robust multistep
#' @export
robust.multistep <- function(fit) {

  result <- information.multistep(fit)

  #### Result ####

  return(result)

}

# Explicit concrete-class method for deterministic S3 dispatch with S4 multiple
# inheritance.
#'
#' @rdname robust.multistep
#' @method robust multistep_lcfa
#' @export
robust.multistep_lcfa <- function(fit) {

  result <- robust.multistep(fit)

  #### Result ####

  return(result)

}
