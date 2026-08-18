# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 18/08/2026
#'
#' Robust Variance-Covariance Matrix for Latent Models
#'
#' Default robust method for latent objects without a class-specific robust
#' estimator. The information covariance matrix is returned.
#'
#' @param fit A fitted object inheriting from class \code{"latent"}.
#'
#' @return The result of \code{information(fit)}.
#'
#' @method robust latent
#' @export
robust.latent <- function(fit) {

  result <- information(fit)

  #### Result ####

  return(result)

}
