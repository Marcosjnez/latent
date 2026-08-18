# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 17/08/2026
#'
#' Information Variance-Covariance Matrix for Confirmatory Factor Models
#'
#' @param fit A fitted object of class \code{"lcfa"}.
#'
#' @return A list containing the CFA Hessian and variance-covariance matrix.
#'
#' @method information lcfa
#' @export
information.lcfa <- function(fit) {

  result <- vcov_lcfa(fit)

  #### Result ####

  return(result)

}
