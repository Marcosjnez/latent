# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 15/08/2026
#'
#' Hessian and Variance-Covariance Matrix matrix for Confirmatory Factor Models
#' using the Information method
#'
#' @export
information.lcfa <- function(fit) {

  return(robust.lcfa(fit))

}
