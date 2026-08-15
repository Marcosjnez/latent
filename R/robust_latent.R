# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 15/08/2026
#'
#' Hessian and Variance-Covariance Matrix matrix for Latent Models using the
#' Information method
#'
#' @export
robust.latent <- function(fit) {

  return(information.latent(fit))

}
