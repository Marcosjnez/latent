# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 15/08/2026
#'
#' Hessian and Variance-Covariance Matrix matrix for Latent Models using the
#' Information method
#'
#' @export
information.latent <- function(fit) {

  fit@modelInfo$control_optimizer$parameters[[1]] <- fit@Optim$parameters
  fit@modelInfo$control_optimizer$transparameters[[1]] <- fit@Optim$transparameters

  #### Compute the Hessian ####

  H <- get_hess(fit@modelInfo$control_manifold,
                fit@modelInfo$control_transform,
                fit@modelInfo$control_estimator,
                fit@modelInfo$control_optimizer)$h
  rownames(H) <- colnames(H) <- fit@modelInfo$parameters_labels
  VCOV <- solve(H)

  #### Return ####

  result <- list(H = H, VCOV = VCOV, se = sqrt(diag(VCOV)))

  return(result)

}
