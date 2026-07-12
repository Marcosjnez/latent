# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 12/07/2026

standard_se <- function(fit) {

  # Compute the variance-covariance matrix of the parameters:

  control_manifold <- fit@modelInfo$control_manifold
  control_transform <- fit@modelInfo$control_transform
  control_estimator <- fit@modelInfo$control_estimator
  control_optimizer <- fit@modelInfo$control_optimizer
  control_optimizer$parameters[[1]] <- fit@Optim$parameters
  control_optimizer$transparameters[[1]] <- fit@Optim$transparameters

  H <- get_hess(control_manifold = control_manifold,
                control_transform = control_transform,
                control_estimator = control_estimator,
                control_optimizer = control_optimizer)$h

  # Get the variance-covariance matrix between the parameters:
  result <- get_vcov(control_manifold = control_manifold,
                     control_transform = control_transform,
                     control_estimator = control_estimator,
                     control_optimizer = control_optimizer,
                     H = H)

  # Name the parameters:
  colnames(H) <- rownames(H) <- fit@modelInfo$parameters_labels
  result$H <- H

  result$se <- as.vector(result$se)
  names(result$se) <- colnames(result$vcov) <- rownames(result$vcov) <-
    fit@modelInfo$parameters_labels

  # Return:
  return(result)

}
