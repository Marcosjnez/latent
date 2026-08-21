# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 21/08/2026

#### Information covariance for saturated multivariate-normal moments ####

compute_se_lmvnorm <- function(dataList, modelInfo, Optim) {

  control_optimizer <- modelInfo$control_optimizer
  control_optimizer$parameters[[1L]] <- Optim$parameters
  control_optimizer$transparameters[[1L]] <- Optim$transparameters

  labels <- modelInfo$parameters_labels

  H <- get_hess(
    modelInfo$control_manifold,
    modelInfo$control_transform,
    modelInfo$control_estimator,
    control_optimizer
  )$h
  H <- validate_covariance_matrix(
    H,
    labels = labels,
    object_name = "saturated-moment Hessian"
  )

  VCOV <- invert_information_matrix(
    H,
    labels = labels,
    object_name = "saturated-moment Hessian"
  )

  se <- standard_errors_from_vcov(
    VCOV,
    object_name = "saturated-moment variance-covariance matrix"
  )

  #### Result ####

  result <- list(
    H = H,
    VCOV = VCOV,
    se = se,
    type = "information"
  )

  return(result)

}
