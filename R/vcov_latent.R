# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 16/08/2026
#'
#' Variance-Covariance Matrix for Latent Objects
#'
#' @export
vcov.latent <- function(fit, v, parameters = NULL) {

  if(!identical(nrow(v), fit@modelInfo$nparam) ||
     !identical(ncol(v), fit@modelInfo$nparam) ||
     !isSymmetric.matrix(v)) {
    stop("v should be an square, symmetric matrix of dimensions ",
         fit@modelInfo$nparam, "x", fit@modelInfo$nparam)
  }

  if(is.null(parameters)) {
    parameters <- fit@modelInfo$parameters_labels
  } else if(!any(unlist(parameters) %in% fit@modelInfo$transparameters_labels)) {
    stop("Unknown parameters.")
  }

  fit@modelInfo$control_optimizer$parameters[[1]] <- fit@Optim$parameters
  fit@modelInfo$control_optimizer$transparameters[[1]] <- fit@Optim$transparameters
  fit@modelInfo$control_optimizer$idx_transforms <-
    trans_depends(fit@modelInfo, parameters)

  VCOV <- get_vcov(fit@modelInfo$control_manifold,
                   fit@modelInfo$control_transform,
                   fit@modelInfo$control_estimator,
                   fit@modelInfo$control_optimizer,
                   vcov = v)

  rownames(VCOV$vcov) <- colnames(VCOV$vcov) <- names(VCOV$se) <-
    fit@modelInfo$transparameters_labels

  selected_parameters <- unique(unlist(parameters))
  selected_idx <- match(selected_parameters,
                        fit@modelInfo$transparameters_labels)

  VCOV$vcov <- VCOV$vcov[selected_idx, selected_idx, drop = FALSE]
  # VCOV$vcov <- (VCOV$vcov + t(VCOV$vcov)) / 2
  VCOV$se <- VCOV$se[selected_idx]
  VCOV$jacob <- VCOV$jacob[selected_idx, selected_idx, drop = FALSE]

  return(VCOV)

}
