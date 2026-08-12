#' @export
jacobian.latent <- function(fit, parameters = NULL) {

  if(is.null(parameters)) {
    parameters <- fit@modelInfo$trans[names(fit@modelInfo$param)]
  } else if(!any(unlist(parameters) %in% fit@modelInfo$transparameters_labels)) {
    stop("Unknown parameters.")
  }

  fit@modelInfo$control_optimizer$parameters[[1]] <- fit@Optim$parameters
  fit@modelInfo$control_optimizer$transparameters[[1]] <- fit@Optim$transparameters
  fit@modelInfo$control_optimizer$idx_transforms <- trans_depends(fit, parameters)

  jacob <- get_jacob(fit@modelInfo$control_manifold,
                     fit@modelInfo$control_transform,
                     fit@modelInfo$control_estimator,
                     fit@modelInfo$control_optimizer)$jacob
  rownames(jacob) <- colnames(jacob) <- fit@modelInfo$transparameters_labels

  selected_parameters <- unique(unlist(parameters))
  selected_idx <- match(selected_parameters, fit@modelInfo$transparameters_labels)
  jacob <- jacob[selected_idx, selected_idx]

  return(jacob)

}
