#' @export
constraints_derivs.latent <- function(fit, parameters = NULL) {

  if(is.null(parameters)) {
    parameters <- fit@modelInfo$trans[names(fit@modelInfo$param)]
  } else if(!any(unlist(parameters) %in% fit@modelInfo$transparameters_labels)) {
    stop("Unknown parameters.")
  }

  fit@modelInfo$control_optimizer$parameters[[1]] <- fit@Optim$parameters
  fit@modelInfo$control_optimizer$transparameters[[1]] <- fit@Optim$transparameters

  fit@modelInfo$control_optimizer$idx_transforms <- trans_depends(fit, parameters)

  dconstr <- get_dconstr(fit@modelInfo$control_manifold,
                         fit@modelInfo$control_transform,
                         fit@modelInfo$control_estimator,
                         fit@modelInfo$control_optimizer)$dconstr
  rownames(dconstr) <- fit@modelInfo$transparameters_labels
  idx_transforms <- fit@modelInfo$control_optimizer$idx_transforms
  colnames(dconstr) <- sapply(fit@modelInfo$control_transform[idx_transforms+1L],
                              FUN = \(x) x$transform)

  selected_parameters <- unique(unlist(parameters))
  selected_idx <- match(selected_parameters, fit@modelInfo$transparameters_labels)
  dconstr <- dconstr[selected_idx, ]

  return(dconstr)

}
