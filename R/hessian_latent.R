#' @export
hessian.latent <- function(fit) {

  fit@modelInfo$control_optimizer$parameters[[1]] <- fit@Optim$parameters
  fit@modelInfo$control_optimizer$transparameters[[1]] <- fit@Optim$transparameters

  H <- get_hess(fit@modelInfo$control_manifold,
                fit@modelInfo$control_transform,
                fit@modelInfo$control_estimator,
                fit@modelInfo$control_optimizer)$h

  rownames(H) <- colnames(H) <- fit@modelInfo$parameters_labels

  return(H)

}
