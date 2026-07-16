# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 16/07/2026

standard_se <- function(fit, parameters = NULL) {

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

  # Get the variance-covariance matrix between the parameters of interest:
  if(is.null(parameters)) {
    parameters <- fit@modelInfo$trans[names(fit@parameters)]
  } else if(!any(unlist(parameters) %in% fit@modelInfo$transparameters_labels)) {
    stop("Unknown parameters.")
  }

  control_optimizer$idx_transforms <- trans_depends(fit, parameters)
  result <- get_vcov(control_manifold = control_manifold,
                     control_transform = control_transform,
                     control_estimator = control_estimator,
                     control_optimizer = control_optimizer,
                     H = H)

  # Name the parameters in the Hessian matrix:
  colnames(H) <- rownames(H) <- fit@modelInfo$parameters_labels
  result$H <- H
  # The Hessian only contains the untransformed parameters: the parameters
  # tha are directly optimized.

  # vcov and se contains all the parameters but their values are mostly 0 if
  # you don't request a custom parameters.
  selected_parameters <- unique(unlist(parameters))
  idx <- match(selected_parameters, fit@modelInfo$transparameters_labels)
  result$se <- as.vector(result$se[idx])
  result$vcov <- result$vcov[idx, idx, drop = FALSE]
  rownames(result$vcov) <- colnames(result$vcov) <- names(result$se) <-
    selected_parameters

  # Return:
  return(result)

}

trans_depends <- function(fit, trans_subset) {

  # We want to compute the standard errors of a subset of parameters. Some
  # parameters may be transformed parameters. In this scenario, we first must
  # identify what transformations do the parameters depend on. Because we will
  # use them to compute the VCOV matrix using the delta method:
  # vcov(output parameters) = jacobian %*% vcov(input parameters) %*% t(jacobian)
  # We only want to extract the subset of transformations that produce the
  # transformed parameters we are interested in (trans_subset).
  # After this, we will identify all the parameters that play a role as inputs.
  # We already have the VCOV of the untransformed parameters, so this is just
  # computing the jacobians and applying the delta method bottom-up.

  control_transform <- fit@modelInfo$control_transform
  parameters_labels <- fit@modelInfo$parameters_labels
  transparameters_labels <- fit@modelInfo$transparameters_labels

  # Unique parameters in the subset of parameters:
  unq_trans_subset <- unique(c(unlist(trans_subset), parameters_labels))

  # Some transformed parameters may depend on other transformed parameters.
  # Recursively add every transformation input until no additional
  # transformations are required:
  idx_transforms <- integer(0L)

  repeat {

    # Get the positions of the current subset parameters in the full vector:
    idx <- match(unq_trans_subset, transparameters_labels)-1L

    # Select transformations producing parameters currently in the subset:
    candidate_transforms <- which(vapply(control_transform, FUN = \(x) {
      any(unlist(x$indices_out) %in% idx)
    }, FUN.VALUE = logical(1L)))-1L

    # Retain only transformations that have not already been processed:
    new_transforms <- setdiff(candidate_transforms, idx_transforms)

    if(length(new_transforms) == 0L) {
      break
    }

    idx_transforms <- sort(unique(c(idx_transforms, new_transforms)))

    # Add every input parameter used by the newly identified transformations:
    input_idx <- unique(unlist(lapply(control_transform[new_transforms+1L],
                                      FUN = \(x) unlist(x$indices_in)),
                               use.names = FALSE))

    if(length(input_idx) > 0L) {
      input_labels <- transparameters_labels[input_idx+1L]
      unq_trans_subset <- unique(c(unq_trans_subset, input_labels))
    }

  }

  #### Result ####

  return(idx_transforms)

}
