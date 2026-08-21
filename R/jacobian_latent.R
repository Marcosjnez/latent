# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 21/08/2026
#'
#' Jacobian Matrix for Latent Models
#'
#' Compute the cumulative Jacobian from the freely estimated parameters to
#' selected transformed parameters.
#'
#' @param fit A fitted object inheriting from class \code{"latent"}.
#' @param parameters Optional parameter specification identifying the transformed
#'   parameters whose derivatives should be returned.
#'
#' @return A numeric matrix whose rows correspond to selected transformed
#'   parameters and whose columns correspond to the freely estimated parameters.
#'
#' @details
#' Only transformations required by the selected outputs are evaluated. Local
#' Jacobians are composed in dependency order, so chains of transformations are
#' represented correctly.
#'
#' @method jacobian latent
#' @export
jacobian.latent <- function(fit, parameters = NULL) {

  if(is.null(parameters)) {

    parameters <- fit@modelInfo$trans[names(fit@modelInfo$param)]

  } else {

    selected_parameters <- unique(unlist(parameters))

    if(!all(selected_parameters %in%
            fit@modelInfo$transparameters_labels)) {
      stop("Unknown parameters.")
    }

  }

  fit@modelInfo$control_optimizer$parameters[[1L]] <-
    fit@Optim$parameters
  fit@modelInfo$control_optimizer$transparameters[[1L]] <-
    fit@Optim$transparameters
  fit@modelInfo$control_optimizer$idx_transforms <-
    trans_depends(fit@modelInfo, parameters)

  jacob <- get_jacob(fit@modelInfo$control_manifold,
                     fit@modelInfo$control_transform,
                     fit@modelInfo$control_estimator,
                     fit@modelInfo$control_optimizer)$jacob

  trans_labels <- fit@modelInfo$transparameters_labels
  free_labels <- fit@modelInfo$parameters_labels
  rownames(jacob) <- colnames(jacob) <- trans_labels

  selected_parameters <- unique(unlist(parameters))
  selected_idx <- match(selected_parameters, trans_labels)
  free_idx <- match(free_labels, trans_labels)

  if(anyNA(selected_idx) || anyNA(free_idx)) {
    stop("The parameter labels could not be matched to the transformed ",
         "parameter coordinates.")
  }

  result <- as.matrix(
    jacob[selected_idx, free_idx, drop = FALSE]
  )
  rownames(result) <- selected_parameters
  colnames(result) <- free_labels

  #### Result ####

  return(result)

}
