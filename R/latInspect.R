# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# Modification date: 10/09/2026

#' Inspect Latent Model Results
#'
#' @param x A fitted latent model or collection of latent class models.
#' @param ... Arguments passed to the inspection method, including \code{what}.
#' @param sort Logical. Apply factor ordering/sign orientation or class ordering
#'   to the requested output. Defaults to TRUE. FALSE returns native outputs.
#'
#' @details
#' Sorting never modifies the supplied object, its optimizer coordinates, or
#' its model information. All uncertainty calculations use the native fit.
#' \code{what = "structures"} returns the estimate and template structures;
#' \code{"param"} and \code{"trans"} return the corresponding label templates.
#' \code{"labels"} has the layout of \code{"est"} with original labels.
#' The actual transformed-estimate slot is named \code{transformed_pars}.
#' See \code{sort_latent()} for the ordering and sign conventions.
#' Unrotated requests keep native source coordinates. Standard errors are
#' never negated; inspected covariance matrices receive both axis signs.
#' A shared label with opposite display signs can occur twice in an inspected
#' covariance matrix, without inventing new labels. Previously saved objects
#' using constructor-side sorting should be refitted.
#'
#' @return The requested model component, without altering \code{x}.
#' @export
latInspect <- function(x, ..., sort = TRUE) {

  UseMethod("latInspect")

}

#' @rdname latInspect
#' @param fit A fitted latent model. This method also handles lrotate results,
#'   which inherit from \code{"latent"} or \code{"multistep"}.
#' @param what Character string identifying the requested component.
#' @method latInspect latent
#' @export
latInspect.latent <- function(fit, what = "est", sort = TRUE) {

  if(!inherits(fit, "latent")) stop("fit must inherit from class 'latent'.")
  if(!is.character(what) || length(what) != 1L || is.na(what)) {
    stop("what must be a single character string.")
  }
  if(!is.logical(sort) || length(sort) != 1L || is.na(sort)) {
    stop("sort must be TRUE or FALSE.")
  }
  what <- tolower(what)
  output <- if(sort) sort_latent(fit) else {
    list(parameters = fit@parameters, transformed_pars = fit@transformed_pars,
         param = fit@modelInfo$param, trans = fit@modelInfo$trans)
  }
  data_param <- fit@modelInfo$data_param
  rotation <- !is.null(data_param$X_group)
  source <- if(inherits(fit, "lefa")) lefa_source_lcfa(fit) else {
    if(rotation) fit@extra$efa else NULL
  }
  source_blocks <- character(0L)

  if(rotation) {
    blocks <- unique(c(data_param$lambda_group, data_param$psi_group,
                        data_param$alpha_group))
    if(inherits(source, "lcfa")) {
      source_data <- source@dataList$data_param
      source_blocks <- setdiff(names(latInspect(source, what = "est", sort = FALSE)),
                               c(source_data$lambda_group, source_data$psi_group,
                                 source_data$alpha_group))
    }
    blocks <- c(blocks, source_blocks)
  } else if(inherits(fit, "lcfa")) {
    data_param <- fit@dataList$data_param
    blocks <- unique(c(data_param$lambda_group, data_param$alpha_group,
                        data_param$theta_group, data_param$psi_group,
                        data_param$nu_group, data_param$kappa_group,
                        data_param$delta_group))
  } else {
    blocks <- names(output$parameters)
  }

  if(what == "structures") {
    result <- output
  } else if(what == "param") {
    result <- output$param
  } else if(what == "trans") {
    result <- output$trans
  } else if(what %in% c("transparameters", "transformed_pars")) {
    result <- output$transformed_pars
  } else if(what == "labels") {
    result <- output$trans[intersect(blocks, names(output$trans))]
  } else if(what %in% c("est", "estimates", "parameters", "fixed")) {
    result <- output$transformed_pars[intersect(blocks, names(output$transformed_pars))]
  } else if(what %in% c("se", "standard.errors", "standard_errors", "vcov", "covariance")) {

    # Obtain uncertainty from the original object before changing presentation.
    SE <- fit@Optim$SE
    if(is.null(SE) && rotation && inherits(fit, "multistep")) {
      parameters <- fit@modelInfo$trans[unique(c(data_param$lambda_group,
                                                 data_param$psi_group,
                                                 data_param$alpha_group))]
      SE <- se.multistep(fit, parameters = parameters, digits = NULL)
    } else if(is.null(SE) && inherits(fit, "llca")) {
      SE <- se(fit, digits = NULL)
    }
    covariance <- what %in% c("vcov", "covariance")
    result <- if(covariance) SE$VCOV else SE$se
    if(is.null(result) && covariance) result <- SE$vcov
    if(is.null(result)) stop("No requested standard errors/covariance matrix is available.")

    if(sort && !is.null(attr(output, "order", exact = TRUE))) {
      labels <- if(covariance) rownames(result) else names(result)
      if(is.null(labels)) stop("Standard-error outputs must have parameter labels.")
      templates <- output$trans[unique(c(intersect(blocks, names(output$trans)),
                                         names(output$trans)))]
      signs <- attr(output, "signs", exact = TRUE)
      mapping <- vector("list", length(templates))
      for(i in seq_along(templates)) {
        nm <- names(templates)[i]
        direction <- signs[[nm]]
        if(is.null(direction)) direction <- rep(1, length(templates[[i]]))
        mapping[[i]] <- data.frame(label = c(templates[[i]]), sign = c(direction),
                                   stringsAsFactors = FALSE)
      }
      mapping <- unique(do.call(rbind, mapping))
      mapping <- mapping[mapping$label %in% labels, , drop = FALSE]
      remaining <- setdiff(labels, mapping$label)
      mapping <- rbind(mapping, data.frame(label = remaining, sign = rep(1, length(remaining))))
      if(covariance) {
        result <- result[mapping$label, mapping$label, drop = FALSE]*
          outer(mapping$sign, mapping$sign)
      } else {
        result <- result[mapping$label]
      }
    }

  } else if(rotation && what %in% c("lambda", "loadings", "psi", "alpha", "latent.means",
                                     "latent_means", "x", "rotation", "rotation.matrix",
                                     "rotation_matrix", "xinv")) {
    selected <- switch(what, lambda =, loadings = data_param$lambda_group,
                        psi = data_param$psi_group, alpha =, latent.means =,
                        latent_means = data_param$alpha_group, xinv = data_param$Xinv_group,
                        data_param$X_group)
    result <- output$transformed_pars[intersect(selected, names(output$transformed_pars))]
  } else if(rotation && what %in% c("efa", "unrotated.fit", "unrotated_fit")) {
    if(is.null(source)) stop("No fitted source model was supplied to lrotate().")
    result <- source
  } else if(rotation && what %in% c("unrotated", "unrotated.est", "unrotated_est")) {
    result <- if(inherits(source, "lcfa")) latInspect(source, what = "est", sort = FALSE) else {
      fit@transformed_pars[c(data_param$ulambda_group, data_param$upsi_group,
                             data_param$ualpha_group)]
    }
  } else if(rotation && what %in% c("unrotated.lambda", "unrotated_lambda", "unrotated.psi",
                                     "unrotated_psi", "unrotated.alpha", "unrotated_alpha")) {
    selected <- if(grepl("lambda$", what)) data_param$ulambda_group else {
      if(grepl("psi$", what)) data_param$upsi_group else data_param$ualpha_group
    }
    result <- fit@transformed_pars[selected]
  } else if(rotation && what %in% c("rotation.loss", "rotation_loss")) {
    result <- fit@Optim$f
  } else if(rotation && inherits(source, "lcfa")) {
    result <- latInspect(source, what = what, sort = FALSE)
  } else {
    stop("Unknown request: ", what)
  }

  #### Result ####

  return(result)

}
