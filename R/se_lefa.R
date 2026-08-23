# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 23/08/2026
#'
#' Standard Errors for Exploratory Factor Analysis
#'
#' @param fit A fitted object inheriting from class \code{"lefa"}.
#' @param type Optional compatibility argument.
#' @param parameters Optional parameter specification. By default, standard
#'   errors are returned for the rotated loadings, factor covariance matrices,
#'   and factor means.
#' @param digits Non-negative integer used to format parameter tables. If
#'   \code{NULL}, values are not rounded.
#' @param ... Additional arguments reserved for future methods.
#'
#' @return Standard errors propagated from the fitted unrotated \code{lcfa}
#'   model to the rotated EFA parameters.
#'
#' @method se lefa
#' @export
se.lefa <- function(fit, type = NULL, parameters = NULL,
                    digits = 4L, ...) {

  if(!is.null(type)) {
    type <- match.arg(tolower(type),
                      c("standard", "information", "robust"))
  }

  if(is.null(parameters)) {
    parameters <- lefa_rotated_parameters(fit)
  }

  result <- se_lrotate_multistep(fit = fit,
                                 parameters = parameters,
                                 digits = digits)

  #### Result ####

  return(result)

}

#'
#' @rdname information.multistep
#' @method information lefa
#' @export
information.lefa <- function(fit) {

  SE <- se_lrotate_multistep(
    fit = fit,
    parameters = fit@modelInfo$parameters_labels,
    digits = NULL
  )

  VCOV <- SE$free_VCOV
  labels <- fit@modelInfo$parameters_labels

  rownames(VCOV) <- colnames(VCOV) <- labels

  se <- sqrt(diag(VCOV))
  names(se) <- labels

  #### Result ####

  result <- list(H = SE$H,
                 B = SE$B,
                 C = SE$C,
                 VCOV = VCOV,
                 se = se,
                 joint_vcov = SE$joint_vcov)

  return(result)

}

#'
#' @rdname robust.multistep
#' @method robust lefa
#' @export
robust.lefa <- function(fit) {

  result <- information.lefa(fit)

  #### Result ####

  return(result)

}
