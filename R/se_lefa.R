# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 23/08/2026
#'
#' Standard Errors for Exploratory Factor Analysis
#'
#' @param fit A fitted object inheriting from class \code{"lefa"}.
#' @param type Optional compatibility argument passed to
#'   \code{se.multistep()}.
#' @param parameters Optional parameter specification. By default, standard
#'   errors are returned for the rotated loadings, factor covariance matrices,
#'   and factor means.
#' @param digits Non-negative integer used to format parameter tables. If
#'   \code{NULL}, values are not rounded.
#' @param ... Additional arguments passed to \code{se.multistep()}.
#'
#' @return The result of \code{se.multistep()}.
#'
#' @method se lefa
#' @export
se.lefa <- function(fit, type = NULL, parameters = NULL,
                    digits = 4L, ...) {

  if(is.null(parameters)) {
    parameters <- lefa_rotated_parameters(fit)
  }

  result <- se.multistep(fit = fit,
                         type = type,
                         parameters = parameters,
                         digits = digits,
                         ...)

  #### Result ####

  return(result)

}

#'
#' @rdname information.multistep
#' @method information lefa
#' @export
information.lefa <- function(fit) {

  result <- information.multistep(fit)

  #### Result ####

  return(result)

}

#'
#' @rdname robust.multistep
#' @method robust lefa
#' @export
robust.lefa <- function(fit) {

  result <- robust.multistep(fit)

  #### Result ####

  return(result)

}
