# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 10/09/2026
#'
#' Inspect Fitted Exploratory Factor Analysis Objects
#'
#' @param fit A fitted object inheriting from class \code{"lefa"}.
#' @param what Character string identifying the requested component.
#' @param sort Logical. Sort/orient only the requested output; defaults to TRUE.
#'
#' @return A rotated parameter list, the unrotated \code{lcfa} object, or a
#'   component delegated to \code{latInspect.lcfa()}.
#'
#' @param x Fitted model or collection of fitted models passed to the method.
#' @param ... Additional arguments to the generic.
#' @method latInspect lefa
#' @export
latInspect.lefa <- function(x, what = "est", sort = TRUE, ..., fit = x) {
  if(!missing(fit) && !missing(x)) return(.match_legacy_model("fit"))

  if(!inherits(fit, "lefa")) {
    stop("fit must inherit from class 'lefa'.")
  }

  result <- latInspect.latent(fit, what = what, sort = sort)

  #### Result ####

  return(result)

}
