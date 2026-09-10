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
#' @method latInspect lefa
#' @export
latInspect.lefa <- function(fit, what = "est", sort = TRUE) {

  if(!inherits(fit, "lefa")) {
    stop("fit must inherit from class 'lefa'.")
  }

  result <- latInspect.latent(fit, what = what, sort = sort)

  #### Result ####

  return(result)

}
