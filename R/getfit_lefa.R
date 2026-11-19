# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 23/08/2026
#'
#' Fit Indices for Exploratory Factor Analysis
#'
#' @param model A fitted object inheriting from class \code{"lefa"}.
#' @param digits Number of decimal places used in the returned values.
#' @param fit_matrix Optional precomputed internal fit matrix passed to
#'   \code{getfit.lcfa()}.
#'
#' @return A named numeric vector containing model dimensions and fit
#'   statistics from the unrotated model. Rotation does not change model fit.
#'
#' @param x Fitted model or collection of fitted models passed to the method.
#' @param ... Additional arguments to the generic.
#' @method getfit lefa
#' @export
getfit.lefa <- function(x, digits = 3L, fit_matrix = NULL, ..., model = x) {
  if(!missing(model) && !missing(x)) return(.match_legacy_model("model"))

  #### Check inputs ####

  if(!inherits(model, "lefa")) {
    stop("model must inherit from class 'lefa'.")
  }

  source <- lefa_source_lcfa(model)
  result <- getfit.lcfa(model = source,
                        digits = digits,
                        fit_matrix = fit_matrix)
  class(result) <- c("getfit.lefa", "numeric")

  #### Result ####

  return(result)

}
