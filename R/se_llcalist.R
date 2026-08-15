# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 15//08/2026
#'
#' Standard errors for a collection of latent class models
#'
#' Applies \code{se()} to every fitted model in an \code{"llcalist"} object.
#'
#' @param model An object of class \code{"llcalist"} containing fitted
#'   \code{"llca"} objects.
#' @param type Character string indicating the standard-error estimator. See
#'   \code{se.llca()}.
#' @param digits Non-negative integer indicating the number of decimal places used
#'   in the formatted tables, or \code{NULL} to avoid rounding.
#' @param ... Additional arguments passed to \code{se.llca()}.
#'
#' @details
#' Existing names are preserved. Consequently, class-enumeration results retain
#' names such as \code{"nclasses=2"}, whereas multiple-step models retain names
#' such as \code{"measurement"} and \code{"structural"}. Unnamed elements are
#' labelled according to their number of latent classes.
#'
#' @return A named list with one standard-error result per fitted model and class
#'   \code{"se.llcalist"}.
#'
#' @examples
#' \dontrun{
#' fits <- lca(data = empathy, nclasses = 2:4,
#'             gaussian = c("ec1", "ec2", "ec3"))
#' se(fits)
#' }
#'
#' @method se llcalist
#' @export
se.llcalist <- function(model, type = "information", parameters = NULL,
                        digits = 4L, ...) {

  if(!inherits(model, "llcalist")) {
    stop("model must inherit from class 'llcalist'.")
  }

  if(length(model) == 0L) {
    stop("model must contain at least one fitted llca object.")
  }

  nmodels <- length(model)
  result <- vector("list", length = nmodels)

  for(i in seq_len(nmodels)) {
    result[[i]] <- se(model[[i]], type = type, parameters = parameters,
                      digits = digits, ...)
  }

  result_names <- names(model)
  if(is.null(result_names)) {
    result_names <- rep("", nmodels)
  }

  unnamed <- is.na(result_names) | result_names == ""
  if(any(unnamed)) {
    nclasses <- vapply(model[unnamed], FUN = function(x) {
      ncol(x@modelInfo$trans$class)
    }, FUN.VALUE = integer(1L))
    result_names[unnamed] <- paste0("nclasses=", nclasses)
  }

  names(result) <- make.unique(result_names)
  class(result) <- "se.llcalist"

  #### Result ####

  return(result)

}
