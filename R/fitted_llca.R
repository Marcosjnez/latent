# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 11/07/2026
#'
#' Fitted Latent Class Membership Probabilities
#'
#' Extract fitted prior latent class membership probabilities for the subjects
#' used to estimate a latent class model.
#'
#' @description
#' \code{fitted.llca()} returns the fitted latent class membership probabilities
#' for the observations retained in an \code{"llca"} model. These probabilities
#' depend on the class-membership regression model and its covariates, but not on
#' the observed indicator responses.
#'
#' @usage
#' \method{fitted}{llca}(object, ...)
#'
#' @param object A fitted object of class \code{"llca"}.
#' @param ... Additional arguments. They are currently ignored.
#'
#' @details
#' The returned values are the prior class-membership probabilities
#' \eqn{P(C_i = c \mid x_i)}. They should not be confused with posterior
#' probabilities, which additionally condition on the observed indicator
#' responses and can be obtained with
#' \code{latInspect(object, what = "posterior")}.
#'
#' Class probabilities are stored internally for the unique response and
#' covariate patterns. This method expands them back to the subject level using
#' the mapping stored in the fitted object. Consequently, the number and order
#' of rows correspond to the observations retained in \code{object@dataList$data}.
#' Observations removed before estimation because of missing covariates are not
#' included.
#'
#' If the model has no covariates, every subject receives the same vector of
#' estimated class proportions.
#'
#' @return
#' A numeric matrix with one row per retained subject and one column per latent
#' class. Row names correspond to the retained observations and column names
#' identify the latent classes.
#'
#' @seealso
#' \code{\link{predict.llca}}, \code{\link{latInspect}}
#'
#' @references
#' None yet.
#'
#' @method fitted llca
#' @export
fitted.llca <- function(object, ...) {

  if(!inherits(object, "llca")) {
    stop("object must inherit from class 'llca'.")
  }

  if(length(object@transformed_pars) == 0L ||
     is.null(object@transformed_pars$class)) {
    stop("The llca object must be fitted before fitted values can be extracted.")
  }

  class_patterns <- object@transformed_pars$class
  short2full <- object@dataList$short2full

  if(is.null(short2full)) {
    short2full <- seq_len(nrow(class_patterns))
  }

  if(any(short2full < 1L) || any(short2full > nrow(class_patterns))) {
    stop("The stored pattern-to-subject mapping is invalid.")
  }

  result <- class_patterns[short2full, , drop = FALSE]
  rownames(result) <- rownames(object@dataList$data)

  if(is.null(colnames(result))) {
    colnames(result) <- paste("Class", seq_len(ncol(result)), sep = "")
  }

  #### Result ####

  return(result)

}

#' Fitted Probabilities for Lists of Latent Class Models
#'
#' Extract fitted prior class-membership probabilities from an
#' \code{"llcalist"} object.
#'
#' @description
#' \code{fitted.llcalist()} applies \code{fitted()} to latent class models stored
#' in an \code{"llcalist"}. For Bakk--Kuha or ML adjustment results containing
#' elements named \code{measurement} and \code{structural}, only the structural
#' model is returned because it contains the class-membership regression.
#'
#' @usage
#' \method{fitted}{llcalist}(object, ...)
#'
#' @param object An object of class \code{"llcalist"} containing fitted
#'   \code{"llca"} objects.
#' @param ... Additional arguments passed to \code{fitted.llca()}.
#'
#' @return
#' For an adjusted model with a named \code{structural} component, a numeric
#' matrix containing the fitted class-membership probabilities from that model.
#' Otherwise, a list containing one such matrix for each fitted \code{"llca"}
#' model. The latter result has class \code{"fitted.llcalist"}.
#'
#' @seealso
#' \code{\link{fitted.llca}}, \code{\link{predict.llcalist}}
#'
#' @references
#' None yet.
#'
#' @method fitted llcalist
#' @export
fitted.llcalist <- function(object, ...) {

  if(!inherits(object, "llcalist")) {
    stop("object must inherit from class 'llcalist'.")
  }

  if(length(object) == 0L) {
    stop("object does not contain any fitted models.")
  }

  if("structural" %in% names(object) &&
     inherits(object[["structural"]], "llca")) {

    result <- fitted.llca(object[["structural"]], ...)

    #### Result ####

    return(result)

  }

  is_llca <- vapply(object, inherits, logical(1L), what = "llca")

  if(!all(is_llca)) {
    stop("All elements of object must inherit from class 'llca'.")
  }

  result <- lapply(object, fitted.llca, ...)

  if(is.null(names(result)) || any(names(result) == "")) {
    names(result) <- vapply(object, FUN = function(x) {
      paste("nclasses=", ncol(x@transformed_pars$class), sep = "")
    }, FUN.VALUE = character(1L))
  }

  class(result) <- "fitted.llcalist"

  #### Result ####

  return(result)

}
