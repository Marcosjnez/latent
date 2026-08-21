# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 21/08/2026
#'
#' Confidence Intervals for Latent Class Models
#'
#' Compute normal-approximation confidence intervals for parameters of a fitted
#' latent class model.
#'
#' @param fit A fitted object inheriting from class \code{"llca"}.
#' @param type Character string selecting the standard-error estimator.
#'   Available options are \code{"information"} and \code{"robust"}.
#'   \code{"standard"} is retained as an alias of \code{"information"}.
#' @param confidence Numeric scalar strictly between zero and one specifying the
#'   confidence level.
#' @param parameters Optional parameter specification identifying the parameters
#'   or transformed parameters for which intervals should be returned.
#' @param digits Non-negative integer indicating the number of decimal places
#'   used in the formatted confidence-interval table.
#' @param ... Additional arguments passed to \code{se()}.
#'
#' @return A list containing formatted and numeric confidence limits and the
#'   standard-error result used to construct them.
#'
#' @method ci llca
#' @export
ci.llca <- function(fit, type = "information", confidence = 0.95,
                    parameters = NULL, digits = 3L, ...) {

  #### Check inputs ####

  if(!inherits(fit, "llca")) {
    stop("fit must inherit from class 'llca'.")
  }

  if(length(fit@Optim$parameters) == 0L) {
    stop("The llca object has not been fitted.")
  }

  if(is.null(parameters)) {

    parameters <- fit@modelInfo$trans[names(fit@modelInfo$param)]

  } else {

    selected_parameters <- unique(unlist(parameters))

    if(!all(selected_parameters %in%
            fit@modelInfo$transparameters_labels)) {
      stop("Unknown parameters.")
    }

  }

  type <- match.arg(tolower(type),
                    c("information", "standard", "robust"))

  if(type == "standard") {
    type <- "information"
  }

  if(!is.numeric(confidence) ||
     length(confidence) != 1L ||
     !is.finite(confidence) ||
     confidence <= 0 ||
     confidence >= 1) {
    stop("confidence must be a finite numeric value strictly between 0 and 1.")
  }

  if(!is.numeric(digits) ||
     length(digits) != 1L ||
     !is.finite(digits) ||
     digits < 0L ||
     digits != as.integer(digits)) {
    stop("digits must be a non-negative integer.")
  }

  #### Standard errors ####

  SE <- se(fit,
           type = type,
           parameters = parameters,
           digits = NULL,
           ...)

  critical <- stats::qnorm(0.5+confidence/2)
  selected_parameters <- names(SE$se)
  estimates <- fit@Optim$transparameters[selected_parameters]

  lower <- estimates-critical*SE$se
  upper <- estimates+critical*SE$se
  names(lower) <- names(upper) <- selected_parameters

  #### Parameter table ####

  est <- fill_in(parameters,
                 fit@Optim$transparameters,
                 miss = NA)
  lower_table <- fill_in(parameters,
                         lower,
                         miss = NA)
  upper_table <- fill_in(parameters,
                         upper,
                         miss = NA)

  table <- combine_est_ci(lower = lower_table,
                          est = est,
                          upper = upper_table,
                          digits = digits)

  #### Result ####

  result <- list(table = table,
                 lower_table = lower_table,
                 upper_table = upper_table,
                 lower = lower,
                 upper = upper,
                 se = SE$se,
                 vcov = SE$vcov,
                 VCOV = SE,
                 SE = SE)

  return(result)

}

#'
#' @rdname ci.llca
#' @param object A legacy structural-after-measurement object containing a
#'   fitted structural \code{"llca"} component.
#' @method ci llca_sam
#' @export
ci.llca_sam <- function(object, ...) {

  result <- ci(object$structural, ...)

  #### Result ####

  return(result)

}

#'
#' @rdname ci.llca
#' @param model For the \code{"llcalist"} method, a collection of fitted
#'   \code{"llca"} models.
#' @method ci llcalist
#' @export
ci.llcalist <- function(model, type = "information", confidence = 0.95,
                        parameters = NULL, digits = 3L, ...) {

  #### Check inputs ####

  if(!inherits(model, "llcalist")) {
    stop("model must inherit from class 'llcalist'.")
  }

  if(length(model) == 0L) {
    stop("model must contain at least one fitted llca object.")
  }

  if(!all(vapply(model, inherits, logical(1L), what = "llca"))) {
    stop("All elements of model must inherit from class 'llca'.")
  }

  #### Confidence intervals ####

  result <- lapply(
    model,
    FUN = \(x) ci(x,
                  type = type,
                  confidence = confidence,
                  parameters = parameters,
                  digits = digits,
                  ...)
  )

  result_names <- names(model)

  if(is.null(result_names)) {
    result_names <- rep("", length(model))
  }

  unnamed <- is.na(result_names) | result_names == ""

  if(any(unnamed)) {

    nclasses <- vapply(
      model[unnamed],
      FUN = \(x) ncol(x@modelInfo$trans$class),
      FUN.VALUE = integer(1L)
    )

    result_names[unnamed] <- paste0("nclasses=", nclasses)

  }

  names(result) <- make.unique(result_names)
  class(result) <- "ci.llcalist"

  #### Result ####

  return(result)

}
