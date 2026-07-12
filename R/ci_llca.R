# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 12/07/2026
#'
#' Confidence intervals for latent class models
#'
#' Computes confidence intervals for the parameters of a fitted latent class
#' model.
#'
#' @param fit A fitted object of class \code{"llca"}.
#' @param type Character string indicating the standard-error estimator used to
#'   construct the intervals. Available options are \code{"standard"} and
#'   \code{"robust"}. See \code{se.llca()}.
#' @param confidence Numeric scalar strictly between zero and one specifying the
#'   confidence level.
#' @param digits Non-negative integer indicating the number of decimal places used
#'   in the formatted confidence-interval table.
#' @param ... Additional arguments passed to other methods.
#'
#' @details
#' Confidence limits are computed using the asymptotic normal approximation. The
#' critical value is obtained as the square root of the corresponding one-degree-
#' of-freedom chi-squared quantile. Standard errors are obtained through
#' \code{se()}, so two-step and robust adjustments are used when requested.
#'
#' The \code{digits} argument affects only the formatted \code{table}. The numeric
#' confidence limits, standard errors, and covariance matrices are returned
#' without rounding.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{table}}{A list of formatted parameter tables showing the estimate
#'   and its confidence interval.}
#'   \item{\code{lower_table}}{The lower confidence limits arranged in the model
#'   parameter structure.}
#'   \item{\code{upper_table}}{The upper confidence limits arranged in the model
#'   parameter structure.}
#'   \item{\code{lower}}{A named numeric vector of lower confidence limits.}
#'   \item{\code{upper}}{A named numeric vector of upper confidence limits.}
#'   \item{\code{se}}{A named numeric vector of standard errors.}
#'   \item{\code{vcov}}{The variance-covariance matrix used to construct the
#'   intervals.}
#'   \item{\code{B}}{The empirical or two-step correction matrix, when available.}
#'   \item{\code{H}}{The Hessian matrix, when available.}
#'   \item{\code{newH}}{The adjusted Hessian used by the robust estimator, when
#'   available.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- lca(data = empathy, nclasses = 3L,
#'            gaussian = c("ec1", "ec2", "ec3"))
#'
#' ci(fit)
#' ci(fit, type = "robust", confidence = 0.90, digits = 4L)
#' }
#'
#' @seealso \code{se()}
#'
#' @method ci llca
#' @export
ci.llca <- function(fit, type = "standard", confidence = 0.95,
                    digits = 3L, ...) {

  if(!inherits(fit, "llca")) {
    stop("fit must inherit from class 'llca'.")
  }

  if(length(fit@Optim$parameters) == 0L) {
    stop("The llca object has not been fitted.")
  }

  type <- match.arg(tolower(type), c("standard", "robust"))

  if(!is.numeric(confidence) || length(confidence) != 1L ||
     !is.finite(confidence) || confidence <= 0 || confidence >= 1) {
    stop("confidence must be a finite numeric value strictly between 0 and 1.")
  }

  if(!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
     digits < 0L || digits != as.integer(digits)) {
    stop("digits must be a non-negative integer.")
  }

  SE <- se(fit, type = type, digits = NULL, ...)

  parameter_names <- fit@modelInfo$parameters_labels
  x <- c(fit@Optim$parameters)
  if(is.null(names(x))) {
    names(x) <- parameter_names
  }
  x <- x[parameter_names]

  se_values <- c(SE$se)
  if(is.null(names(se_values))) {
    names(se_values) <- parameter_names
  }
  se_values <- se_values[parameter_names]

  if(anyNA(x) || anyNA(se_values)) {
    stop("The parameter estimates and standard errors could not be aligned by name.")
  }

  critical <- sqrt(stats::qchisq(confidence, df = 1L))
  lower <- x - critical*se_values
  upper <- x + critical*se_values
  names(lower) <- names(upper) <- parameter_names

  est <- fill_in(fit@modelInfo$trans[names(fit@modelInfo$param)],
                 c(fit@Optim$parameters, fit@Optim$transparameters), miss = NA)
  lower_ci <- fill_in(fit@modelInfo$param, lower)
  upper_ci <- fill_in(fit@modelInfo$param, upper)
  table <- combine_est_ci(lower_ci, est, upper_ci, digits = digits)

  result <- list(table = table, lower_table = lower_ci,
                 upper_table = upper_ci, lower = lower, upper = upper,
                 se = se_values, vcov = SE$vcov, B = SE$B,
                 H = SE$H, newH = SE$newH)

  #### Result ####

  return(result)

}

#' Confidence intervals for a collection of latent class models
#'
#' Applies \code{ci()} to every fitted model in an \code{"llcalist"} object.
#'
#' @param model An object of class \code{"llcalist"} containing fitted
#'   \code{"llca"} objects.
#' @param type Character string indicating the standard-error estimator. See
#'   \code{ci.llca()}.
#' @param confidence Numeric scalar strictly between zero and one specifying the
#'   confidence level.
#' @param digits Non-negative integer indicating the number of decimal places used
#'   in the formatted tables.
#' @param ... Additional arguments passed to \code{ci.llca()}.
#'
#' @details
#' Existing names are preserved. Consequently, class-enumeration results retain
#' names such as \code{"nclasses=2"}, whereas multiple-step models retain names
#' such as \code{"measurement"} and \code{"structural"}. Unnamed elements are
#' labelled according to their number of latent classes.
#'
#' @return A named list with one confidence-interval result per fitted model and
#'   class \code{"ci.llcalist"}.
#'
#' @examples
#' \dontrun{
#' fits <- lca(data = empathy, nclasses = 2:4,
#'             gaussian = c("ec1", "ec2", "ec3"))
#' ci(fits, confidence = 0.90)
#' }
#'
#' @method ci llcalist
#' @export
ci.llcalist <- function(model, type = "standard", confidence = 0.95,
                        digits = 3L, ...) {

  if(!inherits(model, "llcalist")) {
    stop("model must inherit from class 'llcalist'.")
  }

  if(length(model) == 0L) {
    stop("model must contain at least one fitted llca object.")
  }

  if(!all(vapply(model, inherits, logical(1L), what = "llca"))) {
    stop("All elements of model must inherit from class 'llca'.")
  }

  nmodels <- length(model)
  result <- vector("list", length = nmodels)

  for(i in seq_len(nmodels)) {
    result[[i]] <- ci(model[[i]], type = type, confidence = confidence,
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
  class(result) <- "ci.llcalist"

  #### Result ####

  return(result)

}
