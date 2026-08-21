# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 21/08/2026
#'
#' Summary of a Latent Class Model
#'
#' Print optimization information, likelihood-based fit indices, and the
#' estimated latent class profile of a fitted latent class model.
#'
#' @param fit A fitted object inheriting from class \code{"llca"}.
#' @param digits Non-negative integer giving the number of decimal places used
#'   in printed numeric results.
#' @param ... Additional arguments reserved for future summary options.
#'
#' @return The latent class profile returned by
#'   \code{latInspect(fit, what = "profile")}, invisibly.
#'
#' @method summary llca
#' @export
summary.llca <- function(fit, digits = 3L, ...) {

  #### Check inputs ####

  if(!inherits(fit, "llca")) {
    stop("fit must inherit from class 'llca'.")
  }

  if(length(fit@Optim) == 0L ||
     length(fit@Optim$transparameters) == 0L) {
    stop("The llca object has not been fitted.")
  }

  if(!is.numeric(digits) ||
     length(digits) != 1L ||
     !is.finite(digits) ||
     digits < 0L ||
     digits != as.integer(digits)) {
    stop("digits must be a non-negative integer.")
  }

  #### Model information ####

  converged <- isTRUE(fit@Optim$convergence)
  iterations <- fit@Optim$iterations
  control <- fit@modelInfo$control_optimizer
  regularized <- isTRUE(control$reg)

  estimator <- if(regularized) {
    "Penalized maximum likelihood"
  } else {
    "Maximum likelihood"
  }

  status <- if(converged) {
    "converged"
  } else {
    "did not converge"
  }

  cat("latent ", fit@version, " ", status, sep = "")

  if(length(iterations) == 1L && is.finite(iterations)) {
    cat(" after ", iterations, " iterations", sep = "")
  }

  cat("\n\n")

  cat(sprintf("  %-45s %s\n", "Estimator", estimator))
  cat(sprintf("  %-45s %s\n", "Optimization method", control$opt))
  cat(sprintf("  %-45s %d\n", "Number of model parameters",
              fit@modelInfo$nparam))
  cat(sprintf("  %-45s %d\n", "Number of observations",
              fit@dataList$nobs))
  cat(sprintf("  %-45s %d\n", "Number of response patterns",
              fit@dataList$npatterns))
  cat(sprintf("  %-45s %d\n\n", "Number of possible patterns",
              fit@dataList$npossible_patterns))

  #### Fit indices ####

  fit_indices <- getfit(fit, digits = NULL)

  cat("Model fit:\n")
  cat(strrep("-", 72L), "\n", sep = "")

  print_names <- c(
    loglik = "Log-likelihood",
    penalized_loglik = "Penalized log-likelihood",
    L2 = "Likelihood-ratio statistic",
    dof = "Degrees of freedom",
    pvalue = "P-value",
    AIC = "AIC",
    BIC = "BIC",
    SABIC = "Sample-size-adjusted BIC",
    ICL = "Integrated classification likelihood",
    R2_entropy = "Entropy R-squared"
  )

  for(name in names(print_names)) {

    if(!name %in% names(fit_indices)) {
      next
    }

    value <- fit_indices[[name]]

    if(length(value) != 1L || !is.finite(value)) {
      next
    }

    value_digits <- if(name == "dof") 0L else digits

    cat(sprintf("  %-45s %.*f\n",
                print_names[[name]], value_digits, value))

  }

  #### Latent class profile ####

  profile <- latInspect(fit, what = "profile", digits = digits)

  cat("\nLatent class profile:\n")
  cat(strrep("-", 72L), "\n", sep = "")
  print(profile, digits = digits)

  #### Result ####

  return(invisible(profile))

}

#'
#' @rdname summary.llca
#' @param model For the \code{"llcalist"} method, a collection of fitted
#'   \code{"llca"} models.
#' @method summary llcalist
#' @export
summary.llcalist <- function(model, digits = 3L, ...) {

  if(!inherits(model, "llcalist")) {
    stop("model must inherit from class 'llcalist'.")
  }

  if(length(model) == 0L) {
    stop("model must contain at least one fitted llca object.")
  }

  result <- vector("list", length = length(model))

  for(i in seq_along(model)) {
    result[[i]] <- summary(model[[i]],
                           digits = digits,
                           ...)
  }

  names(result) <- names(model)
  class(result) <- "summary.llcalist"

  #### Result ####

  return(invisible(result))

}
