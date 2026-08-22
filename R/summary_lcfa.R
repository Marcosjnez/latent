# Author: Marcos Jimenez
# Author: Mauricio Garnier-Villarreal
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 22/08/2026
#'
#' Summarize a Fitted CFA Model
#'
#' @param object A fitted object inheriting from class \code{"lcfa"}.
#' @param digits Number of decimal places used when printing results.
#' @param fit Alias retained for backward compatibility. If supplied,
#'   \code{object} is ignored.
#' @param ... Additional arguments reserved for future summary options.
#'
#' @return Invisibly, an object of class \code{"summary.lcfa"} containing the
#'   convergence information, fit indices, and free-parameter table.
#'
#' @method summary lcfa
#' @export
summary.lcfa <- function(object, digits = 3L, fit = NULL, ...) {

  if(!is.null(fit)) {
    object <- fit
  }

  #### Check inputs ####

  if(!inherits(object, "lcfa")) {
    stop("object must inherit from class 'lcfa'.")
  }

  if(!is.numeric(digits) ||
     length(digits) != 1L ||
     !is.finite(digits) ||
     digits < 0L ||
     digits != as.integer(digits)) {
    stop("digits must be a non-negative integer.")
  }

  converged <- isTRUE(object@Optim$convergence)
  iterations <- object@Optim$iterations

  if(is.null(iterations) || length(iterations) == 0L) {
    iterations <- NA_integer_
  }

  estimator <- switch(
    object@dataList$estimator,
    ml = "Maximum likelihood",
    fml = "Maximum likelihood",
    uls = "Unweighted least squares",
    dwls = "Diagonally weighted least squares",
    object@dataList$estimator
  )

  optimization <- object@modelInfo$control_optimizer$opt
  nparam <- object@modelInfo$nparam
  npatterns <- sum(unlist(object@dataList$npatterns))
  nobs <- sum(unlist(object@dataList$nobs))

  fit_indices <- getfit(object, digits = digits)
  parameter_table <- lcfa_parameter_table(object)

  #### Print model information ####

  status <- if(converged) {
    "converged"
  } else {
    "did not converge"
  }

  cat("latent ", object@version, " ", status, sep = "")

  if(is.finite(iterations)) {
    cat(" after ", iterations, " iterations", sep = "")
  }

  cat("\n\n")

  cat(sprintf("  %-38s %s\n", "Estimator", estimator))
  cat(sprintf("  %-38s %s\n", "Optimization method", optimization))
  cat(sprintf("  %-38s %d\n", "Number of model parameters", nparam))
  cat(sprintf("  %-38s %d\n", "Number of sample statistics", npatterns))
  cat(sprintf("  %-38s %d\n\n", "Number of observations", nobs))

  #### Print fit indices ####

  cat("Model fit:\n")
  cat(strrep("-", 58L), "\n", sep = "")

  print_fit <- c(
    "Log-likelihood" = fit_indices[["loglik"]],
    "AIC" = fit_indices[["AIC"]],
    "BIC" = fit_indices[["BIC"]],
    "Chi-square" = fit_indices[["chisq"]],
    "Degrees of freedom" = fit_indices[["dof"]],
    "P-value" = fit_indices[["pvalue"]],
    "CFI" = fit_indices[["cfi"]],
    "RMSEA" = fit_indices[["rmsea"]],
    "SRMR" = fit_indices[["srmr"]]
  )

  for(name in names(print_fit)) {

    value <- print_fit[[name]]

    if(is.finite(value)) {
      cat(sprintf("  %-38s %.*f\n",
                  name, digits, value))
    } else {
      cat(sprintf("  %-38s %s\n",
                  name, "NA"))
    }

  }

  #### Print parameter table ####

  cat("\nParameter estimates:\n")
  cat(strrep("-", 78L), "\n", sep = "")
  cat(sprintf("%-28s %11s %11s %11s %11s\n",
              "label", "estimate", "se", "z", "p-value"))
  cat(strrep("-", 78L), "\n", sep = "")

  for(i in seq_len(nrow(parameter_table))) {

    row <- parameter_table[i, ]

    cat(sprintf("%-28s %11.*f %11.*f %11.*f %11.*f\n",
                row$label,
                digits, row$estimate,
                digits, row$se,
                digits, row$z,
                digits, row$pvalue))

  }

  result <- list(
    call = object@call,
    version = object@version,
    converged = converged,
    iterations = iterations,
    estimator = estimator,
    optimization = optimization,
    fit = fit_indices,
    parameters = parameter_table
  )
  class(result) <- "summary.lcfa"

  #### Result ####

  return(invisible(result))

}
