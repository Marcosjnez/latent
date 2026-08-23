# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 23/08/2026
#'
#' Summarize a Fitted Exploratory Factor Analysis
#'
#' @param object A fitted object inheriting from class \code{"lefa"}.
#' @param digits Number of decimal places used when printing results.
#' @param fit Alias retained for compatibility. If supplied, \code{object} is
#'   ignored.
#' @param ... Additional arguments reserved for future summary options.
#'
#' @return Invisibly, an object of class \code{"summary.lefa"} containing the
#'   convergence information, fit indices, and rotated parameter table.
#'
#' @method summary lefa
#' @export
summary.lefa <- function(object, digits = 3L, fit = NULL, ...) {

  if(!is.null(fit)) {
    object <- fit
  }

  #### Check inputs ####

  if(!inherits(object, "lefa")) {
    stop("object must inherit from class 'lefa'.")
  }

  if(!is.numeric(digits) ||
     length(digits) != 1L ||
     !is.finite(digits) ||
     digits < 0L ||
     digits != as.integer(digits)) {
    stop("digits must be a non-negative integer.")
  }

  source <- lefa_source_lcfa(object)
  converged <- isTRUE(object@Optim$convergence)
  iterations <- object@Optim$iterations

  if(is.null(iterations) || length(iterations) == 0L) {
    iterations <- NA_integer_
  }

  estimator <- lefa_estimator_label(object)
  projection <- object@dataList$projection
  rotation <- object@dataList$rotation
  optimization <- object@modelInfo$control_optimizer$opt
  nparam <- source@modelInfo$nparam
  npatterns <- sum(unlist(source@dataList$npatterns))
  nobs <- sum(unlist(source@dataList$nobs))

  fit_indices <- getfit(object, digits = digits)
  parameter_table <- lefa_parameter_table(object)

  #### Print model information ####

  status <- if(converged) {
    "converged"
  } else {
    "did not converge"
  }

  cat("latent ", object@version, " ", status, sep = "")

  if(is.finite(iterations)) {
    cat(" after ", iterations, " rotation iterations", sep = "")
  }

  cat("\n\n")

  cat(sprintf("  %-38s %s\n", "Estimator", estimator))
  cat(sprintf("  %-38s %s\n", "Rotation", rotation))
  cat(sprintf("  %-38s %s\n", "Projection", projection))
  cat(sprintf("  %-38s %s\n", "Optimization method", optimization))
  cat(sprintf("  %-38s %d\n", "Number of model parameters", nparam))
  cat(sprintf("  %-38s %d\n", "Number of sample statistics", npatterns))
  cat(sprintf("  %-38s %d\n\n", "Number of observations", nobs))

  if(length(object@Optim$f) == 1L && is.finite(object@Optim$f)) {
    cat(sprintf("  %-38s %.*f\n\n",
                "Rotation loss", digits, object@Optim$f))
  }

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

  cat("\nRotated parameter estimates:\n")
  cat(strrep("-", 88L), "\n", sep = "")
  cat(sprintf("%-38s %11s %11s %11s %11s\n",
              "label", "estimate", "se", "z", "p-value"))
  cat(strrep("-", 88L), "\n", sep = "")

  for(i in seq_len(nrow(parameter_table))) {

    row <- parameter_table[i, ]

    cat(sprintf("%-38s %11.*f %11.*f %11.*f %11.*f\n",
                row$label,
                digits, row$estimate,
                digits, row$se,
                digits, row$z,
                digits, row$pvalue))

  }

  lefa_call <- object@dataList$lefa_call

  if(is.null(lefa_call)) {
    lefa_call <- object@call
  }

  result <- list(
    call = lefa_call,
    version = object@version,
    converged = converged,
    iterations = iterations,
    estimator = estimator,
    rotation = rotation,
    projection = projection,
    optimization = optimization,
    fit = fit_indices,
    parameters = parameter_table
  )
  class(result) <- "summary.lefa"

  #### Result ####

  return(invisible(result))

}
