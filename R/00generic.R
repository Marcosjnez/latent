# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 22/08/2026

.show_latent_row <- function(label, value) {

  cat(sprintf("  %-45s %s\n", label, value))

  #### Result ####

  return(invisible(NULL))

}

.show_latent_number <- function(label, value, digits = 3L) {

  if(length(value) == 1L && is.finite(value)) {
    .show_latent_row(label, formatC(value, format = "f", digits = digits))
  }

  #### Result ####

  return(invisible(NULL))

}

setMethod("show", "llca", function(object) {

  fitted <- length(object@Optim) > 0L &&
    length(object@Optim$transparameters) > 0L

  if(!fitted) {

    cat("Unfitted latent class model specification\n\n")
    .show_latent_row("Object class", class(object)[1L])
    .show_latent_row("Number of model parameters",
                     object@modelInfo$nparam)

    #### Result ####

    return(invisible(object))

  }

  converged <- isTRUE(object@Optim$convergence)
  iterations <- object@Optim$iterations

  status <- if(converged) {
    "converged"
  } else {
    "did not converge"
  }

  cat("latent ", object@version, " ", status, sep = "")

  if(length(iterations) == 1L && is.finite(iterations)) {
    cat(" after ", iterations, " iterations", sep = "")
  }

  cat("\n\n")

  control <- object@modelInfo$control_optimizer
  regularized <- isTRUE(control$reg)

  estimator <- if(regularized) {
    "Penalized maximum likelihood"
  } else {
    "Maximum likelihood"
  }

  .show_latent_row("Estimator", estimator)
  .show_latent_row("Optimization method", control$opt)
  .show_latent_row("Number of model parameters",
                   object@modelInfo$nparam)
  .show_latent_row("Number of observations",
                   object@dataList$nobs)
  .show_latent_row("Number of response patterns",
                   object@dataList$npatterns)
  .show_latent_row("Number of possible patterns",
                   object@dataList$npossible_patterns)

  fit_indices <- tryCatch(
    getfit(object, digits = NULL),
    error = function(e) NULL
  )

  if(!is.null(fit_indices)) {

    cat("\nModel fit:\n")
    cat(strrep("-", 58L), "\n", sep = "")

    .show_latent_number("Log-likelihood",
                        fit_indices[["loglik"]])

    if("penalized_loglik" %in% names(fit_indices)) {
      .show_latent_number("Penalized log-likelihood",
                          fit_indices[["penalized_loglik"]])
    }

    .show_latent_number("Likelihood-ratio statistic",
                        fit_indices[["L2"]])
    .show_latent_number("Degrees of freedom",
                        fit_indices[["dof"]], digits = 0L)
    .show_latent_number("P-value",
                        fit_indices[["pvalue"]])

  }

  #### Result ####

  return(invisible(object))

})

setMethod("show", "lcfa", function(object) {

  fitted <- length(object@Optim) > 0L &&
    length(object@Optim$transparameters) > 0L

  if(!fitted) {

    cat("Unfitted confirmatory factor model specification\n\n")
    .show_latent_row("Object class", class(object)[1L])
    .show_latent_row("Number of model parameters",
                     object@modelInfo$nparam)

    #### Result ####

    return(invisible(object))

  }

  converged <- isTRUE(object@Optim$convergence)
  iterations <- object@Optim$iterations

  status <- if(converged) {
    "converged"
  } else {
    "did not converge"
  }

  cat("latent ", object@version, " ", status, sep = "")

  if(length(iterations) == 1L && is.finite(iterations)) {
    cat(" after ", iterations, " iterations", sep = "")
  }

  cat("\n\n")

  estimator <- switch(
    object@dataList$estimator,
    ml = "Maximum likelihood",
    fml = "Maximum likelihood",
    uls = "Unweighted least squares",
    dwls = "Diagonally weighted least squares",
    object@dataList$estimator
  )

  if(isTRUE(object@modelInfo$control_optimizer$reg)) {
    estimator <- paste("Penalized", tolower(estimator))
  }

  .show_latent_row("Estimator", estimator)
  .show_latent_row(
    "Optimization method",
    object@modelInfo$control_optimizer$opt
  )
  .show_latent_row("Number of model parameters",
                   object@modelInfo$nparam)
  .show_latent_row("Number of sample statistics",
                   sum(unlist(object@dataList$npatterns)))
  .show_latent_row("Number of observations",
                   sum(unlist(object@dataList$nobs)))

  # Keep show() inexpensive. Full likelihood-ratio and incremental fit indices,
  # including the saturated incomplete-data model for direct FIML, are computed
  # by summary() or getfit().
  fit_matrix <- tryCatch(
    lcfa_fit_matrix(object, compute_h1 = FALSE),
    error = function(e) NULL
  )

  if(!is.null(fit_matrix)) {

    cat("\nModel objective:\n")
    cat(strrep("-", 58L), "\n", sep = "")

    .show_latent_number(
      "Loss",
      fit_matrix["loss", "overall"]
    )
    .show_latent_number(
      "Penalized loss",
      fit_matrix["penalized_loss", "overall"]
    )
    .show_latent_number(
      "Log-likelihood",
      fit_matrix["loglik", "overall"]
    )
    .show_latent_number(
      "Penalized log-likelihood",
      fit_matrix["penalized_loglik", "overall"]
    )

  }

  #### Result ####

  return(invisible(object))

})
