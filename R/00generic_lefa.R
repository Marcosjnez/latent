# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 24/08/2026

setMethod("show", "lefa", function(object) {

  fitted <- length(object@Optim) > 0L &&
    length(object@Optim$transparameters) > 0L

  if(!fitted) {

    cat("Unfitted exploratory factor model specification\n\n")
    .show_latent_row("Object class", class(object)[1L])

    #### Result ####

    return(invisible(object))

  }

  source <- lefa_source_lcfa(object)
  converged <- isTRUE(object@Optim$convergence)
  iterations <- object@Optim$iterations

  status <- if(converged) {
    "converged"
  } else {
    "did not converge"
  }

  cat("latent ", object@version, " ", status, sep = "")

  if(length(iterations) == 1L && is.finite(iterations)) {
    cat(" after ", iterations, " rotation iterations", sep = "")
  }

  cat("\n\n")

  .show_latent_row("Estimator", lefa_estimator_label(object))
  .show_latent_row("Rotation", object@dataList$rotation)
  .show_latent_row("Projection", object@dataList$projection)
  .show_latent_row(
    "Optimization method",
    object@modelInfo$control_optimizer$opt
  )
  .show_latent_row("Number of model parameters",
                   source@modelInfo$nparam)
  .show_latent_row("Number of sample statistics",
                   sum(unlist(source@dataList$npatterns)))
  .show_latent_row("Number of observations",
                   sum(unlist(source@dataList$nobs)))

  if(length(object@Optim$f) == 1L && is.finite(object@Optim$f)) {
    .show_latent_number("Rotation loss", object@Optim$f)
  }

  fit_matrix <- tryCatch(
    lcfa_fit_matrix(source, compute_h1 = FALSE),
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

#'
#' Standard Errors for Exploratory Factor Analysis
#'
#' @param fit A fitted object inheriting from class \code{"lefa"}.
#' @param type Optional compatibility argument passed to
#'   \code{se.multistep()}.
#' @param parameters Optional parameter specification. By default, standard
#'   errors are returned for the rotated loadings, factor covariance matrices,
#'   and factor means.
#' @param digits Non-negative integer used to format parameter tables. If
#'   \code{NULL}, values are not rounded.
#' @param ... Additional arguments passed to \code{se.multistep()}.
#'
#' @return The result of \code{se.multistep()}.
#'
#' @method se lefa
#' @export
se.lefa <- function(fit, type = NULL, parameters = NULL,
                    digits = 4L, ...) {

  if(is.null(parameters)) {
    parameters <- lefa_rotated_parameters(fit)
  }

  result <- se.multistep(fit = fit,
                         type = type,
                         parameters = parameters,
                         digits = digits,
                         ...)

  #### Result ####

  return(result)

}
