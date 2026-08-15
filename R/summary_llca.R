# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 04/05/2026
#'
#' Summary of a Latent Class Model
#'
#' Print a summary of a fitted latent class model.
#'
#' @description
#' \code{summary.llca()} displays optimization information, the estimation
#' method, the number of parameters and observations, available model-fit
#' statistics, and the estimated latent class profile.
#'
#' @param fit A fitted object of class \code{"llca"}.
#'
#' @return
#' The latent class profile returned by
#' \code{latInspect(fit, what = "profile")}, invisibly.
#'
#' @details
#' The printed output reports whether optimization converged, the optimization
#' method, the number of freely estimated parameters, and information about the
#' observed response patterns.
#'
#' For fully multinomial models, the likelihood-ratio statistic, degrees of
#' freedom, and corresponding p-value are also displayed.
#'
#' @seealso
#' \code{\link{getfit.llca}}, \code{\link{latInspect.llca}},
#' \code{\link{se.llca}}
#'
#' @method summary llca
#' @export
summary.llca <- function(fit) {

  #### Print fit ####

  conv <- fit@Optim$convergence
  # Print header with model name and version
  if(conv) {

    cat(sprintf("%s %s converged after %d iterations\n\n",
                "latent", as.character(packageVersion('latent')),
                fit@Optim$iterations))

  } else {

    cat(sprintf("%s %s did not converged after %d iterations\n\n",
                "latent", as.character(packageVersion('latent')),
                fit@Optim$iterations))

  }


  # Print Estimator, Optimization, and Parameters section
  if(isFALSE(fit@Optim$control$penalties)) {
    est <- "ML"
  } else {
    est <- "Penalized-ML"
  }

  cat(sprintf("  %-45s %s\n", "Estimator", est))
  cat(sprintf("  %-45s %s\n", "Optimization method", fit@Optim$control$opt))
  cat(sprintf("  %-45s %d\n\n", "Number of model parameters", fit@modelInfo$nparam))

  # Print Number of Observations
  cat(sprintf("  %-45s %d\n\n", "Number of observations", fit@modelInfo$nobs))
  cat(sprintf("  %-45s %d\n\n", "Number of response patterns (include NA)", fit@modelInfo$npatterns))
  cat(sprintf("  %-45s %d\n\n", "Number of possible patterns", fit@modelInfo$npossible_patterns))

  # Print Model Test Section
  if(sum(fit@modelInfo$item != "multinomial") == 0) {
    summary_table <- latInspect(fit, what = "summary")
    ni <- summary_table$Observed
    mi <- summary_table$Estimated
    dof <- fit@modelInfo$dof
    L2 <- 2*sum(ni*log(ni/mi))
    pv <- 1-pchisq(L2, dof)
  } else {
    L2 <- NA
    pv <- NA
    dof <- NA
  }
  cat("Model Test User Model:\n")
  cat("  ", paste(rep("-", 54), collapse = ""), "\n\n", sep = "")
  cat(sprintf("  %-45s %.3f\n", "Test statistic (L2)", L2))
  cat(sprintf("  %-45s %d\n", "Degrees of freedom", dof))
  cat(sprintf("  %-45s %.3f\n", "P-value (L2)", pv))

  result <- latInspect(fit, what = "profile")

  invisible(result)

}

#' @rdname summary.llca
#' @param model For the \code{"llcalist"} method, a collection of fitted
#'   \code{"llca"} and/or \code{"llca_sam"} models.
#' @details
#' For an \code{"llcalist"}, \code{summary()} is applied to each model in the
#' collection.
#' @method summary llcalist
#' @export
summary.llcalist <- function(model) {

  out <- lapply(model, FUN = summary)

  class(out) <- "summary.llcalist"

  invisible(out)

}
