# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 22/05/2025
#'
#' @title
#' Latent Class Analysis.
#' @description
#'
#' Print the information of latent class models.
#'
#' @usage
#'
#' lca(data, model = rep("multinomial", ncol(data)), nclasses = 2L,
#' control = list(opt = "lbfgs", rstarts = 30L, cores = 1L))
#'
#' @param model Character vector with the model for each item.
#'
#' @details None.
#'
#' @return Stuff:
#' \item{df}{Degrees of freedom}
#'
#' @references
#'
#' None yet.
#'
#' @export
print.lca <- function(model) {

  # Print header with model name and version
  cat(sprintf("%s %s ended normally after %d iterations\n\n",
              "latent", "0.1.0", model$opt$iterations))

  # Print Estimator, Optimization, and Parameters section
  cat(sprintf("  %-45s %s\n", "Estimator", "Multinomial"))
  cat(sprintf("  %-45s %s\n", "Optimization method", model$opt$control$opt))
  cat(sprintf("  %-45s %d\n\n", "Number of model parameters", model$modelInfo$nparam))

  # Print Number of Observations
  cat(sprintf("  %-45s %d\n\n", "Number of observations", model$modelInfo$nobs))
  cat(sprintf("  %-45s %d\n\n", "Number of response patterns", model$modelInfo$npatterns))

  # Print Model Test Section
  chisq <- -2*model$loglik
  pval <- dchisq(chisq, df = model$modelInfo$df)
  cat("Model Test User Model:\n")
  cat("  ", paste(rep("-", 54), collapse = ""), "\n\n", sep = "")
  cat(sprintf("  %-45s %.3f\n", "Test statistic", chisq))
  cat(sprintf("  %-45s %d\n", "Degrees of freedom", model$modelInfo$df))
  cat(sprintf("  %-45s %.3f\n", "P-value (Chi-square)", pval))

  invisible(model)

}
