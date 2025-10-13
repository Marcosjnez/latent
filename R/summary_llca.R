# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 13/10/2025
#'
#' @title
#' Fit indices
#' @description
#'
#' Compute fit indices from any model.
#'
#' @usage
#'
#' getfit(model)
#'
#' @param model data.frame or matrix of response.
#'
#' @details \code{getfit} computes all the fit indices related to a specific model.
#'
#' @return List with the following fit indices:
#' \item{AIC}{.}
#' \item{BIC}{.}
#'
#' @references
#'
#' None yet.
#'
#' @method summary llca
#' @export
summary.llca <- function(fit, digits = 3) {

  #### Print fit ####

  conv <- fit@Optim$opt$convergence
  # Print header with model name and version
  if(conv) {

    cat(sprintf("%s %s converged after %d iterations\n\n",
                "latent", as.character( packageVersion('latent') ),
                fit@Optim$opt$iterations))

  } else {

    cat(sprintf("%s %s did not converged after %d iterations\n\n",
                "latent", as.character( packageVersion('latent') ),
                fit@Optim$opt$iterations))

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
    ni <- fit@summary_table$Observed
    mi <- fit@summary_table$Estimated
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

#' @method summary llcalist
#' @export
summary.llcalist <- function(model, digits = 3) {

  nmodels <- length(model)
  out <- vector("list", length = nmodels)
  for(i in 1:nmodels) {

    out[[i]] <- summary.llca(model[[i]], digits = digits)
    names(out)[i] <- paste("nclasses = ", model[[i]]@modelInfo$nclasses,
                           sep = "")

  }

  class(out) <- "summary.llcalist"

  invisible(out)

}
