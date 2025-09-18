# Author: Marcos Jimenez
# Author: Mauricio Garnier-Villarreal
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 03/09/2025
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
#' @method summary lcfa
#' @export
summary.lcfa <- function(fit, digits = 3) {

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
  if(is.null(fit@loglik)) {
    est <- "ULS"
  } else{
    est <- "ML"
  }

  cat(sprintf("  %-45s %s\n", "Estimator", est))
  cat(sprintf("  %-45s %s\n", "Optimization method", fit@Optim$control$opt))
  cat(sprintf("  %-45s %d\n\n", "Number of model parameters", fit@modelInfo$nparam))
  cat(sprintf("  %-45s %d\n\n", "Number of patterns", fit@modelInfo$npatterns))

  # Print Number of Observations
  N <- fit@modelInfo$nobs
  cat(sprintf("  %-45s %d\n\n", "Number of observations", fit@modelInfo$nobs))

  dof <- fit@modelInfo$dof
  loglik <- fit@loss
  X2 <- loglik*N
  pval <- 1-pchisq(X2, df = dof)
  cat("Model Test User Model:\n")
  cat("  ", paste(rep("-", 54), collapse = ""), "\n\n", sep = "")
  cat(sprintf("  %-45s %.3f\n", "Test statistic (Chi-square)", X2))
  cat(sprintf("  %-45s %d\n", "Degrees of freedom", dof))
  cat(sprintf("  %-45s %.3f\n", "P-value (Chi-square)", pval))

  #### Summary stats ####

  lambda <- fit@transformed_pars[[1]]$lambda
  psi <- fit@transformed_pars[[1]]$psi
  theta <- fit@transformed_pars[[1]]$theta
  rhat <- lambda %*% psi %*% t(lambda) + theta

  select <- match(fit@modelInfo$transparameters_labels,
                  unlist(fit@modelInfo$cfa_param))
  select <- select[!is.na(select)]

  est <- c(lambda, psi, theta)[select]
  se <- unname(se(fit, type = "standard", model = "user", digits = 3)$se)
  z <- abs(est) / se
  pval <- 2*(1-pnorm(z))

  labels <- fit@modelInfo$parameters_labels
  numeric_matrix <- round(cbind(est = est, se = se, z = z, pval = pval),
                          digits = digits)
  result <- data.frame(labels = labels, numeric_matrix, stringsAsFactors = FALSE)

  ## Header
  cat("\nParameter Estimates:\n")
  cat("-----------------------------------------------------\n")

  ## Column names
  cat(sprintf("%-20s %10s %10s %10s %10s\n",
              "label", "est", "se", "z", "pval"))
  cat("-----------------------------------------------------\n")

  ## Rows
  for (i in seq_len(nrow(result))) {
    cat(sprintf("%-20s %10.3f %10.3f %10.3f %10.3f\n",
                result$labels[i],
                result$est[i],
                result$se[i],
                result$z[i],
                result$pval[i]))
  }

  class(result) <- "summary.lcfa"

  invisible(result)

}
