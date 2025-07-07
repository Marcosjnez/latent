

setMethod("show", "llca", function(object) {

  # Print header with model name and version
  cat(sprintf("%s %s ended normally after %d iterations\n\n",
              "latent", as.character( packageVersion('latent') ),
              object@Optim$iterations))

  # Print Estimator, Optimization, and Parameters section
  cat(sprintf("  %-45s %s\n", "Estimator", "Multinomial"))
  cat(sprintf("  %-45s %s\n", "Optimization method", object@Optim$control$opt))
  cat(sprintf("  %-45s %d\n\n", "Number of model parameters", object@modelInfo$nparam))

  # Print Number of Observations
  cat(sprintf("  %-45s %d\n\n", "Number of observations", object@modelInfo$nobs))
  cat(sprintf("  %-45s %d\n\n", "Number of response patterns", object@modelInfo$npatterns))

  # Print Model Test Section
  chisq <- -2*object@loglik
  pval <- 1 - pchisq(chisq, df = object@modelInfo$df)
  cat("Model Test User Model:\n")
  cat("  ", paste(rep("-", 54), collapse = ""), "\n\n", sep = "")
  cat(sprintf("  %-45s %.3f\n", "Test statistic", chisq))
  cat(sprintf("  %-45s %d\n", "Degrees of freedom", object@modelInfo$df))
  cat(sprintf("  %-45s %.3f\n", "P-value (Chi-square)", pval))

  invisible(object)

})




