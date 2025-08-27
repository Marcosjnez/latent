

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
  if(sum(model@modelInfo$item != "multinomial") == 0){
    ni <- model@summary_table$Observed
    mi <- model@summary_table$Estimated
    df <- model@modelInfo$df
    L2 <- 2*sum(ni*log(ni/mi))
    pv <- 1-pchisq(L2, df)
  }else{
    L2 <- NA
    pv <- NA
    df <- NA
  }
  cat("Model Test User Model:\n")
  cat("  ", paste(rep("-", 54), collapse = ""), "\n\n", sep = "")
  cat(sprintf("  %-45s %.3f\n", "Test statistic", L2))
  cat(sprintf("  %-45s %d\n", "Degrees of freedom", df))
  cat(sprintf("  %-45s %.3f\n", "P-value (L2)", pv))

  invisible(object)

})




