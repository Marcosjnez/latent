

setMethod("show", "llca", function(object) {

  # Print header with model name and version
  cat(sprintf("%s %s ended normally after %d iterations\n\n",
              "latent", as.character( packageVersion('latent') ),
              object@Optim$iterations))

  # Print Estimator, Optimization, and Parameters section
  if(isFALSE(object@Optim$control$penalties)){
    est <- "ML"
  }else{
    est <- "Penalized-ML"}

  cat(sprintf("  %-45s %s\n", "Estimator", est))
  cat(sprintf("  %-45s %s\n", "Optimization method", object@Optim$control$opt))
  cat(sprintf("  %-45s %d\n\n", "Number of model parameters", object@modelInfo$nparam))

  # Print Number of Observations
  cat(sprintf("  %-45s %d\n\n", "Number of observations", object@modelInfo$nobs))
  cat(sprintf("  %-45s %d\n\n", "Number of response patterns (include NA)", object@modelInfo$npatterns))
  cat(sprintf("  %-45s %d\n\n", "Number of possible patterns", object@modelInfo$npossible_patterns))

  # Print Model Test Section
  if(sum(object@modelInfo$item != "multinomial") == 0){
    ni <- object@summary_table$Observed
    mi <- object@summary_table$Estimated
    df <- object@modelInfo$df
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




