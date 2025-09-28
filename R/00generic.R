setMethod("show", "llca", function(fit) {

  conv <- fit@Optim$opt$convergence
  # Print header with model name and version
  if(conv){
    cat(sprintf("%s %s converged after %d iterations\n\n",
                "latent", as.character( packageVersion('latent') ),
                fit@Optim$opt$iterations))
  }else{
    cat(sprintf("%s %s did not converged after %d iterations\n\n",
                "latent", as.character( packageVersion('latent') ),
                fit@Optim$opt$iterations))
  }


  # Print Estimator, Optimization, and Parameters section
  if(isFALSE(fit@Optim$control$penalties)){
    est <- "ML"
  }else{
    est <- "Penalized-ML"}

  cat(sprintf("  %-45s %s\n", "Estimator", est))
  cat(sprintf("  %-45s %s\n", "Optimization method", fit@Optim$control$opt))
  cat(sprintf("  %-45s %d\n\n", "Number of model parameters", fit@modelInfo$nparam))

  # Print Number of Observations
  cat(sprintf("  %-45s %d\n\n", "Number of observations", fit@modelInfo$nobs))
  cat(sprintf("  %-45s %d\n\n", "Number of response patterns (include NA)", fit@modelInfo$npatterns))
  cat(sprintf("  %-45s %d\n\n", "Number of possible patterns", fit@modelInfo$npossible_patterns))

  # Print Model Test Section
  fit_ind <- getfit(fit)
  L2 <- fit_ind[["L2"]]
  pv <- fit_ind[["pvalue"]]
  dof <- fit_ind[["dof"]]

  cat("  ", paste(rep("-", 54), collapse = ""), "\n\n", sep = "")

  cat("Model Test User Model:\n")
  cat(sprintf("  %-45s %.3f\n", "Test statistic (L2)", L2))
  cat(sprintf("  %-45s %d\n", "Degrees of freedom", dof))
  cat(sprintf("  %-45s %.3f\n", "P-value (L2)", pv))

  invisible(fit)

})

setMethod("show", "lcfa", function(fit) {

  conv <- fit@Optim$opt$convergence
  # Print header with model name and version
  if(conv){
    cat(sprintf("%s %s converged after %d iterations\n\n",
                "latent", as.character( packageVersion('latent') ),
                fit@Optim$opt$iterations))
  }else{
    cat(sprintf("%s %s did not converged after %d iterations\n\n",
                "latent", as.character( packageVersion('latent') ),
                fit@Optim$opt$iterations))
  }


  N <- fit@modelInfo$nobs
  dof <- fit@modelInfo$dof


  # Print Estimator, Optimization, and Parameters section
  if(length(fit@loglik) == 0) {

    loglik <- NA
    X2 <- NA
    pval <- NA
    est <- "ULS"
    test_message <- sprintf("  %-45s %.3f\n", "Test statistic", fit@loss)
    pval_message <- sprintf("  %-45s %.3f\n", "P-value (Unknown)", pval)

  } else {

    loglik <- abs(fit@loss)
    X2 <- loglik*N
    pval <- 1-pchisq(X2, df = dof)
    est <- "ML"
    test_message <- sprintf("  %-45s %.3f\n", "Test statistic (Chi-square)", X2)
    pval_message <- sprintf("  %-45s %.3f\n", "P-value (Chi-square)", pval)

  }

  cat(sprintf("  %-45s %s\n", "Estimator", est))
  cat(sprintf("  %-45s %s\n", "Optimization method", fit@Optim$control$opt))
  cat(sprintf("  %-45s %d\n\n", "Number of model parameters", fit@modelInfo$nparam))
  cat(sprintf("  %-45s %d\n\n", "Number of patterns", fit@modelInfo$npatterns))

  # Print Number of Observations
  cat(sprintf("  %-45s %d\n\n", "Number of observations", fit@modelInfo$nobs))

  cat("  ", paste(rep("-", 54), collapse = ""), "\n\n", sep = "")

  cat("Model Test User Model:\n")
  cat(test_message)
  cat(sprintf("  %-45s %d\n", "Degrees of freedom", dof))
  cat(pval_message)

  invisible(fit)

})



