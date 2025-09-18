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
  if(sum(fit@modelInfo$item != "multinomial") == 0){
    ni <- fit@summary_table$Observed
    mi <- fit@summary_table$Estimated
    dof <- fit@modelInfo$dof
    L2 <- abs(2*sum(ni*log(ni/mi)))
    pv <- 1-pchisq(L2, dof)
  }else{
    L2 <- NA
    pv <- NA
    dof <- NA
  }
  cat("Model Test User Model:\n")
  cat("  ", paste(rep("-", 54), collapse = ""), "\n\n", sep = "")
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
  loglik <- abs(fit@loss)
  X2 <- loglik*N
  pval <- 1-pchisq(X2, df = dof)
  cat("Model Test User Model:\n")
  cat("  ", paste(rep("-", 54), collapse = ""), "\n\n", sep = "")
  cat(sprintf("  %-45s %.3f\n", "Test statistic (Chi-square)", X2))
  cat(sprintf("  %-45s %d\n", "Degrees of freedom", dof))
  cat(sprintf("  %-45s %.3f\n", "P-value (Chi-square)", pval))

  invisible(fit)

})



