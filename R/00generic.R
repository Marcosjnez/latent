# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 03/11/2025

setMethod("show", "llca", function(object) {

  conv <- object@Optim$opt$convergence
  # Print header with model name and version
  if(conv){
    cat(sprintf("%s %s converged after %d iterations\n\n",
                "latent", as.character( packageVersion('latent') ),
                object@Optim$opt$iterations))
  }else{
    cat(sprintf("%s %s did not converged after %d iterations\n\n",
                "latent", as.character( packageVersion('latent') ),
                object@Optim$opt$iterations))
  }


  # Print Estimator, Optimization, and Parameters section
  if(isFALSE(object@modelInfo$control$penalties)){
    est <- "ML"
  }else{
    est <- "Penalized-ML"}

  cat(sprintf("  %-45s %s\n", "Estimator", est))
  cat(sprintf("  %-45s %s\n", "Optimization method", object@modelInfo$control$opt))
  cat(sprintf("  %-45s %d\n\n", "Number of model parameters", object@modelInfo$nparam))

  # Print Number of Observations
  cat(sprintf("  %-45s %d\n\n", "Number of observations", object@modelInfo$nobs))
  cat(sprintf("  %-45s %d\n\n", "Number of response patterns (include NA)", object@modelInfo$npatterns))
  cat(sprintf("  %-45s %d\n\n", "Number of possible patterns", object@modelInfo$npossible_patterns))

  # Print Model Test Section
  fit_ind <- getfit(object)
  L2 <- fit_ind[["L2"]]
  pv <- fit_ind[["pvalue"]]
  dof <- fit_ind[["dof"]]

  cat("  ", paste(rep("-", 54), collapse = ""), "\n\n", sep = "")

  cat("Model Test User Model:\n")
  cat(sprintf("  %-45s %.3f\n", "Test statistic (L2)", L2))
  cat(sprintf("  %-45s %d\n", "Degrees of freedom", dof))
  cat(sprintf("  %-45s %.3f\n", "P-value (L2)", pv))

  invisible()

})

setMethod("show", "lcfa", function(object) {

  conv <- object@Optim$opt$convergence
  # Print header with model name and version
  if(conv) {
    cat(sprintf("%s %s converged after %d iterations\n\n",
                "latent", as.character( packageVersion('latent') ),
                object@Optim$opt$iterations))
  } else {
    cat(sprintf("%s %s did not converged after %d iterations\n\n",
                "latent", as.character( packageVersion('latent') ),
                object@Optim$opt$iterations))
  }


  N <- object@modelInfo$nobs
  dof <- object@modelInfo$dof

  opt_method <- object@modelInfo$control$opt
  nparam <- object@modelInfo$nparam
  npatterns <- sum(unlist(object@modelInfo$npatterns))
  nobs <- sum(unlist(object@modelInfo$nobs))

  # Print Estimator, Optimization, and Parameters section
  if(length(object@loglik) == 0) {

    loglik <- NA
    X2 <- NA
    pval <- NA
    est <- "ULS"
    test_message <- sprintf("  %-45s %.3f\n", "Test statistic", object@loss)
    pval_message <- sprintf("  %-45s %.3f\n", "P-value (Unknown)", pval)

  } else {

    loglik <- object@loglik
    # X2 <- 2*loglik
    X2 <- object@loss*nobs
    pval <- 1-pchisq(X2, df = dof)
    est <- "ML"
    test_message <- sprintf("  %-45s %.3f\n", "Test statistic (Chi-square)", X2)
    pval_message <- sprintf("  %-45s %.3f\n", "P-value (Chi-square)", pval)

  }

  cat(sprintf("  %-45s %s\n", "Estimator", est))
  cat(sprintf("  %-45s %s\n", "Optimization method", opt_method))
  cat(sprintf("  %-45s %d\n\n", "Number of model parameters", nparam))
  cat(sprintf("  %-45s %d\n\n", "Number of patterns", npatterns))

  # Print Number of Observations
  cat(sprintf("  %-45s %d\n\n", "Number of observations", nobs))

  cat("  ", paste(rep("-", 54), collapse = ""), "\n\n", sep = "")

  cat("Model Test User Model:\n")
  cat(test_message)
  cat(sprintf("  %-45s %d\n", "Degrees of freedom", dof))
  cat(pval_message)

  invisible()

})



