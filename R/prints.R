print.efa <- function(x, nobs = NULL, ...) {

  efa <- x
  # Check if nobs was provided
  if(is.null(nobs)) {
    if(is.null(efa$modelInfo$nobs) | isTRUE(efa$modelInfo$nobs == 0L)) {
      warning("Sample size was not provided. Some Chi-squared-based statistics will not be computed.")
      nobs <- NA
    } else {
      nobs <- efa$modelInfo$nobs
    }
  }
  ### Pattern matrix with communalities, uniqueness, and complexity
  lambda <- efa$rotation$lambda
  phi <- efa$rotation$phi

  ObjFn <- efa$efa$f
  ordering <- order(diag(phi %*% t(lambda) %*% lambda), decreasing=T)
  fit <- suppressWarnings(fit(efa, nobs))
  phi <- phi[ordering, ordering]
  rownames(phi) <- colnames(phi) <- paste("F",sprintf(paste("%0",nchar(efa$modelInfo$nfactors),"d",sep=""),
                                                      ordering),sep="")
  colnames(lambda) <- rownames(phi) <-
    colnames(phi) <- paste("F",sprintf(paste("%0",nchar(ncol(lambda)),"d",sep=""),1:ncol(lambda)),sep="")
  if(is.null(colnames(efa$modelInfo$correlation))) {
    rownames(lambda) <- paste("item_",sprintf(paste("%0",nchar(nrow(lambda)),"d",sep=""),1:nrow(lambda)),sep="")
  } else {
    rownames(lambda) <- colnames(efa$modelInfo$correlation)
  }
  # Print
  cat("Factor Analysis using estimator = ", efa$modelInfo$estimator, "\n", sep="")
  cat("Elapsed time of ", round(efa$elapsed * 1e-9, 3), " seconds", "\n", sep="")
  cat("\n","Goodness-of-fit and model misfit indices", sep="")
  cat("\n","The standardized root mean square residual (SRMR) is ", round(fit$indices["SRMR", 1], 3), sep="")
  cat("\n","The largest absolute value of standardized residual correlation is ", round(fit$indices["Maximum Residual", 1], 3), sep="")
  cat("\n","The degrees of freedom for the model are ", efa$modelInfo$df,
      " and the objective function was ", round(ObjFn, 2), "\n", sep="")
  if(!is.null(nobs)) {
    cat("The total number of observations was ", nobs, " with corrected Chi-squared = ",
        round(fit$indices["Chi-square", 2], 1), " with prob < ", fit$indices["p-value", 2], sep="")
    cat("\n","Corrected Tucker Lewis Index of factoring reliability = ", round(fit$indices["TLI", 2], 3), sep="")
    cat("\n","Corrected RMSEA index = ", round(fit$indices["RMSEA", 2], 3), sep="")
    cat("\n","Corrected AIC = ", round(fit$information["AIC", 2], 1), sep="")
    cat("\n","Corrected BIC = ", round(fit$information["BIC", 2], 1), sep="")
    cat("\n","Corrected HQ = ", round(fit$information["HQ", 2], 1), "\n", sep="")
  }
  # Loadings
  cat("Standardized loadings (pattern matrix)\n", sep = ""); print(round(lambda, 2))
  invisible(NULL)

}

print.bifactor <- function(x, nobs = NULL, ...) {

  efa <- x
  # Check if nobs was provided
  if(is.null(nobs)) {
    if(is.null(efa$modelInfo$nobs) | isTRUE(efa$modelInfo$nobs == 0L)) {
      warning("Sample size was not provided. Some Chi-squared-based statistics will not be computed.")
      nobs <- NA
    } else {
      nobs <- efa$modelInfo$nobs
    }
  }

  ### Pattern matrix with communalities, uniqueness, and complexity
  lambda <- efa$bifactor$lambda
  phi <- efa$bifactor$phi

  ObjFn <- efa$efa$f
  ordering <- order(diag(phi %*% t(lambda) %*% lambda), decreasing=T)
  fit <- suppressWarnings(fit(efa, nobs))
  phi <- phi[ordering, ordering]
  rownames(phi) <- colnames(phi) <- paste("F",sprintf(paste("%0",nchar(efa$modelInfo$nfactors),"d",sep=""),
                                                      ordering),sep="")
  colnames(lambda) <- rownames(phi) <-
    colnames(phi) <- paste("F",sprintf(paste("%0",nchar(ncol(lambda)),"d",sep=""),1:ncol(lambda)),sep="")
  if(is.null(colnames(efa$modelInfo$correlation))) {
    rownames(lambda) <- paste("item_",sprintf(paste("%0",nchar(nrow(lambda)),"d",sep=""),1:nrow(lambda)),sep="")
  } else {
    rownames(lambda) <- colnames(efa$modelInfo$correlation)
  }
  # Print
  cat("Factor Analysis using estimator = ", efa$modelInfo$estimator, "\n", sep="")
  cat("Elapsed time of ", round(efa$elapsed * 1e-9, 3), " seconds", "\n", sep="")
  cat("\n","Goodness-of-fit and model misfit indices", sep="")
  cat("\n","The standardized root mean square residual (SRMR) is ", round(fit$indices["SRMR", 1], 3), sep="")
  cat("\n","The largest absolute value of standardized residual correlation is ", round(fit$indices["Maximum Residual", 1], 3), sep="")
  cat("\n","The degrees of freedom for the model are ", efa$modelInfo$df,
      " and the objective function was ", round(ObjFn, 2), "\n", sep="")
  if(!is.null(nobs)) {
    cat("The total number of observations was ", nobs, " with corrected Chi-squared = ",
        round(fit$indices["Chi-square", 2], 1), " with prob < ", fit$indices["p-value", 2], sep="")
    cat("\n","Corrected Tucker Lewis Index of factoring reliability = ", round(fit$indices["TLI", 2], 3), sep="")
    cat("\n","Corrected RMSEA index = ", round(fit$indices["RMSEA", 2], 3), sep="")
    cat("\n","Corrected AIC = ", round(fit$information["AIC", 2], 1), sep="")
    cat("\n","Corrected BIC = ", round(fit$information["BIC", 2], 1), sep="")
    cat("\n","Corrected HQ = ", round(fit$information["HQ", 2], 1), "\n", sep="")
  }
  # Loadings
  cat("Standardized loadings (pattern matrix)\n", sep = ""); print(round(lambda, 2))
  invisible(NULL)

}
