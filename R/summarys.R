summary.efa <- function(object, nobs = NULL, suppress = 0, order = FALSE, digits = 2, ...) {

  efa <- object
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
  Rhat <- efa$efa$Rhat
  R <- efa$modelInfo$correlation

  uniquenesses <- c(efa$efa$uniquenesses)
  ObjFn <- efa$efa$f
  SSloads <- diag(phi %*% t(lambda) %*% lambda) # Generalizing to oblique rotation
  ordering <- order(SSloads, decreasing=T)
  colnames(lambda) <- paste("F",sprintf(paste("%0",nchar(ncol(lambda)),"d",sep=""),1:ncol(lambda)),sep="")
  lambda <- lambda[,ordering]
  if(is.null(colnames(efa$modelInfo$correlation))) {
    rownames(lambda) <- paste("item_",sprintf(paste("%0",nchar(nrow(lambda)),"d",sep=""),1:nrow(lambda)),sep="")
  } else {
    rownames(lambda) <- colnames(efa$modelInfo$correlation)
  }
  lambda <- lambda * {abs(lambda) > abs(suppress)}
  if(order) {
    orderItems <- apply(abs(lambda), 1, which.max)
    lambda  <- lambda[order(orderItems, decreasing=FALSE),]
    uniquenesses <- uniquenesses[order(orderItems, decreasing=FALSE)]
  }
  h2     <- 1 - uniquenesses
  u2     <- uniquenesses
  com    <- {rowSums(lambda ^ 2)^2}/rowSums(lambda ^ 4)
  names(h2) <- names(u2) <- names(com) <- rownames(lambda)
  loadsM <- round(data.frame(lambda, h2, u2, com), digits)
  loadsM[, 1:ncol(lambda)][loadsM[, 1:ncol(lambda)] == 0] <- ""

  ### Variance accounted for
  SSloads <- sort(SSloads, decreasing=T)
  propVar <- SSloads/nrow(lambda)
  cumsVar <- cumsum(propVar)
  propExp <- propVar/sum(propVar)
  cumsExp <- cumsum(propExp)
  VAF <- round(t(data.frame(SSloads, propVar, cumsVar, propExp, cumsExp)),2)
  rownames(VAF) <- c("SS loadings", "Proportion Var", "Cumulative Var",
                     "Proportion Explained", "Cumulative Proportion")
  colnames(VAF) <- colnames(lambda)

  ### Factor correlations
  phi <- phi[ordering, ordering]
  rownames(phi) <- colnames(phi) <- colnames(lambda)

  ### Factor score indeterminacy
  Reliability <- diag(phi %*% t(lambda) %*% solve(Rhat) %*% lambda %*% phi)
  Indeterminacy <- sqrt(Reliability)
  min_cor <- 2*Reliability-1
  RELIABILITY <- round(t(data.frame(Reliability, Indeterminacy, min_cor)), digits)
  rownames(RELIABILITY)[3] <- "Minimum Correlation"
  colnames(RELIABILITY) <- colnames(lambda)

  ### Fit statistics
  fit  <- suppressWarnings(fit(efa, nobs))

  ### Print
  # Basic info
  cat("Factor Analysis using estimator = ", efa$modelInfo$estimator, "\n", sep="")
  cat("Elapsed time of ", round(efa$elapsed * 1e-9, 3), " seconds", "\n", sep="")
  # Loadings
  cat("Standardized loadings (pattern matrix) based upon correlation matrix\n", sep=""); print(loadsM)
  # Variance accounted for
  cat("\n","Variance accounted for after rotation\n",sep=""); print(VAF)
  # Latent correlations
  cat("\n","Factor correlations after rotation\n",sep=""); print(round(phi, digits))
  # Factor Indeterminacy
  cat("\n","Factor score indeterminacy\n",sep=""); print(round(RELIABILITY, digits))
  # Fit
  cat("\n","Goodness-of-fit and model misfit indices", sep="")
  cat("\n","Mean item complexity = ", round(mean(com),1), sep="")
  cat("\n","The standardized root mean square residual (SRMR) is ", round(fit$indices["SRMR", 2], 3), sep="")
  cat("\n","The largest absolute value of standardized residual correlation is ", round(fit$indices["Maximum Residual", 2], 3), sep="")
  cat("\n","The degrees of freedom for the null model are ", efa$modelInfo$df_null,
      " and the objective function was ", round(efa$modelInfo$f_null,2), sep="")
  cat("\n","The degrees of freedom for the model are ", efa$modelInfo$df,
      " and the objective function was ", round(ObjFn, 2), "\n", sep="")
  if(!is.null(nobs)) {
    cat("The total number of observations was ", nobs, " with Corrected Chi-squared = ",
        round(fit$indices["Chi-square", 2], 1), " with prob < ", fit$indices["p-value", 2] , sep="")
    cat("\n","Corrected Tucker Lewis Index of factoring reliability = ", round(fit$indices["TLI", 2], 3), sep="")
    cat("\n","Corrected RMSEA index = ", round(fit$indices["RMSEA", 2], 3), sep="")
    cat("\n","Corrected AIC = ", round(fit$information["AIC", 2], 1), sep="")
    cat("\n","Corrected BIC = ", round(fit$information["BIC", 2], 1), sep="")
    cat("\n","Corrected HQ = ", round(fit$information["HQ", 2], 1), "\n", sep="")
  }

  Results <- list("loadings"=lambda, "communalities"=h2, "uniqueness"=u2, "complexity"=com,
                  "VAF"=VAF, "phi"=phi, "fit"=fit)
  invisible(Results)

}

summary.bifactor <- function(object, nobs = NULL, suppress = 0, order = FALSE, digits = 2, ...) {

  efa <- object
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
  Rhat <- efa$efa$Rhat
  R <- efa$modelInfo$correlation

  uniquenesses <- c(efa$efa$uniquenesses)
  ObjFn <- efa$efa$f
  SSloads <- diag(phi %*% t(lambda) %*% lambda) # Generalizing to oblique rotation
  ordering <- order(SSloads, decreasing=T)
  colnames(lambda) <- paste("F",sprintf(paste("%0",nchar(ncol(lambda)),"d",sep=""),1:ncol(lambda)),sep="")
  lambda <- lambda[,ordering]
  if(is.null(colnames(efa$modelInfo$correlation))) {
    rownames(lambda) <- paste("item_",sprintf(paste("%0",nchar(nrow(lambda)),"d",sep=""),1:nrow(lambda)),sep="")
  } else {
    rownames(lambda) <- colnames(efa$modelInfo$correlation)
  }
  lambda <- lambda * {abs(lambda) > abs(suppress)}
  if(order) {
    orderItems <- apply(abs(lambda), 1, which.max)
    lambda  <- lambda[order(orderItems, decreasing=FALSE),]
    uniquenesses <- uniquenesses[order(orderItems, decreasing=FALSE)]
  }
  h2     <- 1 - uniquenesses
  u2     <- uniquenesses
  com    <- {rowSums(lambda ^ 2)^2}/rowSums(lambda ^ 4)
  names(h2) <- names(u2) <- names(com) <- rownames(lambda)
  loadsM <- round(data.frame(lambda, h2, u2, com), digits)
  loadsM[, 1:ncol(lambda)][loadsM[, 1:ncol(lambda)] == 0] <- ""

  ### Variance accounted for
  SSloads <- sort(SSloads, decreasing=T)
  propVar <- SSloads/nrow(lambda)
  cumsVar <- cumsum(propVar)
  propExp <- propVar/sum(propVar)
  cumsExp <- cumsum(propExp)
  VAF <- round(t(data.frame(SSloads, propVar, cumsVar, propExp, cumsExp)),2)
  rownames(VAF) <- c("SS loadings", "Proportion Var", "Cumulative Var",
                     "Proportion Explained", "Cumulative Proportion")
  colnames(VAF) <- colnames(lambda)

  ### Factor correlations
  phi <- phi[ordering, ordering]
  rownames(phi) <- colnames(phi) <- colnames(lambda)

  ### Factor score indeterminacy
  Reliability <- diag(phi %*% t(lambda) %*% solve(Rhat) %*% lambda %*% phi)
  Indeterminacy <- sqrt(Reliability)
  min_cor <- 2*Reliability-1
  RELIABILITY <- round(t(data.frame(Reliability, Indeterminacy, min_cor)), digits)
  rownames(RELIABILITY)[3] <- "Minimum Correlation"
  colnames(RELIABILITY) <- colnames(lambda)

  ### Fit statistics
  fit  <- suppressWarnings(fit(efa, nobs))

  ### Print
  # Basic info
  cat("Factor Analysis using estimator = ", efa$modelInfo$estimator, "\n", sep="")
  cat("Elapsed time of ", round(efa$elapsed * 1e-9, 3), " seconds", "\n", sep="")
  # Loadings
  cat("Standardized loadings (pattern matrix) based upon correlation matrix\n", sep=""); print(loadsM)
  # Variance accounted for
  cat("\n","Variance accounted for after bifactor rotation\n",sep=""); print(VAF)
  # Latent correlations
  cat("\n","Factor correlations after bifactor rotation\n",sep=""); print(round(phi, digits))
  # Factor Indeterminacy
  cat("\n","Factor score indeterminacy\n",sep=""); print(round(RELIABILITY, digits))
  # Fit
  cat("\n","Goodness-of-fit and model misfit indices", sep="")
  cat("\n","Mean item complexity = ", round(mean(com),1), sep="")
  cat("\n","The standardized root mean square residual (SRMR) is ", round(fit$indices["SRMR", 2], 3), sep="")
  cat("\n","The largest absolute value of standardized residual correlation is ", round(fit$indices["Maximum Residual", 2], 3), sep="")
  cat("\n","The degrees of freedom for the null model are ", efa$modelInfo$df_null,
      " and the objective function was ", round(efa$modelInfo$f_null,2), sep="")
  cat("\n","The degrees of freedom for the model are ", efa$modelInfo$df,
      " and the objective function was ", round(ObjFn, 2), "\n", sep="")
  if(!is.null(nobs)) {
    cat("The total number of observations was ", nobs, " with Corrected Chi-squared = ",
        round(fit$indices["Chi-square", 2], 1), " with prob < ", fit$indices["p-value", 2] , sep="")
    cat("\n","Corrected Tucker Lewis Index of factoring reliability = ", round(fit$indices["TLI", 2], 3), sep="")
    cat("\n","Corrected RMSEA index = ", round(fit$indices["RMSEA", 2], 3), sep="")
    cat("\n","Corrected AIC = ", round(fit$information["AIC", 2], 1), sep="")
    cat("\n","Corrected BIC = ", round(fit$information["BIC", 2], 1), sep="")
    cat("\n","Corrected HQ = ", round(fit$information["HQ", 2], 1), "\n", sep="")
  }

  Results <- list("loadings"=lambda, "communalities"=h2, "uniqueness"=u2, "complexity"=com,
                  "VAF"=VAF, "phi"=phi, "fit"=fit)
  invisible(Results)

}
