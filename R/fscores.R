#' @title
#' Compute factor scores
#' @description
#'
#' Compute the factor scores from an exploratory factor model
#'
#' @usage
#'
#' fscores(fit, scores = NULL, method = "regression")
#'
#' @param fit object of class `efa` or `bifactor`.
#' @param scores Matrix of raw scores.
#' @param method Method to compute the factor scores.
#'
#' @details ...
#'
#' @return List...
#'
#' @author
#'
#' Marcos Jiménez & Vithor R. Franco
#'
#' @export
fscores <- function(fit, scores = NULL, method = "regression") {

  if(is.null(scores)) stop("Please, provide the matrix of observed scores")

  R <- fit$modelInfo$correlation
  invR <- solve(R)
  n <- nrow(scores)
  p <- fit$modelInfo$nvars
  q <- fit$modelInfo$nfactors
  z <- scale(scores)

  if(fit$modelInfo$rotation == "none"){

    Lambda <- fit$efa$lambda
    phi <- diag(q)

  } else {

    if(inherits(fit, "efa")) {
      Lambda <- fit$rotation$lambda
      phi <- fit$rotation$phi
    } else if(inherits(fit, "bifactor")){
      Lambda <- fit$bifactor$lambda
      phi <- fit$bifactor$phi
    }

  }

  S <- Lambda %*% phi # Correlations between factors and items

  if(method == "regression") {

    weights <- solve(R, S)

  } else if(method == "tenBerge") {

    SVD <- svd(phi)
    phi12 <- SVD$u %*% diag(sqrt(SVD$d)) %*% t(SVD$v)
    SVD <- svd(R)
    R12 <- SVD$u %*% diag(1/sqrt(SVD$d)) %*% t(SVD$v)
    L <- Lambda %*% phi12
    SVD <- svd(t(L) %*% invR %*% L)
    LRL12 <- SVD$u %*% diag(1/sqrt(SVD$d)) %*% t(SVD$v)
    C <- R12 %*% L %*% LRL12
    weights <- R12 %*% C %*% phi12

  } else if(method == "Bartlett") {

    U <- c(fit$efa$uniquenesses)
    U2 <- diag(1/(U*U))
    weights <- U2 %*% Lambda %*% solve(t(Lambda) %*% U2 %*% Lambda)

  } else if(method == "Harman") {

    weights <- solve(fit$efa$Rhat) %*% Lambda

  }

  fs <- z %*% weights # Factor scores

  # Validity coefficients:
  # invL <- diag(1/apply(fs, MARGIN = 2, FUN = sd))
  C <- t(weights) %*% R %*% weights
  invL <- diag(sqrt(diag(C))) # Standard deviations of the factor scores
  validity_univocality <- t(S) %*% weights %*% invL
  validity <- matrix(diag(validity_univocality), nrow = 1)
  rownames(validity) <- ""
  univocality <- validity_univocality
  diag(univocality) <- NA

  # Accuracy:
  accuracy <- stats::cor(fs)

  # standard errors for factor scores:
  r <- matrix(diag(invR), ncol = 1)
  Rj <- matrix(1-c(validity^2), nrow = 1)
  se <- sqrt(r %*% Rj / (n-p-1))
  se <- matrix(se, nrow = p, ncol = q)

  colnames(fs) <- colnames(weights) <- colnames(validity) <-
    colnames(univocality) <- rownames(univocality) <-
    colnames(accuracy) <- rownames(accuracy) <- colnames(se) <-
    paste("F", sprintf(paste("%0",nchar(q),"d",sep=""), 1:q),sep="")

  result <- list(fscores = fs, weights = weights, validity = validity,
                 univocality = univocality, accuracy = accuracy, se = se)

  class(result) <- "scores"

  return(result)

}
