# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 07/09/2025

correlation <- function(data, item_names = colnames(data),
                        cor = "pearson", estimator = "uls",
                        group = NULL, missing = "pairwise.complete.obs") {

  result <- list()

  X <- as.matrix(data[, item_names])
  p <- nrow(X) # Number of rows
  q <- ncol(X) # Number of columns

  # Get the correlation matrix:
  if(cor == "pearson") {

    if(p == q) {

      R <- X
      result$R <- R

      if(estimator == "uls") {
        W <- matrix(1, nrow = q, ncol = q)
        # diag(W) <- 0
        result$W <- W
      } else if(estimator == "dwls") {
        asymp <- asymptotic_normal(R)
        W <- matrix(diag(asymp), nrow = q, ncol = q)
        W <- 1 / W
        # diag(W) <- 0
        result$W <- W
      }

    } else {

      R <- stats::cor(X, use = missing)
      result$R <- R

      if(estimator == "uls") {
        W <- matrix(1, nrow = q, ncol = q)
        # diag(W) <- 0
        result$W <- W
      } else if(estimator == "dwls") {
        asymp <- asymptotic_general(X)
        W <- matrix(diag(asymp), nrow = q, ncol = q)
        W <- 1 / W
        # diag(W) <- 0
        result$W <- W
      }

    }

  } else if(cor == "poly") {

    if(p == q) {

      R <- X
      result$R <- R

      if(estimator == "uls") {
        W <- matrix(1, nrow = q, ncol = q)
        # diag(W) <- 0
        result$W <- W
      } else if(estimator == "dwls") {
        warning("The full data was not provided. The variance of the polychoric
              correlations will be approximated")
        asymp <- asymptotic_normal(R)
        W <- matrix(diag(asymp), nrow = q, ncol = q)
        W <- 1 / W
        # diag(W) <- 0
        result$W <- W
      }

    } else {

      polychorics <- polyfast(X)
      R <- polychorics$correlation
      taus <- polychorics$thresholds
      cumprop <- polychorics$cumulative_freqs
      n <- polychorics$contingency_tables
      result$R <- R
      result$thresholds <- taus
      result$cumprop <- cumprop
      results$contingency_tables <- n
      if(any(eigen(R)$values < 0)) {
        warning("The matrix of polychoric correlations was non-positive
                definite with the two-step method. \n
                The closest positive semidefinite solution was estimated
                instead.")
        newfit <- poly_positive(R = R, taus = taus, n = n,
                                control = list(eps = 1e-04))
        R <- matrix(newfit$fit$outputs$estimators$matrices[[1]][[1]], p, p)
        taus <- newfit$fit$outputs$estimators$list_vectors[[1]][[1]]
        cumprop <- newfit$fit$outputs$estimators$list_vectors[[1]][[2]]
        result$R <- R
        result$thresholds <- taus
        result$cumprop <- cumprop
      }

      if(estimator == "uls") {
        W <- matrix(1, nrow = q, ncol = q)
        # diag(W) <- 0
        result$W <- W
      } else if(estimator == "dwls") {
        W <- 1 / DACOV2(p, R, polychorics$contingency_tables,
                        polychorics$thresholds, polychorics$cumulative_freqs)
        # diag(W) <- 0
        result$W <- W
      }
    }

  }

  return(result)

}
