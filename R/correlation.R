# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 11/09/2025

asymptotic_poly <- function(X, taus) {

  fit_poly <- lpoly(data = X, do.fit = FALSE,
                    model = list(taus = taus),
                    method = "crossprodn",
                    penalties = FALSE,
                    control = NULL)

  control_manifold <- fit_poly@modelInfo$control_manifold
  control_transform <- fit_poly@modelInfo$control_transform
  control_estimator <- fit_poly@modelInfo$control_estimator
  control_optimizer <- fit_poly@modelInfo$control

  x <- get_hess(control_manifold, control_transform,
                control_estimator, control_optimizer)
  ACOV <- solve(x$h)

  return(ACOV)

}

correlation <- function(data, item_names = colnames(data),
                        cor = "pearson", estimator = "ml",
                        acov = "standard", nobs = NULL,
                        missing = "pairwise.complete.obs") {

  result <- list()

  cor <- tolower(cor)
  estimator <- tolower(estimator)
  acov <- tolower(acov)
  X <- as.matrix(data[, item_names])
  p <- nrow(X) # Number of rows
  q <- ncol(X) # Number of columns

  #### Compute the covariance and asymptotic covariance matrix ####

  if(p < q) {

    stop("Please provide either a full-rank matrix of scores of a covariance matrix")

  } else if(p == q) {

    if(is.null(nobs)) {
      stop("A covariance matrix was provided but nobs is missing")
    }

    result$R <- X
    result$ACOV <- asymptotic_normal(result$R)

  } else if(cor == "poly" || cor == "polys" ||
            cor == "polychoric" || cor == "polychorics") {

    polychorics <- polyfast(X)
    result$R <- polychorics$correlation
    result$thresholds <- polychorics$thresholds
    result$cumprop <- polychorics$cumulative_freqs
    result$contingency_tables <- polychorics$contingency_tables

    result$ACOV <- diag(1/c(polychorics$hess))#*nrow(X)
    # result$ACOV <- asymptotic_poly(X, taus = result$thresholds)
    # result$ACOV <- diag(c(DACOV2(p, result$R,
    #                              polychorics$contingency_tables,
    #                              polychorics$thresholds,
    #                              polychorics$cumulative_freqs)))

  } else if(cor == "pearson") {

    result$R <- stats::cor(X, use = missing)

    if(acov == "standard") {

      result$ACOV <- asymptotic_normal(result$R)

    } else if(acov == "robust") {

      result$ACOV <- asymptotic_general(X)

    } else {
      stop("Unknown acov argument")
    }

  } else {

    stop("Unknown covariance matrix method")

  }

  #### Compute the weight matrix ####

  if(estimator == "uls" || estimator == "ulsr" ||
     estimator == "ml" || estimator == "mlr") {

    result$W <- matrix(1, nrow = q, ncol = q)
    diag(result$W) <- 1

  } else if(estimator == "dwls" || estimator == "dwlsr") {

    result$W <- matrix(NA, nrow = q, ncol = q)
    result$W[lower.tri(result$W, diag = FALSE)] <- diag(result$ACOV)
    result$W[upper.tri(result$W)] <- t(result$W)[upper.tri(result$W)]
    result$W <- 1 / result$W
    diag(result$W) <- 1

  } else {

    stop("Unknown esitmator")

  }

  #### Return ####

  return(result)

}
