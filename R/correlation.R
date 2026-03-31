# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 29/03/2026

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

split_by_missing_pattern <- function(data) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data.frame or matrix.")
  }

  p <- ncol(data)

  if (nrow(data) == 0L) {
    return(list())
  }

  miss <- is.na(data)

  # One key per row, based on its missing-data pattern
  pattern_key <- apply(miss, 1L, function(x) paste(as.integer(x), collapse = ""))

  # Keep patterns in order of first appearance
  pattern_levels <- unique(pattern_key)

  out <- vector("list", length(pattern_levels))

  for (k in seq_along(pattern_levels)) {
    idx_rows <- which(pattern_key == pattern_levels[k])
    idx_vars <- !miss[idx_rows[1L], ]   # TRUE = observed in this pattern

    out[[k]] <- list(
      data = data[idx_rows, idx_vars, drop = FALSE],
      vars = idx_vars,
      nobs = length(idx_rows)
    )
  }

  names(out) <- NULL
  out
}

correlation <- function(data, item_names = colnames(data),
                        cor = "pearson", estimator = "ml",
                        acov = "standard", nobs = NULL,
                        missing = "pairwise.complete.obs",
                        likelihood = NULL) {

  result <- list()

  cor <- tolower(cor)
  estimator <- tolower(estimator)
  acov <- tolower(acov)
  missing <- tolower(missing)
  X <- as.matrix(data[, item_names])
  p <- nrow(X) # Number of rows
  q <- ncol(X) # Number of columns

  if(is.null(likelihood)) {
    if(estimator == "ml") {
      likelihood <- "normal"
    } else {
      likelihood <- "wishart"
    }
  }

  #### Compute the covariance and asymptotic covariance matrix ####

  if(p < q) {

    stop("Please provide either a full-rank matrix of scores or a covariance matrix")

  } else if(p == q) {

    if(is.null(nobs)) {
      stop("A covariance matrix was provided but nobs is missing. Standard errors and some statistics will not be computed")
    }

    result$R <- X
    result$ACOV <- asymptotic_normal(result$R)

  } else if(cor == "poly" || cor == "polys" ||
            cor == "polychoric" || cor == "polychorics") {

    polychorics <- polyfast(X)
    result$R <- polychorics$correlation
    result$thresholds <- polychorics$thresholds
    result$thresholds <- lapply(result$thresholds, FUN = \(x) x[-c(1, length(x))])
    result$cumprop <- polychorics$cumulative_freqs
    result$contingency_tables <- polychorics$contingency_tables

    result$ACOV <- diag(1/c(polychorics$hess))
    # result$ACOV <- asymptotic_poly(X, taus = thresholds)

  } else if(cor == "pearson") {

    if(is.null(nobs)) nobs <- p
    if(missing == "fiml") {
      result$R <- stats::cor(X, use = "pairwise.complete.obs")
    } else {
      result$R <- stats::cor(X, use = missing)
    }

    if(acov == "standard") {

      result$ACOV <- asymptotic_normal(result$R)

    } else if(acov == "robust") {

      result$ACOV <- asymptotic_general(X)

    } else {
      stop("Unknown acov argument")
    }

    if(likelihood == "normal") {

      result$R <- result$R*(nobs-1L)/nobs # SQRT??
      result$ACOV <- result$ACOV*(nobs-1L)/nobs

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

    stop("Unknown estimator")

  }

  #### FIML ####

  if(missing == "fiml") {

    result$missing <- split_by_missing_pattern(X)
    estimator <- "ml"
    missing <- "pairwise.complete.obs"

    for(i in 1:length(result$missing)) {

      nobsij <- result$missing[[i]]$nobs
      item_namesij <- colnames(result$missing[[i]]$data)
      result$missing[[i]] <- correlation(data = result$missing[[i]]$data,
                                         item_names = colnames(result$missing[[i]]$data),
                                         cor = cor, estimator = estimator,
                                         acov = acov, nobs = result$missing[[i]]$nobs,
                                         missing = missing,
                                         likelihood = likelihood)

      result$missing[[i]]$nobs <- nobsij
      result$missing[[i]]$item_names <- item_namesij

    }

  }

  #### Return ####

  return(result)

}
