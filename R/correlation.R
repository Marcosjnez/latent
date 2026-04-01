# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 01/04/2026

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

  if (nrow(data) == 0L) {
    return(list())
  }

  miss <- is.na(data)

  # One key per row, based on its missing-data pattern
  pattern_key <- apply(miss, 1L, function(x) paste(as.integer(x), collapse = ""))

  # Unique patterns in first-appearance order
  pattern_levels <- unique(pattern_key)

  # Number of rows per pattern
  counts <- tabulate(match(pattern_key, pattern_levels), nbins = length(pattern_levels))

  # Order by decreasing number of rows; ties keep first appearance order
  ord <- order(-counts, seq_along(pattern_levels))
  pattern_levels <- pattern_levels[ord]
  counts <- counts[ord]

  out <- vector("list", length(pattern_levels))

  for (k in seq_along(pattern_levels)) {
    idx_rows <- which(pattern_key == pattern_levels[k])
    idx_vars <- !miss[idx_rows[1L], ]  # TRUE = observed in this pattern

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

  cor <- tolower(cor)
  estimator <- tolower(estimator)
  acov <- tolower(acov)
  missing <- tolower(missing)

  X <- as.matrix(data[, item_names, drop = FALSE])

  compute_one_pattern <- function(X, item_names, cor, estimator, acov, nobs,
                                  missing, likelihood, vars = rep(TRUE, length(item_names))) {

    out <- list()
    out$item_names <- item_names
    out$vars <- vars
    out$nobs <- if (is.null(nobs)) nrow(X) else nobs

    p <- nrow(X)  # rows
    q <- ncol(X)  # cols

    if (is.null(likelihood)) {
      likelihood <- if (estimator == "ml") "normal" else "wishart"
    }

    #### Compute covariance/correlation and ACOV ####

    if (p < q) {
      stop("Please provide either a full-rank matrix of scores or a covariance matrix.")

    } else if (p == q) {

      if (is.null(nobs)) {
        stop("A covariance matrix was provided but `nobs` is missing.")
      }

      out$R <- X
      out$ACOV <- asymptotic_normal(out$R)

    } else if (cor %in% c("poly", "polys", "polychoric", "polychorics")) {

      polychorics <- polyfast(X)
      out$R <- polychorics$correlation
      out$thresholds <- lapply(polychorics$thresholds, function(x) x[-c(1, length(x))])
      out$cumprop <- polychorics$cumulative_freqs
      out$contingency_tables <- polychorics$contingency_tables

      out$ACOV <- diag(1 / c(polychorics$hess))
      # out$ACOV <- asymptotic_poly(X, taus = out$thresholds)

    } else if (cor == "pearson") {

      if (is.null(nobs)) nobs <- p
      out$nobs <- nobs

      out$R <- stats::cor(X, use = missing)

      if (acov == "standard") {
        out$ACOV <- asymptotic_normal(out$R)
      } else if (acov == "robust") {
        out$ACOV <- asymptotic_general(X)
      } else {
        stop("Unknown `acov` argument.")
      }

      if (likelihood == "normal") {
        out$R <- out$R * (nobs - 1L) / nobs
        out$ACOV <- out$ACOV * (nobs - 1L) / nobs
      }

    } else {
      stop("Unknown covariance/correlation method.")
    }

    #### Compute weight matrix ####

    if (estimator %in% c("uls", "ulsr", "ml", "mlr")) {

      out$W <- matrix(1, nrow = q, ncol = q)
      diag(out$W) <- 1

    } else if (estimator %in% c("dwls", "dwlsr")) {

      out$W <- matrix(NA_real_, nrow = q, ncol = q)
      out$W[lower.tri(out$W, diag = FALSE)] <- diag(out$ACOV)
      out$W[upper.tri(out$W)] <- t(out$W)[upper.tri(out$W)]
      out$W <- 1 / out$W
      diag(out$W) <- 1

    } else {
      stop("Unknown `estimator`.")
    }

    out
  }

  #### FIML-style missing-pattern split ####

  if (missing == "fiml") {

    patterns <- split_by_missing_pattern(X)

    out <- lapply(seq_along(patterns), function(i) {
      pat <- patterns[[i]]

      Xi <- if (!is.null(pat$data)) pat$data else pat[[1L]]
      vars_i <- if (!is.null(pat$vars)) pat$vars else pat[[2L]]
      item_names_i <- colnames(Xi)
      nobs_i <- nrow(Xi)

      compute_one_pattern(
        X = Xi,
        item_names = item_names_i,
        cor = cor,
        estimator = "ml",
        acov = acov,
        nobs = nobs_i,
        missing = "pairwise.complete.obs",
        likelihood = likelihood,
        vars = vars_i
      )
    })

    names(out) <- paste0("pattern", seq_along(out))
    out$npatterns <- length(patterns)
    # out$nobs <- sum(unlist(lapply(out, FUN = \(x) x$nobs)))

    return(out)

  }

  #### Non-FIML ####

  out <- list(
    pattern1 = compute_one_pattern(
      X = X,
      item_names = item_names,
      cor = cor,
      estimator = estimator,
      acov = acov,
      nobs = nobs,
      missing = missing,
      likelihood = likelihood,
      vars = rep(TRUE, length(item_names))
    )
  )

  out$npatterns <- 1L
  # out$nobs <- nobs

  return(out)

}
