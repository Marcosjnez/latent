# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 20/04/2026

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

lcov <- function(data, item_names = colnames(data),
                 cor = "pearson", estimator = "ml",
                 acov = "standard", nobs = NULL,
                 missing = "pairwise.complete.obs",
                 std.ov = FALSE, likelihood = NULL,
                 meanstructure = TRUE) {

  cor <- tolower(cor)
  estimator <- tolower(estimator)
  acov <- tolower(acov)
  missing <- tolower(missing)

  X <- as.matrix(data[, item_names, drop = FALSE])

  compute_one_pattern <- function(X, item_names, cor, estimator, acov, nobs,
                                  missing, std.ov, likelihood, meanstructure,
                                  vars = rep(TRUE, length(item_names))) {

    out <- list()
    out$item_names <- item_names
    out$vars <- vars
    out$nobs <- if (is.null(nobs)) nrow(X) else nobs

    p <- nrow(X)
    q <- ncol(X)

    if(is.null(likelihood)) {
      likelihood <- if (estimator %in% c("ml", "fml", "means_fml")) "normal" else "wishart"
    }

    #### Compute covariance/correlation and ACOV ####

    if(cor %in% c("poly", "polys", "polychoric", "polychorics")) {

      polychorics <- polyfast(X)
      out$S <- polychorics$correlation
      out$thresholds <- lapply(polychorics$thresholds,
                               FUN = function(x) x[-c(1, length(x))])
      rownames(out$S) <- colnames(out$S) <- names(out$thresholds) <- out$item_names
      out$cumprop <- polychorics$cumulative_freqs
      out$contingency_tables <- polychorics$contingency_tables
      out$ACOV <- diag(1 / c(polychorics$hess))

      # fit_poly <- lpoly(data = X,
      #                   model = NULL,
      #                   method = "none", penalties = FALSE,
      #                   do.fit = TRUE)
      # out$ACOV <- asymptotic_poly(fit_poly, model = NULL)

    } else if(cor == "pearson") {

      if(nobs < 2) {
        out$S <- t(X) %*% X
        acov <- "standard"
        likelihood <- "wishart"
      } else {
        out$S <- stats::cov(X, use = missing)
      }
      rownames(out$S) <- colnames(out$S) <- out$item_names

      if(std.ov) {
        inv_sqrtdiagS <- diag(1/sqrt(diag(out$S)))
        out$S <- inv_sqrtdiagS %*% out$S %*% inv_sqrtdiagS
      }

      if(likelihood == "normal") {
        out$S <- out$S * (nobs - 1L) / nobs
        # out$ACOV <- out$ACOV * (nobs - 1L) / nobs
      }

      if(acov == "standard") {
        out$ACOV <- asymptotic_normal(out$S, cov = !std.ov)
      } else if (acov == "robust") {
        out$ACOV <- asymptotic_general(X, cov = !std.ov)
      } else {
        stop("Unknown `acov` argument")
      }

    } else {
      stop("Unknown covariance/correlation method")
    }

    #### Compute weight matrices ####

    if (estimator %in% c("uls", "means_uls", "ml", "fml", "means_fml")) {

      out$W <- matrix(1, nrow = q, ncol = q)
      out$w_means <- rep(1, times = q)

    } else if (estimator %in% c("dwls", "means_dwls")) {

      if(cor %in% c("poly", "polys", "polychoric", "polychorics")) {
        fix_diag <- TRUE
      } else {
        fix_diag <- FALSE
      }

      out$W <- matrix(NA_real_, nrow = q, ncol = q)
      out$W[lower.tri(out$W, diag = !fix_diag)] <- diag(out$ACOV)
      out$W[upper.tri(out$W)] <- t(out$W)[upper.tri(out$W)]
      out$W <- 1 / out$W
      if(std.ov) diag(out$W) <- 1

      var_means <- apply(X, MARGIN = 2, FUN = var, na.rm = TRUE)
      out$w_means <- 1/var_means

    } else {
      stop("Unknown `estimator`")
    }

    #### Mean structure ####

    if(meanstructure) {

      if(std.ov) {
        out$means <- rep(0, times = q)
      } else {
        out$means <- colMeans(X, na.rm = TRUE)
      }

      if(out$nobs < 2) {
        acov_means <- rep(0, times = length(out$vars))
      } else {
        acov_means <- apply(X, MARGIN = 2, FUN = var, na.rm = TRUE)
      }
      out$ACOV <- block_diag(list(diag(acov_means), out$ACOV))

    }

    out$NACOV <- out$ACOV * nobs

    #### Return ####

    return(out)

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
        std.ov = std.ov,
        likelihood = likelihood,
        meanstructure = meanstructure,
        vars = vars_i
      )
    })

    names(out) <- paste0("pattern", seq_along(out))
    out$npatterns <- length(patterns)
    out$patterns_names <- names(patterns)
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
      std.ov = std.ov,
      meanstructure = meanstructure,
      vars = rep(TRUE, length(item_names))
    )
  )

  out$npatterns <- 1L
  out$patterns_names <- "pattern1"
  # out$nobs <- nobs

  #### Return ####

  return(out)

}

# if (p < q) {
#
#   stop("Please provide either a full-rank matrix of scores or a covariance matrix.")
#
# } else if (p == q) {
#
#   if (is.null(nobs)) {
#     stop("A covariance matrix was provided but `nobs` is missing.")
#   }
#
#   out$S <- X
#   if(std.ov) {
#     inv_sqrtdiagS <- diag(1/sqrt(diag(out$S)))
#     out$S <- inv_sqrtdiagS %*% out$S %*% inv_sqrtdiagS
#   }
#   out$ACOV <- asymptotic_normal(out$S, cov = !std.ov)
#
# }
