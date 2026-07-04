# Author: Mauricio Garnier-Villarreal
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 03/07/2026 by Marcos Jimenez
#'
#' Local Bivariate Residuals for Latent Class Analysis
#'
#' Computes bivariate residuals (BVR) and related diagnostics for a fitted latent
#' class model with mixed indicators (continuous and categorical). The function
#' returns a matrix of residuals, p-values, a correlation-like effect size, and a
#' summary table. The root mean square residual (RMSR) provides a global fit measure.
#'
#' @param model An object of class `"lca"` from the \pkg{latent} package.
#' @param digits Integer; number of decimal places to use in rounding the output
#'   matrices and the summary table. Default is 4.
#'
#' @return A list of class `"lbvr"` with components:
#'   \item{residual_matrix}{A symmetric matrix of residual measures:
#'     \itemize{
#'       \item For continuous–continuous pairs: Pearson correlation of residuals.
#'       \item For categorical–categorical pairs: maximum absolute standardized residual.
#'       \item For continuous–categorical pairs: maximum absolute Z‑score of mean residuals per category.
#'     }}
#'   \item{pvalue_matrix}{Symmetric matrix of p-values corresponding to the
#'     residual measures (using correlation test, chi‑square test, or ANOVA F‑test).}
#'   \item{r_mat}{Symmetric matrix of effect sizes comparable to a correlation:
#'     \itemize{
#'       \item Continuous–continuous: Pearson \eqn{r}.
#'       \item Categorical–categorical: difference of Cramér's V (\eqn{V_{obs} - V_{exp}}).
#'       \item Continuous–categorical: \eqn{\eta} (square root of eta‑squared from ANOVA).
#'     }}
#'   \item{rmsr}{Root mean square of the upper triangle of `r_mat`; a global fit index.}
#'   \item{res_tab}{A data frame with one row per variable pair, containing:
#'     \code{pair}, \code{statistic} (the residual measure), \code{p_value},
#'     \code{r} (the correlation‑like effect size), \code{var_types}, and
#'     \code{test} (a descriptive label). The table is sorted by absolute \code{r}.}
#'
#' @details
#' The function extracts data, variable types, posterior probabilities, class‑specific
#' means, and response probabilities from a fitted \pkg{latent} model. It then computes
#' residuals for each pair of indicators using only pairwise complete observations.
#' Missing data are handled by scaling expected frequencies to the number of complete
#' cases for the pair.
#'
#' For categorical variables, the original data (which may be zero‑based) are
#' incremented by 1 to ensure proper table construction.
#'
#' @examples
#' \dontrun{
#' library(latent)
#' # Fit a 3‑class model with multinomial indicators
#' fit <- lca(data = gss82, nclasses = 3,item = rep("multinomial", ncol(gss82)))
#' # Compute bivariate residuals
#' lbvr(fit)
#' }
#'
#' @importFrom car recode
#' @importFrom stats complete.cases cor.test pchisq sd aov
#' @export
lbvr <- function(model, digits = 4) {

  types <- model@dataList$variable_type
  posterior <- latInspect(model, "posterior")

  n <- model@dataList$nobs
  J <- model@dataList$nitems
  item_names <- names(types)
  K <- ncol(posterior)
  short2full <- model@dataList$short2full

  cont_vars <- model@dataList$gaussian$gaussian_names
  cat_vars  <- model@dataList$multinomial$multinomial_names
  data_gaussian <- model@dataList$gaussian$patterns_gaussian[short2full, , drop = FALSE]
  data_multinom <- model@dataList$multinomial$patterns_multinomial[short2full, , drop = FALSE] + 1L
  n_cont <- model@dataList$gaussian$ngaussian
  n_cat <- model@dataList$multinomial$nmultinomial

  profile <- latInspect(model, what = "item")
  class_means <- sapply(profile[cont_vars], FUN = function(x){ x[1, ] })
  class_vars <- sapply(profile[cont_vars], FUN = function(x){ x[2, ]^2 })

  class_probs <- latInspect(model, what = "item")[cat_vars]

  pi <- latInspect(model, what = "class")

  # Check and align class_means
  if (n_cont > 0) {
    if (is.null(class_means)) stop("class_means required for continuous variables")
    if (nrow(class_means) != K || ncol(class_means) != n_cont)
      stop("class_means must be K x n_cont matrix")
    if (is.null(colnames(class_means))) colnames(class_means) <- cont_vars
    if (!all(cont_vars %in% colnames(class_means)))
      stop("class_means columns must match continuous variable names")
  }

  # Check and align class_probs
  if (n_cat > 0) {
    if (is.null(class_probs)) stop("class_probs required for categorical variables")
    if (length(class_probs) != n_cat) stop("class_probs must have length n_cat")
    if (is.null(names(class_probs))) names(class_probs) <- cat_vars
    if (!all(cat_vars %in% names(class_probs)))
      stop("class_probs names must match categorical variable names")
  }

  # Output matrices
  resid_mat <- matrix(NA, J, J, dimnames = list(item_names, item_names))
  raw_resid_mat <- matrix(NA, J, J, dimnames = list(item_names, item_names))
  pval_mat <- matrix(NA, J, J, dimnames = list(item_names, item_names))
  r_mat <- matrix(NA, J, J, dimnames = list(item_names, item_names))
  names_mat <- matrix(NA, J, J, dimnames = list(item_names, item_names))
  types_mat <- matrix(NA, J, J, dimnames = list(item_names, item_names))

  details <- list()

  # ----- Helper functions -----

  cont_cont_resid_pval <- function(v1, v2) {
    e1 <- data_gaussian[, v1] - posterior %*% class_means[, v1, drop = FALSE]
    e2 <- data_gaussian[, v2] - posterior %*% class_means[, v2, drop = FALSE]
    comp <- complete.cases(e1, e2)
    if (sum(comp) < 3) return(list(resid = NA, pval = NA))
    ct <- cor.test(e1[comp], e2[comp], use = "complete.obs")
    raw_resid <- cov(e1[comp], e2[comp], use = "complete.obs")

    raw_resid <- latentgold_cont_cont_bvr(data_gaussian, v1, v2, posterior,
                                          class_means, class_vars,
                                          covariance = "class_specific",
                                          weights = NULL, tol = 1e-10)$resid


    list(resid = ct$estimate, pval = ct$p.value, rlike = ct$estimate,
         raw_resid = raw_resid)
  }

  cat_cat_resid_pval <- function(v1, v2) {
    comp <- complete.cases(data_multinom[, v1], data_multinom[, v2])
    if (sum(comp) == 0) return(list(resid = NA, pval = NA))
    obs_tab <- table(data_multinom[comp, v1], data_multinom[comp, v2])

    N_pair <- sum(comp)

    prob1 <- class_probs[[v1]]   # K x L1
    prob2 <- class_probs[[v2]]   # K x L2
    L1 <- nrow(prob1)
    L2 <- nrow(prob2)
    exp_tab <- matrix(0, L1, L2)
    for (r in 1:L1) {
      for (c in 1:L2) {
        exp_tab[r, c] <- sum(pi * prob1[r, ] * prob2[c, ]) * N_pair
      }
    }
    obs_cat1 <- as.numeric(rownames(obs_tab))
    obs_cat2 <- as.numeric(colnames(obs_tab))
    exp_tab_aligned <- exp_tab[obs_cat1, obs_cat2, drop = FALSE]

    # Standardized residuals and chi-square
    std_resid <- (obs_tab - exp_tab_aligned) / sqrt(exp_tab_aligned)
    std_resid[exp_tab_aligned == 0] <- NA
    chi2 <- sum(std_resid^2, na.rm = TRUE)
    df <- (nrow(obs_tab)-1) * (ncol(obs_tab)-1)
    pval <- pchisq(chi2, df, lower.tail = FALSE)

    resid_val <- max(abs(std_resid), na.rm = TRUE)

    cramer_v <- function(tab) {
      if (any(dim(tab) < 2)) return(NA)
      n <- sum(tab)
      chi2 <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
      df <- min(nrow(tab)-1, ncol(tab)-1)
      sqrt(chi2 / (n * df))
    }
    V_obs <- cramer_v(obs_tab)
    V_exp <- cramer_v(exp_tab_aligned)
    resid_cv <- V_obs - V_exp

    # raw residual LG:
    P <- (L1-1L) * (L2-1L)
    raw_resid <- sum((obs_tab - exp_tab_aligned)^2 / (exp_tab_aligned))/P

    list(resid = resid_val, pval = pval, rlike = resid_cv,
         raw_resid = raw_resid)
  }

  cont_cat_resid_pval <- function(v_cont, v_cat) {
    e <- data_gaussian[, v_cont] - posterior %*% class_means[, v_cont, drop = FALSE]
    cat_var <- data_multinom[, v_cat]
    # Keep only complete cases for both
    comp <- complete.cases(e, cat_var)
    if (sum(comp) < 3) return(list(resid = NA, pval = NA))
    e_comp <- e[comp]
    cat_comp <- cat_var[comp]

    # Residual metric: max absolute Z per category
    levs <- levels(cat_comp)
    if (is.null(levs)) levs <- sort(unique(cat_comp))
    sd_e <- sd(e_comp, na.rm = TRUE)
    if (is.na(sd_e) || sd_e == 0) return(list(resid = NA, pval = NA))

    z_vals <- numeric(length(levs))
    for (l in seq_along(levs)) {
      idx <- which(cat_comp == levs[l])
      n_l <- length(idx)
      if (n_l < 2) {
        z_vals[l] <- NA
        next
      }
      mean_resid <- mean(e_comp[idx])
      se <- sd_e / sqrt(n_l)
      z_vals[l] <- mean_resid / se
    }
    resid_val <- max(abs(z_vals), na.rm = TRUE)

    # P-value from ANOVA F-test
    df <- data.frame(e = e_comp, cat = cat_comp)
    aov_fit <- aov(e ~ cat, data = df)
    pval <- summary(aov_fit)[[1]][["Pr(>F)"]][1]

    # Eta (sqrt of eta-squared)
    ss_between <- summary(aov_fit)[[1]]$`Sum Sq`[1]
    ss_total <- sum(ss_between, summary(aov_fit)[[1]]$`Sum Sq`[2])
    eta <- sqrt(ss_between / ss_total)

    list(resid = resid_val, pval = pval, rlike = eta,
         raw_resid = NA)
  }

  # ----- Main loop -----

  var_names <- item_names
  for (i in 1:(J-1)) {
    for (j in (i+1):J) {
      v1 <- var_names[i]
      v2 <- var_names[j]
      t1 <- types[v1]
      t2 <- types[v2]

      res <- list(resid = NA, pval = NA)

      if (t1 == "gaussian" && t2 == "gaussian") {
        res <- cont_cont_resid_pval(v1, v2)
      } else if (t1 == "multinomial" && t2 == "multinomial") {
        res <- cat_cat_resid_pval(v1, v2)
      } else if (t1 == "gaussian" && t2 == "multinomial") {
        res <- cont_cat_resid_pval(v1, v2)
      } else if (t1 == "multinomial" && t2 == "gaussian") {
        res <- cont_cat_resid_pval(v2, v1)
      }

      resid_mat[v1, v2] <- resid_mat[v2, v1] <- res$resid
      raw_resid_mat[v1, v2] <- raw_resid_mat[v2, v1] <- res$raw_resid
      pval_mat[v1, v2] <- pval_mat[v2, v1] <- res$pval
      r_mat[v1, v2] <- r_mat[v2, v1] <- res$rlike

      names_mat[v1, v2] <- names_mat[v2, v1] <- paste(v1, v2, sep = ".vs.")
      types_mat[v1, v2] <- types_mat[v2, v1] <- paste(t1, t2, sep = "-")

      details[[paste(v1, v2, sep = ".vs.")]] <- list(
        type = paste(t1, t2, sep = "-"),
        residual = res$resid,
        raw_resid = res$raw_resid,
        p_value = res$pval
      )
    }
  }

  ## root mean square residual correlation (RMSR)
  upper_tris <- r_mat[upper.tri(r_mat)]
  rmsr <- sqrt(mean(upper_tris^2))

  ## make it a user friendly table
  res_tab <- data.frame(
    pair = names_mat[upper.tri(names_mat)],
    residual = round(raw_resid_mat[upper.tri(raw_resid_mat)], digits),
    statistic = round(resid_mat[upper.tri(resid_mat)], digits),
    p_value = round(pval_mat[upper.tri(pval_mat)], digits),
    r = round(r_mat[upper.tri(r_mat)], digits),
    var_types = types_mat[upper.tri(types_mat)]
  )

  res_tab$test <- car::recode(res_tab$var_types,
                              " 'gaussian-gaussian' = 'Pearson r';
      'multinomial-multinomial' = 'Pearson chi2';
      'gaussian-multinomial' = 'Z-score';
      'multinomial-gaussian' = 'Z-score' ")
  res_tab <- res_tab[order(abs(res_tab$r), decreasing = TRUE), ]

  out <- list(
    residual_matrix = round(resid_mat, digits),
    raw_residual_matrix = round(raw_resid_mat, digits),
    pvalue_matrix = round(pval_mat, digits),
    r_mat = round(r_mat, digits),
    rmsr = round(rmsr, digits),
    res_tab = res_tab
  )

  class(out) <- "lbvr"

  return(out)
}

#' Print method for lbvr objects
#'
#' @param x An object of class \code{"lbvr"}.
#' @param ... Additional arguments passed to \code{print.data.frame}.
#'
#' @method print lbvr
#' @export
print.lbvr <- function(x, ...) {

  cat("\n——— Residual diagnostics ———\n")
  cat(rep("─", 50), "\n\n", sep = "")

  # Print RMSR
  cat("◆ Root Mean Square Residual (RMSR):", x$rmsr, "\n\n")

  # Print residual table
  cat("◆  Bivariate Residuals (BVR) Table:\n")
  print(x$res_tab, row.names = FALSE, ...)

  cat("\n")

  invisible(x)
}

# Bivariate residuals in LatentGold's style (created with GPT):
latentgold_cont_cont_bvr <- function(data,
                                     v1,
                                     v2,
                                     posterior,
                                     class_means,
                                     class_vars,
                                     weights = NULL,
                                     covariance = c("class_specific"),
                                     tol = 1e-10) {
  covariance <- match.arg(covariance)

  y1 <- data[, v1]
  y2 <- data[, v2]

  if (is.null(weights)) {
    weights <- rep(1, nrow(data))
  }

  get_class_param <- function(x, v, K) {
    if (is.matrix(x) || is.data.frame(x)) {
      out <- as.numeric(x[, v])
    } else if (is.list(x)) {
      out <- as.numeric(x[[v]])
    } else {
      stop("class_means and class_vars should be matrices/data.frames with named columns, or lists.")
    }

    if (length(out) != K) {
      stop("Class-specific parameter length does not match number of classes.")
    }

    out
  }

  K <- ncol(posterior)

  mu1 <- get_class_param(class_means, v1, K)
  mu2 <- get_class_param(class_means, v2, K)

  var1 <- get_class_param(class_vars, v1, K)
  var2 <- get_class_param(class_vars, v2, K)

  if (any(!is.finite(var1)) || any(!is.finite(var2)) ||
      any(var1 <= 0) || any(var2 <= 0)) {
    stop("All class-specific variances must be positive and finite.")
  }

  comp <- complete.cases(y1, y2, posterior, weights)

  if (sum(comp) < 3) {
    return(list(resid = NA_real_, bvr = NA_real_, pval = NA_real_))
  }

  y1 <- y1[comp]
  y2 <- y2[comp]
  tau <- posterior[comp, , drop = FALSE]
  w <- weights[comp]

  rs <- rowSums(tau)
  good <- is.finite(rs) & rs > 0 & is.finite(w) & w > 0

  y1 <- y1[good]
  y2 <- y2[good]
  tau <- tau[good, , drop = FALSE]
  w <- w[good]

  if (length(y1) < 3) {
    return(list(resid = NA_real_, bvr = NA_real_, pval = NA_real_))
  }

  # Normalize posterior probabilities just in case of tiny numerical deviations.
  tau <- tau / rowSums(tau)

  # n x K residual matrices
  e1 <- outer(y1, mu1, "-")
  e2 <- outer(y2, mu2, "-")

  v12 <- var1 * var2

  # First derivative of log class density wrt covariance sigma_12 at sigma_12 = 0
  u <- e1 * e2 / matrix(v12, nrow = length(y1), ncol = K, byrow = TRUE)

  # Second derivative of log class density wrt covariance sigma_12 at sigma_12 = 0
  h <- matrix(1 / v12, nrow = length(y1), ncol = K, byrow = TRUE) -
    e1^2 / matrix(var1^2 * var2, nrow = length(y1), ncol = K, byrow = TRUE) -
    e2^2 / matrix(var1 * var2^2, nrow = length(y1), ncol = K, byrow = TRUE)

  safe_solve <- function(A, b) {
    out <- tryCatch(
      solve(A, b),
      error = function(e) NULL
    )

    if (!is.null(out)) return(out)

    A <- (A + t(A)) / 2
    ee <- eigen(A, symmetric = TRUE)
    keep <- ee$values > tol * max(abs(ee$values))

    if (!any(keep)) {
      return(rep(NA_real_, length(b)))
    }

    ee$vectors[, keep, drop = FALSE] %*%
      ((t(ee$vectors[, keep, drop = FALSE]) %*% b) / ee$values[keep])
  }

  if (covariance == "common") {
    # One covariance parameter shared across classes.
    ui <- rowSums(tau * u)

    # Observed second derivative of the mixture log-likelihood.
    hi <- rowSums(tau * (h + u^2)) - ui^2

    score <- sum(w * ui)
    info <- -sum(w * hi)

    if (!is.finite(info) || info <= tol) {
      return(list(
        resid = NA_real_,
        bvr = NA_real_,
        pval = NA_real_,
        score = score,
        information = info,
        df = 1L
      ))
    }

    stat <- score^2 / info
    bvr <- stat
    pval <- pchisq(stat, df = 1, lower.tail = FALSE)

    return(list(
      resid = bvr,
      bvr = bvr,
      pval = pval,
      score_stat = stat,
      score = score,
      information = info,
      df = 1L
    ))
  }

  if (covariance == "class_specific") {
    # One covariance parameter per class.
    U <- colSums((tau * u) * w)

    A <- tau * u
    H <- diag(colSums((tau * (h + u^2)) * w), nrow = K) -
      t(A) %*% (A * w)

    info <- -H
    sol <- safe_solve(info, U)

    if (anyNA(sol)) {
      return(list(
        resid = NA_real_,
        bvr = NA_real_,
        pval = NA_real_,
        score = U,
        information = info,
        df = K
      ))
    }

    score_stat <- as.numeric(t(U) %*% sol)

    # LatentGOLD reports BVR as improvement per added parameter.
    bvr <- score_stat / K

    # Approximate only. LatentGOLD itself recommends/bootstrap-reports p-values
    # because the asymptotic distribution of BVRs is not generally known.
    pval <- pchisq(score_stat, df = K, lower.tail = FALSE)

    return(list(
      resid = bvr,
      bvr = bvr,
      pval = pval,
      score_stat = score_stat,
      score = U,
      information = info,
      df = K
    ))
  }
}
