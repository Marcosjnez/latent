# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 12/07/2026
#'
#' @title Local bivariate residuals for latent class analysis
#'
#' @description
#' Computes pairwise residual diagnostics for the indicators of a fitted latent
#' class model. Different diagnostics are used for pairs of Gaussian indicators,
#' pairs of multinomial indicators, and mixed Gaussian-multinomial pairs.
#'
#' @param x A fitted object of class \code{"llca"}.
#' @param ... Additional arguments reserved for other methods.
#' @param digits A non-negative integer indicating the number of decimal places
#'   used to round the returned matrices and summary table. The default is
#'   \code{4L}.
#'
#' @details
#' The function evaluates every pair of measurement indicators. Covariates and
#' distal outcomes are not included. Pairwise complete observations are used for
#' each diagnostic.
#'
#' For two Gaussian indicators, the standardized diagnostic is the Pearson
#' correlation between their posterior-weighted residuals. The raw residual is
#' obtained with the LatentGold-style score calculation implemented in
#' \code{latentgold_cont_cont_bvr()}.
#'
#' For two multinomial indicators, the function compares the observed
#' cross-classification table with the model-implied table. It reports the
#' maximum absolute standardized residual, its chi-squared p-value, the
#' difference between the observed and expected Cramer's V values, and the raw
#' residual score used by the package.
#'
#' For a Gaussian and a multinomial indicator, the function reports the maximum
#' absolute category-specific Z-score of the Gaussian residuals, the p-value
#' from the ANOVA calculation, and eta as a correlation-like effect size. A raw
#' residual is not defined for these mixed pairs and is returned as
#' \code{NA_real_}.
#'
#' Indicator pairs already included as residual dependencies in the fitted model
#' are assigned a residual of zero, a p-value of one, and an effect size of zero.
#'
#' @return An object of class \code{"lbvr"}, consisting of a list with the
#'   following components:
#'   \item{residual_matrix}{A symmetric matrix containing the standardized
#'     pairwise diagnostic.}
#'   \item{raw_residual_matrix}{A symmetric matrix containing the raw residual
#'     score. Mixed Gaussian-multinomial pairs contain \code{NA_real_}.}
#'   \item{pvalue_matrix}{A symmetric matrix containing the corresponding
#'     p-values.}
#'   \item{r_mat}{A symmetric matrix containing the correlation-like effect
#'     sizes.}
#'   \item{rmsr}{The root mean square of the upper-triangular values of
#'     \code{r_mat}.}
#'   \item{res_tab}{A data frame with one row per indicator pair and columns
#'     \code{pair}, \code{residual}, \code{statistic}, \code{p_value},
#'     \code{r}, \code{var_types}, and \code{test}. Rows are ordered by the
#'     absolute value of \code{r}.}
#'
#' @examples
#' \dontrun{
#' fit <- lca(data = gss82, nclasses = 3L,
#'            multinomial = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"))
#'
#' residuals <- lbvr(fit)
#' residuals$res_tab
#' }
#'
#' @importFrom stats aov chisq.test complete.cases cor.test pchisq sd
#' @method lbvr llca
#' @export
lbvr.llca <- function(x, digits = 4L, ...) {

  model <- x

  if(!inherits(model, "llca")) {
    stop("model must be a fitted object of class 'llca'")
  }
  if(length(model@transformed_pars) == 0L) {
    stop("model must be fitted before computing bivariate residuals")
  }
  if(length(digits) != 1L || is.na(digits) || digits < 0L ||
     digits != as.integer(digits)) {
    stop("digits must be a single non-negative integer")
  }
  digits <- as.integer(digits)

  indicators_names <- model@dataList$indicators_names
  variable_type <- model@dataList$variable_type
  if(is.null(names(variable_type))) {
    stop("model@dataList$variable_type must be a named vector")
  }

  types <- variable_type[indicators_names]
  if(anyNA(types)) {
    stop("The variable type is missing for at least one indicator")
  }

  item_names <- indicators_names
  J <- length(item_names)
  if(J < 2L) {
    stop("At least two indicators are required to compute bivariate residuals")
  }

  posterior <- latInspect(model, what = "posterior")
  K <- ncol(posterior)
  short2full <- model@dataList$short2full

  continuous_names <- c(model@dataList$gaussian$gaussian_names,
                        model@dataList$mvgaussian$mvgaussian_names)
  cont_vars <- intersect(indicators_names, continuous_names)
  cat_vars <- intersect(indicators_names,
                        model@dataList$multinomial$multinomial_names)
  unsupported <- setdiff(indicators_names, c(cont_vars, cat_vars))
  if(length(unsupported) > 0L) {
    stop("Unsupported indicator type for: ", paste(unsupported, collapse = ", "))
  }

  n_cont <- length(cont_vars)
  n_cat <- length(cat_vars)

  if(n_cont > 0L) {
    continuous_patterns <- cbind(model@dataList$gaussian$patterns_gaussian,
                                 model@dataList$mvgaussian$patterns_mvgaussian)
    data_cont <- continuous_patterns[short2full, cont_vars, drop = FALSE]
  } else {
    data_cont <- matrix(numeric(0L), nrow = length(short2full), ncol = 0L,
                        dimnames = list(NULL, character(0L)))
  }

  if(n_cat > 0L) {
    data_multinom <- model@dataList$multinomial$patterns_multinomial[
      short2full, cat_vars, drop = FALSE] + 1L
  } else {
    data_multinom <- matrix(integer(0L), nrow = length(short2full), ncol = 0L,
                            dimnames = list(NULL, character(0L)))
  }

  profile <- latInspect(model, what = "item")

  class_means <- matrix(numeric(0L), nrow = K, ncol = 0L,
                        dimnames = list(colnames(posterior), character(0L)))
  if(n_cont > 0L) {
    class_means <- matrix(NA_real_, nrow = K, ncol = n_cont,
                          dimnames = list(colnames(posterior), cont_vars))
    for(v in cont_vars) {
      indicator_profile <- profile[[v]]
      if(is.null(indicator_profile) || is.null(rownames(indicator_profile)) ||
         !"mean" %in% rownames(indicator_profile)) {
        stop("Class-specific means are missing for continuous indicator '", v, "'")
      }
      if(ncol(indicator_profile) != K) {
        stop("The class-specific profile for '", v,
             "' does not match the number of latent classes")
      }
      class_means[, v] <- as.numeric(indicator_profile["mean", , drop = TRUE])
    }
  }

  class_probs <- profile[cat_vars]
  if(n_cat > 0L) {
    if(length(class_probs) != n_cat || !all(cat_vars %in% names(class_probs))) {
      stop("Class-specific probabilities are missing for categorical indicators")
    }
    for(v in cat_vars) {
      if(is.null(class_probs[[v]]) || ncol(class_probs[[v]]) != K) {
        stop("The class-specific probabilities for '", v,
             "' do not match the number of latent classes")
      }
    }
  }

  pi <- as.numeric(latInspect(model, what = "class"))
  if(length(pi) != K) {
    stop("The estimated class proportions do not match the number of latent classes")
  }

  resid_mat <- matrix(NA_real_, J, J, dimnames = list(item_names, item_names))
  raw_resid_mat <- matrix(NA_real_, J, J,
                          dimnames = list(item_names, item_names))
  pval_mat <- matrix(NA_real_, J, J, dimnames = list(item_names, item_names))
  r_mat <- matrix(NA_real_, J, J, dimnames = list(item_names, item_names))
  names_mat <- matrix(NA_character_, J, J,
                      dimnames = list(item_names, item_names))
  types_mat <- matrix(NA_character_, J, J,
                      dimnames = list(item_names, item_names))

  pair_is_modeled <- function(pairs, v1, v2) {

    result <- FALSE
    if(!is.null(pairs) && length(pairs) > 0L) {
      pairs <- as.matrix(pairs)
      if(ncol(pairs) < 2L && nrow(pairs) >= 2L) {
        pairs <- t(pairs)
      }
      if(ncol(pairs) >= 2L) {
        for(i in seq_len(nrow(pairs))) {
          current_pair <- as.character(pairs[i, , drop = TRUE])
          if(v1 %in% current_pair && v2 %in% current_pair) {
            result <- TRUE
            break
          }
        }
      }
    }

    #### Result ####

    return(result)

  }

  cont_cont_resid_pval <- function(v1, v2) {

    result <- list(resid = NA_real_, pval = NA_real_, rlike = NA_real_,
                   raw_resid = NA_real_)
    modeled <- pair_is_modeled(model@dataList$dependencies$gaussian_pairs,
                               v1, v2)

    if(modeled) {
      result <- list(resid = 0, pval = 1, rlike = 0, raw_resid = 0)
    } else {
      e1 <- data_cont[, v1] - posterior %*% class_means[, v1, drop = FALSE]
      e2 <- data_cont[, v2] - posterior %*% class_means[, v2, drop = FALSE]
      comp <- complete.cases(e1, e2)

      if(sum(comp) >= 3L) {
        e1_complete <- as.numeric(e1[comp])
        e2_complete <- as.numeric(e2[comp])
        sd1 <- sd(e1_complete)
        sd2 <- sd(e2_complete)

        if(is.finite(sd1) && sd1 > 0 && is.finite(sd2) && sd2 > 0) {
          ct <- cor.test(e1_complete, e2_complete, use = "complete.obs")
          raw_resid <- latentgold_cont_cont_bvr(model, v1, v2)$resid
          result <- list(resid = as.numeric(ct$estimate), pval = ct$p.value,
                         rlike = as.numeric(ct$estimate), raw_resid = raw_resid)
        }
      }
    }

    #### Result ####

    return(result)

  }

  cat_cat_resid_pval <- function(v1, v2) {

    cramer_v <- function(tab) {

      result <- NA_real_
      if(all(dim(tab) >= 2L) && sum(tab) > 0) {
        n <- sum(tab)
        chi2 <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
        df <- min(nrow(tab) - 1L, ncol(tab) - 1L)
        if(df > 0L) {
          result <- sqrt(as.numeric(chi2) / (n * df))
        }
      }

      #### Result ####

      return(result)

    }

    result <- list(resid = NA_real_, pval = NA_real_, rlike = NA_real_,
                   raw_resid = NA_real_)
    modeled <- pair_is_modeled(model@dataList$dependencies$multinomial_pairs,
                               v1, v2)

    if(modeled) {
      result <- list(resid = 0, pval = 1, rlike = 0, raw_resid = 0)
    } else {
      comp <- complete.cases(data_multinom[, v1], data_multinom[, v2])

      if(sum(comp) > 0L) {
        obs_tab <- table(data_multinom[comp, v1], data_multinom[comp, v2])
        N_pair <- sum(comp)

        prob1 <- class_probs[[v1]]
        prob2 <- class_probs[[v2]]
        L1 <- nrow(prob1)
        L2 <- nrow(prob2)
        exp_tab <- matrix(0, L1, L2)

        for(r in seq_len(L1)) {
          for(c in seq_len(L2)) {
            exp_tab[r, c] <- sum(pi * prob1[r, ] * prob2[c, ]) * N_pair
          }
        }

        obs_cat1 <- as.numeric(rownames(obs_tab))
        obs_cat2 <- as.numeric(colnames(obs_tab))
        exp_tab_aligned <- exp_tab[obs_cat1, obs_cat2, drop = FALSE]

        std_resid <- (obs_tab - exp_tab_aligned) / sqrt(exp_tab_aligned)
        std_resid[exp_tab_aligned == 0] <- NA_real_
        chi2 <- sum(std_resid^2, na.rm = TRUE)
        df <- (nrow(obs_tab) - 1L) * (ncol(obs_tab) - 1L)
        if(df > 0L) {
          pval <- pchisq(chi2, df, lower.tail = FALSE)
        } else {
          pval <- NA_real_
        }

        finite_residuals <- abs(std_resid[is.finite(std_resid)])
        if(length(finite_residuals) > 0L) {
          resid_val <- max(finite_residuals)
        } else {
          resid_val <- NA_real_
        }

        V_obs <- cramer_v(obs_tab)
        V_exp <- cramer_v(exp_tab_aligned)
        resid_cv <- V_obs - V_exp

        P <- (L1 - 1L) * (L2 - 1L)
        if(P > 0L) {
          raw_resid <- sum((obs_tab - exp_tab_aligned)^2 / exp_tab_aligned) / P
        } else {
          raw_resid <- NA_real_
        }

        result <- list(resid = resid_val, pval = pval, rlike = resid_cv,
                       raw_resid = raw_resid)
      }
    }

    #### Result ####

    return(result)

  }

  cont_cat_resid_pval <- function(v_cont, v_cat) {

    result <- list(resid = NA_real_, pval = NA_real_, rlike = NA_real_,
                   raw_resid = NA_real_)
    e <- data_cont[, v_cont] - posterior %*% class_means[, v_cont, drop = FALSE]
    cat_var <- data_multinom[, v_cat]
    comp <- complete.cases(e, cat_var)

    if(sum(comp) >= 3L) {
      e_comp <- as.numeric(e[comp])
      cat_comp <- cat_var[comp]
      levs <- levels(cat_comp)
      if(is.null(levs)) {
        levs <- sort(unique(cat_comp))
      }

      sd_e <- sd(e_comp, na.rm = TRUE)
      if(is.finite(sd_e) && sd_e > 0 && length(unique(cat_comp)) >= 2L) {
        z_vals <- numeric(length(levs))
        for(l in seq_along(levs)) {
          idx <- which(cat_comp == levs[l])
          n_l <- length(idx)
          if(n_l < 2L) {
            z_vals[l] <- NA_real_
          } else {
            mean_resid <- mean(e_comp[idx])
            se <- sd_e / sqrt(n_l)
            z_vals[l] <- mean_resid / se
          }
        }

        finite_z <- abs(z_vals[is.finite(z_vals)])
        if(length(finite_z) > 0L) {
          resid_val <- max(finite_z)
        } else {
          resid_val <- NA_real_
        }

        anova_data <- data.frame(e = e_comp, cat = cat_comp)
        aov_fit <- aov(e ~ cat, data = anova_data)
        aov_summary <- summary(aov_fit)[[1]]
        pval <- aov_summary[["Pr(>F)"]][1]
        ss_between <- aov_summary$`Sum Sq`[1]
        ss_total <- sum(ss_between, aov_summary$`Sum Sq`[2])
        if(is.finite(ss_total) && ss_total > 0) {
          eta <- sqrt(ss_between / ss_total)
        } else {
          eta <- NA_real_
        }

        result <- list(resid = resid_val, pval = pval, rlike = eta,
                       raw_resid = NA_real_)
      }
    }

    #### Result ####

    return(result)

  }

  continuous_types <- c("gaussian", "mvgaussian")
  for(i in seq_len(J - 1L)) {
    for(j in seq.int(i + 1L, J)) {
      v1 <- item_names[i]
      v2 <- item_names[j]
      t1 <- types[[v1]]
      t2 <- types[[v2]]
      result_pair <- list(resid = NA_real_, pval = NA_real_, rlike = NA_real_,
                          raw_resid = NA_real_)

      if(t1 %in% continuous_types && t2 %in% continuous_types) {
        result_pair <- cont_cont_resid_pval(v1, v2)
      } else if(t1 == "multinomial" && t2 == "multinomial") {
        result_pair <- cat_cat_resid_pval(v1, v2)
      } else if(t1 %in% continuous_types && t2 == "multinomial") {
        result_pair <- cont_cat_resid_pval(v1, v2)
      } else if(t1 == "multinomial" && t2 %in% continuous_types) {
        result_pair <- cont_cat_resid_pval(v2, v1)
      }

      resid_mat[v1, v2] <- resid_mat[v2, v1] <- result_pair$resid
      raw_resid_mat[v1, v2] <- raw_resid_mat[v2, v1] <- result_pair$raw_resid
      pval_mat[v1, v2] <- pval_mat[v2, v1] <- result_pair$pval
      r_mat[v1, v2] <- r_mat[v2, v1] <- result_pair$rlike

      if(t1 %in% continuous_types) {
        display_t1 <- "gaussian"
      } else {
        display_t1 <- t1
      }
      if(t2 %in% continuous_types) {
        display_t2 <- "gaussian"
      } else {
        display_t2 <- t2
      }

      pair_name <- paste(v1, v2, sep = ".vs.")
      pair_type <- paste(display_t1, display_t2, sep = "-")
      names_mat[v1, v2] <- names_mat[v2, v1] <- pair_name
      types_mat[v1, v2] <- types_mat[v2, v1] <- pair_type
    }
  }

  upper_tris <- r_mat[upper.tri(r_mat)]
  rmsr <- sqrt(mean(upper_tris^2))

  res_tab <- data.frame(pair = names_mat[upper.tri(names_mat)],
                        residual = round(raw_resid_mat[upper.tri(raw_resid_mat)],
                                         digits),
                        statistic = round(resid_mat[upper.tri(resid_mat)], digits),
                        p_value = round(pval_mat[upper.tri(pval_mat)], digits),
                        r = round(r_mat[upper.tri(r_mat)], digits),
                        var_types = types_mat[upper.tri(types_mat)])

  test_labels <- c("gaussian-gaussian" = "Pearson r",
                   "multinomial-multinomial" = "Pearson chi2",
                   "gaussian-multinomial" = "Z-score",
                   "multinomial-gaussian" = "Z-score")
  res_tab$test <- unname(test_labels[res_tab$var_types])
  res_tab <- res_tab[order(abs(res_tab$r), decreasing = TRUE), , drop = FALSE]
  rownames(res_tab) <- NULL

  result <- list(residual_matrix = round(resid_mat, digits),
                 raw_residual_matrix = round(raw_resid_mat, digits),
                 pvalue_matrix = round(pval_mat, digits),
                 r_mat = round(r_mat, digits), rmsr = round(rmsr, digits),
                 res_tab = res_tab)
  class(result) <- "lbvr"

  #### Result ####

  return(result)

}


#' @title Local bivariate residuals for a collection of latent class models
#'
#' @description
#' Applies \code{lbvr()} to every fitted \code{"llca"} object contained in an
#' object of class \code{"llcalist"}. Model names are preserved. Unnamed models
#' are labelled according to their number of latent classes.
#'
#' @param x An object of class \code{"llcalist"} containing fitted
#'   \code{"llca"} models.
#' @param digits A non-negative integer indicating the number of decimal places
#'   used to round the returned diagnostics. The default is \code{4L}.
#' @param ... Additional arguments passed to the \code{lbvr.llca()} method.
#'
#' @return A named list containing one object of class \code{"lbvr"} for each
#'   fitted model. The returned list has class \code{"lbvr.llcalist"}.
#'
#' @examples
#' \dontrun{
#' fits <- lca(data = gss82, nclasses = 2:4,
#'             multinomial = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"))
#' residuals <- lbvr(fits)
#' }
#'
#' @method lbvr llcalist
#' @export
lbvr.llcalist <- function(x, digits = 4L, ...) {

  if(!inherits(x, "llcalist")) {
    stop("x must be an object of class 'llcalist'")
  }
  if(length(x) == 0L) {
    stop("x must contain at least one fitted 'llca' model")
  }
  valid_models <- vapply(x, FUN = inherits, FUN.VALUE = logical(1L), what = "llca")
  if(!all(valid_models)) {
    stop("Every element of x must be a fitted object of class 'llca'")
  }

  result <- lapply(x, FUN = lbvr, digits = digits, ...)
  model_names <- names(x)
  if(is.null(model_names)) {
    model_names <- rep("", length(x))
  }
  unnamed <- !nzchar(model_names)
  if(any(unnamed)) {
    nclasses <- vapply(x[unnamed], FUN = function(model) {
      ncol(model@transformed_pars$class)
    }, FUN.VALUE = integer(1L))
    model_names[unnamed] <- paste0("nclasses=", nclasses)
  }
  names(result) <- model_names
  class(result) <- c("lbvr.llcalist", "list")

  #### Result ####

  return(result)

}

#' @title Print local bivariate residual diagnostics
#'
#' @description
#' Prints the root mean square residual and the pairwise residual summary table
#' stored in an object returned by \code{lbvr()}.
#'
#' @param x An object of class \code{"lbvr"}.
#' @param ... Additional arguments passed to \code{print.data.frame()} when
#'   printing the residual summary table.
#'
#' @return The input object \code{x}, returned invisibly.
#'
#' @method print lbvr
#' @export
print.lbvr <- function(x, ...) {

  cat("\n——— Residual diagnostics ———\n")
  cat(rep("─", 50L), "\n\n", sep = "")
  cat("◆ Root Mean Square Residual (RMSR):", x$rmsr, "\n\n")
  cat("◆ Bivariate Residuals (BVR) Table:\n")
  print(x$res_tab, row.names = FALSE, ...)
  cat("\n")

  #### Result ####

  return(invisible(x))

}

#' @title LatentGold-style residual for two continuous indicators
#'
#' @description
#' Computes the score-based raw bivariate residual used by \code{lbvr()} for a
#' pair of continuous indicators. The fitted model is recreated with the selected
#' residual dependency while all previously estimated parameters remain fixed.
#'
#' @param fit A fitted object of class \code{"llca"}.
#' @param v1 Character string naming the first continuous indicator.
#' @param v2 Character string naming the second continuous indicator.
#'
#' @return A list with components \code{resid}, \code{pval},
#'   \code{score_stat}, and \code{df}.
#'
#' @keywords internal
latentgold_cont_cont_bvr <- function(fit, v1, v2) {

  args <- fit@dataList$args
  gaussian_pairs <- fit@dataList$dependencies$gaussian_pairs
  if(!is.null(gaussian_pairs) && length(gaussian_pairs) > 0L) {
    gaussian_pairs <- as.matrix(gaussian_pairs)
    if(ncol(gaussian_pairs) < 2L && nrow(gaussian_pairs) >= 2L) {
      gaussian_pairs <- t(gaussian_pairs)
    }
    model_dependencies <- apply(gaussian_pairs, MARGIN = 1L, FUN = paste,
                                collapse = "~~")
    args$model <- c(paste(c(v1, v2), collapse = "~~"), model_dependencies, fit)
  } else {
    args$model <- c(paste(c(v1, v2), collapse = "~~"), fit)
  }

  args$do.fit <- FALSE
  args$adjustment <- "none"
  args$control <- list(rstarts = 1L, cores = 1L, free_beta = FALSE)
  fit2 <- do.call(lca, args)

  fit2@modelInfo$control_optimizer$parameters[[1]][] <- 0
  control_manifold <- fit2@modelInfo$control_manifold
  control_transform <- fit2@modelInfo$control_transform
  control_estimator <- fit2@modelInfo$control_estimator
  control_optimizer <- fit2@modelInfo$control_optimizer

  U <- get_grad(control_manifold, control_transform, control_estimator,
                control_optimizer)$g
  I <- get_hess(control_manifold, control_transform, control_estimator,
                control_optimizer)$h
  K <- length(U)

  score_stat <- as.numeric(t(U) %*% solve(I) %*% U)
  resid <- score_stat / K
  pval <- 1 - pchisq(resid, df = K)

  result <- list(resid = resid, pval = pval, score_stat = score_stat, df = K)

  #### Result ####

  return(result)

}
