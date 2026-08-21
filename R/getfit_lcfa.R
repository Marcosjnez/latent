# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 21/08/2026
#'
#' Fit Indices for CFA Models
#'
#' @param model A fitted object inheriting from class \code{"lcfa"}.
#' @param digits Number of decimal places used in the returned values.
#' @param fit_matrix Optional precomputed internal fit matrix. This argument is
#'   used by package post-processing methods to avoid refitting the direct-FIML
#'   saturated model more than once.
#'
#' @return A named numeric vector containing model dimensions and available fit
#'   statistics.
#'
#' @method getfit lcfa
#' @export
getfit.lcfa <- function(model, digits = 3L, fit_matrix = NULL) {

  #### Check inputs ####

  if(!inherits(model, "lcfa")) {
    stop("model must inherit from class 'lcfa'.")
  }

  if(!is.numeric(digits) ||
     length(digits) != 1L ||
     !is.finite(digits) ||
     digits < 0L ||
     digits != as.integer(digits)) {
    stop("digits must be a non-negative integer.")
  }

  nobs <- sum(unlist(model@dataList$nobs))
  nitems <- unlist(model@dataList$nitems)
  nfactors <- unlist(model@dataList$nfactors)
  nparam <- model@modelInfo$nparam
  dof <- model@modelInfo$dof
  dof_base <- sum(nitems*(nitems-1L)/2)

  if(is.null(fit_matrix)) {
    fit_matrix <- lcfa_fit_matrix(model)
  }

  likelihood_model <- model@dataList$estimator %in%
    c("ml", "fml", "means_fml")

  loglik <- NA_real_
  AIC <- NA_real_
  BIC <- NA_real_
  X2 <- NA_real_
  X2_base <- NA_real_

  if(likelihood_model) {

    loglik <- fit_matrix["loglik", "overall"]
    loglik_sat <- fit_matrix["loglik_sat", "overall"]
    loglik_base <- fit_matrix["loglik_base", "overall"]

    if(all(is.finite(c(loglik, loglik_sat)))) {
      X2 <- 2*(loglik_sat-loglik)
    }

    if(all(is.finite(c(loglik_base, loglik_sat)))) {
      X2_base <- 2*(loglik_sat-loglik_base)
    }

    if(is.finite(loglik)) {
      AIC <- -2*loglik+2*nparam
      BIC <- -2*loglik+log(nobs)*nparam
    }

  } else {

    loss <- fit_matrix["loss", "overall"]
    loss_base <- fit_matrix["loss_base", "overall"]

    if(is.finite(loss)) {
      X2 <- 2*nobs*loss
    }

    if(is.finite(loss_base)) {
      X2_base <- 2*nobs*loss_base
    }

  }

  pvalue <- if(is.finite(X2) && dof > 0L) {
    stats::pchisq(X2, df = dof, lower.tail = FALSE)
  } else {
    NA_real_
  }

  CFI <- NA_real_

  if(all(is.finite(c(X2, X2_base)))) {

    numerator <- max(X2-dof, 0)
    denominator <- max(X2-dof,
                       X2_base-dof_base,
                       0)

    CFI <- if(denominator > 0) {
      1-numerator/denominator
    } else {
      1
    }

    CFI <- min(max(CFI, 0), 1)

  }

  RMSEA <- if(is.finite(X2) &&
              dof > 0L &&
              nobs > 1L) {
    sqrt(max((X2-dof)/(dof*(nobs-1L)), 0))
  } else {
    NA_real_
  }

  #### Standardized residual mean square ####

  residuals <- tryCatch(
    lcfa_residuals(model),
    error = function(e) NULL
  )
  standardized_residuals <- numeric(0L)

  if(!is.null(residuals)) {

    for(i in seq_along(residuals)) {

      covariance <- residuals[[i]]$covariance
      sample_covariance <- model@dataList$sample.cov[[i]]
      variances <- diag(sample_covariance)

      if(all(is.finite(variances)) &&
         all(variances > 0)) {

        scale <- outer(sqrt(variances),
                       sqrt(variances))
        residual_correlation <- covariance/scale

        standardized_residuals <- c(
          standardized_residuals,
          residual_correlation[
            lower.tri(residual_correlation, diag = TRUE)
          ]
        )

      }

    }

  }

  SRMR <- if(length(standardized_residuals) > 0L &&
             all(is.finite(standardized_residuals))) {
    sqrt(mean(standardized_residuals^2))
  } else {
    NA_real_
  }

  result <- c(
    ngroups = model@dataList$ngroups,
    nfactors = sum(nfactors),
    npar = nparam,
    nobs = nobs,
    loglik = loglik,
    AIC = AIC,
    BIC = BIC,
    chisq = X2,
    dof = dof,
    pvalue = pvalue,
    cfi = CFI,
    rmsea = RMSEA,
    srmr = SRMR
  )

  result <- round(result, digits = digits)
  class(result) <- c("getfit.lcfa", "numeric")

  #### Result ####

  return(result)

}
