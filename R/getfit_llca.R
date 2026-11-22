# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 10/09/2026
#'
#' Fit indices for latent class models
#'
#' @description
#' Computes likelihood-based fit indices, classification measures, and entropy
#' for a fitted latent class model or for a list of fitted latent class models.
#'
#' @param model An object of class \code{"llca"} fitted with \code{lca()}, or an
#'   object of class \code{"llcalist"} containing several fitted
#'   \code{"llca"} models.
#' @param digits Integer giving the number of decimal places used to round the
#'   output. Use \code{NULL} to return the unrounded values.
#'
#' @details
#' For a single \code{"llca"} model, \code{getfit()} returns the number of
#' classes, model parameters, and observations; the log-likelihood; information
#' criteria; the integrated classification likelihood; and entropy-based
#' pseudo-\eqn{R^2}.
#'
#' The likelihood-ratio statistic \eqn{L^2}, its degrees of freedom, and its
#' p-value use complete multinomial responses, comparing frequencies within
#' each observed covariate pattern. The response-pattern test is unavailable
#' with missing responses or Gaussian variables. Posterior entropy and the
#' classification contribution to ICL use the fitted response-pattern weights.
#'
#' For hidden Markov models, information criteria use the trajectory
#' log-likelihoods and the total number of observations across trajectories.
#' The response-pattern test, its degrees of freedom and p-value, the
#' classification likelihood, and entropy-based pseudo-\eqn{R^2} are
#' returned as \code{NA}.
#'
#' When penalization is active, the output also contains the penalized
#' log-likelihood and the corresponding penalized information criteria. The
#' attribute \code{"penalized"} indicates whether these additional indices are
#' present.
#'
#' For an \code{"llcalist"} object, the indices are collected in a named
#' list with one element per fitted model.
#'
#' @return
#' For an object of class \code{"llca"}, a named numeric vector of class
#' \code{"getfit.llca"} containing:
#' \describe{
#'   \item{nclasses}{Number of latent classes.}
#'   \item{npar}{Number of freely estimated model parameters.}
#'   \item{nobs}{Number of observations used in estimation.}
#'   \item{loglik}{Model log-likelihood.}
#'   \item{penalized_loglik}{Penalized log-likelihood, when penalization is active.}
#'   \item{L2}{Likelihood-ratio statistic for a fully multinomial model.}
#'   \item{dof}{Degrees of freedom associated with \code{L2}.}
#'   \item{pvalue}{P-value associated with \code{L2}.}
#'   \item{AIC}{Akaike information criterion.}
#'   \item{BIC}{Bayesian information criterion.}
#'   \item{AIC3}{AIC using a penalty of three per parameter.}
#'   \item{CAIC}{Consistent Akaike information criterion.}
#'   \item{KIC}{Kullback information criterion.}
#'   \item{SABIC}{Sample-size-adjusted BIC.}
#'   \item{ICL}{Integrated classification likelihood.}
#'   \item{AICp, BICp, AIC3p, CAICp, KICp, SABICp, ICLp}{Penalized versions of
#'     the corresponding indices, returned only when penalization is active.}
#'   \item{R2_entropy}{Entropy-based pseudo-\eqn{R^2}.}
#' }
#'
#' For an object of class \code{"llcalist"}, a list of class
#' \code{"getfit.llcalist"} with one vector of fit indices per model.
#'
#' @references
#' Akaike, H. (1974). A New Look at the Statistical Model Identification.
#' \emph{In: Parzen, E., Tanabe, K., Kitagawa, G. (eds) Selected Papers of Hirotugu Akaike.
#' Springer Series in Statistics. Springer, New York, NY.}
#' https://doi.org/10.1007/978-1-4612-1694-0_16
#'
#' @examples
#' \dontrun{
#' fit <- lca(data = gss82, nclasses = 3L,
#'            multinomial = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"))
#'
#' getfit(fit)
#' getfit(fit, digits = NULL)
#'
#' fits <- lca(data = gss82, nclasses = 1:4,
#'             multinomial = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"))
#'
#' getfit(fits)
#' }
#'
#' @aliases getfit
#' @param x Fitted model or collection of fitted models passed to the method.
#' @param ... Additional arguments to the generic.
#' @method getfit llca
#' @export
getfit.llca <- function(x, digits = 4L, ..., model = x) {
  if(!missing(model) && !missing(x)) return(.match_legacy_model("model"))

  #### Check inputs ####

  if(!inherits(model, "llca")) {
    stop("model must inherit from class 'llca'.")
  }

  if(length(model@Optim) == 0L) {
    stop("getfit() requires a fitted llca object.")
  }

  if(!is.null(digits) &&
     (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
      digits < 0L || digits != as.integer(digits))) {
    stop("digits must be NULL or a non-negative integer.")
  }

  #### Extract model information ####

  penalized <- model@modelInfo$control_optimizer$reg
  temporal <- !is.null(model@dataList$lengths)
  if(!temporal) {
    posterior <- .llca_posterior(model)
    weights <- model@dataList$pattern_weights
  }
  fit_loglik <- latInspect(model, what = "loglik", sort = FALSE)

  nclasses <- ncol(model@transformed_pars$class)
  k <- model@modelInfo$nparam
  nobs <- model@dataList$nobs
  loglik <- fit_loglik["loglik", "overall"]
  penalized_loglik <- fit_loglik["penalized_loglik", "overall"]

  #### Likelihood-ratio test ####

  d <- model@dataList
  if(!temporal && all(d$variable_type == "multinomial") &&
     !anyNA(d$multinomial$patterns_multinomial)) {

    # Expected counts use the total frequency at the same covariate values.
    strata <- .llca_strata(d)
    total <- strata$total
    dof <- model@modelInfo$dof
    if(length(d$covariates_names)) {
      cells <- prod(d$multinomial$multinomial_factor_lengths)
      dof <- min(nobs, strata$npatterns*(cells-1)) - k
    }
    ll <- as.numeric(model@Optim$outputs$estimators$vectors[[1L]][[1L]])
    L2 <- max(0, 2*sum(weights*(log(weights/total)-ll)))
    pv <- if(dof > 0) pchisq(L2, dof, lower.tail = FALSE) else NA_real_

  } else {

    L2 <- NA_real_
    pv <- NA_real_
    dof <- if(temporal || all(d$variable_type == "multinomial")) NA_real_ else model@modelInfo$dof

  }

  #### Information criteria ####

  AIC <- (-2*loglik) + (k*2)
  AIC3 <- (-2*loglik) + (k*3)
  BIC <- k*log(nobs) - 2*loglik
  CAIC <- -2*loglik + (k*(log(nobs) + 1))
  KIC <- -2*loglik + (3*(k + 1))
  SABIC <- -2*loglik + (k*log((nobs + 2)/24))
  if(!temporal) {
    modal <- max.col(posterior, ties.method = "first")
    modal_loglik <- sum(weights * log(posterior[cbind(seq_along(modal), modal)]))
  }
  ICL <- if(temporal) NA_real_ else -BIC + 2*modal_loglik

  #### Entropy ####

  if(temporal) {
    entropyR2 <- NA_real_
  } else if(nclasses < 2L) {

    entropyR2 <- 1.00

  } else {

    classes <- colSums(model@transformed_pars$class * weights) / sum(weights)
    positive <- classes > 0
    prior_entropy <- -sum(classes[positive]*log(classes[positive]))
    log_post <- matrix(0, nrow(posterior), ncol(posterior))
    positive <- posterior > 0
    log_post[positive] <- log(posterior[positive])
    post_entropy <- -sum(weights * rowSums(posterior * log_post))/sum(weights)
    entropyR2 <- if(prior_entropy > 0) 1 - post_entropy/prior_entropy else NA_real_

  }

  #### Collect fit indices ####

  if(penalized) {

    AICp <- (-2*penalized_loglik) + (k*2)
    AIC3p <- (-2*penalized_loglik) + (k*3)
    BICp <- k*log(nobs) - 2*penalized_loglik
    CAICp <- -2*penalized_loglik + (k*(log(nobs) + 1))
    KICp <- -2*penalized_loglik + (3*(k + 1))
    SABICp <- -2*penalized_loglik + (k*log((nobs + 2)/24))
    ICLp <- if(temporal) NA_real_ else -BICp + 2*modal_loglik

    result <- c(nclasses = nclasses, npar = k, nobs = nobs,
                loglik = loglik, penalized_loglik = penalized_loglik,
                L2 = L2, dof = dof, pvalue = pv,
                AIC = AIC, BIC = BIC, AIC3 = AIC3,
                CAIC = CAIC, KIC = KIC, SABIC = SABIC, ICL = ICL,
                AICp = AICp, BICp = BICp, AIC3p = AIC3p,
                CAICp = CAICp, KICp = KICp, SABICp = SABICp, ICLp = ICLp,
                R2_entropy = entropyR2)

  } else {

    result <- c(nclasses = nclasses, npar = k, nobs = nobs,
                loglik = loglik, L2 = L2, dof = dof, pvalue = pv,
                AIC = AIC, BIC = BIC, AIC3 = AIC3,
                CAIC = CAIC, KIC = KIC, SABIC = SABIC, ICL = ICL,
                R2_entropy = entropyR2)

  }

  if(!is.null(digits)) {
    result <- round(result, digits = digits)
  }

  class(result) <- "getfit.llca"
  attr(result, "penalized") <- penalized

  #### Result ####

  return(result)

}

# Frequencies of the conditioning covariate patterns for multinomial responses.
.llca_strata <- function(d) {
  if(!length(d$covariates_names)) return(list(npatterns = 1L, total = sum(d$pattern_weights)))
  covariates <- d$data[d$full2short, d$covariates_names, drop = FALSE]
  strata <- make_patterns(covariates, d$pattern_weights, any_continuous = FALSE)
  list(npatterns = strata$npatterns, total = strata$pattern_weights[strata$short2full])
}

#' @rdname getfit.llca
#' @param model An object of class \code{"llcalist"} containing fitted
#'   \code{"llca"} or \code{"llca_sam"} objects.
#' @method getfit llcalist
#' @export
getfit.llcalist <- function(x, digits = 4L, ..., model = x) {
  if(!missing(model) && !missing(x)) return(.match_legacy_model("model"))

  #### Check inputs ####

  if(!inherits(model, "llcalist")) {
    stop("model must inherit from class 'llcalist'.")
  }

  if(length(model) == 0L) {
    stop("model must contain at least one llca object.")
  }

  if(!is.null(digits) &&
     (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
      digits < 0L || digits != as.integer(digits))) {
    stop("digits must be NULL or a non-negative integer.")
  }

  #### Extract fit indices ####

  result <- lapply(model, FUN = getfit, digits = digits)

  class(result) <- "getfit.llcalist"

  #### Result ####

  return(result)

}
