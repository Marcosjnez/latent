# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 11/07/2026
#'
#' Fit indices for latent class models
#'
#' @description
#' Computes likelihood-based fit indices, classification measures, and entropy
#' for a fitted latent class model or for a list of fitted latent class models.
#'
#' @usage
#' \method{getfit}{llca}(model, digits = 4)
#'
#' \method{getfit}{llcalist}(model, digits = 4)
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
#' p-value are returned only when all modeled variables use a multinomial
#' likelihood. They are returned as \code{NA} for models containing Gaussian
#' variables.
#'
#' When penalization is active, the output also contains the penalized
#' log-likelihood and the corresponding penalized information criteria. The
#' attribute \code{"penalized"} indicates whether these additional indices are
#' present.
#'
#' For an \code{"llcalist"} object, the indices are collected into a matrix
#' with one row per fitted model.
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
#' For an object of class \code{"llcalist"}, a numeric matrix of class
#' \code{"getfit.llcalist"} with one row per model and one column per available
#' fit index.
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
#' @method getfit llca
#' @export
getfit.llca <- function(model, digits = 4L) {

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

  #### Auxiliary functions ####

  return_na <- function(e) {

    #### Result ####

    return(NA_real_)

  }

  icl_default <- function(post_prob, BIC) {

    result <- tryCatch({

      if(!is.null(dim(post_prob))) {
        max_prob <- apply(post_prob, MARGIN = 1L, FUN = max)
        classification <- post_prob == max_prob
        value <- -BIC + 2*sum(classification*log(max_prob + 1e-12))
      } else {
        value <- -BIC + 2*sum(post_prob*log(post_prob))
      }

      value

    }, error = return_na)

    #### Result ####

    return(result)

  }

  entropy <- function(p) {

    p <- p[p > sqrt(.Machine$double.eps)]
    result <- sum(-p*log(p))

    #### Result ####

    return(result)

  }

  entropy.R2 <- function(prop, post) {

    error_prior <- entropy(prop)
    error_post <- mean(apply(post, MARGIN = 1L, FUN = entropy))
    result <- (error_prior-error_post)/error_prior

    #### Result ####

    return(result)

  }

  #### Extract model information ####

  penalized <- model@modelInfo$control_optimizer$reg
  posterior <- latInspect(model, what = "posterior")
  summary_table <- latInspect(model, what = "summary")
  fit_loglik <- latInspect(model, what = "loglik")

  nclasses <- ncol(posterior)
  k <- model@modelInfo$nparam
  nobs <- model@dataList$nobs
  loglik <- fit_loglik["loglik", "overall"]
  penalized_loglik <- fit_loglik["penalized_loglik", "overall"]

  #### Likelihood-ratio test ####

  if(all(model@dataList$variable_type == "multinomial")) {

    ni <- summary_table$Observed
    mi <- summary_table$Estimated
    dof <- model@modelInfo$dof
    term <- ni*log(ni/mi)
    term[is.infinite(term)] <- 0
    L2 <- 2*sum(term)
    pv <- 1-pchisq(L2, dof)

  } else {

    L2 <- NA_real_
    pv <- NA_real_
    dof <- model@modelInfo$dof

  }

  #### Information criteria ####

  AIC <- (-2*loglik) + (k*2)
  AIC3 <- (-2*loglik) + (k*3)
  BIC <- k*log(nobs) - 2*loglik
  CAIC <- -2*loglik + (k*(log(nobs) + 1))
  KIC <- -2*loglik + (3*(k + 1))
  SABIC <- -2*loglik + (k*log((nobs + 2)/24))
  ICL <- icl_default(posterior, BIC)

  #### Entropy ####

  if(nclasses < 2L) {

    entropyR2 <- 1.00

  } else {

    pattern_weights <- model@dataList$pattern_weights
    classes <- colSums(model@transformed_pars$class*pattern_weights) /
      sum(pattern_weights)
    entropyR2 <- entropy.R2(classes, posterior)

  }

  #### Collect fit indices ####

  if(penalized) {

    AICp <- (-2*penalized_loglik) + (k*2)
    AIC3p <- (-2*penalized_loglik) + (k*3)
    BICp <- k*log(nobs) - 2*penalized_loglik
    CAICp <- -2*penalized_loglik + (k*(log(nobs) + 1))
    KICp <- -2*penalized_loglik + (3*(k + 1))
    SABICp <- -2*penalized_loglik + (k*log((nobs + 2)/24))
    ICLp <- icl_default(posterior, BICp)

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

#' @rdname getfit.llca
#' @method getfit llcalist
#' @export
getfit.llcalist <- function(model, digits = 4L) {

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

  out <- t(vapply(model, FUN = getfit.llca,
                  FUN.VALUE = getfit.llca(model[[1]], digits = digits),
                  digits = digits))

  nclasses <- vapply(model, FUN = function(x) {

    result <- ncol(x@modelInfo$trans$class)

    #### Result ####

    return(result)

  }, FUN.VALUE = integer(1L))

  rownames(out) <- paste0("nclasses=", nclasses)

  class(out) <- "getfit.llcalist"
  attr(out, "penalized") <- model[[1]]@modelInfo$control_optimizer$reg

  #### Result ####

  return(out)

}
