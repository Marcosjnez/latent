# Author: Marcos Jimenez
# Author: Mauricio Garnier-Villarreal
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 03/09/2025
#'
#' @title
#' Fit indices
#' @description
#'
#' Compute fit indices from any model.
#'
#' @usage
#'
#' getfit(model)
#'
#' @param model data.frame or matrix of response.
#'
#' @details \code{getfit} computes all the fit indices related to a specific model.
#'
#' @return List with the following fit indices:
#' \item{AIC}{.}
#' \item{BIC}{.}
#'
#' @references
#'
#' None yet.
#'
#' @method getfit llca
#' @export
getfit.llca <- function(object, digits = 3) {
  if(is(object) == "llca" ){
    out <- getfit0(object, digits = digits)
  }
  if(is(object) == "llca.list"){
    out <- t(sapply(fit, getfit0, digits=digits))
  }

  class(out) <- "getfit.llca"

  return(out)

}

getfit0 <- function(model, digits = 3) {

  icl_default <- function (post_prob, BIC){
    tryCatch({
      if (!is.null(dim(post_prob))) {
        C <- post_prob == apply(post_prob, 1, max)
        (-1 * BIC) + 2 * sum(C * log(apply(post_prob, 1,
                                           function(x) {
                                             x[which.max(x)]
                                           }) + 1e-12))
      }
      else {
        (-1 * BIC) + 2 * sum(post_prob * log(post_prob))
      }
    }, error = function(e) {
      NA
    })
  }

  entropy <- function(p) {
    p <- p[p > sqrt(.Machine$double.eps)]
    sum(-p * log(p))
  }

  entropy.R2 <- function(prop, post) {
    error_prior <- entropy(prop)
    error_post <- mean(apply(post, 1, entropy))
    R2_entropy <- (error_prior - error_post)/error_prior
    R2_entropy
  }

  ##
  penalized <- isFALSE(fit@Optim$control$penalties) == FALSE

  ##
  ##
  nclasses <- ncol(model@posterior)
  k <- model@modelInfo$nparam
  loglik <- model@loglik


  ##
  # if(sum(model@modelInfo$item != "multinomial") == 0){
  if(all(model@modelInfo$item == "multinomial")){
    ni <- model@summary_table$Observed
    mi <- model@summary_table$Estimated
    dof <- model@modelInfo$dof
    term <- ni*log(ni/mi)
    # Set to zero Inf terms due to 0 values in Estimated
    term[is.infinite(term)] <- 0
    L2 <- 2*sum(term)
    pv <- 1-pchisq(L2, dof)
  }else{
    L2 <- NA
    pv <- NA
    dof <- NA
  }

  nobs <- model@modelInfo$nobs
  ###
  AIC  <- (-2 * loglik) + (k * 2)
  AIC3 <- (-2 * loglik) + (k * 3)
  BIC <- k * log(sum(nobs)) - 2 * loglik
  CAIC <- -2 * loglik + (k * (log(nobs) + 1))
  KIC <- -2 * loglik + (3 * (k + 1))
  SABIC <- -2 * loglik + (k * log(((nobs + 2)/24)))
  ICL <- icl_default(model@posterior, BIC )

  if(nclasses < 2) {
    entropyR2 <- 1.00 # To match LG output
  } else {
    entropyR2 <- entropy.R2(model@transformed_pars$class,
                            model@posterior)
  }

  if(penalized){
    penalized_loglik <- model@penalized_loglik

    AICp  <- (-2 * penalized_loglik) + (k * 2)
    AIC3p <- (-2 * penalized_loglik) + (k * 3)
    BICp <- k * log(sum(nobs)) - 2 * penalized_loglik
    CAICp <- -2 * penalized_loglik + (k * (log(nobs) + 1))
    KICp <- -2 * penalized_loglik + (3 * (k + 1))
    SABICp <- -2 * penalized_loglik + (k * log(((nobs + 2)/24)))
    ICLp <- icl_default(model@posterior, BICp )

    result <- c(nclasses = nclasses,
                npar = k, nobs = nobs,
                loglik = loglik,
                penalized_loglik = penalized_loglik,
                L2 = L2,
                dof = dof, pvalue = pv,
                AIC = AIC, BIC = BIC, AIC3 = AIC3,
                CAIC = CAIC, KIC = KIC, SABIC = SABIC,
                ICL = ICL,
                AICp = AICp, BICp = BICp, AIC3p = AIC3p,
                CAICp = CAICp, KICp = KICp, SABICp = SABICp,
                ICLp = ICLp,
                R2_entropy = entropyR2)

  }else{
    result <- c(nclasses = nclasses,
                npar = k, nobs = nobs,
                loglik = loglik,
                L2 = L2,
                dof = dof, pvalue = pv,
                AIC = AIC, BIC = BIC, AIC3 = AIC3,
                CAIC = CAIC, KIC = KIC, SABIC = SABIC,
                ICL = ICL,
                R2_entropy = entropyR2)
  }


  return(round(result, digits = digits))

}
