# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 22/05/2025
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
#' @export
getfit0 <- function(model, digits=3){

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
  ##
  nclasses <- ncol(model@posterior)
  k <- model@modelInfo$nparam
  loglik <- model@loglik
  penalized_loglik <- model@penalized_loglik
  if(loglik == penalized_loglik){
    penalized_loglik <- NA
  }

  ##
  if(sum(model@modelInfo$item != "multinomial") == 0){
    ni <- model@summary_table$Observed
    mi <- model@summary_table$Estimated
    df <- model@modelInfo$df
    L2 <- 2*sum(ni*log(ni/mi))
    pv <- 1-pchisq(L2, df)
  }else{
    L2 <- NA
    pv <- NA
    df <- NA
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


  entropyR2 <- entropy.R2(model@transformed_pars$classes,
                          model@posterior)

  result <- c(nclasses = nclasses,
              npar = k, nobs = nobs,
              loglik = loglik,
              penalized_loglik = penalized_loglik,
              L2 = L2,
              df = df, pvalue = pv,
              AIC = AIC, BIC = BIC, AIC3 = AIC3,
              CAIC = CAIC, KIC = KIC, SABIC = SABIC,
              ICL = ICL,
              R2_entropy = entropyR2)

  return(round(result, digits) )

}

getfit <- function(object, digits=3){
  if(is(object) == "llca" ){
    out <- getfit0(object, digits = digits)
  }
  if(is(object) == "llca.list"){
    out <- t(sapply(fit, getfit0, digits=digits))
  }
  return(out)
}
