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
getfit <- function(model) {

  entropy <- function(p) {
    p <- p[p > sqrt(.Machine$double.eps)] # since Lim_{p->0} p log(p) = 0
    sum(-p * log(p))
  }
  entropy.R2 <- function(prop, post) {
    error_prior <- entropy(prop) # Class proportions
    error_post <- mean(apply(post, 1, entropy))
    R2_entropy <- (error_prior - error_post) / error_prior
    R2_entropy
  }

  df <- fit$modelInfo$df
  k <- length(fit$modelInfo$parameters)
  loglik <- sum(model$logliks)
  AIC <- 2*k - 2*loglik
  BIC <- k*log(sum(fit$modelInfo$nobs)) - 2*loglik
  entropyR2 <- entropy.R2(model$classes, model$posterior[, 1:length(model$classes)])

  result <- list(loglik = loglik, df = df, AIC = AIC, BIC = BIC,
                 "entropy R squared" = entropyR2)

  return(result)

}
