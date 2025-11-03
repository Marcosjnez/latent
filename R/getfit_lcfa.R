# Author: Marcos Jimenez
# Author: Mauricio Garnier-Villarreal
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 03/11/2025
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
#' @method getfit lcfa
#' @export
getfit.lcfa <- function(model, digits = 3) {

  nobs <- sum(unlist(model@modelInfo$nobs))
  nitems <- model@Optim$data_list$nitems
  nfactors <- model@Optim$data_list$nfactors
  nparam <- model@modelInfo$nparam
  dof <- model@modelInfo$dof
  S <- model@Optim$data_list$correl[[1]]$R

  # Compute fit statistics if using ML:
  if(length(model@loglik) == 0) {

    est <- "ULS"
    # loglik <- NULL
    # X2 <- NULL
    # pval <- NULL
    # CFI <- NULL
    # RMSEA <- NULL

    loglik <- model@loss
    # X2 <- loglik*(nobs-1L)
    X2 <- model@loss*(nobs-1L)
    pval <- NULL
    # Compute loss value
    Sigma <- diag(nitems) # Identity model
    F_id <- 0.5*sum((S - Sigma)^2)
    X2_id <- F_id*(nobs-1L)
    dof_id <- nitems
    t1 <- max(c(X2 - dof, 0))
    t2 <- max(c(X2 - dof, X2_id - dof_id, 0))
    CFI <- 1-t1/t2
    RMSEA <- sqrt(max(c(loglik/dof - 1/nobs, 0)))

  } else {

    est <- "ML"
    loglik <- model@loss
    X2 <- loglik*nobs
    pval <- 1-pchisq(X2, df = dof)
    # Compute loss value
    Sigma <- diag(nitems) # Identity model
    F_id <- log(det(Sigma)) - log(det(S))
    X2_id <- F_id*nobs
    dof_id <- nitems
    t1 <- max(c(X2 - dof, 0))
    t2 <- max(c(X2 - dof, X2_id - dof_id, 0))
    CFI <- 1-t1/t2
    RMSEA <- sqrt(max(c(loglik/dof - 1/nobs, 0)))

  }

  lambda <- model@transformed_pars[[1]]$lambda
  psi <- model@transformed_pars[[1]]$psi
  theta <- model@transformed_pars[[1]]$theta
  rhat <- lambda %*% psi %*% t(lambda) + theta
  residuals <- (S-rhat)[lower.tri(S, diag = TRUE)]
  SRMR <- sqrt(mean(residuals^2))

  result <- c(nfactors = nfactors,
              npar = nparam, nobs = nobs,
              loglik = loglik,
              chisq = X2,
              dof = dof,
              pvalue = pval,
              cfi = CFI,
              rmsea = RMSEA,
              srmr = SRMR)

  class(result) <- "getfit.lcfa"

  return(round(result, digits = digits))

}
