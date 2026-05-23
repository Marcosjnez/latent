# Author: Marcos Jimenez
# Author: Mauricio Garnier-Villarreal
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 05/05/2026
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

  nobs <- sum(unlist(model@dataList$nobs))
  nitems <- model@dataList$nitems[[1]]
  nfactors <- model@dataList$nfactors[[1]]
  nparam <- model@modelInfo$nparam
  dof <- model@modelInfo$dof

  # Compute fit statistics if using ML:
  if(model@modelInfo$control_optimizer$reg) {

    fit_mat <- latInspect(model, what = "loss")

    loglik <- NA
    X2 <- -fit_mat["loss", "overall"]*(nobs-1L)
    pval <- NULL
    F_id <- fit_mat["loss_base", "overall"]
    X2_id <- -F_id*(nobs-1L)
    dof_id <- nitems
    t1 <- max(c(X2 - dof, 0))
    t2 <- max(c(X2 - dof, X2_id - dof_id, 0))
    CFI <- 1-t1/t2
    RMSEA <- sqrt(max(c((X2-dof)/(dof*(nobs-1)), 0)))

  } else {

    fit_mat <- latInspect(model, what = "loglik")

    llsat <- fit_mat["loglik_sat", "overall"]
    ll <- fit_mat["loglik", "overall"]
    llbas <- fit_mat["loglik_base", "overall"]
    X2 <- 2*(llsat - ll)
    pval <- 1-pchisq(X2, df = dof)
    X2_id <- 2*(llsat - llbas)
    dof_id <- nitems
    t1 <- max(c(X2 - dof, 0))
    t2 <- max(c(X2 - dof, X2_id - dof_id, 0))
    CFI <- 1-t1/t2
    RMSEA <- sqrt(max(c((X2-dof)/(dof*(nobs-1)), 0)))

  }

  # resids <- latInspect(model, what = "residuals")
  # SRMR_list <- lapply(resids, FUN = \(x) {
  #   x[lower.tri(x, diag = TRUE)]
  # })
  # SRMR <- sqrt(mean(unlist(SRMR_list)^2))
  SRMR <- NA

  result <- c(nfactors = nfactors,
              npar = nparam,
              nobs = nobs,
              loglik = model@loglik,
              chisq = X2,
              dof = dof,
              pvalue = pval,
              cfi = CFI,
              rmsea = RMSEA,
              srmr = SRMR)

  class(result) <- "getfit.lcfa"

  return(round(result, digits = digits))

}
