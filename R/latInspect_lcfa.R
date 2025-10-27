# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 29/08/2025
#'
#' @title
#' Standard Errors
#' @description
#'
#' Compute standard errors.
#'
#' @usage
#'
#' se(fit)
#'
#' @param fit model fitted with lca.
#' @param confidence Coverage of the confidence interval.
#'
#' @details Compute standard errors.
#'
#' @return List with the following objects:
#' \item{vcov}{Variance-covariance matrix between the parameters.}
#' \item{se}{Standard errors.}
#' \item{SE}{Standard errors in the model list.}
#'
#' @references
#'
#' None yet.
#'
#' @method latInspect lcfa
#' @export
latInspect.lcfa <- function(fit,
                            what = "est",
                            digits = 3) {

  # fit must inherit from class llca
  stopifnot(inherits(fit, "lcfa"))

  # be case insensitive
  what <- tolower(what)

  lambda <- lapply(fit@transformed_pars, FUN = \(x) round(x$lambda, digits = digits))
  psi <- lapply(fit@transformed_pars, FUN = \(x) round(x$psi, digits = digits))
  theta <- lapply(fit@transformed_pars, FUN = \(x) round(x$theta, digits = digits))

  ngroups <- fit@Optim$data_list$ngroups
  rhat <- resids <- vector("list", length = ngroups)
  for(i in 1:ngroups) {

    p <- sqrt(length(fit@Optim$opt$outputs$estimators$matrices[[i]][[4]]))
    temp <- fit@Optim$opt$outputs$estimators$matrices[[i]][[4]]
    rhat[[i]] <- matrix(round(temp, digits = digits), nrow = p, ncol = p)
    temp <- fit@Optim$opt$outputs$estimators$matrices[[i]][[5]]
    resids[[i]] <- matrix(round(temp, digits = digits), nrow = p, ncol = p)

  }

  if(what == "est" ||
     what == "estimates" ||
     what == "parameters" ||
     what == "fixed" ||
     what == "items") {

    groups <- vector("list", length = ngroups)
    for(i in 1:ngroups) {
      groups[[i]]$lambda <- round(lambda[[i]], digits = digits)
      groups[[i]]$psi <- round(psi[[i]], digits = digits)
      groups[[i]]$theta <- round(theta[[i]], digits = digits)
    }

    return(groups)

  } else if(what == "rhat" ||
            what == "model") {

    return(rhat)

  } else if(what == "resid" ||
            what == "residuals") {

    return(resids)

  } else if(what == "lambda" ||
            what == "loadings") {

    return(lambda)

  } else if(what == "psi") {

    return(psi)

  } else if(what == "theta") {

    return(theta)

  } else if(what == "uniquenesses") {

    u <- lapply(theta, FUN = \(x) diag(x))
    return(u)

  } else {

    stop("Unknown request")

  }

}
