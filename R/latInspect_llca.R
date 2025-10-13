# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 13/10/2025
#'
#' @title
#' Inspect objects from fitted lca models.
#' @description
#'
#' Inspect objects.
#'
#' @usage
#'
#' latInspect(fit)
#'
#' @param fit model fitted with lca.
#' @param digits Number of digits to print.
#'
#' @details Extract objects from fitted lca object.
#'
#' @return List with one of the following objects:
#' \item{class}{.}
#' \item{posterior}{.}
#' \item{profile}{.}
#'
#' @references
#'
#' None yet.
#'
#' @method latInspect llca
#' @export
latInspect.llca <- function(fit,
                            what = "classconditional",
                            digits = 3) {

  # fit must inherit from class llca
  stopifnot(inherits(fit, "llca"))

  # be case insensitive
  what <- tolower(what)

  ### result fits

  if (what == "profile") {

    weights <- fit@Optim$data_list$weights
    classes <- colSums(fit@transformed_pars$class * weights) / sum(weights)
    temp <- fit@ClassConditional

    if(!is.null(digits)) {
      for(j in 1:length(temp)) {
        temp[[j]] <- round(temp[[j]], digits = digits)
      }
      classes <- round(classes, digits = digits)
    }

    return(list(class = classes, item = temp))

  } else if (what == "classconditional" ||
             what == "items" ||
             what == "item") {

    temp <- fit@ClassConditional

    if(!is.null(digits)) {
      for(j in 1:length(temp)) {
        temp[[j]] <- round(temp[[j]], digits = digits)
      }
    }

    return(temp)

  } else if (what == "class" ||
             what == "classes" ||
             what == "cluster" ||
             what == "clusters") {

    weights <- fit@Optim$data_list$weights
    temp <- colSums(fit@transformed_pars$class * weights) / sum(weights)
    names(temp) <- paste0("Class", 1:length(temp))

    if(!is.null(digits)) {
      temp <- round(temp, digits = digits)
    }

    return(temp)

  } else if (what == "fullclasses") {

    temp <- fit@transformed_pars$class
    colnames(temp) <- paste0("Class", 1:ncol(temp))
    rownames(temp) <- rownames(fit@Optim$data)

    if(!is.null(digits)) {
      temp <- round(temp, digits = digits)
    }

    return(temp)

  } else if (what == "beta"  ||
             what == "betas" ||
             what == "coef"  ||
             what == "coefs" ||
             what == "coef"  ||
             what == "coefficient" ||
             what == "coefficients") {

    temp <- fit@transformed_pars$beta
    colnames(temp) <- paste0("Class", 1:ncol(temp))
    rownames(temp) <- colnames(fit@Optim$data_list$cov_patterns)

    if(!is.null(digits)) {
      temp <- round(temp, digits = digits)
    }

    return(temp)

  } else if (what == "respconditional") {

    temp <- fit@RespConditional

    for(j in 1:length(temp)){
      temp[[j]] <- round(temp[[j]], digits)
    }

    return(temp)

  } else if (what == "convergence") {

    return(fit@Optim$convergence)

  } else if (what == "data") {

    return(fit@Optim$data)

  } else if (what == "pattern") {

    temp <- cbind(fit@Optim$data_list$patterns,
                  times = fit@Optim$data_list$weights)

    return(temp)

  } else if (what == "table") {

    return(fit@summary_table)

  } else if (what == "posterior") {

    temp <- fit@posterior

    if(!is.null(digits)) {
      temp <- round(temp, digits = digits)
    }

    return(temp)

  } else if (what == "state") {

    fit@state

  }else if (what == "loglik_case" ||
            what == "casell") {

    return(fit@loglik_case)

  } else if (what == "loglik" ||
             what == "ll" ||
             what == "LL") {

    return(fit@loglik)

  } else if (what == "probcat") {

    temp <- fit@probCat

    if(!is.null(digits)) {
      for(j in 1:length(temp)) {
        temp[[j]] <- round(temp[[j]], digits = digits)
      }
    }

  } else if (what == "timing") {

    return(fit@timing)

  }

}

#' @method latInspect llcalist
#' @export
latInspect.llcalist <- function(model, what = "classconditional",
                                digits = 3) {

  nmodels <- length(model)
  out <- vector("list", length = nmodels)
  for(i in 1:nmodels) {

    out[[i]] <- latInspect.llca(model[[i]], what = what, digits = digits)
    names(out)[i] <- paste("nclasses = ", model[[i]]@modelInfo$nclasses,
                           sep = "")

  }

  class(out) <- "latInspect.llcalist"

  return(out)

}
