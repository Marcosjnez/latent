# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 06/10/2025
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

    for(j in 1:length(temp)){
      temp[[j]] <- round(temp[[j]], digits)
    }

    return(list(classes, temp))

  } else if (what == "classconditional" ||
             what == "items" ||
             what == "item") {

    temp <- fit@ClassConditional

    for(j in 1:length(temp)){
      temp[[j]] <- round(temp[[j]], digits)
    }

    return(temp)

  } else if (what == "class" ||
             what == "classes" ||
             what == "cluster" ||
             what == "clusters") {

    weights <- fit@Optim$data_list$weights
    classes <- colSums(fit@transformed_pars$class * weights) / sum(weights)
    temp <- round(classes, digits)
    names(temp) <- paste0("Class", 1:length(temp))

    return(temp)

  } else if (what == "fullclasses") {

    temp <- round(fit@transformed_pars$class, digits)
    colnames(temp) <- paste0("Class", 1:ncol(temp))
    rownames(temp) <- rownames(fit@Optim$data)

    return(temp)

  } else if (what == "beta"  ||
             what == "betas" ||
             what == "coef"  ||
             what == "coefs" ||
             what == "coef"  ||
             what == "coefficient" ||
             what == "coefficients") {

    temp <- round(fit@transformed_pars$beta, digits)
    colnames(temp) <- paste0("Class", 1:ncol(temp))
    rownames(temp) <- colnames(fit@Optim$data_list$cov_patterns)

    return(temp)

  } else if (what == "respconditional") {

    temp <- fit@RespConditional

    for(j in 1:length(temp)){
      temp[[j]] <- round(temp[[j]], digits)
    }

    return(temp)

  } else if (what == "convergence") {

    fit@Optim$convergence

  } else if (what == "data") {

    fit@Optim$data

  } else if (what == "pattern") {

    cbind(fit@Optim$data_list$patterns,
          times = fit@Optim$data_list$weights)

  } else if (what == "table") {

    fit@summary_table

  } else if (what == "posterior") {

    round(fit@posterior, digits)

  } else if (what == "state") {

    fit@state

  }else if (what == "loglik_case" ||
            what == "casell") {

    fit@loglik_case

  } else if (what == "loglik" ||
             what == "ll" ||
             what == "LL") {

    fit@loglik

  } else if (what == "probcat") {

    temp <- fit@probCat

    for(j in 1:length(temp)){
      temp[[j]] <- round(temp[[j]], digits)
    }

  } else if (what == "timing") {

    fit@timing

  }


}

