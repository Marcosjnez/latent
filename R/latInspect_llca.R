# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 04/05/2026
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
                            what = "profile") {

  # fit must inherit from class llca
  stopifnot(inherits(fit, "llca"))

  # be case insensitive
  what <- tolower(what)

  list2env(fit@dataList, envir = environment())
  list2env(fit@modelInfo, envir = environment())

  # Logarithm likelihood of each response pattern:
  loglik_case <- fit@Optim$outputs$estimators$vectors[[1]][[1]]
  # Sum of logarithm likelihoods by response pattern:
  loglik_pattern <- weights * loglik_case

  # Estimated "counts" for each response pattern:
  estimated <- exp(loglik_case) * nobs
  # Posterior:
  posterior <- exp(matrix(fit@Optim$outputs$estimators$matrices[[1]][[2]],
                          nrow = npatterns, ncol = nclasses))
  colnames(posterior) <- paste("P(", "Class", 1:nclasses, "|data)", sep = "")
  # Posterior classification:
  state <- apply(posterior, MARGIN = 1, FUN = which.max)
  # Data table of response patterns:
  summary_table <- cbind(Pattern = patterns + 1,
                         Observed = weights,
                         Estimated = estimated,
                         Posterior = posterior,
                         State = state,
                         loglik_case = loglik_case,
                         loglik_pattern = loglik_pattern)
  summary_table <- as.data.frame(summary_table)
  # Sort the patterns by increasing order:
  summary_table <- summary_table[do.call(order, as.data.frame(summary_table)), ]
  rownames(summary_table) <- paste("pattern", 1:nrow(summary_table), sep = "")

  # Check the existence of gaussian items:
  gauss <- "gaussian" %in% item
  # Check the existence of multinomial items:
  multin <- "multinomial" %in% item

  classes <- colSums(fit@transformed_pars$class * weights) / sum(weights)
  ClassConditional <- fit@transformed_pars[item_names]
  RespConditional <- probCat <- list() # Only for full multinomial models

  # Additional outputs for full multinomial models:
  if(all(item == "multinomial")) {

    classes <- colMeans(fit@transformed_pars$class)
    conditionals <- ClassConditional

    probCat <- lapply(conditionals, FUN = \(mat) {
      # Calculate P(y|X)*P(X), the joint probability:
      jointp <- t(mat) * classes
      # Calculate P(y), the denominator of the posterior:
      probCat <- colSums(jointp)
      return(probCat)
    })

    RespConditional <- lapply(conditionals, FUN = \(mat) {
      # Calculate P(y|X)*P(X), the joint probability:
      jointp <- t(mat) * classes
      # Calculate P(y), the denominator of the posterior:
      probCat <- colSums(jointp)
      # Calculate P(X|y) = P(y|X)*P(X)/P(y), the posterior:
      posterior <- t(jointp) / probCat
      return(posterior)
    })

    names(RespConditional) <- colnames(data)

  }

  loglik_case <- loglik_case[short2full]
  posterior <- posterior[short2full, , drop = FALSE]
  state <- state[short2full]
  rownames(posterior) <- names(state) <-
    names(loglik_case) <- rownames(data)

  profile <- list(class = classes, item = ClassConditional)

  #### Result ####

  if (what == "profile") {

    return(profile)

  } else if (what == "classconditional" ||
             what == "items" ||
             what == "item") {

    return(ClassConditional)

  } else if (what == "class" ||
             what == "classes" ||
             what == "cluster" ||
             what == "clusters") {

    return(classes)

  } else if (what == "fullclasses") {

    return(fit@transformed_pars$class)

  } else if (what == "beta"  ||
             what == "betas" ||
             what == "coef"  ||
             what == "coefs" ||
             what == "coef"  ||
             what == "coefficient" ||
             what == "coefficients") {

    return(fit@transformed_pars$beta)

  } else if (what == "respconditional") {

    return(RespConditional)

  } else if (what == "convergence") {

    return(fit@Optim$convergence)

  } else if (what == "data") {

    return(fit@dataList$data)

  } else if (what == "pattern") {

    return(cbind(patterns, times = weights))

  } else if (what == "table" ||
             what == "summary") {

    return(summary_table)

  } else if (what == "posterior") {

    return(posterior)

  } else if (what == "state") {

    return(state)

  } else if (what == "loglik_case" ||
             what == "loglik.case" ||
             what == "case") {

    return(loglik_case)

  } else if (what == "loglik_pattern" ||
             what == "loglik.pattern" ||
             what == "pattern") {

    return(loglik_pattern)

  } else if (what == "loglik" ||
             what == "ll" ||
             what == "LL") {

    return(fit@loglik)

  } else if (what == "probcat") {

    return(probCat)

  } else if (what == "timing" ||
             what == "elapsed") {

    return(fit@Optim$elapsed)

  }

}

#' @method latInspect llcalist
#' @export
latInspect.llcalist <- function(model, what = "profile",
                                digits = 3) {

  nmodels <- length(model)
  out <- vector("list", length = nmodels)
  for(i in 1:nmodels) {

    out[[i]] <- latInspect.llca(model[[i]], what = what, digits = digits)
    names(out)[i] <- paste("nclasses = ", model[[i]]@datalist$nclasses,
                           sep = "")

  }

  class(out) <- "latInspect.llcalist"

  return(out)

}
