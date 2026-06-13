# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 13/06/2026
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
                            what = "profile",
                            digits = 4L) {

  # fit must inherit from class llca
  stopifnot(inherits(fit, "llca"))

  # be case insensitive
  what <- tolower(what)

  list2env(fit@dataList, envir = environment())
  list2env(fit@modelInfo, envir = environment())
  nclasses <- ncol(trans$class)

  # Logarithm likelihood of each response pattern:
  loglik_case <- fit@Optim$outputs$estimators$vectors[[1]][[1]]
  pattern_weights <- fit@dataList$pattern_weights
  # Sum of logarithm likelihoods by response pattern:
  loglik_patterns <- pattern_weights * loglik_case

  # Estimated "counts" for each response pattern:
  estimated <- exp(loglik_case) * nobs
  # Posterior:
  posterior <- exp(matrix(fit@Optim$outputs$estimators$matrices[[1]][[1]],
                          nrow = npatterns, ncol = nclasses))
  colnames(posterior) <- paste("P(", "Class", 1:nclasses, "|data)", sep = "")
  # Posterior classification:
  state <- apply(posterior, MARGIN = 1, FUN = which.max)
  # Data table of response patterns:
  summary_table <- cbind(patterns_original[, -1],
                         Observed = pattern_weights,
                         Estimated = estimated,
                         Posterior = posterior,
                         State = state,
                         loglik_case = loglik_case,
                         loglik_patterns = loglik_patterns)
  summary_table <- as.data.frame(summary_table)
  # Sort the patterns by increasing order:
  summary_table <- summary_table[do.call(order, summary_table), ]
  rownames(summary_table) <- paste("pattern", 1:nrow(summary_table), sep = "")

  # Check the existence of gaussian items:
  gauss <- "gaussian" %in% item
  # Check the existence of multinomial items:
  multin <- "multinomial" %in% item

  classes <- colSums(fit@transformed_pars$class * pattern_weights) /
    sum(pattern_weights)
  ClassConditional <- fit@transformed_pars[item_names]
  RespConditional <- probCat <- list() # Only for full multinomial models

  # Additional outputs for full multinomial models:
  if(all(item == "multinomial")) {

    classes <- colMeans(fit@transformed_pars$class)
    probCat <- lapply(ClassConditional, FUN = \(mat) {
      # Calculate P(y|X)*P(X), the joint probability:
      jointp <- t(mat) * classes
      # Calculate P(y), the denominator of the posterior:
      probCat <- colSums(jointp)
      return(probCat)
    })

    RespConditional <- lapply(ClassConditional, FUN = \(mat) {
      # Calculate P(y|X)*P(X), the joint probability:
      jointp <- t(mat) * classes
      # Calculate P(y), the denominator of the posterior:
      probCat <- colSums(jointp)
      # Calculate P(X|y) = P(y|X)*P(X)/P(y), the posterior:
      posterior <- t(jointp) / probCat
      return(posterior)
    })

    names(RespConditional) <- colnames(measurement)

  }

  loglik_case <- loglik_case[short2full]
  posterior <- posterior[short2full, , drop = FALSE]
  state <- state[short2full]
  rownames(posterior) <- names(state) <-
    names(loglik_case) <- rownames(measurement)

  if(fit@dataList$any_gaussian) {
    gaussian_names <- fit@modelInfo$control_optimizer$gaussian$gaussian_names
    ClassConditional[gaussian_names] <- lapply(ClassConditional[gaussian_names],
                                               FUN = \(x) {
                                                 x[c(1, 4), ]
                                               })
  }
  profile <- list(class_size = colMeans(posterior), item = ClassConditional)

  #### Extract the fit ####

  doubles <- fit@Optim$outputs$estimators$doubles

  fit_matrix <- vapply(
    doubles,
    FUN = function(x) {
      c(
        loss             = x[[1]],
        loss_base        = x[[2]],
        loss_sat         = x[[3]],
        loglik           = x[[4]],
        loglik_base      = x[[5]],
        loglik_sat       = x[[6]],
        penalty          = x[[7]],
        penalized_loss   = x[[1]] + x[[7]],
        penalized_loglik = x[[4]] - x[[7]]
      )
    },
    FUN.VALUE = numeric(9)
  )

  colnames(fit_matrix) <- unlist(lapply(fit@modelInfo$control_estimator,
                                        FUN = \(x) x$double_names))
  fit_matrix <- cbind(fit_matrix, overall = rowSums(fit_matrix))

  #### Result ####

  if (what == "profile") {

    result <- profile

  } else if (what == "classconditional" ||
             what == "items" ||
             what == "item") {

    result <- ClassConditional

  } else if (what == "class" ||
             what == "classes" ||
             what == "cluster" ||
             what == "clusters") {

    result <- classes

  } else if (what == "fullclasses") {

    result <- fit@transformed_pars$class

  } else if (what == "beta"  ||
             what == "betas" ||
             what == "coef"  ||
             what == "coefs" ||
             what == "coef"  ||
             what == "coefficient" ||
             what == "coefficients") {

    result <- list(beta        = fit@transformed_pars$beta,
                   distal_beta = fit@transformed_pars$distal_beta)

  } else if (what == "respconditional") {

    result <- RespConditional

  } else if (what == "convergence") {

    conv <- data.frame(Iterations = fit@Optim$iterations,
                       convergence = fit@Optim$convergence,
                       grad.norm = fit@Optim$ng)

    return(conv)

  } else if (what == "gradient") {

    gradient.info <- data.frame(name = fit@modelInfo$parameters_labels,
                                gradient = fit@Optim$g,
                                rgradient = fit@Optim$rg,
                                dir = fit@Optim$dir)

    return(list(grad.norm = fit@Optim$ng, gradient.info = gradient.info))

  } else if (what == "data") {

    result <- fit@dataList$data

  } else if (what == "measurement") {

    result <- fit@dataList$measurement

  } else if (what == "pattern" ||
             what == "patterns") {

    patterns <- data.frame(patterns_original[, -1], Observed = pattern_weights)
    # Sort the patterns by increasing order:
    patterns <- patterns[do.call(order, patterns), ]
    rownames(patterns) <- paste("pattern", 1:nrow(patterns), sep = "")

    result <- patterns

  } else if (what == "table" ||
             what == "summary") {

    result <- summary_table

  } else if (what == "posterior") {

    result <- posterior

  } else if (what == "state") {

    return(state)

  } else if (what == "loglik_case" ||
             what == "loglik.case" ||
             what == "case") {

    result <- loglik_case

  } else if (what == "loglik_pattern" ||
             what == "loglik.pattern" ||
             what == "pattern") {

    result <- loglik_patterns

  } else if (what == "loss" ||
             what == "losses") {

    result <- fit_matrix[c("loss", "penalized_loss"), "overall", drop = FALSE]

  } else if (what == "loglik" ||
             what == "ll" ||
             what == "LL") {

    result <- fit_matrix[c("loglik", "penalized_loglik"), "overall", drop = FALSE]

  } else if (what == "fit_matrix" ||
             what == "fit.matrix") {

    result <- fit_matrix

  } else if (what == "probcat") {

    result <- probCat

  } else if (what == "timing" ||
             what == "elapsed") {

    return(fit@Optim$elapsed)

  }

  return(result)

}

#' @method latInspect llcalist
#' @export
latInspect.llcalist <- function(model, what = "profile") {

  nmodels <- length(model)
  out <- vector("list", length = nmodels)
  for(i in 1:nmodels) {

    out[[i]] <- latInspect.llca(model[[i]], what = what)
    names(out)[i] <- paste("nclasses=",
                           ncol(model[[i]]@modelInfo$trans$class),
                           sep = "")

  }

  class(out) <- "latInspect.llcalist"

  return(out)

}
