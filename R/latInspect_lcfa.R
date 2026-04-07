# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 07/04/2026
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
                            what = "est") {

  # fit must inherit from class lcfa
  stopifnot(inherits(fit, "lcfa"))

  # be case insensitive
  what <- tolower(what)

  groups <- vector("list", length = fit@data_list$ngroups)
  names(groups) <- fit@data_list$group_label

  #### Extract the fit ####

  losses_vector <- unlist(lapply(fit@Optim$outputs$estimators$doubles,
                                 FUN = \(x) x[[1]]))
  logliks_vector <- unlist(lapply(fit@Optim$outputs$estimators$doubles,
                                  FUN = \(x) x[[2]]))
  logliks_indep_vector <- unlist(lapply(fit@Optim$outputs$estimators$doubles,
                                        FUN = \(x) x[[3]]))
  logliks_sat_vector <- unlist(lapply(fit@Optim$outputs$estimators$doubles,
                                      FUN = \(x) x[[4]]))
  penalty_vector <- unlist(lapply(fit@Optim$outputs$estimators$doubles,
                                  FUN = \(x) x[[5]]))

  k <- 1L

  # Fit by group and pattern:
  # For each group...
  for(i in 1:fit@data_list$ngroups) {

    # For each pattern within each group...
    for(j in 1:fit@data_list$correl[[i]]$npatterns) {

      groups[[i]]$loss[[j]] <- losses_vector[k]
      groups[[i]]$loglik[[j]] <- logliks_vector[k]
      groups[[i]]$loglik_indep[[j]] <- logliks_indep_vector[k]
      groups[[i]]$loglik_sat[[j]] <- logliks_sat_vector[k]
      groups[[i]]$penalty[[j]] <- penalty_vector[k]
      k <- k+1L

    }

    names(groups[[i]]$loss) <- names(groups[[i]]$loglik) <-
      names(groups[[i]]$loglik_indep) <- names(groups[[i]]$loglik_sat) <-
      names(groups[[i]]$penalty) <- fit@data_list$correl[[i]]$patterns_names

  }

  # Fit by group:
  loss <- lapply(groups, FUN = \(x) sum(unlist(x$loss)))
  penalized_loss <- lapply(groups, FUN = \(x) sum(unlist(x$loss) +
                                                    unlist(x$penalty)))

  loglik <- lapply(groups, FUN = \(x) sum(unlist(x$loglik)))
  penalized_loglik <- lapply(groups, FUN = \(x) sum(unlist(x$loglik) +
                                                      unlist(x$penalty)))

  loglik_indep <- lapply(groups, FUN = \(x) sum(unlist(x$loglik_indep)))
  loglik_sat <- lapply(groups, FUN = \(x) sum(unlist(x$loglik_sat)))

  penalty <- lapply(groups, FUN = \(x) sum(unlist(x$penalty)))

  group_by_fit <- list(loss = loss,
                       penalized_loss = penalized_loss,
                       loglik = loglik,
                       penalized_loglik = penalized_loglik,
                       loglik_indep = loglik_indep,
                       loglik_sat = loglik_sat,
                       penalty = penalty)

  fit_by_group <- group_lists_by_sublists(loss = loss,
                                          penalized_loss = penalized_loss,
                                          loglik = loglik,
                                          penalized_loglik = penalized_loglik,
                                          loglik_indep = loglik_indep,
                                          loglik_sat = loglik_sat,
                                          penalty = penalty)

  #### Extract the parameters ####

  lambda <- theta <- psi <- model <- xtheta <- xpsi <- S <- resids <-
    vector("list", length = fit@data_list$ngroups)
  names(lambda) <- names(theta) <- names(psi) <- names(xtheta) <-
    names(xpsi) <- names(model) <- names(S) <- names(resids) <-
    fit@data_list$group_label

  lambda_group <- paste("lambda.", fit@data_list$group_label, sep = "")
  psi_group <- paste("psi.", fit@data_list$group_label, sep = "")
  theta_group <- paste("theta.", fit@data_list$group_label, sep = "")
  xpsi_group <- paste("xpsi.", fit@data_list$group_label, sep = "")
  xtheta_group <- paste("xtheta.", fit@data_list$group_label, sep = "")
  model_group <- paste("model.", fit@data_list$group_label, sep = "")
  S_group <- vector("list", length = fit@data_list$ngroups)
  for(i in 1:fit@data_list$ngroups) {
    for(j in 1:fit@data_list$correl[[i]]$npatterns) {
      S_group[[i]][[j]] <- paste("S.", fit@data_list$group_label[i],
                                 ".pattern", j, sep = "")
    }
  }

  # For each group...
  for(i in 1:fit@data_list$ngroups) {

    lambda[[i]] <- fit@transformed_pars[lambda_group[i]]
    theta[[i]] <- fit@transformed_pars[theta_group[i]]
    psi[[i]] <- fit@transformed_pars[psi_group[i]]
    model[[i]] <- fit@transformed_pars[model_group[i]]
    xtheta[[i]] <- fit@transformed_pars[xtheta_group[i]]
    xpsi[[i]] <- fit@transformed_pars[xpsi_group[i]]

    for(j in 1:fit@data_list$correl[[i]]$npatterns) {
      S[[i]][j] <- fit@transformed_pars[S_group[[i]][[j]]]
      resids[[i]][[j]] <- S[[i]][[j]][[1]] - model[[i]][[1]]
    }
    names(S[[i]]) <- names(resids[[i]]) <- fit@data_list$correl[[i]]$patterns_names

  }

  #### Return ####

  if(what == "est" ||
     what == "estimates" ||
     what == "parameters" ||
     what == "fixed" ||
     what == "items") {

    for(i in 1:fit@data_list$ngroups) {
      groups[[i]]$lambda <- lambda[[i]][[1]]
      groups[[i]]$psi <- psi[[i]][[1]]
      groups[[i]]$theta <- theta[[i]][[1]]
    }

    return(lapply(groups, FUN = \(x) x[c("lambda", "theta", "psi")]))

  } else if(what == "rhat" ||
            what == "model") {

    return(model)

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

    uniquenesses <- lapply(theta, FUN = \(x) lapply(x, FUN = diag))
    return(uniquenesses)

  } else if(what == "W") {

    W <- lapply(fit@modelInfo$control_estimator, FUN = \(x) x$W)
    return(W)

  } else if(what == "weights") {

    w <- lapply(fit@modelInfo$control_estimator, FUN = \(x) x$w)
    return(w)

  } else if(what == "loss" ||
            what == "f") {

    # return(group_by_fit[c("loss", "penalized_loss")])

    return(lapply(fit_by_group, FUN = \(x) x[c("loss", "penalized_loss")]))

  } else if(what == "loglik") {

    # return(group_by_fit[c("loglik", "penalized_loglik",
    #                       "loglik_indep", "loglik_sat")])

    return(lapply(fit_by_group, FUN = \(x) x[c("loglik", "penalized_loglik",
                                               "loglik_indep", "loglik_sat")]))

  } else if(what == "fit") {

    # return(group_by_fit)

    return(fit_by_group)

  } else {

    stop("Unknown request")

  }

}
