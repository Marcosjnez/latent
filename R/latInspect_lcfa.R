# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 02/11/2025
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

  # Number of groups:
  ngroups <- fit@Optim$data_list$ngroups
  nitems <- fit@Optim$data_list$nitems
  nfactors <- fit@Optim$data_list$nfactors
  group_label <- fit@Optim$data_list$group_label
  item_label <- fit@Optim$data_list$item_label
  factor_label <- fit@Optim$data_list$factor_label

  # Get the indices of the transformation structures "factor_cor". These are the
  # structures that contain the model parameters:
  all_transforms <- unlist(lapply(fit@modelInfo$control_transform, FUN = \(x) x$transform))
  indices_factor_cor <- which(all_transforms == "factor_cor")

  # Get the indices of the estimator structures "cfa_dwls" and "cfa_ml":
  all_estimators <- unlist(lapply(fit@modelInfo$control_estimator, FUN = \(x) x$estimator))
  indices_cfa <- which(all_estimators == "cfa_dwls" | all_estimators == "cfa_ml")

  # Get the indices of the estimator structures "logdetmat" (penalties):
  indices_logdetmat <- which(all_estimators == "logdetmat")

  # Initialize the objects to be returned:
  jacob <- lambda <- psi <- theta <- uniquenesses <- model <- resids <- W <- w <-
    loss <- penalized_loss <- loglik <- penalized_loglik <- penalty <-
    vector("list", length = ngroups)
  names(jacob) <- names(lambda) <- names(psi) <- names(theta) <- names(uniquenesses) <-
    names(model) <- names(resids) <- names(W) <- names(w) <- names(loss) <-
    names(penalized_loss) <- names(loglik) <- names(penalized_loglik) <-
    names(penalty) <- group_label
  # jacob is the jacobian of the model matrix wrt the parameters

  x <- fit@Optim$opt

  for(i in 1:ngroups) {

    p <- nitems[[i]]
    q <- nfactors[[i]]
    j <- indices_factor_cor[i]
    k <- indices_cfa[i]

    loss[[i]] <- c(x$outputs$estimators$doubles[[k]][[1]])
    loglik[[i]] <- c(x$outputs$estimators$doubles[[k]][[2]])
    w[[i]] <- c(x$outputs$estimators$doubles[[k]][[3]])

    # If there are penalties, add the penalties to the loss or loglik:
    if(length(indices_logdetmat) > 0) {

      l <- indices_logdetmat[i]
      penalty[[i]] <- c(x$outputs$estimators$doubles[[l]][[1]])
      penalized_loss[[i]] <- loss[[i]] + penalty[[i]]
      penalized_loglik[[i]] <- loglik[[i]] + penalty[[i]]

    } else {

      penalized_loss[[i]] <- loss[[i]]
      penalized_loglik[[i]] <- loglik[[i]]

    }

    # Extract model parameters:
    jacob[[i]] <- matrix(x$outputs$transformations$matrices[[j]][[1]], p, q)
    lambda[[i]] <- matrix(x$outputs$transformations$matrices[[j]][[2]], p, q)
    psi[[i]] <- matrix(x$outputs$transformations$matrices[[j]][[3]], q, q)
    theta[[i]] <- matrix(x$outputs$transformations$matrices[[j]][[4]], p, p)
    model[[i]] <- matrix(x$outputs$transformations$matrices[[j]][[5]], p, p)
    uniquenesses[[i]] <- c(x$outputs$transformations$vectors[[j]][[1]])
    resids[[i]] <- matrix(x$outputs$estimators$matrices[[k]][[1]], p, p)
    W[[i]] <- matrix(x$outputs$estimators$matrices[[k]][[2]], p, p)

    # Name model parameters:
    colnames(lambda[[i]]) <- colnames(psi[[i]]) <- rownames(psi[[i]]) <-
      factor_label[[i]]
    rownames(lambda[[i]]) <- rownames(theta[[i]]) <- colnames(theta[[i]]) <-
      rownames(model[[i]]) <- colnames(model[[i]]) <- rownames(resids[[i]]) <-
      colnames(resids[[i]]) <- rownames(W[[i]]) <- colnames(W[[i]]) <-
      names(uniquenesses[[i]]) <- item_label[[i]]

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

    return(uniquenesses)

  } else if(what == "W") {

    return(W)

  } else if(what == "weights") {

    return(w)

  } else if(what == "jacob" ||
            what == "jacobs" ||
            what == "jacobian" ||
            what == "jacobians") {

    return(jacob)

  } else if(what == "loss" ||
            what == "f") {

    return(list(loss = loss, penalized_loss = penalized_loss))

  } else if(what == "loglik") {

    return(list(loglik = loglik, penalized_loglik = penalized_loglik))

  } else if(what == "fit") {

    return(list(loss = loss, penalized_loss = penalized_loss,
                loglik = loglik, penalized_loglik = penalized_loglik))

  } else {

    stop("Unknown request")

  }

}
