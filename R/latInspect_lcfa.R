# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 05/05/2026
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
latInspect.lcfa <- function(fit, what = "est") {

  # fit must inherit from class lcfa
  stopifnot(inherits(fit, "lcfa"))

  # be case insensitive
  what <- tolower(what)

  groups <- vector("list", length = fit@dataList$ngroups)
  names(groups) <- fit@dataList$group_label

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

  colnames(fit_matrix) <- unlist(fit@dataList$data_param$S_group)

  # group_by_fit <- matrix(rowSums(fit_matrix), ncol = 1L,
  #                        dimnames = list(rownames(fit_matrix), "value"))

  if(fit@dataList$ngroups > 1) {

    Ind <- sapply(fit@dataList$group_label, FUN = \(x) {
      idx <- grepl(paste0(x, "($|\\.)"), colnames(fit_matrix)) + 0.00
    })
    Ind <- cbind(Ind, 1)
    fit_by_group <- fit_matrix %*% Ind
    colnames(fit_by_group) <- c(as.character(fit@dataList$group_label),
                                "overall")

  } else {

    fit_by_group <- matrix(rowSums(fit_matrix), ncol = 1L,
                           dimnames = list(rownames(fit_matrix), "overall"))

  }

  #### Compute residuals ####

  matrices <- fit@Optim$outputs$estimators$matrices
  idx_S <- grepl(paste(c("S"), collapse = "|"),
                 names(fit@transformed_pars))
  idx_means <- grepl(paste(c("means"), collapse = "|"),
                     names(fit@transformed_pars))
  resids_S <- vector("list", length = sum(idx_S))
  for(i in 1:sum(idx_S)) {
    p <- sqrt(length(matrices[[1]][[i]]))
    resids_S[[i]] <- matrix(matrices[[1]][[i]], nrow = p, ncol = p)
    rownames(resids_S[[i]]) <- colnames(resids_S[[i]]) <-
      rownames(fit@transformed_pars[idx_S][[i]])
    names(resids_S)[i] <- names(fit@transformed_pars[idx_S])[i]
  }

  if(fit@modelInfo$control_optimizer$meanstructure) {

    resids_means <- vector("list", length = sum(idx_means))
    for(i in 1:sum(idx_S)) {
      p <- length(matrices[[2]][[i]])
      resids_means[[i]] <- matrices[[2]][[i]]
    }
    names(resids_means)[i] <- names(fit@transformed_pars[idx_means])[i]

    resids <- c(resids_S, resids_means)

  } else {
    resids <- resids_S
  }

  #### Return ####

  if(what == "est" ||
     what == "estimates" ||
     what == "parameters" ||
     what == "fixed" ||
     what == "items") {

    idx <- grepl(paste(c("lambda", "theta", "psi"), collapse = "|"),
                 names(fit@transformed_pars))

    return(fit@transformed_pars[idx])

  } else if(what == "rhat" ||
            what == "model") {

    idx <- grepl(paste(c("model"), collapse = "|"),
                 names(fit@transformed_pars))

    return(fit@transformed_pars[idx])

    return(model)

  } else if(what == "resid" ||
            what == "residuals") {

    return(resids)

  } else if(what == "lambda" ||
            what == "loadings") {

    idx <- grepl(paste(c("lambda"), collapse = "|"),
                 names(fit@transformed_pars))

    return(fit@transformed_pars[idx])

  } else if(what == "psi") {

    idx <- grepl(paste(c("psi"), collapse = "|"),
                 names(fit@transformed_pars))

    return(fit@transformed_pars[idx])

  } else if(what == "theta") {

    idx <- grepl(paste(c("theta"), collapse = "|"),
                 names(fit@transformed_pars))

    return(fit@transformed_pars[idx])

  } else if(what == "uniquenesses") {

    idx <- grepl(paste(c("theta"), collapse = "|"),
                 names(fit@transformed_pars))

    return(lapply(fit@transformed_pars[idx], FUN = diag))

  } else if(what == "W") {

    W <- lapply(fit@modelInfo$control_estimator, FUN = \(x) x$W)
    return(W)

  } else if(what == "weights") {

    w <- lapply(fit@modelInfo$control_estimator, FUN = \(x) x$w)
    return(w)

  } else if(what == "loss" ||
            what == "f") {

    return(fit_by_group[c("loss", "penalized_loss",
                          "loss_base", "loss_sat"), , drop = FALSE])

  } else if(what == "loglik") {

    return(fit_by_group[c("loglik", "penalized_loglik",
                          "loglik_base", "loglik_sat"), , drop = FALSE])

  } else if(what == "fit") {

    # return(group_by_fit)

    return(fit_by_group)

  } else if(what == "fit.matrix" ||
            what == "fit_matrix") {

    return(fit_matrix)

  } else {

    stop("Unknown request")

  }

}
