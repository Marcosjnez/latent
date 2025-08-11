# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 22/05/2025
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
#' @export
se <- function(fit, confidence = 0.95, digits = 2) {

  data <- fit@Optim$data
  nclasses <- length(fit@parameters$classes)
  control_manifold <- fit@Optim$manifold_setup
  control_transform <- fit@Optim$transform_setup
  control_estimator <- fit@Optim$estimator_setup
  control_optimizer <- fit@Optim$control_setup

  item <- fit@modelInfo$item
  gauss <- "gaussian" %in% item
  multin <- "multinomial" %in% item

  ntransparam <- length(fit@Optim$transparameters)
  indices <- 1:ntransparam
  X <- vector(length = ntransparam)
  control_manifold <- control_transform <- list()
  control_manifold[[1]] <- list(manifold = "euclidean",
                                indices = indices-1L)
  control_transform[[1]] <- list(transform = "identity",
                                 indices = indices-1L,
                                 target_indices = indices-1L,
                                 vector_indices = indices-1L,
                                 X = X)
  x <- fit@Optim$transparameters
  control_optimizer$parameters <- list()
  control_optimizer$parameters[[1]] <- x

  grad_computation <- function(parameters, control_manifold,
                               control_transform, control_estimator,
                               control_optimizer) {

    grad <- latent::grad_comp(parameters = parameters,
                              control_manifold = control_manifold,
                              control_transform = control_transform,
                              control_estimator = control_estimator,
                              control_optimizer = control_optimizer)$grad

    return(grad)

  }

  # latent:::grad_comp(x, control_manifold,
  #                    control_transform, control_estimator,
  #                    control_optimizer)
  H <- numDeriv::jacobian(func = grad_computation,
                          x = x,
                          control_manifold = control_manifold,
                          control_transform = control_transform,
                          control_estimator = control_estimator,
                          control_optimizer = control_optimizer)

  param_charvector <- unique(unlist(fit@modelInfo$prob_model))
  non_alnum_indices <- grep("^(?!-?\\d+(\\.\\d+)?$)",
                            param_charvector, perl = TRUE)
  rownames(H) <- colnames(H) <- param_charvector[non_alnum_indices]

  nconstraints <- length(fit@Optim$outputs$transformations$vectors)

  # Total number of rows is the number of parameters
  nparam <- length(x)

  # Initialize matrix with zeros
  constraints <- matrix(0, nrow = nparam, ncol = nconstraints)

  for(i in 1:nconstraints) {

    constraints0 <- fit@Optim$outputs$transformations$vectors
    if(!is.null(constraints0[[i]][1][[1]])) {
      v <- constraints0[[i]][[1]]
      ind <- constraints0[[i]][[2]]+1
      constraints[ind, i] <- v
    }

  }

  constraints <- constraints[non_alnum_indices, ]

  ind_target <- fit@Optim$args$indices_full_transparam_vector
  ind_select <- fit@Optim$args$indices_transparam_vector2
  se_type <- fit@modelInfo$prob_model
  se_type$classes[] <- "prob"
  indices <- which(item == "multinomial")
  for(i in indices) {
    se_type$conditionals[[i]][] <- "prob"
  }
  indices <- which(item == "gaussian")
  for(i in indices) {
    se_type$conditionals[[i]][] <- "linear"
  }

  type <- unlist(se_type)[ind_select[!is.na(ind_select)]]

  # REMOVE PROBABILITIES CLOSE TO ZERO in H and constraints:
  param_vector <- x
  param_vector <- param_vector[non_alnum_indices]
  names(param_vector) <- param_charvector[non_alnum_indices]
  type0 <- type[non_alnum_indices]
  condition <- (param_vector < 0.005 | param_vector > 0.995) & type0 == "prob"
  keep <- which(!condition)
  remove <- which(condition)
  if(length(remove) > 0) {
    H <- H[-remove, ][, -remove]
    constraints <- constraints[-remove, ]
  }
  zeros <- colSums(constraints) < 1 | duplicated(t(constraints))
  if(any(zeros)) {
    constraints <- constraints[, -which(zeros)]
  }

  # Visser, I., Raijmakers, M.E.J. and Molenaar, P.C.M. (2000),
  # Confidence intervals for hidden Markov model parameters.
  # British Journal of Mathematical and Statistical Psychology, 53: 317-327.
  # https://doi.org/10.1348/000711000159240

  K <- t(constraints)
  D <- H + t(K) %*% K
  D_inv <- solve(D)
  C0 <- D_inv - D_inv %*% t(K) %*% solve(K %*% D_inv %*% t(K)) %*% K %*% D_inv
  se0 <- sqrt(diag(C0))

  se <- vector(length = nparam)
  se[keep] <- se0
  se[remove] <- NA
  names(se) <- fit@Optim$args$transparameter_vector

  Se <- vector(length = length(unlist(fit@parameters)))
  Se[ind_target] <- se[ind_select[!is.na(ind_select)]]
  SE <- fill_list_with_vector(fit@transformed_pars, Se)

  C <- matrix(NA, nrow = nparam, ncol = nparam)
  C[keep, keep] <- C0
  rownames(C) <- colnames(C) <- fit@Optim$args$transparameter_vector

  conf <- function(x, se, confidence, type) {

    z <- stats::qnorm((1+confidence)/2)

    if(type == "prob") {

      logit_prob <- stats::qlogis(x)
      logit_se <- se / (x * (1 - x))  # Delta method

      logit_lower <- logit_prob - z * logit_se
      logit_upper <- logit_prob + z * logit_se

      # Back-transform to probability scale
      lower <- stats::plogis(logit_lower)
      upper <- stats::plogis(logit_upper)

    } else if(type == "linear") {

      lower <- x - z * se
      upper <- x + z * se

    } else if(type == "sd") {

      lower <- exp(log(x) - z * se)
      upper <- exp(log(x) + z * se)

    }


    return(c(lower, upper))

  }

  ci <- sapply(1:nparam, FUN = \(i) conf(x[i], se[i], confidence = confidence,
                                         type = type[i]))
  rownames(ci) <- c("lower", "upper")
  colnames(ci) <- fit@Optim$args$transparameter_vector

  ci0 <- format(round(rbind(x, ci), digits), nsmall = 2)
  Ci <- apply(ci0, MARGIN = 2, FUN = \(x) paste(x[1], " (", x[2], " - ", x[3], ")",
                                                sep = "", collapse = ""))
  CI <- vector(length = length(unlist(fit@parameters)))
  CI[ind_target] <- Ci[ind_select[!is.na(ind_select)]]
  CI <- fill_list_with_vector(fit@transformed_pars, CI)

  result <- list()
  result$se <- se
  result$SE <- SE
  result$vcov <- C
  result$ci <- ci
  result$CI <- CI

  return(result)

}

