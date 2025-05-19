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
se <- function(fit) {

  data <- fit$opt$data
  nclasses <- length(fit$parameters$classes)
  control_manifold <- fit$opt$manifold_setup
  control_transform <- fit$opt$transform_setup
  control_estimator <- fit$opt$estimator_setup
  control_optimizer <- fit$opt$control

  item_model <- fit$modelInfo$item_model
  gauss <- "gaussian" %in% item_model
  multin <- "multinomial" %in% item_model

  ntransparam <- length(fit$opt$transparameters)
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
  x <- fit$opt$transparameters
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
  # round(cbind(x, H-t(H)), 3)

  param_charvector <- unlist(fit$modelInfo$prob_model)
  non_alnum_indices <- grep("^(?!-?\\d+(\\.\\d+)?$)",
                            param_charvector, perl = TRUE)
  rownames(H) <- colnames(H) <- param_charvector[non_alnum_indices]

  nconstraints <- length(fit$opt$outputs$transformations$vectors)

  # Total number of rows is the number of parameters
  nparam <- length(x)

  # Initialize matrix with zeros
  constraints <- matrix(0, nrow = nparam, ncol = nconstraints)

  for(i in 1:nconstraints) {

    v <- fit$opt$outputs$transformations$vectors[[i]][[1]]
    ind <- fit$opt$outputs$transformations$vectors[[i]][[2]]+1
    constraints[ind, i] <- v

  }

  constraints <- constraints[non_alnum_indices, ]

  # REMOVE PROBABILITIES CLOSE TO ZERO in H and constraints:
  param_vector <- x
  param_vector <- param_vector[non_alnum_indices]
  names(param_vector) <- param_charvector[non_alnum_indices]
  condition <- param_vector < 0.005 | param_vector > 0.995
  keep <- which(!condition)
  remove <- which(condition)
  H <- H[-remove, ][, -remove]
  constraints <- constraints[-remove, ]
  zeros <- colSums(constraints) < 1
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
  names(se) <- param_charvector
  SE <- fill_list_with_vector(fit$transformed_parameters, se)

  C <- matrix(NA, nrow = nparam, ncol = nparam)
  C[keep, keep] <- C0
  rownames(C) <- colnames(C) <- param_charvector

  confidence <- function(x, se) {

    z <- qnorm(0.975)
    logit_prob <- qlogis(x)
    logit_se <- se / (x * (1 - x))  # Delta method

    logit_lower <- logit_prob - z * logit_se
    logit_upper <- logit_prob + z * logit_se

    # Back-transform to probability scale
    lower_prob <- plogis(logit_lower)
    upper_prob <- plogis(logit_upper)

    return(c(lower_prob, upper_prob))

  }

  ci <- sapply(1:nparam, FUN = \(i) confidence(x[i], se[i]))
  rownames(ci) <- c("lower", "upper")
  colnames(ci) <- param_charvector
  ci0 <- format(round(rbind(x, ci), 2), nsmall = 2)
  CI <- apply(ci0, MARGIN = 2, FUN = \(x) paste(x[1], " (", x[2], " - ", x[3], ")",
                                                sep = "", collapse = ""))
  CI <- fill_list_with_vector(fit$transformed_parameters, CI)

  result <- list()
  result$se <- se
  result$SE <- SE
  result$vcov <- C
  result$ci <- ci
  result$CI <- CI

  return(result)

}

