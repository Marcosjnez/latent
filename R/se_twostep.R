# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 28/03/2026

se_twostep <- function(fit2, type = "standard") {

  fit1 <- fit2@modelInfo$original_model

  # Variance covariance of step 1:
  # fit1@modelInfo$param
  VCOV1 <- se(fit1, type = type)
  # Variance covariance of step 2:
  # fit2@modelInfo$param
  fit2@modelInfo$original_model <- NULL # Avoid infinite recursion
  VCOV2 <- se(fit2, type = type)

  # Get the full model structure without constraints:
  args <- fit1@data_list$args
  args$do.fit <- FALSE
  args$X <- fit2@data_list$original_X
  fit <- do.call(lca, args)
  # fit@modelInfo$param

  # Get the hessian matrix:
  control_manifold <- fit@modelInfo$control_manifold
  control_transform <- fit@modelInfo$control_transform
  control_estimator <- fit@modelInfo$control_estimator
  control_optimizer <- fit@modelInfo$control_optimizer
  parameters <- fit2@Optim$transparameters[fit@modelInfo$parameters_labels]
  transparameters <- fit2@Optim$transparameters[fit@modelInfo$transparameters_labels]
  control_optimizer$parameters[[1]] <- parameters
  control_optimizer$transparameters[[1]] <- transparameters
  x <- get_hess(control_manifold, control_transform,
                control_estimator, control_optimizer)
  colnames(x$h) <- rownames(x$h) <- fit@modelInfo$parameters_labels

  # Get the second-order devarivates between fixed and estimated parameters:
  model_pars <- fit@modelInfo$parameters_labels %in% fit2@modelInfo$parameters_labels
  nuisance_pars <- !model_pars
  df2_dparamdR <- x$h[nuisance_pars, model_pars]

  # Pick the right coefficients from the VCOV of the first step:
  ACOV <- VCOV1$vcov[fit@modelInfo$parameters_labels[nuisance_pars],
                     fit@modelInfo$parameters_labels[nuisance_pars]]
  # rownames(ACOV)
  # rownames(df2_dparamdR)
  # Ham of sandwich estimator:
  B <- t(df2_dparamdR) %*% ACOV %*% df2_dparamdR
  VCOV2$B <- B

  # Get the hessian matrix of second-step model:
  H_inv <- solve(VCOV2$H)
  VCOV2$vcov <- H_inv %*% B %*% H_inv

  # Update the standard errors of the model parameters in the VCOV2 object:
  VCOV2$se <- sqrt(diag(VCOV2$vcov))
  names(VCOV2$se) <- fit2@modelInfo$parameters_labels

  # # Create the tables of parameters with standard errors:
  # indices <- match(unlist(fit2@modelInfo$param),
  #                  fit2@modelInfo$parameters_labels)
  # values <- VCOV2$se[indices]
  # values[is.na(values)] <- 0
  # VCOV2$table_se <- fill_list_with_vector(fit2@modelInfo$param, values)
  # VCOV2$table_se <- allnumeric(VCOV2$table_se)

  # se(fit2, type = type)$table_se$beta
  # VCOV2$table_se$beta

  # effects_vcov <- effects_coding(fit2@parameters$beta, VCOV2$vcov)
  # rownames(effects_vcov$vcov_new) <- colnames(effects_vcov$vcov_new) <-
  #   names(effects_vcov$se_new) <- c(fit2@modelInfo$lca_all$beta)

  # round(effects_vcov$se_new, 4)
  # round(effects_vcov$vcov_new, 4)

  return(VCOV2)

}
