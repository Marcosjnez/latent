#' @export
vcov.latent <- function(fit, parameters = NULL, H = NULL) {

  #### Parameters ####

  if(is.null(parameters)) {

    parameters <- fit@modelInfo$trans[names(fit@modelInfo$param)]

  } else if(!any(unlist(parameters) %in%
                 fit@modelInfo$transparameters_labels)) {

    stop("Unknown parameters.")

  }

  #### Identify nested latent model ####

  model <- fit@modelInfo$control_optimizer$model
  latent_idx <- integer(0L)

  if(inherits(model, "latent")) {

    model <- list(model)
    latent_idx <- 1L

  } else if(is.list(model)) {

    latent_idx <- which(vapply(model,
                               FUN = \(x) inherits(x, "latent"),
                               FUN.VALUE = logical(1L)))

  }

  #### Ordinary model ####

  # If no nested latent model is found in fit, get directly the vcov:
  if(length(latent_idx) == 0L) {

    fit@modelInfo$control_optimizer$parameters[[1]] <-
      fit@Optim$parameters

    fit@modelInfo$control_optimizer$transparameters[[1]] <-
      fit@Optim$transparameters

    if(is.null(H)) H <- hessian(fit)

    fit@modelInfo$control_optimizer$idx_transforms <-
      trans_depends(fit, parameters)

    VCOV <- get_vcov(fit@modelInfo$control_manifold,
                     fit@modelInfo$control_transform,
                     fit@modelInfo$control_estimator,
                     fit@modelInfo$control_optimizer,
                     H = H)

    rownames(VCOV$vcov) <- colnames(VCOV$vcov) <- names(VCOV$se) <-
      fit@modelInfo$transparameters_labels

    selected_parameters <- unique(unlist(parameters))
    selected_idx <- match(selected_parameters,
                          fit@modelInfo$transparameters_labels)

    VCOV$vcov <- VCOV$vcov[selected_idx, selected_idx, drop = FALSE]
    VCOV$se <- VCOV$se[selected_idx]
    VCOV$jacob <- VCOV$jacob[selected_idx, selected_idx, drop = FALSE]
    VCOV$H <- H
    VCOV$B <- matrix(numeric(0L), nrow = 0L, ncol = 0L)

    #### Result ####

    return(VCOV)

  }

  #### Structural after measurement model ####

  if(length(latent_idx) > 1L) {

    stop("Only one nested latent model at each model level is currently supported.")

  }

  latent_idx <- latent_idx[1L]

  # fit0 is the measurement model:
  fit0 <- model[[latent_idx]]

  # Retain every other model specification:
  remaining_model <- model[-latent_idx]

  if(length(remaining_model) == 0L) {
    remaining_model <- NULL
  }

  #### Measurement-model uncertainty ####

  # Recursive call:
  #
  # If fit0 itself contains another latent object, vcov.latent()
  # automatically applies the same correction to fit0 first.
  VCOV0 <- vcov(fit0, parameters = fit0@modelInfo$parameters_labels, H = NULL)

  #### Structural-model covariance ####

  # Hessian of the structural model:
  H2 <- hessian(fit)

  #### Full unrestricted model ####

  args <- fit@dataList$args
  # Remove the parameter constraints from previous fitted objects:
  args$model <- remaining_model  # Paste only the model dependencies
  args$do.fit <- FALSE           # Don't fit, just get the model
  args$adjustment <- "none"      # One-step fit
  fit_full <- do.call(lca, args) # Unrestricted model

  # Paste the parameter estimates of the structural model into the
  # corresponding positions of the unrestricted model:
  fit_full@Optim$parameters <-
    fit@Optim$transparameters[fit_full@modelInfo$parameters_labels]

  fit_full@Optim$transparameters <-
    fit@Optim$transparameters[fit_full@modelInfo$transparameters_labels]

  #### Full-model Hessian ####

  # Hessian matrix of the unrestricted model:
  H_full <- hessian(fit_full)

  # Parameters estimated in the structural model:
  structural_pars <- fit_full@modelInfo$parameters_labels %in%
    fit@modelInfo$parameters_labels

  # Measurement parameters fixed during structural estimation:
  measurement_pars <- !structural_pars

  if(!any(structural_pars) || !any(measurement_pars)) {
    stop("The model does not contain both structural and measurement parameters.")
  }

  #### Measurement-structural second-order derivatives ####

  C <- H_full[measurement_pars, structural_pars, drop = FALSE]

  measurement_pars_names <-
    fit_full@modelInfo$parameters_labels[measurement_pars]

  structural_pars_names <-
    fit_full@modelInfo$parameters_labels[structural_pars]

  A <- as.matrix(VCOV0$vcov[measurement_pars_names,
                            measurement_pars_names, drop = FALSE])

  # Additional structural uncertainty propagated from the measurement model:
  B <- t(C) %*% A %*% C

  # Remove small numerical asymmetries:
  B <- (B + t(B))/2

  #### Corrected Hessian ####

  # The joint asymptotic covariance of measurement and structural parameters is:
  #
  #   [ A                     -A C H2^-1                 ]
  #   [ -H2^-1 C' A           H2^-1 + H2^-1 C' A C H2^-1 ]
  #
  # Rather than constructing and inverting this covariance matrix, construct
  # its inverse directly:
  #
  #   [ A^-1 + C H2^-1 C'     C  ]
  #   [ C'                    H2 ]
  #
  # This retains the non-zero covariance between measurement and structural
  # estimates and allows vcov() to use its ordinary Hessian inversion machinery.

  newH11 <- solve(A) + C %*% solve(H2, t(C))
  newH <- rbind(cbind(newH11, C), cbind(t(C), H2))
  newH_names <- c(measurement_pars_names, structural_pars_names)
  rownames(newH) <- colnames(newH) <- newH_names

  newH <- newH[fit_full@modelInfo$parameters_labels,
               fit_full@modelInfo$parameters_labels,
               drop = FALSE]

  #### Variance-covariance matrix ####

  result <- vcov(fit_full, parameters = parameters, H = newH)
  result$B <- B
  colnames(result$B) <- rownames(result$B) <- structural_pars_names

  result$newH <- newH

  #### Result ####

  return(result)

}
