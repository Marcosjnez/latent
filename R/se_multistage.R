# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 15/08/2026
#'
#' Standard Errors for Latent Objects with multistage estimation.
#'
#' @export
se.multistage <- function(fit, type = "information", parameters = NULL,
                          digits = 4L, ...) {

  #### Parameters ####

  if(is.null(parameters)) {
    parameters <- fit@modelInfo$trans[names(fit@modelInfo$param)]
    # parameters <- fit@modelInfo$parameters_labels
  } else if(!any(unlist(parameters) %in%
                 fit@modelInfo$transparameters_labels)) {

    stop("Unknown parameters.")

  }

  #### Identify the nested latent models ####

  # Identifying more than one model is possible

  model <- fit@extra
  latent_idx <- integer(0L)

  latent_idx <- which(vapply(model,
                             FUN = \(x) inherits(x, "latent"),
                             FUN.VALUE = logical(1L)))

  #### Structural after measurement model ####

  nmodels <- length(latent_idx)
  SE <- vector("list", length = nmodels)
  for(i in seq_len(nmodels)) {

    # se.multistage is applied recursively if a nested model is also multistage:
    SE[[i]] <- se(model[[latent_idx[i]]], type = type, parameters = NULL, ...)

  }

  #### Free all the constraints in the fit model ####

  # Full unrestricted model:
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

  # Parameters estimated in the first model:
  second_pars <- fit_full@modelInfo$parameters_labels %in%
    fit@modelInfo$parameters_labels

  # Measurement parameters fixed during structural estimation:
  first_pars <- !second_pars

  if(!any(second_pars) || !any(first_pars)) {
    stop("The model does not contain both structural and measurement parameters.")
  }

  #### Second derivatives between first-model and second-model parameters ####

  C <- H_full[first_pars, second_pars, drop = FALSE]

  first_pars_names <- fit_full@modelInfo$parameters_labels[first_pars]
  second_pars_names <- fit_full@modelInfo$parameters_labels[second_pars]

  A <- as.matrix(VCOV0$vcov[first_pars_names,
                            first_pars_names, drop = FALSE])

  # Additional structural uncertainty propagated from the measurement model:
  B <- t(C) %*% A %*% C

  # Remove small numerical asymmetries:
  B <- (B + t(B))/2

  H2_inv <- solve(H2)

  newH11 <- solve(A) +
    C %*% H2_inv %*% P2 %*% H2_inv %*% t(C)

  newH12 <- C %*% H2_inv %*% P2

  newH <- rbind(
    cbind(newH11, newH12),
    cbind(t(newH12), P2)
  )

  newH_names <- c(first_pars_names, second_pars_names)
  rownames(newH) <- colnames(newH) <- newH_names

  newH <- newH[fit_full@modelInfo$parameters_labels,
               fit_full@modelInfo$parameters_labels,
               drop = FALSE]

  #### Variance-covariance matrix ####

  result <- vcov(fit_full, parameters = parameters, vcov = solve(newH))
  result$B <- B
  colnames(result$B) <- rownames(result$B) <- second_pars_names

  #### Result ####

  return(result)

}
