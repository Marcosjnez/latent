# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 17/08/2026
#'
#' Variance-Covariance Matrix for Confirmatory Factor Models
#'
#' Compute the covariance matrix of CFA parameters by propagating the sampling
#' covariance matrix of the means, covariances, correlations, and thresholds
#' used as sample statistics.
#'
#' @param fit A fitted object of class \code{"lcfa"}.
#'
#' @return A list containing the Hessian, cross-derivative matrix, sandwich
#'   middle matrix, parameter variance-covariance matrix, and standard errors.
#'
#' @keywords internal
robust.lcfa <- function(fit) {

  result <- vcov_lcfa(fit)

  #### Result ####

  return(result)

}

vcov_lcfa <- function(fit) {

  dataList <- fit@dataList
  modelInfo <- fit@modelInfo
  Optim <- fit@Optim
  data_param <- dataList$data_param

  #### Hessian of the fitted CFA model ####

  H <- hessian(fit)
  H <- (H+t(H))/2

  eigenvalues <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  tolerance <- sqrt(.Machine$double.eps)*max(1, max(abs(eigenvalues)))

  if(any(eigenvalues <= tolerance)) {
    warning("The CFA Hessian is not positive definite; an approximate inverse was used.")
  }

  Hinv <- approx_Hinv(H)

  if(any(!is.finite(Hinv))) {
    stop("The CFA Hessian could not be inverted.")
  }

  Hinv <- (Hinv+t(Hinv))/2
  rownames(Hinv) <- colnames(Hinv) <- modelInfo$parameters_labels

  #### Variance-covariance matrix of the sample statistics ####

  VCOV_cov <- data_param$VCOV_cov
  VCOV_means <- data_param$VCOV_means

  if(modelInfo$control_optimizer$meanstructure) {
    sample_VCOV <- block_diag(c(VCOV_means, VCOV_cov))
    # Means and covariance statistics are currently treated as independent
    # blocks. Their within-block covariances are retained.
  } else {
    sample_VCOV <- block_diag(VCOV_cov)
  }

  #### Unrestricted sample-statistic model ####

  control_free <- modelInfo$control_optimizer
  control_free$free_taus <- TRUE
  control_free$free_S <- TRUE
  control_free$free_M <- TRUE
  control_free$rstarts <- 1L
  control_free$cores <- 1L
  control_free$start <- NULL

  full_model_free <- create_lcfa_model(dataList = dataList,
                                       model = dataList$model,
                                       control = control_free)

  modelInfo_free <- create_lcfa_modelInfo(dataList = dataList,
                                          full_model = full_model_free,
                                          control = control_free)

  parameters_free <-
    Optim$transparameters[modelInfo_free$parameters_labels]
  transparameters_free <-
    Optim$transparameters[modelInfo_free$transparameters_labels]

  if(anyNA(parameters_free) || anyNA(transparameters_free)) {
    stop("Could not match all unrestricted-model parameters to the fitted model.")
  }

  control_optimizer_free <- modelInfo_free$control_optimizer
  control_optimizer_free$parameters[[1L]] <- parameters_free
  control_optimizer_free$transparameters[[1L]] <- transparameters_free

  H_full <- get_hess(control_manifold = modelInfo_free$control_manifold,
                     control_transform = modelInfo_free$control_transform,
                     control_estimator = modelInfo_free$control_estimator,
                     control_optimizer = control_optimizer_free,
                     cores = 1L)$h

  rownames(H_full) <- colnames(H_full) <-
    modelInfo_free$parameters_labels

  #### Sample-statistic/model cross derivatives ####

  model_pars <- modelInfo_free$parameters_labels %in%
    modelInfo$parameters_labels
  nuisance_pars <- !model_pars

  if(!any(model_pars)) {
    stop("No CFA model parameters were identified in the unrestricted model.")
  }

  if(!any(nuisance_pars)) {
    stop("No sample-statistic nuisance parameters were identified.")
  }

  C <- H_full[nuisance_pars, model_pars, drop = FALSE]

  missing_model_pars <- setdiff(modelInfo$parameters_labels, colnames(C))

  if(length(missing_model_pars) > 0L) {
    stop("Could not identify all fitted CFA parameters in the unrestricted model: ",
         paste(missing_model_pars, collapse = ", "))
  }

  C <- C[, modelInfo$parameters_labels, drop = FALSE]

  nuisance_labels <- rownames(C)
  missing_VCOV <- setdiff(nuisance_labels, rownames(sample_VCOV))

  if(length(missing_VCOV) > 0L) {
    stop("The sample-statistic variance-covariance matrix is missing parameter(s): ",
         paste(missing_VCOV, collapse = ", "))
  }

  sample_VCOV <- sample_VCOV[nuisance_labels,
                              nuisance_labels, drop = FALSE]

  B <- t(C) %*% sample_VCOV %*% C
  B <- (B+t(B))/2

  #### Sandwich covariance matrix ####

  VCOV <- Hinv %*% B %*% Hinv
  VCOV <- (VCOV+t(VCOV))/2
  rownames(VCOV) <- colnames(VCOV) <- modelInfo$parameters_labels

  se <- sqrt(diag(VCOV))
  names(se) <- modelInfo$parameters_labels

  #### Result ####

  result <- list(H = H,
                 C = C,
                 B = B,
                 sample_VCOV = sample_VCOV,
                 VCOV = VCOV,
                 se = se)

  return(result)

}
