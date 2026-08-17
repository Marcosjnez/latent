# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 17/08/2026
#'
#' Robust Variance-Covariance Matrix
#'
#' @keywords internal
robust.lcfa <- function(fit) {

  dataList <- fit@dataList
  modelInfo <- fit@modelInfo
  Optim <- fit@Optim

  #### Compute the Hessian ####

  VCOV_fit <- information.latent(fit)
  H <- VCOV_fit$H
  Hinv <- VCOV_fit$VCOV

  #### Asymptotic covariance of the sample statistics ####

  data_param <- dataList$data_param

  scale_acov_lcfa <- function(acov, nobs, object_name) {

    if(!is.list(acov)) {
      acov <- list(acov)
    }

    nobs <- unlist(nobs, use.names = FALSE)

    if(length(nobs) == 1L && length(acov) > 1L) {
      nobs <- rep(nobs, length(acov))
    }

    if(length(nobs) != length(acov)) {
      stop("The number of ", object_name,
           " ACOV matrices does not match the number of sample sizes.")
    }

    result <- lapply(seq_along(acov), FUN = \(j) {
      acov[[j]]/nobs[j]
    })

    #### Result ####

    return(result)

  }

  ACOV_covij <- vector("list", length = dataList$ngroups)
  ACOV_meansij <- vector("list", length = dataList$ngroups)

  for(i in seq_len(dataList$ngroups)) {

    ACOV_covij[[i]] <- scale_acov_lcfa(
      acov = data_param$acov_cov[[i]],
      nobs = data_param$nobs_ij[[i]],
      object_name = "covariance"
    )

    ACOV_meansij[[i]] <- scale_acov_lcfa(
      acov = data_param$acov_means[[i]],
      nobs = data_param$nobs_ij[[i]],
      object_name = "mean"
    )

  }

  if(modelInfo$control_optimizer$meanstructure) {
    ACOV <- block_diag(c(ACOV_meansij, ACOV_covij))
    # This assumes that means and covariances are estimated independently
  } else {
    ACOV <- block_diag(ACOV_covij)
  }

  #### Unrestricted sample-statistic model ####

  # Recreate only the model and optimization structures. The sample means,
  # covariances, thresholds, and their ACOV matrices already stored in dataList
  # are reused, so this does not recompute the correlations or means.
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

  common <- intersect(rownames(ACOV), rownames(C))

  if(length(common) == 0L) {
    stop("The sample-statistic covariance matrix and cross-Hessian have no matching parameters.")
  }

  C <- C[common, , drop = FALSE]
  ACOV <- ACOV[common, common, drop = FALSE]

  B <- t(C) %*% ACOV %*% C
  B <- (B+t(B))/2

  #### Sandwich covariance matrix ####

  VCOV <- Hinv %*% B %*% Hinv
  VCOV <- (VCOV+t(VCOV))/2

  #### Result ####

  result <- list(H = H, C = C, VCOV = VCOV, se = sqrt(diag(VCOV)))

  return(result)

}
