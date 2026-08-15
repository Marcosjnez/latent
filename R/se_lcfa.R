# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 14/08/2026
#'
#' Standard Errors for Confirmatory Factor Analysis
#'
#' Extract or compute standard errors for a fitted \code{"lcfa"} object.
#'
#' @param fit A fitted object of class \code{"lcfa"}.
#' @param type Character string identifying the requested standard-error type.
#' @param digits Optional number of decimal places used for formatted output.
#' @param ... Additional arguments.
#'
#' @return A list containing the standard errors, variance-covariance matrix,
#'   Hessian, sandwich middle matrix, and the parameter-shaped standard-error
#'   table.
#'
#' @method se lcfa
#' @export
se.lcfa <- function(fit, type = "standard", digits = 5L, ...) {

  #### Check inputs ####

  if(!inherits(fit, "lcfa")) {
    stop("fit must inherit from class 'lcfa'.")
  }

  if(length(fit@Optim$parameters) == 0L) {
    stop("The lcfa object has not been fitted.")
  }

  type <- match.arg(tolower(type), c("standard", "robust"))

  if(!is.null(digits) &&
     (!is.numeric(digits) || length(digits) != 1L ||
      !is.finite(digits) || digits < 0L ||
      digits != as.integer(digits))) {
    stop("digits must be NULL or a non-negative integer.")
  }

  #### Standard errors ####

  # Standard errors are computed before the final lcfa object is created and
  # stored in Optim$SE. If the model was fitted with se = FALSE, compute them
  # here on demand using the same internal routine.
  if(!is.null(fit@Optim$SE) &&
     !is.null(fit@Optim$SE$se) &&
     !is.null(fit@Optim$SE$vcov)) {

    result <- fit@Optim$SE

  } else {

    result <- compute_se_lcfa(dataList = fit@dataList,
                              modelInfo = fit@modelInfo,
                              Optim = fit@Optim)

  }

  # The asymptotic covariance of the sample statistics is selected when lcfa()
  # is fitted through the acov argument. The historical type argument is kept
  # for compatibility with the se() interface.
  result$type <- type

  #### Result ####

  return(result)

}

#### Function to compute lcfa standard errors before object construction ####

compute_se_lcfa <- function(dataList, modelInfo, Optim) {

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

  #### Hessian of the fitted CFA model ####

  control_optimizer <- modelInfo$control_optimizer
  control_optimizer$parameters[[1L]] <- Optim$parameters
  control_optimizer$transparameters[[1L]] <- Optim$transparameters

  H <- get_hess(control_manifold = modelInfo$control_manifold,
                control_transform = modelInfo$control_transform,
                control_estimator = modelInfo$control_estimator,
                control_optimizer = control_optimizer,
                cores = 1L)$h

  rownames(H) <- colnames(H) <- modelInfo$parameters_labels

  #### Sandwich covariance matrix ####

  H_inv <- solve(H)
  VCOV <- H_inv %*% B %*% H_inv
  VCOV <- (VCOV+t(VCOV))/2
  rownames(VCOV) <- colnames(VCOV) <- modelInfo$parameters_labels

  SE <- sqrt(diag(VCOV))
  names(SE) <- modelInfo$parameters_labels

  table_se <- fill_in(modelInfo$param, SE, miss = 0)

  #### Result ####

  result <- list(table_se = table_se,
                 se = SE,
                 vcov = VCOV,
                 H = H,
                 B = B,
                 ACOV = ACOV,
                 C = C)

  return(result)

}

#### Auxiliary functions ####

block_diag <- function(mats) {

  if(!is.list(mats) || length(mats) == 0L) {
    stop("mats must be a non-empty list of square matrices or nested lists.")
  }

  flatten_mats <- function(x) {

    if(is.matrix(x)) {

      result <- list(x)

    } else if(is.list(x)) {

      result <- unlist(lapply(x, flatten_mats), recursive = FALSE)

    } else {

      stop("All elements must be matrices or lists containing matrices.")

    }

    #### Result ####

    return(result)

  }

  mats <- flatten_mats(mats)

  if(length(mats) == 0L) {
    stop("No matrices found in mats.")
  }

  dims <- vapply(mats, FUN = function(M) {

    nr <- nrow(M)
    nc <- ncol(M)

    if(nr != nc) {
      stop("All matrices must be square.")
    }

    #### Result ####

    return(nr)

  }, FUN.VALUE = integer(1L))

  n_tot <- sum(dims)

  rn <- unlist(Map(function(M, d) {

    if(is.null(rownames(M))) {
      result <- rep(NA_character_, d)
    } else {
      result <- rownames(M)
    }

    #### Result ####

    return(result)

  }, mats, dims), use.names = FALSE)

  cn <- unlist(Map(function(M, d) {

    if(is.null(colnames(M))) {
      result <- rep(NA_character_, d)
    } else {
      result <- colnames(M)
    }

    #### Result ####

    return(result)

  }, mats, dims), use.names = FALSE)

  result <- matrix(0,
                   nrow = n_tot,
                   ncol = n_tot,
                   dimnames = list(rn, cn))

  idx <- 0L

  for(k in seq_along(mats)) {

    d <- dims[k]
    r <- (idx+1L):(idx+d)
    result[r, r] <- mats[[k]]
    idx <- idx+d

  }

  #### Result ####

  return(result)

}

# Backward-compatible name used by older code.
general_se <- function(fit, type = "standard") {

  result <- compute_se_lcfa(dataList = fit@dataList,
                            modelInfo = fit@modelInfo,
                            Optim = fit@Optim)

  #### Result ####

  return(result)

}
