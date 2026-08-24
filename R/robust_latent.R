# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 24/08/2026
#'
#' Robust Variance-Covariance Matrix for Latent Models
#'
#' Compute or retrieve the robust variance-covariance matrix of the freely
#' estimated parameters of a fitted latent-variable model.
#'
#' @param fit A fitted object inheriting from class \code{"latent"}.
#'
#' @return A list containing the Hessian when available, the empirical score
#'   covariance when available, the variance-covariance matrix, standard errors,
#'   and covariance-method metadata.
#'
#' @details
#' Deterministic multistep models use the same propagated covariance for
#' \code{information()} and \code{robust()} because the covariance method of
#' every preceding statistic is fixed when the multistep model is fitted.
#'
#' Pearson covariance and correlation estimators compute their robust covariance
#' directly from the observed data. Other latent models reuse a stored robust
#' covariance when one is available. If no robust covariance has been defined,
#' the information covariance is returned with a warning.
#'
#' Direct robust likelihood scores are not yet available for ordinary
#' \code{"lcfa"} likelihood models.
#'
#' @method robust latent
#' @export
robust.latent <- function(fit) {

  #### Check inputs ####

  if(!inherits(fit, "latent")) {
    stop("fit must inherit from class 'latent'.")
  }

  if(length(fit@Optim$transparameters) == 0L) {
    stop("The latent object has not been fitted.")
  }

  #### Deterministic multistep covariance ####

  if(inherits(fit, "multistep")) {

    result <- information.latent(fit)

    #### Result ####

    return(result)

  }

  #### Direct CFA likelihood ####

  if(inherits(fit, "lcfa")) {
    stop("Robust standard errors are not yet implemented for ordinary lcfa ",
         "likelihood models. Use information standard errors, or a multistep ",
         "analysis whose parent statistic supplies a robust VCOV matrix.")
  }

  labels <- fit@modelInfo$parameters_labels
  function_name <- tail(as.character(fit@call[[1L]]), 1L)

  #### Pearson robust covariance ####

  if(function_name == "lpearson") {

    control <- fit@modelInfo$control_optimizer
    fit@dataList$VCOV_type <- "robust"

    robust_SE <- compute_se_lpearson(
      dataList = fit@dataList,
      modelInfo = fit@modelInfo,
      Optim = fit@Optim,
      control = control
    )

    VCOV <- validate_covariance_matrix(
      robust_SE$VCOV,
      labels = labels,
      object_name = "Pearson robust variance-covariance matrix"
    )

    se <- standard_errors_from_vcov(
      VCOV,
      object_name = "Pearson robust variance-covariance matrix"
    )

    #### Result ####

    result <- list(H = NULL,
                   B = NULL,
                   VCOV = VCOV,
                   se = se,
                   type = "robust")

    return(result)

  }

  #### Stored robust covariance ####

  stored_SE <- tryCatch(fit@Optim$SE,
                        error = function(e) NULL)
  stored_type <- tryCatch(tolower(stored_SE$type),
                          error = function(e) character(0L))
  stored_VCOV <- tryCatch(stored_SE$VCOV,
                          error = function(e) NULL)

  if(length(stored_type) == 1L &&
     identical(stored_type, "robust") &&
     !is.null(stored_VCOV)) {

    VCOV <- validate_covariance_matrix(
      stored_VCOV,
      labels = labels,
      object_name = "stored robust variance-covariance matrix"
    )

    se <- standard_errors_from_vcov(
      VCOV,
      object_name = "stored robust variance-covariance matrix"
    )

    H <- tryCatch(stored_SE$H,
                  error = function(e) NULL)
    B <- tryCatch(stored_SE$B,
                  error = function(e) NULL)

    #### Result ####

    result <- list(H = H,
                   B = B,
                   VCOV = VCOV,
                   se = se,
                   type = "robust")

    return(result)

  }

  warning("No robust covariance has been defined for ",
          paste(class(fit), collapse = "/"),
          "; the information covariance matrix was returned.")

  result <- information.latent(fit)
  result$type <- "information_fallback"

  #### Result ####

  return(result)

}
