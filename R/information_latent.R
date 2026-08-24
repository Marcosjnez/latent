# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 24/08/2026
#'
#' Information Variance-Covariance Matrix for Latent Models
#'
#' Compute the variance-covariance matrix of the freely estimated parameters of
#' a fitted latent-variable model.
#'
#' @param fit A fitted object inheriting from class \code{"latent"}.
#'
#' @return A list containing the Hessian, variance-covariance matrix, standard
#'   errors, and covariance-method metadata. Multistep objects additionally
#'   return the propagation matrices and joint covariance matrix.
#'
#' @details
#' Three covariance calculations are handled directly. Deterministic multistep
#' models use the covariance propagated from preceding estimation steps. Source
#' estimators such as \code{lmean()}, \code{lpearson()}, \code{lpoly()},
#' \code{lyule()}, and \code{lmvnorm()} use their analytic information
#' covariance. Remaining latent models use the inverse Hessian.
#'
#' For ordinary CFA maximum-likelihood models, the Hessian is evaluated with the
#' corresponding \code{cfa_ml} estimator because \code{cfa_fml} represents the
#' fitted discrepancy function rather than the information contribution used
#' for the covariance matrix.
#'
#' @method information latent
#' @export
information.latent <- function(fit) {

  #### Check inputs ####

  if(!inherits(fit, "latent")) {
    stop("fit must inherit from class 'latent'.")
  }

  if(length(fit@Optim$transparameters) == 0L) {
    stop("The latent object has not been fitted.")
  }

  labels <- fit@modelInfo$parameters_labels

  #### Deterministic multistep covariance ####

  if(inherits(fit, "multistep")) {

    SE <- se.multistep(fit = fit,
                       parameters = labels,
                       digits = NULL)

    VCOV <- validate_covariance_matrix(
      SE$free_VCOV,
      labels = labels,
      object_name = "multistep variance-covariance matrix"
    )

    se <- standard_errors_from_vcov(
      VCOV,
      object_name = "multistep variance-covariance matrix"
    )

    #### Result ####

    result <- list(H = SE$H,
                   B = SE$B,
                   C = SE$C,
                   VCOV = VCOV,
                   se = se,
                   joint_vcov = SE$joint_vcov,
                   steps = SE$steps,
                   type = "multistep",
                   sample_se = SE$sample_se)

    return(result)

  }

  if(length(labels) == 0L) {

    empty <- matrix(numeric(0L), nrow = 0L, ncol = 0L,
                    dimnames = list(labels, labels))

    #### Result ####

    result <- list(H = empty,
                   VCOV = empty,
                   se = numeric(0L),
                   type = "information")

    return(result)

  }

  #### Analytic covariance of source estimators ####

  function_name <- tail(as.character(fit@call[[1L]]), 1L)
  control <- fit@modelInfo$control_optimizer

  source_SE <- switch(
    function_name,
    lmean = compute_se_lmean(
      dataList = fit@dataList,
      modelInfo = fit@modelInfo,
      Optim = fit@Optim,
      control = control
    ),
    lpearson = {

      fit@dataList$VCOV_type <- "standard"

      compute_se_lpearson(
        dataList = fit@dataList,
        modelInfo = fit@modelInfo,
        Optim = fit@Optim,
        control = control
      )

    },
    lpoly = compute_se_lpoly(
      dataList = fit@dataList,
      modelInfo = fit@modelInfo,
      Optim = fit@Optim
    ),
    lyule = {

      fit_output <- fit_lyule(
        dataList = fit@dataList,
        modelInfo = fit@modelInfo,
        control = control
      )

      compute_se_lyule(
        dataList = fit@dataList,
        modelInfo = fit@modelInfo,
        yule_output = fit_output$yule_output
      )

    },
    lmvnorm = compute_se_lmvnorm(
      dataList = fit@dataList,
      modelInfo = fit@modelInfo,
      Optim = fit@Optim
    ),
    NULL
  )

  if(!is.null(source_SE)) {

    VCOV <- validate_covariance_matrix(
      source_SE$VCOV,
      labels = labels,
      object_name = paste0(function_name,
                           " information variance-covariance matrix")
    )

    H <- tryCatch(source_SE$H,
                  error = function(e) NULL)

    if(!is.null(H)) {
      H <- validate_covariance_matrix(
        H,
        labels = labels,
        object_name = paste0(function_name, " Hessian")
      )
    }

    se <- standard_errors_from_vcov(
      VCOV,
      object_name = paste0(function_name,
                           " information variance-covariance matrix")
    )

    #### Result ####

    result <- list(H = H,
                   VCOV = VCOV,
                   se = se,
                   type = "information")

    return(result)

  }

  #### Information Hessian ####

  fit_information <- fit

  for(i in seq_along(fit_information@modelInfo$control_estimator)) {

    estimator <- fit_information@modelInfo$control_estimator[[i]]$estimator

    if(identical(estimator, "cfa_fml")) {
      fit_information@modelInfo$control_estimator[[i]]$estimator <- "cfa_ml"
    }

  }

  H <- hessian(fit_information)
  H <- validate_covariance_matrix(H,
                                  labels = labels,
                                  object_name = "Hessian")

  VCOV <- invert_information_matrix(H,
                                    labels = labels,
                                    object_name = "Hessian")

  se <- standard_errors_from_vcov(
    VCOV,
    object_name = "information variance-covariance matrix"
  )

  #### Result ####

  result <- list(H = H,
                 VCOV = VCOV,
                 se = se,
                 type = "information")

  return(result)

}
