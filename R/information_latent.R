# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 21/08/2026
#'
#' Information Variance-Covariance Matrix for Latent Models
#'
#' @param fit A fitted object inheriting from class \code{"latent"}.
#'
#' @return A list containing the Hessian, variance-covariance matrix, standard
#'   errors, and covariance-method metadata.
#'
#' @details
#' The fitted object is never modified. For ordinary CFA likelihood models, a
#' temporary derivative copy replaces \code{cfa_fml} by \code{cfa_ml} and
#' \code{cfa_means_fml} by \code{cfa_means_ml}. Stored covariance matrices are
#' reused only when they were explicitly produced by an information/standard
#' method.
#'
#' @method information latent
#' @export
information.latent <- function(fit) {

  if(inherits(fit, "multistep")) {

    result <- information.multistep(fit)

    #### Result ####

    return(result)

  }

  labels <- fit@modelInfo$parameters_labels

  if(length(labels) == 0L) {

    empty <- matrix(numeric(0L), nrow = 0L, ncol = 0L,
                    dimnames = list(labels, labels))

    #### Result ####

    return(list(H = empty, VCOV = empty, se = numeric(0L),
                type = "information"))

  }

  stored_SE <- tryCatch(fit@Optim$SE,
                        error = function(e) NULL)
  stored_type <- tryCatch(tolower(stored_SE$type),
                          error = function(e) character(0L))
  stored_VCOV <- tryCatch(stored_SE$VCOV,
                          error = function(e) NULL)

  if(length(stored_type) == 0L &&
     !is.null(stored_VCOV)) {

    function_name <- tail(as.character(fit@call[[1L]]), 1L)

    if(function_name == "lpearson") {
      stored_type <- tolower(fit@dataList$VCOV_type)
    } else if(function_name %in%
              c("lmean", "lpoly", "lyule", "lmvnorm")) {
      stored_type <- "information"
    }

  }

  function_name <- tail(as.character(fit@call[[1L]]), 1L)

  if(function_name == "lpearson" &&
     length(stored_type) == 1L &&
     identical(stored_type, "robust")) {

    control <- fit@modelInfo$control_optimizer
    std_ov <- isTRUE(control$std.ov)

    VCOV <- asymptotic_normal(
      fit@dataList$S,
      cov = !std_ov,
      diag = FALSE
    )/fit@dataList$nobs
    VCOV <- validate_covariance_matrix(
      VCOV,
      labels = labels,
      object_name = "Pearson information variance-covariance matrix"
    )
    se <- standard_errors_from_vcov(
      VCOV,
      object_name = "Pearson information variance-covariance matrix"
    )

    #### Result ####

    return(list(H = NULL,
                VCOV = VCOV,
                se = se,
                type = "information"))

  }

  reuse_stored <- length(stored_type) == 1L &&
    stored_type %in% c("standard", "information") &&
    !is.null(stored_VCOV)

  if(!reuse_stored &&
     is.null(stored_VCOV)) {

    lazy_SE <- lazy_information_se_latent(
      fit = fit,
      function_name = function_name
    )

    if(!is.null(lazy_SE)) {
      stored_SE <- lazy_SE
      stored_VCOV <- lazy_SE$VCOV
      stored_type <- "information"
      reuse_stored <- TRUE
    }

  }

  if(reuse_stored) {

    VCOV <- validate_covariance_matrix(
      stored_VCOV,
      labels = labels,
      object_name = "stored information variance-covariance matrix"
    )

    H <- tryCatch(stored_SE$H,
                  error = function(e) NULL)

    if(!is.null(H)) {
      H <- validate_covariance_matrix(H, labels = labels,
                                      object_name = "stored Hessian")
    }

  } else {

    fit_information <- information_model_latent(fit)

    H <- hessian(fit_information)
    H <- validate_covariance_matrix(H, labels = labels,
                                    object_name = "Hessian")

    VCOV <- invert_information_matrix(
      H,
      labels = labels,
      object_name = "Hessian"
    )

  }

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

#### Temporary information model ####

information_model_latent <- function(fit) {

  result <- fit

  # if(isTRUE(result@modelInfo$control_optimizer$reg)) {
  #   stop("Information standard errors for penalized latent models require a ",
  #        "penalty-specific covariance derivation and are not currently ",
  #        "implemented.")
  # }

  if(inherits(result, "lcfa") &&
     !is.null(result@dataList$likelihood) &&
     !identical(tolower(result@dataList$likelihood), "normal")) {
    stop("Information standard errors for lcfa currently require ",
         "likelihood = 'normal'; Wishart mean/covariance scaling has not ",
         "yet been implemented in the temporary information estimator.")
  }

  if(length(result@modelInfo$control_estimator) == 0L) {

    #### Result ####

    return(result)

  }

  for(i in seq_along(result@modelInfo$control_estimator)) {

    estimator <- result@modelInfo$control_estimator[[i]]$estimator

    if(identical(estimator, "cfa_fml")) {
      result@modelInfo$control_estimator[[i]]$estimator <- "cfa_ml"
    } else if(identical(estimator, "cfa_means_fml")) {
      result@modelInfo$control_estimator[[i]]$estimator <- "cfa_means_ml"
    }

  }

  #### Result ####

  return(result)

}


#### Deferred information covariance for source estimators ####

lazy_information_se_latent <- function(fit, function_name = NULL) {

  if(is.null(function_name)) {
    function_name <- tail(as.character(fit@call[[1L]]), 1L)
  }

  control <- fit@modelInfo$control_optimizer

  result <- switch(
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

  if(is.null(result)) {

    #### Result ####

    return(NULL)

  }

  result$VCOV <- validate_covariance_matrix(
    result$VCOV,
    labels = fit@modelInfo$parameters_labels,
    object_name = paste0(function_name,
                         " information variance-covariance matrix")
  )
  result$se <- standard_errors_from_vcov(
    result$VCOV,
    object_name = paste0(function_name,
                         " information variance-covariance matrix")
  )
  result$type <- "information"

  #### Result ####

  return(result)

}
