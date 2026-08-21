# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 21/08/2026
#'
#' Robust Variance-Covariance Matrix for General Latent Models
#'
#' @param fit A fitted object inheriting from class \code{"latent"}.
#'
#' @return A stored robust covariance result when available. If no
#'   class-specific robust covariance exists, the information covariance is
#'   returned with a warning and an explicit fallback label.
#'
#' @method robust latent
#' @export
robust.latent <- function(fit) {

  labels <- fit@modelInfo$parameters_labels
  function_name <- tail(as.character(fit@call[[1L]]), 1L)

  stored_SE <- tryCatch(
    fit@Optim$SE,
    error = function(e) NULL
  )
  stored_type <- tryCatch(
    tolower(stored_SE$type),
    error = function(e) character(0L)
  )

  if(length(stored_type) == 0L) {

    if(function_name == "lpearson") {
      stored_type <- tolower(fit@dataList$VCOV_type)
    }

  }

  stored_robust <- length(stored_type) == 1L &&
    identical(stored_type, "robust") &&
    !is.null(stored_SE$VCOV)

  if(function_name == "lpearson" &&
     !stored_robust) {

    control <- fit@modelInfo$control_optimizer
    fit@dataList$VCOV_type <- "robust"
    stored_SE <- compute_se_lpearson(
      dataList = fit@dataList,
      modelInfo = fit@modelInfo,
      Optim = fit@Optim,
      control = control
    )
    stored_SE$type <- "robust"
    stored_type <- "robust"
    stored_robust <- TRUE

  }

  if(length(stored_type) == 1L &&
     identical(stored_type, "robust") &&
     !is.null(stored_SE$VCOV)) {

    VCOV <- validate_covariance_matrix(
      stored_SE$VCOV,
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

    return(list(H = H,
                B = B,
                VCOV = VCOV,
                se = se,
                type = "robust"))

  }

  warning("No class-specific robust covariance method is implemented for ",
          paste(class(fit), collapse = "/"),
          "; the information covariance matrix was returned.")

  result <- information(fit)
  result$type <- "information_fallback"

  #### Result ####

  return(result)

}
