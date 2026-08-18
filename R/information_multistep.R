# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 18/08/2026
#'
#' Information Covariance for Multistep Models
#'
#' Return the multistep covariance matrix implied by the sample-statistic
#' covariance estimators selected when the model was fitted.
#'
#' @param fit A fitted object inheriting from class \code{"multistep"}.
#'
#' @return A list containing the top-level Hessian, covariance matrix, and
#'   standard errors of the freely estimated parameters.
#'
#' @method information multistep
#' @export
information.multistep <- function(fit) {

  SE <- se.multistep(fit = fit,
                     parameters = fit@modelInfo$parameters_labels)

  VCOV <- SE$free_VCOV
  labels <- fit@modelInfo$parameters_labels

  rownames(VCOV) <- colnames(VCOV) <- labels

  se <- sqrt(diag(VCOV))
  names(se) <- labels

  #### Result ####

  result <- list(H = SE$H,
                 B = SE$B,
                 C = SE$C,
                 VCOV = VCOV,
                 se = se,
                 joint_vcov = SE$joint_vcov)

  return(result)

}

# Explicit concrete-class method for deterministic S3 dispatch with S4 multiple
# inheritance.
#'
#' @rdname information.multistep
#' @method information multistep_lcfa
#' @export
information.multistep_lcfa <- function(fit) {

  result <- information.multistep(fit)

  #### Result ####

  return(result)

}
