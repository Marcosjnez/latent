# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 17/08/2026
#'
#' Standard Errors for Latent Objects
#'
#' @param fit A fitted object inheriting from class \code{"latent"}.
#' @param type Character string selecting \code{"information"} or
#'   \code{"robust"} covariance estimation.
#' @param parameters Optional parameter specification.
#' @param digits Non-negative integer used to format parameter tables, or
#'   \code{NULL} to avoid rounding.
#' @param ... Additional arguments.
#'
#' @return A list containing parameter tables, standard errors, covariance
#'   matrices, and derivative information.
#'
#' @method se latent
#' @export
se.latent <- function(fit, type = "information", parameters = NULL,
                      digits = 4L, ...) {

  #### Check inputs ####

  if(!inherits(fit, "latent")) {
    stop("fit must inherit from class 'latent'.")
  }

  if(length(fit@Optim$transparameters) == 0L) {
    stop("The latent object has not been fitted.")
  }

  if(is.null(parameters)) {

    parameters <- fit@modelInfo$trans[names(fit@modelInfo$param)]

  } else {

    selected_parameters <- unique(unlist(parameters))

    if(!all(selected_parameters %in%
            fit@modelInfo$transparameters_labels)) {
      stop("Unknown parameters.")
    }

  }

  type <- match.arg(tolower(type), c("information", "robust"))

  if(!is.null(digits) &&
     (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits < 0L || digits != as.integer(digits))) {
    stop("digits must be NULL or a non-negative integer.")
  }

  #### Covariance of the untransformed parameters ####

  if(type == "information") {
    base_VCOV <- information(fit)
  } else {
    base_VCOV <- robust(fit)
  }

  if(is.null(base_VCOV$VCOV)) {
    stop("The selected standard-error method did not return a VCOV matrix.")
  }

  #### Covariance of selected transformed parameters ####

  SE <- vcov(fit,
             parameters = parameters,
             v = base_VCOV$VCOV)

  est <- fill_in(parameters,
                 fit@Optim$transparameters,
                 miss = NA)
  table_se <- fill_in(parameters,
                      SE$se,
                      miss = NA)
  table <- combine_est_se(est, table_se, digits = digits)

  B <- base_VCOV$B
  if(is.null(B)) {
    B <- matrix(numeric(0L), nrow = 0L, ncol = 0L)
  }

  #### Result ####

  result <- list(table = table,
                 table_se = table_se,
                 se = c(SE$se),
                 vcov = SE$vcov,
                 VCOV = SE$vcov,
                 H = base_VCOV$H,
                 B = B,
                 C = base_VCOV$C,
                 jacob = SE$jacob,
                 type = type)

  return(result)

}
