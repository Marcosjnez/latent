# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 15/08/2026
#'
#' Variances-Covariance Matrix for Latent Objects
#'
#' @export
se.latent <- function(fit, type = "information", parameters = NULL,
                      digits = 4L, ...) {

  if(is.null(parameters)) {
    parameters <- fit@modelInfo$trans[names(fit@modelInfo$param)]
    # parameters <- fit@modelInfo$parameters_labels
  } else if(!any(unlist(parameters) %in%
                 fit@modelInfo$transparameters_labels)) {
    stop("Unknown parameters.")
  }

  if(!inherits(fit, "latent")) {
    stop("fit must inherit from class 'latent'.")
  }

  if(length(fit@Optim$parameters) == 0L) {
    stop("The llca object has not been fitted.")
  }

  type <- match.arg(tolower(type), c("information", "robust"))

  if(!is.null(digits) &&
     (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits < 0L || digits != as.integer(digits))) {
    stop("digits must be NULL or a non-negative integer.")
  }

  if(type == "information") {
    se_method <- information
  } else if(type == "robust") {
    se_method <- robust
  } else {
    stop("Unknown type method. Available types: 'information' and 'robust'")
  }

  VCOV <- se_method(fit)
  SE <- vcov(fit, parameters = parameters, v = VCOV$VCOV)

  est <- fill_in(parameters, fit@Optim$transparameters, miss = NA)
  table_se <- fill_in(parameters, SE$se, miss = NA)
  table <- combine_est_se(est, table_se, digits = digits)

  result <- list(table = table, table_se = table_se, se = c(SE$se),
                 vcov = SE$vcov, H = SE$H, B = SE$B, jacob = SE$jacob)

  #### Result ####

  return(result)

}
