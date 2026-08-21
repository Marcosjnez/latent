# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 21/08/2026
#'
#' Inspect Fitted CFA Objects
#'
#' @param fit A fitted object inheriting from class \code{"lcfa"}.
#' @param what Character string identifying the requested component.
#'
#' @return A parameter list, residual list, fit matrix, or estimator-specific
#'   control component, depending on \code{what}.
#'
#' @method latInspect lcfa
#' @export
latInspect.lcfa <- function(fit, what = "est") {

  #### Check inputs ####

  if(!inherits(fit, "lcfa")) {
    stop("fit must inherit from class 'lcfa'.")
  }

  if(!is.character(what) ||
     length(what) != 1L ||
     is.na(what)) {
    stop("what must be a single character string.")
  }

  what <- tolower(what)
  data_param <- fit@dataList$data_param

  #### Parameter blocks ####

  block_names <- unique(c(
    data_param$lambda_group,
    data_param$theta_group,
    data_param$psi_group,
    data_param$nu_group,
    data_param$delta_group
  ))
  block_names <- intersect(block_names,
                           names(fit@transformed_pars))

  if(what %in% c("est", "estimates", "parameters",
                 "fixed")) {

    result <- fit@transformed_pars[block_names]

  } else if(what == "items") {

    result <- fit@dataList$item_label

  } else if(what %in% c("rhat", "model",
                        "implied", "implied.cov")) {

    names_model <- intersect(
      data_param$model_group,
      names(fit@transformed_pars)
    )
    result <- fit@transformed_pars[names_model]

  } else if(what %in% c("resid", "residuals")) {

    result <- lcfa_residuals(fit)

  } else if(what %in% c("lambda", "loadings")) {

    names_lambda <- intersect(
      data_param$lambda_group,
      names(fit@transformed_pars)
    )
    result <- fit@transformed_pars[names_lambda]

  } else if(what == "psi") {

    names_psi <- intersect(
      data_param$psi_group,
      names(fit@transformed_pars)
    )
    result <- fit@transformed_pars[names_psi]

  } else if(what == "theta") {

    names_theta <- intersect(
      data_param$theta_group,
      names(fit@transformed_pars)
    )
    result <- fit@transformed_pars[names_theta]

  } else if(what %in% c("nu", "means")) {

    names_nu <- intersect(
      data_param$nu_group,
      names(fit@transformed_pars)
    )
    result <- fit@transformed_pars[names_nu]

  } else if(what == "uniquenesses") {

    names_theta <- intersect(
      data_param$theta_group,
      names(fit@transformed_pars)
    )
    result <- lapply(
      fit@transformed_pars[names_theta],
      FUN = diag
    )

  } else if(what == "w") {

    result <- lapply(
      fit@modelInfo$control_estimator,
      FUN = \(x) x$W
    )
    result <- result[
      !vapply(result, is.null, logical(1L))
    ]

  } else if(what == "weights") {

    result <- lapply(
      fit@modelInfo$control_estimator,
      FUN = \(x) x$w
    )
    result <- result[
      !vapply(result, is.null, logical(1L))
    ]

  } else if(what %in% c("loss", "f")) {

    fit_matrix <- lcfa_fit_matrix(fit)
    result <- fit_matrix[
      c("loss", "penalized_loss",
        "loss_base", "loss_sat"),
      ,
      drop = FALSE
    ]

  } else if(what %in% c("loglik", "logl")) {

    fit_matrix <- lcfa_fit_matrix(fit)
    result <- fit_matrix[
      c("loglik", "penalized_loglik",
        "loglik_base", "loglik_sat"),
      ,
      drop = FALSE
    ]

  } else if(what %in% c("fit", "fit.matrix",
                        "fit_matrix")) {

    result <- lcfa_fit_matrix(fit)

  } else if(what %in% c("fit.components",
                        "fit_components")) {

    result <- lcfa_fit_components(fit)

  } else if(what %in% c("vcov", "covariance")) {

    result <- tryCatch(
      fit@Optim$SE$VCOV,
      error = function(e) NULL
    )

    if(is.null(result)) {
      stop("No variance-covariance matrix is stored in the fitted object.")
    }

  } else if(what %in% c("se", "standard.errors",
                        "standard_errors")) {

    result <- tryCatch(
      fit@Optim$SE$se,
      error = function(e) NULL
    )

    if(is.null(result)) {
      stop("No standard errors are stored in the fitted object.")
    }

  } else {

    stop("Unknown request: ", what)

  }

  #### Result ####

  return(result)

}
