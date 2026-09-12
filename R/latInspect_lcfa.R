# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 13/09/2026
#'
#' Inspect Fitted CFA Objects
#'
#' @param fit A fitted object inheriting from class \code{"lcfa"}.
#' @param what Character string identifying the requested component.
#' @param sort Logical. Sort/orient only the requested output; defaults to TRUE.
#'
#' @return A parameter list, residual list, fit matrix, or estimator-specific
#'   control component, depending on \code{what}.
#'
#' @method latInspect lcfa
#' @export
latInspect.lcfa <- function(fit, what = "est", sort = TRUE) {

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
  output <- latInspect.latent(fit, what = "structures", sort = sort)
  transformed_pars <- output$transformed_pars

  #### Parameter blocks ####

  block_names <- unique(c(
    data_param$lambda_group,
    data_param$alpha_group,
    data_param$theta_group,
    data_param$psi_group,
    data_param$nu_group,
    data_param$kappa_group,
    data_param$delta_group
  ))
  block_names <- intersect(block_names,
                           names(transformed_pars))

  if(what %in% c("est", "estimates", "parameters",
                 "fixed")) {

    result <- transformed_pars[block_names]

  } else if(what == "items") {

    result <- fit@dataList$item_label

  } else if(what %in% c("rhat", "model",
                        "implied", "implied.cov",
                        "implied_cov")) {

    names_model <- intersect(
      data_param$model_group,
      names(transformed_pars)
    )
    result <- transformed_pars[names_model]

  } else if(what %in% c("resid", "residuals")) {

    result <- lcfa_residuals(fit)

  } else if(what %in% c("lambda", "loadings")) {

    names_lambda <- intersect(
      data_param$lambda_group,
      names(transformed_pars)
    )
    result <- transformed_pars[names_lambda]

  } else if(what == "psi") {

    names_psi <- intersect(
      data_param$psi_group,
      names(transformed_pars)
    )
    result <- transformed_pars[names_psi]

  } else if(what == "theta") {

    names_theta <- intersect(
      data_param$theta_group,
      names(transformed_pars)
    )
    result <- transformed_pars[names_theta]

  } else if(what == "alpha") {

    names_alpha <- intersect(
      data_param$alpha_group,
      names(transformed_pars)
    )
    result <- transformed_pars[names_alpha]

  } else if(what %in% c("mu", "latent.means",
                        "latent_means")) {

    mu_group <- data_param$mu_group
    if(is.null(mu_group)) mu_group <- data_param$alpha_group
    names_mu <- intersect(
      mu_group,
      names(transformed_pars)
    )
    result <- transformed_pars[names_mu]

  } else if(what %in% c("nu", "intercepts")) {

    names_nu <- intersect(
      data_param$nu_group,
      names(transformed_pars)
    )
    result <- transformed_pars[names_nu]

  } else if(what %in% c("means", "meanshat",
                        "implied.means", "implied_means")) {

    names_meanshat <- intersect(
      data_param$meanshat_group,
      names(transformed_pars)
    )
    result <- transformed_pars[names_meanshat]

  } else if(what %in% c("kappa", "kappas",
                        "unstandardized.thresholds",
                        "unstandardized_thresholds")) {

    names_kappa <- intersect(
      data_param$kappa_group,
      names(transformed_pars)
    )
    result <- transformed_pars[names_kappa]

  } else if(what %in% c("tauhat", "implied.thresholds",
                        "implied_thresholds")) {

    names_tauhat <- intersect(
      data_param$tauhat_group,
      names(transformed_pars)
    )
    result <- transformed_pars[names_tauhat]

  } else if(what == "uniquenesses") {

    names_theta <- intersect(
      data_param$theta_group,
      names(transformed_pars)
    )
    result <- lapply(
      transformed_pars[names_theta],
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

  } else if(what == "convergence") {

    result <- data.frame(Iterations = fit@Optim$iterations,
                         convergence = fit@Optim$convergence,
                         grad.norm = fit@Optim$ng)

  } else if(what == "gradient") {

    gradient_info <- data.frame(name = fit@modelInfo$parameters_labels,
                                gradient = fit@Optim$g,
                                rgradient = fit@Optim$rg,
                                dir = fit@Optim$dir)

    result <- list(grad.norm = fit@Optim$ng, gradient.info = gradient_info)

  } else {

    result <- latInspect.latent(fit, what = what, sort = sort)

  }

  #### Result ####

  return(result)

}
