# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 13/09/2026
#'
#' Inspect Fitted CFA and SEM Objects
#'
#' @param fit A fitted object inheriting from class \code{"lcfa"}.
#' @param what Character string identifying the requested component. Use
#'   \code{"est"} or \code{"labels"} for the model parameters, or request a
#'   block such as \code{"lambda"}, \code{"theta"}, \code{"psi"} or
#'   \code{"intercepts"}. Use \code{"implied"}, \code{"implied.means"},
#'   \code{"implied.thresholds"} or \code{"residuals"} for model moments.
#'   \code{"loglik"} and \code{"loss"} return fit components by group and
#'   overall. Log-likelihoods are available for likelihood estimators.
#'   Use \code{"convergence"}, \code{"gradient"} or \code{"elapsed"} for
#'   optimizer diagnostics. Elapsed time is in seconds.
#' @param sort Logical. Sort and orient factor-specific outputs. Defaults to TRUE.
#'
#' @return A parameter list, residual list, fit matrix, or estimator-specific
#'   control component, depending on \code{what}.
#'
#' @param x Fitted model or collection of fitted models passed to the method.
#' @param ... Additional arguments to the generic.
#' @method latInspect lcfa
#' @export
latInspect.lcfa <- function(x, what = "est", sort = TRUE, ..., fit = x) {
  if(!missing(fit) && !missing(x)) return(.match_legacy_model("fit"))

  #### Check inputs ####

  if(!inherits(fit, "lcfa")) {
    stop("fit must inherit from class 'lcfa'.")
  }

  if(!is.character(what) ||
     length(what) != 1L ||
     is.na(what)) {
    stop("what must be a single character string.")
  }

  if(!is.logical(sort) || length(sort) != 1L || is.na(sort)) {
    stop("sort must be TRUE or FALSE.")
  }

  what <- tolower(what)
  if(!length(fit@Optim) && what %in% c("resid", "residuals", "gradient", "convergence")) {
    stop("Fit the model to inspect ", what, ".")
  }
  data_param <- fit@dataList$data_param

  #### Parameter blocks ####

  groups <- switch(what,
    labels =, est =, estimates =, parameters =, fixed =
      c("lambda_group", "theta_group", "psi_group", "alpha_group", "nu_group",
        "B_group", "A_group", "latent_cov_group", "latent_means_group",
        "kappa_group", "delta_group"),
    rhat =, model =, implied =, implied.cov =, implied_cov = "model_group",
    lambda =, loadings = "lambda_group",
    psi = "psi_group",
    b =, beta =, structural.coefficients =, structural_coefficients = "B_group",
    a =, structural.inverse =, structural_inverse = "A_group",
    latent.cov =, latent_cov =, cov.lv =, cov_lv =
      if(is.null(data_param$latent_cov_group)) "psi_group" else "latent_cov_group",
    theta =, uniquenesses = "theta_group",
    alpha = "alpha_group",
    latent.means =, latent_means =, mean.lv =, mean_lv =
      if(is.null(data_param$latent_means_group)) "alpha_group" else "latent_means_group",
    nu =, intercepts = "nu_group",
    means =, meanshat =, implied.means =, implied_means = "meanshat_group",
    kappa =, kappas =, unstandardized.thresholds =, unstandardized_thresholds = "kappa_group",
    tauhat =, implied.thresholds =, implied_thresholds = "tauhat_group"
  )

  if(!is.null(groups)) {

    labels <- what == "labels"
    structures <- if(sort) {
      output <- sort_latent(fit)
      if(labels) output$trans else output$transformed_pars
    } else {
      if(labels) fit@modelInfo$trans else fit@transformed_pars
    }
    block_names <- unlist(data_param[groups], use.names = FALSE)
    result <- structures[intersect(block_names, names(structures))]
    if(what == "uniquenesses") result <- lapply(result, diag)

  } else if(what == "items") {

    result <- fit@dataList$item_label

  } else if(what %in% c("resid", "residuals")) {

    result <- lcfa_residuals(fit)

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

  } else if(what %in% c("elapsed", "timing")) {

    result <- fit@timing

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
