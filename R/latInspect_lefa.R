# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 23/08/2026
#'
#' Inspect Fitted Exploratory Factor Analysis Objects
#'
#' @param fit A fitted object inheriting from class \code{"lefa"}.
#' @param what Character string identifying the requested component.
#'
#' @return A rotated parameter list, the unrotated \code{lcfa} object, or a
#'   component delegated to \code{latInspect.lcfa()}.
#'
#' @method latInspect lefa
#' @export
latInspect.lefa <- function(fit, what = "est") {

  #### Check inputs ####

  if(!inherits(fit, "lefa")) {
    stop("fit must inherit from class 'lefa'.")
  }

  if(!is.character(what) ||
     length(what) != 1L ||
     is.na(what)) {
    stop("what must be a single character string.")
  }

  what <- tolower(what)
  data_param <- fit@modelInfo$data_param
  source <- lefa_source_lcfa(fit)

  if(what %in% c("est", "estimates", "parameters",
                 "fixed")) {

    source_est <- latInspect(source, what = "est")
    source_data_param <- source@dataList$data_param
    unrotated_blocks <- unique(c(source_data_param$lambda_group,
                                 source_data_param$psi_group,
                                 source_data_param$alpha_group))
    source_est <- source_est[
      setdiff(names(source_est), unrotated_blocks)
    ]

    rotated_blocks <- lefa_rotated_block_names(fit)
    result <- c(fit@transformed_pars[rotated_blocks],
                source_est)

  } else if(what %in% c("lambda", "loadings")) {

    names_lambda <- intersect(data_param$lambda_group,
                              names(fit@transformed_pars))
    result <- fit@transformed_pars[names_lambda]

  } else if(what == "psi") {

    names_psi <- intersect(data_param$psi_group,
                           names(fit@transformed_pars))
    result <- fit@transformed_pars[names_psi]

  } else if(what %in% c("alpha", "latent.means",
                        "latent_means")) {

    names_alpha <- intersect(data_param$alpha_group,
                             names(fit@transformed_pars))
    result <- fit@transformed_pars[names_alpha]

  } else if(what %in% c("x", "rotation", "rotation.matrix",
                        "rotation_matrix")) {

    names_X <- lefa_rotation_block_names(fit)
    result <- fit@transformed_pars[names_X]

  } else if(what %in% c("efa", "unrotated.fit",
                        "unrotated_fit")) {

    result <- source

  } else if(what %in% c("unrotated", "unrotated.est",
                        "unrotated_est")) {

    result <- latInspect(source, what = "est")

  } else if(what %in% c("unrotated.lambda",
                        "unrotated_lambda")) {

    result <- latInspect(source, what = "lambda")

  } else if(what %in% c("unrotated.psi",
                        "unrotated_psi")) {

    result <- latInspect(source, what = "psi")

  } else if(what %in% c("unrotated.alpha",
                        "unrotated_alpha")) {

    result <- latInspect(source, what = "alpha")

  } else if(what %in% c("rotation.loss",
                        "rotation_loss")) {

    result <- fit@Optim$f

  } else if(what %in% c("vcov", "covariance")) {

    result <- se(fit, digits = NULL)$VCOV

  } else if(what %in% c("se", "standard.errors",
                        "standard_errors")) {

    result <- se(fit, digits = NULL)$se

  } else {

    result <- latInspect(source, what = what)

  }

  #### Result ####

  return(result)

}
