# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 15/08/2026
#'
#' Variance-Covariance Matrix for Latent Models
#'
#' Compute covariance matrices and standard errors for fitted latent variable
#' models, including uncertainty propagated across sequential estimation steps.
#'
#' @description
#' \code{vcov.latent()} computes the variance-covariance matrix of freely
#' estimated or transformed parameters of a fitted object inheriting from
#' class \code{"latent"}.
#'
#' For ordinary one-step models, covariance matrices are obtained from the
#' inverse Hessian, or from an alternative information matrix supplied through
#' \code{H}, and transformed to the requested parameterization using the delta
#' method.
#'
#' When the fitted model contains a previous \code{"latent"} object among the
#' model specifications stored in
#' \code{fit@modelInfo$control_optimizer$model}, uncertainty from the previous
#' estimation step is propagated to the current parameter estimates. Nested
#' fitted models are handled recursively.
#'
#' @param fit A fitted object inheriting from class \code{"latent"}.
#' @param parameters Optional parameter specification identifying the parameters
#'   for which the covariance matrix should be returned. Labels must occur in
#'   \code{fit@modelInfo$transparameters_labels}. If \code{NULL}, the parameter
#'   blocks corresponding to the fitted model parameters are used.
#' @param H Optional information specification for the freely estimated
#'   parameters. If \code{NULL}, the ordinary Hessian is used. A numeric matrix
#'   can be supplied as the information matrix for the current model. A function
#'   can also be supplied; it must take a fitted model as its first argument and
#'   return an information matrix. Function-valued specifications are propagated
#'   recursively so that the information matrix is recomputed separately for
#'   each nested model.
#'
#' @details
#' For a one-step model, let \eqn{P} denote the information matrix used for
#' inference. With ordinary Hessian-based inference, \eqn{P = H}, where \eqn{H}
#' is the Hessian of the objective function. The covariance matrix of the freely
#' estimated parameters is based on \eqn{P^{-1}}. Covariances of transformed
#' parameters are subsequently obtained with the delta method.
#'
#' For a sequential model, suppose that measurement or nuisance parameters
#' \eqn{\theta_M} are estimated in an earlier step and subsequently treated as
#' fixed while structural parameters \eqn{\theta_S} are estimated. Let
#' \eqn{A} denote the covariance matrix of the earlier parameter estimates,
#' \eqn{H_2} the ordinary Hessian for the structural parameters, \eqn{P_2} the
#' information matrix whose inverse gives the structural covariance when the
#' earlier parameters are treated as known, and \eqn{C} the
#' measurement-structural cross-derivative matrix.
#'
#' The propagated structural uncertainty contains the additional term
#' \deqn{H_2^{-1} C^\top A C H_2^{-1}.}
#'
#' The corresponding joint covariance matrix is
#' \deqn{
#' \begin{pmatrix}
#' A & -A C H_2^{-1} \\
#' -H_2^{-1} C^\top A &
#' P_2^{-1} + H_2^{-1} C^\top A C H_2^{-1}
#' \end{pmatrix}.
#' }
#'
#' Rather than constructing and inverting this covariance matrix directly, the
#' method constructs its inverse and passes it through the ordinary covariance
#' machinery. When \eqn{P_2 = H_2}, this reduces to the standard sequential
#' pseudo-likelihood correction.
#'
#' If \code{H} is supplied as a function and the earlier fitted model itself
#' contains a previous \code{"latent"} object, the same function is passed
#' recursively and evaluated independently for each nested fit. This allows, for
#' example, robust information matrices to be used at every estimation stage.
#'
#' At most one nested \code{"latent"} object is currently supported at each
#' estimation level.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{\code{vcov}}{Variance-covariance matrix for the selected parameters.}
#'   \item{\code{se}}{Standard errors for the selected parameters.}
#'   \item{\code{jacob}}{Jacobian matrix used to propagate covariance through
#'     parameter transformations.}
#'   \item{\code{H}}{Information matrix used in the final covariance
#'     calculation.}
#'   \item{\code{B}}{For sequential models, the additional structural
#'     uncertainty component \eqn{C^\top A C}. For ordinary one-step models,
#'     an empty matrix.}
#'   \item{\code{newH}}{For sequential models, the joint corrected precision
#'     matrix used to obtain the final covariance matrix.}
#' }
#'
#' @seealso
#' \code{\link{hessian.latent}}, \code{\link{jacobian.latent}},
#' \code{\link{constraints_derivs.latent}}, \code{\link{se.llca}}
#'
#' @references
#' Bakk, Z., and Kuha, J. Two-step estimation of models between latent classes
#' and external variables.
#'
#' @method vcov latent
#' @export
vcov.latent <- function(fit, parameters = NULL, H = NULL) {

  #### Auxiliary functions ####

  information_matrix <- function(fit, H = NULL) {

    if(is.null(H)) {

      result <- hessian(fit)

    } else if(is.function(H)) {

      result <- H(fit)

    } else {

      result <- H

    }

    #### Result ####

    return(result)

  }

  #### Parameters ####

  if(is.null(parameters)) {

    parameters <- fit@modelInfo$trans[names(fit@modelInfo$param)]

  } else if(!any(unlist(parameters) %in%
                 fit@modelInfo$transparameters_labels)) {

    stop("Unknown parameters.")

  }

  #### Identify nested latent model ####

  model <- fit@modelInfo$control_optimizer$model
  latent_idx <- integer(0L)

  if(inherits(model, "latent")) {

    model <- list(model)
    latent_idx <- 1L

  } else if(is.list(model)) {

    latent_idx <- which(vapply(model,
                               FUN = \(x) inherits(x, "latent"),
                               FUN.VALUE = logical(1L)))

  }

  #### Ordinary model ####

  # If no nested latent model is found in fit, get directly the vcov:
  if(length(latent_idx) == 0L) {

    fit@modelInfo$control_optimizer$parameters[[1]] <-
      fit@Optim$parameters

    fit@modelInfo$control_optimizer$transparameters[[1]] <-
      fit@Optim$transparameters

    P <- information_matrix(fit, H)

    fit@modelInfo$control_optimizer$idx_transforms <-
      trans_depends(fit@modelInfo, parameters)

    VCOV <- get_vcov(fit@modelInfo$control_manifold,
                     fit@modelInfo$control_transform,
                     fit@modelInfo$control_estimator,
                     fit@modelInfo$control_optimizer,
                     H = P)

    rownames(VCOV$vcov) <- colnames(VCOV$vcov) <- names(VCOV$se) <-
      fit@modelInfo$transparameters_labels

    selected_parameters <- unique(unlist(parameters))
    selected_idx <- match(selected_parameters,
                          fit@modelInfo$transparameters_labels)

    VCOV$vcov <- VCOV$vcov[selected_idx, selected_idx, drop = FALSE]
    VCOV$se <- VCOV$se[selected_idx]
    VCOV$jacob <- VCOV$jacob[selected_idx, selected_idx, drop = FALSE]
    VCOV$H <- P
    VCOV$B <- matrix(numeric(0L), nrow = 0L, ncol = 0L)

    #### Result ####

    return(VCOV)

  }

  #### Structural after measurement model ####

  if(length(latent_idx) > 1L) {

    stop("Only one nested latent model at each model level is currently supported.")

  }

  latent_idx <- latent_idx[1L]

  # fit0 is the measurement model:
  fit0 <- model[[latent_idx]]

  # Retain every other model specification:
  remaining_model <- model[-latent_idx]

  if(length(remaining_model) == 0L) {
    remaining_model <- NULL
  }

  #### Measurement-model uncertainty ####

  # Only function-valued information specifications are propagated recursively.
  # A supplied matrix belongs to the current model and must not be reused for a
  # nested model with a different parameterization.
  nested_H <- if(is.function(H)) H else NULL

  # Recursive call:
  #
  # If fit0 itself contains another latent object, vcov.latent() automatically
  # applies the same correction to fit0 first. When H is a function, it is
  # evaluated independently for fit0 and every earlier nested model.
  VCOV0 <- vcov(fit0, parameters = fit0@modelInfo$parameters_labels,
                 H = nested_H)

  #### Structural-model covariance ####

  # H2 is always the ordinary Hessian because it determines how uncertainty
  # from the previous estimation stage propagates into the structural estimates.
  H2 <- hessian(fit)

  # P2 determines the structural covariance when the measurement parameters are
  # treated as known. For standard inference P2 = H2; for robust inference P2
  # can be, for example, robust_LG(fit).
  P2 <- information_matrix(fit, H)

  #### Full unrestricted model ####

  args <- fit@dataList$args

  # Remove the parameter constraints from previous fitted objects:
  args$model <- remaining_model  # Paste only the model dependencies
  args$do.fit <- FALSE           # Don't fit, just get the model
  args$adjustment <- "none"      # One-step fit
  fit_full <- do.call(lca, args) # Unrestricted model

  # Paste the parameter estimates of the structural model into the
  # corresponding positions of the unrestricted model:
  fit_full@Optim$parameters <-
    fit@Optim$transparameters[fit_full@modelInfo$parameters_labels]

  fit_full@Optim$transparameters <-
    fit@Optim$transparameters[fit_full@modelInfo$transparameters_labels]

  #### Full-model Hessian ####

  # Hessian matrix of the unrestricted model:
  H_full <- hessian(fit_full)

  # Parameters estimated in the structural model:
  structural_pars <- fit_full@modelInfo$parameters_labels %in%
    fit@modelInfo$parameters_labels

  # Measurement parameters fixed during structural estimation:
  measurement_pars <- !structural_pars

  if(!any(structural_pars) || !any(measurement_pars)) {
    stop("The model does not contain both structural and measurement parameters.")
  }

  #### Measurement-structural second-order derivatives ####

  C <- H_full[measurement_pars, structural_pars, drop = FALSE]

  measurement_pars_names <-
    fit_full@modelInfo$parameters_labels[measurement_pars]

  structural_pars_names <-
    fit_full@modelInfo$parameters_labels[structural_pars]

  A <- as.matrix(VCOV0$vcov[measurement_pars_names,
                            measurement_pars_names, drop = FALSE])

  # Additional structural uncertainty propagated from the measurement model:
  B <- t(C) %*% A %*% C

  # Remove small numerical asymmetries:
  B <- (B + t(B))/2

  #### Corrected information matrix ####

  # The joint covariance of measurement and structural parameters is:
  #
  #   [ A                    -A C H2^-1                         ]
  #   [ -H2^-1 C' A          P2^-1 + H2^-1 C' A C H2^-1       ]
  #
  # Its inverse is:
  #
  #   [ A^-1 + C H2^-1 P2 H2^-1 C'     C H2^-1 P2 ]
  #   [ P2 H2^-1 C'                       P2         ]
  #
  # When P2 = H2 this reduces to the ordinary sequential correction:
  #
  #   [ A^-1 + C H2^-1 C'     C  ]
  #   [ C'                    H2 ]

  H2_inv <- solve(H2)

  newH11 <- solve(A) +
    C %*% H2_inv %*% P2 %*% H2_inv %*% t(C)

  newH12 <- C %*% H2_inv %*% P2

  newH <- rbind(
    cbind(newH11, newH12),
    cbind(t(newH12), P2)
  )

  newH_names <- c(measurement_pars_names, structural_pars_names)
  rownames(newH) <- colnames(newH) <- newH_names

  newH <- newH[fit_full@modelInfo$parameters_labels,
               fit_full@modelInfo$parameters_labels,
               drop = FALSE]

  #### Variance-covariance matrix ####

  result <- vcov(fit_full, parameters = parameters, H = newH)

  result$B <- B
  colnames(result$B) <- rownames(result$B) <- structural_pars_names

  result$newH <- newH

  #### Result ####

  return(result)

}
