# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 12/08/2026
#'
#' Standard Errors for Latent Class Models
#'
#' Compute standard errors, covariance matrices, Hessians, and transformation
#' Jacobians for fitted latent class models.
#'
#' @description
#' \code{se.llca()} computes standard or robust standard errors for selected
#' parameters of a fitted latent class model.
#'
#' Covariance matrices are obtained through \code{\link{vcov.latent}}, which
#' also propagates uncertainty from earlier estimation stages when the fitted
#' model contains a nested \code{"latent"} object.
#'
#' @param fit A fitted object of class \code{"llca"}.
#' @param type Character string specifying the covariance estimator. Available
#'   options are \code{"standard"} for Hessian-based standard errors and
#'   \code{"robust"} for a LatentGold-style sandwich estimator.
#' @param parameters Optional parameter specification identifying the parameters
#'   for which standard errors should be returned. Labels must occur in
#'   \code{fit@modelInfo$transparameters_labels}. If \code{NULL}, the parameter
#'   blocks corresponding to the fitted model parameters are used.
#' @param digits Non-negative integer specifying the number of decimal places
#'   used in the formatted parameter tables. If \code{NULL}, values are not
#'   rounded.
#' @param ... Additional arguments passed to other methods.
#'
#' @details
#' With \code{type = "standard"}, covariance matrices are based on the inverse
#' Hessian of the fitted objective function.
#'
#' With \code{type = "robust"}, the empirical covariance of case or
#' response-pattern score contributions is combined with the Hessian to produce
#' a sandwich covariance estimator.
#'
#' When a previous fitted \code{"latent"} object is stored in the model
#' specification of \code{fit}, \code{\link{vcov.latent}} propagates uncertainty
#' from that earlier estimation step. Nested sequential models are processed
#' recursively.
#'
#' Standard errors for transformed parameters are obtained using the delta
#' method and the Jacobians of the corresponding parameter transformations.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{\code{table}}{Parameter estimates and standard errors arranged
#'     according to the model parameter blocks.}
#'   \item{\code{table_se}}{Standard errors in the same parameter structure as
#'     the fitted model.}
#'   \item{\code{se}}{Named numeric vector of standard errors.}
#'   \item{\code{vcov}}{Variance-covariance matrix of the selected parameters.}
#'   \item{\code{B}}{Additional uncertainty component for sequential models, or
#'     an empty matrix for ordinary one-step models.}
#'   \item{\code{H}}{Hessian or information matrix used in the calculation.}
#'   \item{\code{newH}}{Corrected joint precision matrix when applicable.}
#'   \item{\code{jacob}}{Jacobian matrix used for covariance propagation to
#'     transformed parameters.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- lca(data = empathy, nclasses = 3L,
#'            gaussian = c("ec1", "ec2", "ec3"))
#'
#' se(fit)
#' se(fit, type = "robust")
#'
#' # Standard errors for selected transformed parameters
#' se(fit, parameters = fit@modelInfo$trans[c("class", "beta")])
#' }
#'
#' @seealso
#' \code{\link{vcov.latent}}, \code{\link{hessian.latent}},
#' \code{\link{jacobian.latent}}, \code{\link{ci}}
#'
#' @method se llca
#' @export
se.llca <- function(fit, type = "standard", parameters = NULL, digits = 4L,
                    ...) {

  if(is.null(parameters)) {
    parameters <- fit@modelInfo$trans[names(fit@modelInfo$param)]
  } else if(!any(unlist(parameters) %in%
                 fit@modelInfo$transparameters_labels)) {
    stop("Unknown parameters.")
  }

  if(!inherits(fit, "llca")) {
    stop("fit must inherit from class 'llca'.")
  }

  if(length(fit@Optim$parameters) == 0L) {
    stop("The llca object has not been fitted.")
  }

  type <- match.arg(tolower(type), c("standard", "robust"))

  if(!is.null(digits) &&
     (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits < 0L || digits != as.integer(digits))) {
    stop("digits must be NULL or a non-negative integer.")
  }

  if(type == "standard") {
    H <- NULL
  } else if(type == "robust") {
    H <- robust_LG(fit = fit)
  } else {
    stop("Unknown type method. Available types: 'standard' and 'robust' ")
  }

  SE <- vcov(fit, parameters = parameters, H = H)

  est <- fill_in(parameters, fit@Optim$transparameters, miss = NA)
  table_se <- fill_in(parameters, SE$se, miss = NA)
  table <- combine_est_se(est, table_se, digits = digits)

  result <- list(table = table, table_se = table_se, se = c(SE$se),
                 vcov = SE$vcov, B = SE$B, H = SE$H, newH = SE$newH,
                 jacob = SE$jacob)

  #### Result ####

  return(result)

}

#' Standard Errors for Structural-After-Measurement Models
#'
#' Compute standard errors for the structural component of a multi-step latent
#' class model.
#'
#' @param fit An object of class \code{"llca_sam"} containing fitted
#'   \code{measurement} and \code{structural} components.
#' @param type Character string specifying the covariance estimator. See
#'   \code{\link{se.llca}}.
#' @param parameters Optional parameter specification passed to
#'   \code{\link{se.llca}}.
#' @param digits Non-negative integer specifying the number of decimal places
#'   used in formatted parameter tables. Use \code{NULL} to avoid rounding.
#' @param ... Additional arguments passed to \code{\link{se.llca}}.
#'
#' @details
#' The method delegates inference to the structural \code{"llca"} model.
#' When that structural model stores the fitted measurement model among its
#' model specifications, \code{\link{vcov.latent}} automatically propagates
#' measurement-model uncertainty to the structural estimates.
#'
#' @return
#' The result of \code{se(fit$structural, ...)}.
#'
#' @seealso
#' \code{\link{se.llca}}, \code{\link{vcov.latent}}, \code{\link{lca}}
#'
#' @method se llca_sam
#' @export
se.llca_sam <- function(fit, type = "standard", parameters = NULL, digits = 4L,
                        ...) {

  return(se.llca(fit$structural, type = type, parameters = parameters,
                 digits = digits, ...))

}

#' Standard errors for a collection of latent class models
#'
#' Applies \code{se()} to every fitted model in an \code{"llcalist"} object.
#'
#' @param model An object of class \code{"llcalist"} containing fitted
#'   \code{"llca"} objects.
#' @param type Character string indicating the standard-error estimator. See
#'   \code{se.llca()}.
#' @param digits Non-negative integer indicating the number of decimal places used
#'   in the formatted tables, or \code{NULL} to avoid rounding.
#' @param ... Additional arguments passed to \code{se.llca()}.
#'
#' @details
#' Existing names are preserved. Consequently, class-enumeration results retain
#' names such as \code{"nclasses=2"}, whereas multiple-step models retain names
#' such as \code{"measurement"} and \code{"structural"}. Unnamed elements are
#' labelled according to their number of latent classes.
#'
#' @return A named list with one standard-error result per fitted model and class
#'   \code{"se.llcalist"}.
#'
#' @examples
#' \dontrun{
#' fits <- lca(data = empathy, nclasses = 2:4,
#'             gaussian = c("ec1", "ec2", "ec3"))
#' se(fits)
#' }
#'
#' @method se llcalist
#' @export
se.llcalist <- function(model, type = "standard", parameters = NULL,
                        digits = 4L, ...) {

  if(!inherits(model, "llcalist")) {
    stop("model must inherit from class 'llcalist'.")
  }

  if(length(model) == 0L) {
    stop("model must contain at least one fitted llca object.")
  }

  nmodels <- length(model)
  result <- vector("list", length = nmodels)

  for(i in seq_len(nmodels)) {
    result[[i]] <- se(model[[i]], type = type, parameters = parameters,
                      digits = digits, ...)
  }

  result_names <- names(model)
  if(is.null(result_names)) {
    result_names <- rep("", nmodels)
  }

  unnamed <- is.na(result_names) | result_names == ""
  if(any(unnamed)) {
    nclasses <- vapply(model[unnamed], FUN = function(x) {
      ncol(x@modelInfo$trans$class)
    }, FUN.VALUE = integer(1L))
    result_names[unnamed] <- paste0("nclasses=", nclasses)
  }

  names(result) <- make.unique(result_names)
  class(result) <- "se.llcalist"

  #### Result ####

  return(result)

}

#' LatentGold-Style Robust Information Matrix
#'
#' Construct an information matrix corresponding to a LatentGold-style
#' sandwich covariance estimator for a fitted latent class model.
#'
#' @param fit A fitted object of class \code{"llca"}.
#'
#' @details
#' Let \eqn{H} denote the Hessian and \eqn{B} the empirical covariance matrix
#' of score contributions. The robust covariance matrix is
#' \deqn{H^{-1} B H^{-1}.}
#'
#' Because \code{\link{vcov.latent}} expects an information matrix that is
#' subsequently inverted, this function returns the equivalent matrix
#' \deqn{H B^{-1} H.}
#'
#' @return
#' A symmetric numeric matrix representing the robust information matrix.
#'
#' @keywords internal
robust_LG <- function(fit) {

  if(fit@dataList$nobs <= 1L) {
    stop("Robust standard errors require more than one observation.")
  }

  #### Compute the Hessian ####

  fit@modelInfo$control_optimizer$parameters[[1]] <-
    fit@Optim$parameters

  fit@modelInfo$control_optimizer$transparameters[[1]] <-
    fit@Optim$transparameters

  H <- hessian(fit)

  #### Collect the gradient by response pattern ####

  transparameters_labels <- fit@modelInfo$transparameters_labels
  pattern_weights <- fit@dataList$pattern_weights
  npatterns <- fit@dataList$npatterns
  nparam <- fit@modelInfo$nparam
  nobs <- fit@dataList$nobs
  nclasses <- fit@dataList$nclasses

  # Remove the lca estimator but keep everything else:
  control_estimator <- fit@modelInfo$control_estimator[-1L]
  K <- length(control_estimator)

  pattern_struct <- vector("list", length = npatterns)
  for(s in seq_len(npatterns)) {

    pattern_struct[[s]] <- list(estimator = "lca",
                                  parameters = list(fit@modelInfo$trans$class[s, ],
                                                    fit@modelInfo$trans$loglik[s, ]),
                                  extra = list(S = 1L,
                                               I = nclasses,
                                               weights = pattern_weights[s]))
  }

  pattern_estimators <- create_estimators(estimators = pattern_struct,
                                          structures = fit@modelInfo$trans)
  control_estimator <- c(control_estimator, pattern_estimators)

  #### Compute the B matrix ####

  control_manifold <- fit@modelInfo$control_manifold
  control_transform <- fit@modelInfo$control_transform
  control_optimizer <- fit@modelInfo$control_optimizer
  B <- matrix(0, nrow = nparam, ncol = nparam)
  for(s in seq_len(npatterns)) {
    idx <- c(seq_len(K), K+s)
    computations <- get_grad(control_manifold = control_manifold,
                             control_transform = control_transform,
                             control_estimator = control_estimator[idx],
                             control_optimizer = control_optimizer)
    gradient <- computations$g/pattern_weights[s]
    B <- B + pattern_weights[s] * gradient %*% t(gradient)
  }

  B <- B*nobs/(nobs-1L)

  #### Update the Hessian ####

  newH <- H %*% solve(B) %*% H
  newH <- (newH + t(newH))/2

  #### Return ####

  return(newH)

}
