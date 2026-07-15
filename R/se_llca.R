# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 14/07/2026
#'
#' Standard errors for latent class models
#'
#' Computes standard errors and variance-covariance matrices for fitted latent
#' class models.
#'
#' @param fit A fitted object of class \code{"llca"}.
#' @param type Character string indicating the standard-error estimator. Available
#'   options are \code{"standard"} for Hessian-based standard errors and
#'   \code{"robust"} for the LatentGold-style sandwich estimator.
#' @param digits Non-negative integer indicating the number of decimal places used
#'   in the formatted table. If \code{NULL}, the table is returned without
#'   rounding.
#' @param ... Additional arguments passed to other methods.
#'
#' @details
#' For a regular one-step model, \code{type = "standard"} obtains the covariance
#' matrix from the inverse Hessian. With \code{type = "robust"}, the covariance
#' matrix is computed from a sandwich estimator using the score contribution of
#' each observed response pattern.
#'
#' When the fitted model contains a previous \code{"llca"} object in
#' \code{fit@modelInfo$control_optimizer$model}, the standard errors are adjusted
#' for two-step estimation through \code{se_twostep()}.
#'
#' The \code{digits} argument affects only the formatted \code{table}. The numeric
#' standard errors and covariance matrices are returned without rounding.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{table}}{A list of formatted parameter tables containing estimates
#'   and standard errors.}
#'   \item{\code{table_se}}{A list containing the standard errors arranged in the
#'   same parameter structure as the fitted model.}
#'   \item{\code{se}}{A named numeric vector of standard errors.}
#'   \item{\code{vcov}}{The variance-covariance matrix of the model parameters.}
#'   \item{\code{B}}{The empirical score covariance matrix. It is an empty matrix
#'   for ordinary Hessian-based standard errors and may be \code{NULL} when it is
#'   not applicable.}
#'   \item{\code{H}}{The Hessian matrix, when available.}
#'   \item{\code{newH}}{The adjusted Hessian used by the robust estimator, when
#'   available.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- lca(data = empathy, nclasses = 3L,
#'            gaussian = c("ec1", "ec2", "ec3"))
#'
#' se(fit)
#' se(fit, type = "robust", digits = 4L)
#' }
#'
#' @seealso \code{ci()}, \code{vcov()}
#'
#' @method se llca
#' @export
se.llca <- function(fit, type = "standard", digits = 3L, ...) {

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

  original_model <- fit@modelInfo$control_optimizer$model

  if(inherits(original_model, "llca")) {

    SE <- se_twostep(fit2 = fit, type = type)

  } else if(type == "standard") {

    SE <- standard_se(fit = fit)
    SE$B <- matrix(numeric(0L), nrow = 0L, ncol = 0L)

  } else {

    SE <- robust_se_LG(fit = fit)

  }

  est <- fill_in(fit@modelInfo$trans[names(fit@modelInfo$param)],
                 c(fit@Optim$parameters, fit@Optim$transparameters), miss = NA)
  table_se <- fill_in(fit@modelInfo$trans[names(fit@modelInfo$param)],
                      SE$se, miss = NA)
  table <- combine_est_se(est, table_se, digits = digits)

  result <- list(table = table, table_se = table_se, se = c(SE$se),
                 vcov = SE$vcov, B = SE$B, H = SE$H, newH = SE$newH)

  #### Result ####

  return(result)

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
se.llcalist <- function(model, type = "standard", digits = 3L, ...) {

  if(!inherits(model, "llcalist")) {
    stop("model must inherit from class 'llcalist'.")
  }

  if(length(model) == 0L) {
    stop("model must contain at least one fitted llca object.")
  }

  if(!all(vapply(model, inherits, logical(1L), what = "llca"))) {
    stop("All elements of model must inherit from class 'llca'.")
  }

  nmodels <- length(model)
  result <- vector("list", length = nmodels)

  for(i in seq_len(nmodels)) {
    result[[i]] <- se(model[[i]], type = type, digits = digits, ...)
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

#' LatentGold-style robust standard errors
#'
#' Computes a sandwich covariance estimator from the Hessian and the score
#' contribution of each observed response pattern.
#'
#' @param fit A fitted object of class \code{"llca"}.
#'
#' @return A list containing the standard errors, covariance matrix, empirical
#'   score covariance matrix \code{B}, Hessian \code{H}, and adjusted Hessian
#'   \code{newH}.
#'
#' @keywords internal
robust_se_LG <- function(fit) {

  if(fit@dataList$nobs <= 1L) {
    stop("Robust standard errors require more than one observation.")
  }

  control_manifold <- fit@modelInfo$control_manifold
  control_transform <- fit@modelInfo$control_transform
  control_estimator <- fit@modelInfo$control_estimator
  control_optimizer <- fit@modelInfo$control_optimizer
  control_optimizer$parameters[[1L]] <- fit@Optim$parameters
  control_optimizer$transparameters[[1L]] <- fit@Optim$transparameters

  cores <- control_optimizer$cores
  if(is.null(cores) || !is.finite(cores) || cores < 1L) {
    cores <- 1L
  }

  H <- get_hess(control_manifold = control_manifold,
                control_transform = control_transform,
                control_estimator = control_estimator,
                control_optimizer = control_optimizer,
                cores = as.integer(cores))$h

  #### Collect the gradient by response pattern ####

  transparameters_labels <- fit@modelInfo$transparameters_labels
  trans <- fit@modelInfo$trans
  full_loglik <- trans$loglik
  pattern_weights <- fit@dataList$pattern_weights
  npatterns <- fit@dataList$npatterns
  nitems <- fit@dataList$nitems
  nparam <- fit@modelInfo$nparam
  nobs <- fit@dataList$nobs

  control_estimator <- control_estimator[-1L]
  nclasses <- ncol(fit@modelInfo$trans$class)
  K <- length(control_estimator)

  for(s in seq_len(npatterns)) {
    indices <- list(match(trans$class[s, ], transparameters_labels)-1L,
                    match(full_loglik[s, , ], transparameters_labels)-1L)

    control_estimator[[K+s]] <- list(estimator = "lca", indices = indices,
                                     S = 1L, J = nitems, I = nclasses,
                                     pattern_weights = pattern_weights[s])
  }

  B <- matrix(0, nrow = nparam, ncol = nparam)
  pattern_estimators <- control_estimator[K+seq_len(npatterns)]

  for(s in seq_len(npatterns)) {
    computations <- get_grad(control_manifold = control_manifold,
                             control_transform = control_transform,
                             control_estimator = pattern_estimators[s],
                             control_optimizer = control_optimizer)
    gradient <- computations$g/pattern_weights[s]
    B <- B + pattern_weights[s] * gradient %*% t(gradient)
  }

  B <- B*nobs/(nobs-1L)
  newH <- H %*% solve(B) %*% H

  result <- get_vcov(control_manifold = control_manifold,
                     control_transform = control_transform,
                     control_estimator = control_estimator,
                     control_optimizer = control_optimizer, H = newH)

  result$B <- B

  names(result$se) <- colnames(result$vcov) <- rownames(result$vcov) <-
    fit@modelInfo$control_optimizer$se_names

  colnames(newH) <- rownames(newH) <- colnames(result$B) <- rownames(result$B) <-
    fit@modelInfo$parameters_labels

  result$H <- H
  result$newH <- newH

  #### Result ####

  return(result)

}

#' Two-step standard-error adjustment
#'
#' Adjusts the covariance matrix of a structural model for uncertainty in the
#' measurement-model parameters estimated in a previous step.
#'
#' @param fit2 A fitted structural \code{"llca"} object whose optimizer control
#'   stores the fitted measurement model.
#' @param type Character string indicating whether standard or robust covariance
#'   matrices are used in the two steps.
#'
#' @return A list containing the combined covariance matrix, standard errors, and
#'   the correction matrix \code{B}.
#'
#' @keywords internal
se_twostep <- function(fit2, type = "standard") {

  fit1 <- fit2@modelInfo$control_optimizer$model

  if(!inherits(fit1, "llca")) {
    stop("The stored first-step model must inherit from class 'llca'.")
  }

  VCOV1 <- se.llca(fit1, type = type, digits = NULL)

  fit2@modelInfo$control_optimizer$model <- NULL
  VCOV2 <- se.llca(fit2, type = type, digits = NULL)

  args <- fit2@dataList$args
  args$model <- NULL
  args$do.fit <- FALSE
  args$adjustment <- "none"
  fit <- do.call(lca, args)

  control_manifold <- fit@modelInfo$control_manifold
  control_transform <- fit@modelInfo$control_transform
  control_estimator <- fit@modelInfo$control_estimator
  control_optimizer <- fit@modelInfo$control_optimizer

  parameters <- fit2@Optim$transparameters[fit@modelInfo$parameters_labels]
  transparameters <- fit2@Optim$transparameters[fit@modelInfo$transparameters_labels]
  control_optimizer$parameters[[1L]] <- parameters
  control_optimizer$transparameters[[1L]] <- transparameters

  x <- get_hess(control_manifold, control_transform,
                control_estimator, control_optimizer)
  colnames(x$h) <- rownames(x$h) <- fit@modelInfo$parameters_labels

  model_pars <- fit@modelInfo$parameters_labels %in%
    fit2@modelInfo$parameters_labels
  nuisance_pars <- !model_pars

  if(!any(model_pars) || !any(nuisance_pars)) {
    stop("The two-step model does not contain both estimated and fixed parameters.")
  }

  df2_dparamdR <- x$h[nuisance_pars, model_pars, drop = FALSE]
  nuisance_names <- fit@modelInfo$parameters_labels[nuisance_pars]
  ACOV <- VCOV1$vcov[nuisance_names, nuisance_names, drop = FALSE]
  B <- t(df2_dparamdR) %*% ACOV %*% df2_dparamdR

  H_inv <- solve(VCOV2$H)
  VCOV2$vcov <- H_inv %*% B %*% H_inv
  VCOV2$se <- sqrt(diag(VCOV2$vcov))
  names(VCOV2$se) <- fit2@modelInfo$parameters_labels

  combined_vcov <- block_diag(list(VCOV1$vcov[nuisance_names, nuisance_names,
                                               drop = FALSE], VCOV2$vcov))
  result <- list(vcov = combined_vcov, se = sqrt(diag(combined_vcov)),
                 B = B, H = NULL, newH = NULL)

  #### Result ####

  return(result)

}
