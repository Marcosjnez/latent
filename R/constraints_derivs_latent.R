# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 05/09/2026
#'
#' Constraint Derivatives for Latent Models
#'
#' Compute first and second derivatives of manifold- and
#' transformation-induced constraints for selected parameters of a fitted
#' latent variable model.
#'
#' @description
#' Some manifolds and parameter transformations imply constraints. For example,
#' parameters on a unit manifold have unit norm and probabilities produced by a
#' softmax transformation sum to one.
#' \code{constraints_derivs.latent()} extracts the first and second derivatives
#' associated with such constraints for selected transformed parameters.
#'
#' @param fit A fitted object inheriting from class \code{"latent"}.
#' @param parameters Optional parameter specification identifying the transformed
#'   parameters for which constraint derivatives should be returned. Parameter
#'   labels must occur in \code{fit@modelInfo$transparameters_labels}. If
#'   \code{NULL}, the parameter blocks corresponding to the fitted model
#'   parameters are used.
#'
#' @return
#' A list with two sparse matrices:
#' \describe{
#'   \item{\code{dconstr}}{Rows correspond to selected transformed parameters
#'     and columns correspond to constraints.}
#'   \item{\code{d2constr}}{A block-diagonal matrix containing one square
#'     constraint Hessian per constraint, in the same order as the columns of
#'     \code{dconstr}.}
#' }
#'
#' @details
#' Only transformations on which the requested parameters depend are evaluated.
#' Manifold constraints are always evaluated. Transformations and manifolds that
#' do not impose an explicit constraint do not contribute a constraint column.
#'
#' If \eqn{p} parameters and \eqn{m} constraints are selected,
#' \code{d2constr} has dimension \eqn{pm} by \eqn{pm}. Its \eqn{j}-th
#' \eqn{p} by \eqn{p} diagonal block is the Hessian of constraint \eqn{j}.
#'
#' @seealso
#' \code{\link{constraints_derivs}}, \code{\link{jacobian.latent}},
#' \code{\link{vcov.latent}}
#'
#' @method constraints_derivs latent
#' @export
constraints_derivs.latent <- function(fit, parameters = NULL) {

  if(is.null(parameters)) {
    parameters <- fit@modelInfo$parameters_labels
  } else {

    selected_parameters <- unique(unlist(parameters))

    if(!all(selected_parameters %in%
            fit@modelInfo$transparameters_labels)) {
      stop("Unknown parameters.")
    }

  }

  # fit@modelInfo$control_optimizer$parameters[[1]] <- fit@Optim$parameters
  # fit@modelInfo$control_optimizer$transparameters[[1]] <- fit@Optim$transparameters

  fit@modelInfo$control_optimizer$idx_transforms <-
    trans_depends(fit@modelInfo, parameters)

  # derivatives <- get_dconstr(fit@modelInfo$control_manifold,
  #                            fit@modelInfo$control_transform,
  #                            fit@modelInfo$control_estimator,
  #                            fit@modelInfo$control_optimizer)
  derivatives <- get_dconstr(fit)

  dconstr <- derivatives$dconstr
  d2constr <- derivatives$d2constr

  rownames(dconstr) <- fit@modelInfo$transparameters_labels

  # Not every transformation contributes with just one column of constraints, so
  # don't put names to columns (do not run this code):
  # colnames(dconstr) <- sapply(fit@modelInfo$control_transform[idx_transforms+1L],
  #                             FUN = \(x) x$transform)

  selected_parameters <- unique(unlist(parameters))
  selected_idx <- match(selected_parameters,
                        fit@modelInfo$transparameters_labels)
  dconstr <- dconstr[selected_idx, , drop = FALSE]

  ntrans <- length(fit@modelInfo$transparameters_labels)
  nconstraints <- ncol(dconstr)

  if(nconstraints > 0L) {

    expected_dimension <- ntrans*nconstraints

    if(nrow(d2constr) != expected_dimension ||
       ncol(d2constr) != expected_dimension) {
      stop("The first- and second-constraint derivative objects are not aligned.")
    }

    selected_d2_idx <- unlist(
      lapply(seq_len(nconstraints), FUN = \(j) {
        (j-1L)*ntrans+selected_idx
      }),
      use.names = FALSE
    )

    d2constr <- d2constr[selected_d2_idx,
                         selected_d2_idx,
                         drop = FALSE]

  }

  #### Result ####

  result <- list(dconstr = dconstr,
                 d2constr = d2constr)

  return(result)

}
