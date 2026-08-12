# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 12/08/2026
#'
#' Base Class for Fitted Latent Models
#'
#' The \code{"latent"} class provides the common representation used by fitted
#' latent-variable models that share the generic derivative and covariance
#' infrastructure of \pkg{latent}.
#'
#' @slot version Character string containing the package version used to fit
#'   the model.
#' @slot call Original model-fitting call.
#' @slot timing Numeric vector containing optimization timing information.
#' @slot dataList List containing processed data and model-specific data objects.
#' @slot modelInfo List containing parameter structures, transformations,
#'   estimators, manifolds, and optimizer controls.
#' @slot Optim List containing raw optimization results.
#' @slot parameters List containing fitted model parameters.
#' @slot transformed_pars List containing transformed parameter estimates.
#' @slot extra List reserved for additional model-specific information.
#'
#' @seealso
#' \code{\link{hessian}}, \code{\link{jacobian}},
#' \code{\link{constraints_derivs}}, \code{\link{vcov}}
#'
#' @name latent-class
#' @rdname latent-class
setClass("latent",
         slots = c(
           version          = "character",
           call             = "call",
           timing           = "numeric",
           dataList         = "list",
           modelInfo        = "list",
           Optim            = "list",
           parameters       = "list",
           transformed_pars = "list",
           extra            = "list"
         )
)

setClass("llca",
         slots = c(
           version          = "character",
           call             = "call",
           timing           = "numeric",
           dataList         = "list",
           modelInfo        = "list",
           Optim            = "list",
           parameters       = "list",
           transformed_pars = "list",
           extra            = "list"
         ),
         contains = "latent"
)

setClass("lcfa",
         slots = c(
           version          = "character",
           call             = "call",
           timing           = "numeric",
           dataList         = "list",
           modelInfo        = "list",
           Optim            = "list",
           parameters       = "list",
           transformed_pars = "list",
           loglik           = "numeric",
           penalized_loglik = "numeric",
           loss             = "numeric",
           penalized_loss   = "numeric",
           extra            = "list"
         )
)
