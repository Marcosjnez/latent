# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 18/08/2026
#'
#' Base Class for Fitted Latent Models
#'
#' The \code{"latent"} class provides the common representation used by fitted
#' latent-variable models that share the generic derivative and covariance
#' infrastructure of \pkg{latent}.
#'
#' The model-specific subclasses inherit the slots of this class. Multistep models
#' propagate uncertainty from earlier fits through deterministic steps. Multistage
#' models also account for sampling variation from the new estimation stage.
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
#' @slot extra List containing fitted models from previous estimation steps or
#'   other model-specific information.
#'
#' @seealso
#' \code{\link{hessian}}, \code{\link{jacobian}},
#' \code{\link{constraints_derivs}}, \code{\link[stats]{vcov}}
#'
#' @name latent-class
#' @aliases lcfa-class llca-class llhmm-class latscores-class multistep-class multistep_lcfa-class multistage-class multistage_lcfa-class multistage_llca-class
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
         contains = "latent")

setClass("lcfa",
         contains = "latent")

# A multistage estimator retains conditional sampling variability in every new
# estimation stage, as in structural-after-measurement LCA.
setClass("multistage",
         contains = "latent")

setClass("multistage_llca",
         contains = c("multistage", "llca"))

setClass("multistage_lcfa",
         contains = c("multistage", "lcfa"))

# A multistep estimator is a deterministic estimating step conditional on the
# outputs of previous fitted models. Its sampling covariance is obtained only
# by propagating uncertainty from those previous steps.
setClass("multistep",
         contains = "latent")

setClass("multistep_lcfa",
         contains = c("multistep", "lcfa"))
