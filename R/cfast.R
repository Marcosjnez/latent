# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 27/10/2025
#'
#' @title
#' Wrapper to the lcfa function to fit a Confirmatory Factor Analysis (CFA) model with lavaan syntax.
#'
#' @usage
#'
#' lcfa(data, model = NULL, cor = "pearson",
#' estimator = "ml", group = NULL,
#' sample.cov = NULL, nobs = NULL, W = NULL,
#' positive = FALSE, penalties = TRUE,
#' missing = "pairwise.complete.obs",
#' std.lv = FALSE, do.fit = TRUE, control = NULL)
#'
#' @param data data frame or matrix.
#' @param model lavaan's model syntax.
#' @param cor Correlation types: "pearson" and "poly". Defaults to "pearson".
#' @param estimator Available estimators: "ml", "uls", and "dwls". Defaults to "ml".
#' @param group .
#' @param sample.cov Covariance matrix between the items. Defaults to NULL.
#' @param nobs Number of observations. Defaults to NULL.
#' @param W Custom weight matrix for "dwls". Defaults to NULL.
#' @param positive Force at least positive-semidefinite solutions. Defaults to FALSE
#' @param penalties list of penalty terms for the parameters.
#' @param missing Method to handle missing data.
#' @param std.lv Provide the parameters of the standardized model.
#' @param do.fit TRUE to fit the model and FALSE to return only the model setup. Defaults to TRUE.
#' @param control List of control parameters for the optimization algorithm. See 'details' for more information.
#'
#' @details \code{cfast} estimates confirmatory factor models.
#'
#' @return List with the following objects:
#' \item{version}{Version number of 'latent' when the model was estimated.}
#' \item{call}{Code used to estimate the model.}
#' \item{ModelInfo}{Model information.}
#' \item{Optim}{Output of the optimizer.}
#' \item{parameters}{Structure with all model parameters.}
#' \item{transparameters}{Structure with all transformed model parameters.}
#' \item{loglik}{Logarithm likelihood of the model.}
#' \item{penalized_loglik}{Logarithm likelihood + logarithm priors of the model.}
#'
#' @examples
#'
#' \dontrun{
#' # The famous Holzinger and Swineford (1939) example
#' HS.model <- ' visual  =~ x1 + x2 + x3
#'               textual =~ x4 + x5 + x6
#'               speed   =~ x7 + x8 + x9 '
#'
#' fit <- lcfa(model = HS.model, data = HolzingerSwineford1939)
#' summary(fit, digits = 3L)
#'}
#'
#' @export
cfast <- function(data, model = NULL, cor = "pearson",
                  estimator = "ml", ...) {

  lcfa(data, model = model, cor = cor,
       estimator = estimator, ...)

}

