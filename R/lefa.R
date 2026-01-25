# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 26/11/2025
#'
#' @title
#' Fit an Exploratory Factor Analysis (EFA) model.
#'
#' @usage
#'
#' lefa(data, nfactors = 1L, model = NULL, estimator = "ml",
#' ordered = FALSE, group = NULL,
#' sample.cov = NULL, nobs = NULL,
#' positive = FALSE, penalties = TRUE,
#' missing = "pairwise.complete.obs",
#' std.lv = FALSE, do.fit = TRUE, mimic = 'latent',
#' control = NULL, ...)
#'
#' @param data data frame or matrix.
#' @param nfactors integer. Number of latent variables.
#' @param model lavaan's model syntax.
#' @param estimator Available estimators: "ml", "uls", and "dwls". Defaults to "ml".
#' @param ordered Logical. Defaults to TRUE.
#' @param group String. Name of the variable that splits the data in different groups.
#' @param sample.cov Covariance matrix between the items. Defaults to NULL.
#' @param nobs Number of observations. Defaults to NULL.
#' @param positive Force a positive-definite solution. Defaults to FALSE.
#' @param penalties list of penalty terms for the parameters.
#' @param missing Method to handle missing data.
#' @param std.lv Provide the parameters of the standardized model.
#' @param do.fit TRUE to fit the model and FALSE to return only the model setup. Defaults to TRUE.
#' @param mimic String. Choose the output you want to obtain. Defaults to 'latent'.
#' @param control List of control parameters for the optimization algorithm. See 'details' for more information.
#' @param ... Additional lavaan arguments. See ?lavaan for more information.
#'
#' @details \code{lefa} estimates confirmatory factor models.
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
#'
#' fit <- lefa(data = HolzingerSwineford1939, nfactors = 3L)
#' summary(fit, digits = 3L)
#'}
#'
#' @export
lefa <- function(data, nfactors = 1L, estimator = "ml",
                 projection = "oblq", rotation = "oblimin",
                 model = NULL, ordered = FALSE, group = NULL,
                 sample.cov = NULL, nobs = NULL,
                 positive = FALSE, penalties = TRUE,
                 missing = "pairwise.complete.obs",
                 std.lv = TRUE, do.fit = TRUE,
                 mimic = "latent", control = NULL,
                 ...) {


  #### Fit the FA model ####

  if(is.null(model)) {
    model <- make_lowerdiag_lavaan(data = data, nfactors = nfactors)
  } else {
    # CHECK DATA
    # THIS USES ALL THE COLUMNS IN DATA
  }

  fit1 <- lcfa(data = data, model = model, estimator = estimator,
               ordered = ordered, group = group,
               sample.cov = sample.cov, nobs = nobs,
               positive = positive, penalties = penalties,
               missing = missing,
               std.lv = std.lv, do.fit = do.fit,
               mimic = mimic, control = control,
               orthogonal = TRUE, ...)

  lambda <- latInspect(fit1, what = "lambda")

  if(nfactors > 1L && do.fit) {

    fit2 <- lrotate(lambda, projection = projection, rotation = rotation,
                     group = group, positive = positive, penalties = penalties,
                     do.fit = do.fit, control = list(opt = "newton"), ...)

  } else {

    return(fit1)
  }

  result <- list(efa = fit1,
                 rotation = fit2)

  return(result)

}

