# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 08/04/2026
#'
#' @title
#' Fit a Confirmatory Factor Analysis (CFA) model with lavaan syntax.
#'
#' @usage
#'
#' lcfa(data, model = NULL, estimator = "ml",
#' ordered = FALSE, group = NULL,
#' sample.cov = NULL, nobs = NULL,
#' positive = FALSE, penalties = TRUE,
#' missing = "pairwise.complete.obs",
#' std.lv = FALSE, do.fit = TRUE,
#' message = FALSE, mimic = 'latent',
#' control = NULL, ...)
#'
#' @param data data frame or matrix.
#' @param model lavaan's model syntax.
#' @param estimator Available estimators: "ml", "uls", and "dwls". Defaults to "ml".
#' @param ordered Logical. Defaults to TRUE.
#' @param group String. Name of the variable that splits the data in different groups.
#' @param sample.cov Covariance matrix between the items. Defaults to NULL.
#' @param nobs Number of observations. Defaults to NULL.
#' @param positive Force a positive-definite solution. Defaults to FALSE.
#' @param penalties list of penalty terms for the parameters.
#' @param missing Method to handle missing data.
#' @param std.lv Logical. Provide the parameters of the standardized model. Default is TRUE.
#' @param std.ov Logical. Standardize the observed variables before fitting. Default is FALSE.
#' @param acov String. "standard" or "robust". Default is "standard".
#' @param meanstructure Logical. Estimate the means of the variables. Default is FALSE.
#' @param do.fit TRUE to fit the model and FALSE to return only the model setup. Defaults to TRUE.
#' @param message Logical. Defaults to FALSE.
#' @param se Logical. Compute standard errors. Defaults to TRUE.
#' @param likelihood String. Use N (normal) or N-1 (wishart) in the denominator. Defaults to "normal" for ML and "wishart" otherwise.
#' @param mimic String. Choose the output you want to obtain. Defaults to 'latent'.
#' @param control List of control parameters for the optimization algorithm. See 'details' for more information.
#' @param ... Additional lavaan arguments. See ?lavaan for more information.
#'
#' @details \code{lcfa} estimates confirmatory factor models.
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
lcfa <- function(data, model = NULL, estimator = "ml",
                 ordered = FALSE, group = NULL,
                 sample.cov = NULL, nobs = NULL,
                 positive = FALSE, penalties = TRUE,
                 missing = "pairwise.complete.obs",
                 std.lv = FALSE, std.ov = FALSE,
                 acov = "standard",
                 meanstructure = NULL,
                 do.fit = TRUE, message = FALSE,
                 likelihood = NULL, se = TRUE,
                 mimic = "latent", control = NULL,
                 ...) {

  ## store original call
  mc  <- match.call()

  if(ordered) {
    cor <- "poly"
    control$deltaparam <- TRUE
    std.ov <- TRUE
  } else {
    cor <- "pearson"
  }

  estimator <- tolower(estimator)
  missing <- tolower(missing)

  if(is.null(meanstructure)) {
    if(std.ov) meanstructure <- FALSE
    if(missing == "fiml") {
      meanstructure <- TRUE
    } else {
      meanstructure <- FALSE
    }
  }

  if(meanstructure) {
    if(estimator == "ml" || estimator == "fml") estimator <- "means_fml"
    if(estimator == "uls") estimator <- "means_uls"
    if(estimator == "dwls") estimator <- "means_dwls"
  }

  control$std.ov <- std.ov
  control$positive <- positive
  control$penalties <- penalties
  control$estimator <- tolower(estimator)
  control$meanstructure <- meanstructure
  control <- lcfa_control(control)

  args <- as.list(match.call(expand.dots = TRUE))[-1]

  #### Create the datalist ####

  dataList <- create_cfa_datalist(
    data = data,
    model = model,
    cor = cor,
    estimator = estimator,
    ordered = ordered,
    group = group,
    sample.cov = sample.cov,
    nobs = nobs,
    positive = positive,
    penalties = penalties,
    missing = missing,
    std.lv = std.lv,
    std.ov = std.ov,
    acov = acov,
    message = message,
    likelihood = likelihood,
    meanstructure = meanstructure,
    args = args,
    ...
  )

  #### Create the model ####

  full_model <- create_cfa_model(dataList = dataList,
                                 model = model,
                                 control = control)
  list2env(full_model, envir = environment())

  #### Create the modelInfo ####

  modelInfo <- create_cfa_modelInfo(dataList = dataList,
                                    full_model = full_model,
                                    control = control)

  #### Fit the model ####

  if(!do.fit) {

    lcfa_list <- new("latent",
                     version            = as.character( packageVersion('latent') ),
                     call               = mc,
                     timing             = numeric(),
                     dataList           = dataList,
                     modelInfo          = modelInfo,
                     Optim              = list(),
                     parameters         = list(),
                     transformed_pars   = list(),
                     loglik             = numeric(),
                     penalized_loglik   = numeric(),
                     loss               = numeric(),
                     penalized_loss     = numeric()
    )

    return(lcfa_list)

  }

  if(message) {
    msg <- "Fitting the model"
    w <- nchar(msg) + 4
    cat("\n", "+", strrep("-", w), "+\n",
        "|  ", msg, "  |\n",
        "+", strrep("-", w), "+\n\n", sep = "")
  }

  modelInfo$control_optimizer$cores <- min(modelInfo$control_optimizer$rstarts,
                                           modelInfo$control_optimizer$cores)
  # Fit the model:
  Optim <- optimizer(control_manifold = modelInfo$control_manifold,
                     control_transform = modelInfo$control_transform,
                     control_estimator = modelInfo$control_estimator,
                     control_optimizer = modelInfo$control_optimizer)
  names(Optim$parameters) <- modelInfo$parameters_labels
  names(Optim$transparameters) <- modelInfo$transparameters_labels

  # Collect all the information about the optimization:

  elapsed <- Optim$elapsed

  #### Estimated model structures ####

  # Create the structures of transformed parameters:
  transformed_pars <- fill_in(modelInfo$trans, Optim$transparameters)

  # Create the structures of untransformed parameters:
  parameters <- transformed_pars[names(modelInfo$param)]

  #### Process the fit information ####

  loss <- sum(unlist(lapply(Optim$outputs$estimators$doubles,
                            FUN = \(x) x[[1]])))
  loglik <- sum(unlist(lapply(Optim$outputs$estimators$doubles,
                              FUN = \(x) x[[2]])))
  penalty <- sum(unlist(lapply(Optim$outputs$estimators$doubles,
                               FUN = \(x) x[[5]])))
  penalized_loss <- loss + penalty
  penalized_loglik <- loglik + penalty

  #### latent object ####

  result <- new("lcfa",
                version            = as.character( packageVersion('latent') ),
                call               = mc,
                timing             = elapsed,
                dataList           = dataList,
                modelInfo          = modelInfo,
                Optim              = Optim,
                parameters         = parameters,
                transformed_pars   = transformed_pars,
                loglik             = loglik,
                penalized_loglik   = penalized_loglik,
                loss               = loss,
                penalized_loss     = penalized_loss
  )

  if(message) {
    msg <- "Computing standard errors"
    w <- nchar(msg) + 4
    cat("\n", "+", strrep("-", w), "+\n",
        "|  ", msg, "  |\n",
        "+", strrep("-", w), "+\n\n", sep = "")
  }

  # Standard errors:
  if(se != "none" || isFALSE(se)) {
    Optim$SE <- se(result, type = "standard", digits = 9)
  }

  # Fit by group:
  fit_by_group <- latInspect(result, what = "fit")
  loss.group <- unlist(lapply(fit_by_group, FUN = \(x) x["loss"]))
  penalized_loss.group <- unlist(lapply(fit_by_group, FUN = \(x) x["penalized_loss"]))
  loglik.group <- unlist(lapply(fit_by_group, FUN = \(x) x["loglik"]))
  penalized_loglik.group <- unlist(lapply(fit_by_group, FUN = \(x) x["penalized_loglik"]))
  penalty.group <- unlist(lapply(fit_by_group, FUN = \(x) x["penalty"]))

  result@Optim <- Optim

  #### Return ####

  return(result)

}


