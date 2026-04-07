# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 02/04/2026
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
                 std.lv = TRUE, std.ov = FALSE,
                 acov = "standard",
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

  control$std.ov <- std.ov
  control$positive <- positive
  control$penalties <- penalties
  control$estimator <- tolower(estimator)
  control <- lcfa_control(control)

  args <- as.list(match.call(expand.dots = TRUE))[-1]

  #### Create the datalist ####

  data_list <- create_cfa_datalist(
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
    args = args,
    ...
  )

  #### Create the model ####

  # Get the model specification:
  full_model <- create_cfa_model(data_list = data_list,
                                 model = model,
                                 control = control)
  list2env(full_model, envir = environment())

  #### Create the modelInfo ####

  modelInfo <- create_cfa_modelInfo(data_list = data_list,
                                    full_model = full_model,
                                    control = control)

  #### Fit the model ####

  if(!do.fit) {

    lcfa_list <- new("lcfa",
                     version            = as.character( packageVersion('latent') ),
                     call               = mc, # matched call
                     timing             = numeric(), # timing information
                     data_list          = data_list,
                     modelInfo          = modelInfo,
                     Optim              = list(),
                     parameters         = list(),
                     transformed_pars   = list(),
                     loglik             = numeric(), # loglik values
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

  loss <- Optim$f
  penalized_loss <- loss
  loglik <- sum(unlist(lapply(Optim$outputs$estimators$doubles,
                              FUN = \(x) x[[2]])))
  penalized_loglik <- loglik

  #### latent object ####

  result <- new("lcfa",
                version            = as.character( packageVersion('latent') ),
                call               = mc, # matched call
                timing             = elapsed, # timing information
                data_list          = data_list,
                modelInfo          = modelInfo,
                Optim              = Optim,
                parameters         = parameters,
                transformed_pars   = transformed_pars,
                loglik             = loglik, # loglik values
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

  if(mimic == "lavaan") {

    ## for lavaan like object
    ## need to fill in these objects into the lav dummy objects
    ## will take me longer to find how to match these slots


    ###
    ### cecks for latent things that cannot be mimimc == lavaan
    ###

    ## optiosn
    lavoptions$se <- "standard" ## se type argument?
    lavoptions$test <- "standard" ## test type argument?
    lavoptions$do.fit <- do.fit

    ## parameters
    parsdf <- data.frame(plabel = names(Optim$parameters),
                         est = Optim$parameters)

    sedf <- data.frame(plabel = modelInfo$parameters_labels,
                       se = Optim$SE$se)

    parsdf <- merge(parsdf, sedf)
    #parsdf

    ## partable
    lavpartable$est <- NULL
    lavpartable <- merge(lavpartable, parsdf, all =T)
    lavpartable$est <- ifelse(lavpartable$free == 0,
                              lavpartable$start,
                              lavpartable$est)
    lavpartable$se <- ifelse(lavpartable$free == 0, 0,
                             lavpartable$se)
    lavpartable <- lavpartable[order(lavpartable$id),]
    lavpartable <- as.list(lavpartable)
    #lp

    ## pta
    pta <- LAV@pta
    pta$names <- names(lavpartable)

    ## lavmodel
    ## only for 1 group for now
    ## how to adjust multiple group labels?

    lavmodel@GLIST$lambda <- transformed_pars$lambda.g1
    lavmodel@GLIST$psi <- transformed_pars$psi.g1
    lavmodel@GLIST$theta <- transformed_pars$theta.g1

    ## implied
    implied <- lav_model_implied(lavmodel)


    ## lavvcov
    ### order vcov from 1  up .. in names
    lavvcov <- LAV@vcov
    lavvcov$se <- lavoptions$se
    lavvcov$vcov <- Optim$SE$vcov

    ## test
    x <- lavpartable$est[lavpartable$free != 0]
    fx <- loss/2
    fx.group <- loss/2
    attr(fx, "fx.group") <- fx.group
    attr(x, "fx") <- fx

    TEST <- lavaan:::lav_model_test(lavobject = NULL,
                                    lavmodel = lavmodel,
                                    lavpartable = lavpartable,
                                    lavsamplestats = lavsamplestats,
                                    lavimplied = implied,
                                    lavh1 = h1,
                                    lavoptions = lavoptions,
                                    x = x,
                                    VCOV = Optim$SE$vcov,
                                    lavcache = lavcache,
                                    lavdata = lavdata,
                                    lavloglik = loglik_lav,
                                    test.UGamma.eigvals = FALSE)

    ### FIT
    attr(x, "iterations") <- Optim$iterations
    attr(x, "converged") <- Optim$convergence
    attr(x, "control") <- modelInfo$control_optimizer

    lavfit <-  lavaan:::lav_model_fit(lavpartable = lavpartable,
                                      lavmodel = lavmodel,
                                      lavimplied = implied,
                                      x = x,
                                      VCOV = Optim$SE$vcov,
                                      TEST = TEST)



    ##
    ## lavoptim
    optnames <- c('x','npar','iterations','converged','fx','fx.group','logl.group',
                  'logl','control')
    lavoptim <- lapply(optnames, function(x) slot(lavfit, x))
    names(lavoptim) <- optnames

    ## h1
    #### baseline model implied
    ## need to adjust with the estimated baseline model
    h1 <- lavaan:::lav_h1_implied_logl(lavdata = lavdata,
                                       lavsamplestats = lavsamplestats,
                                       lavpartable = lavpartable,
                                       lavoptions = lavoptions)

    ######

    ## add baseline
    #####

    ## set NACOC as NULL, so its recauclauted in some functions, to match residuals
    ## the latent NACOV have different dimensions
    lavsamplestats@NACOV <- vector("list", ngroups)

    ###
    loglik_lav <- lavaan:::lav_model_loglik(lavdata = lavdata,
                                            lavsamplestats = lavsamplestats,
                                            lavh1 = h1,
                                            lavimplied = implied,
                                            lavmodel = lavmodel,
                                            lavoptions = lavoptions)

    ####
    ####
    result <- new("lavaan",
                  version      = as.character( packageVersion('latent') ),
                  call         = mc,                  # match.call
                  timing       = timing,              # list
                  Options      = lavoptions,          # list *
                  ParTable     = lavpartable,         # list *
                  pta          = pta,             # list
                  Data         = lavdata,             # S4 class
                  SampleStats  = lavsamplestats,      # S4 class
                  Model        = lavmodel,            # S4 class *
                  Cache        = lavcache,            # list
                  Fit          = lavfit,              # S4 class * blav_model_fit
                  boot         = list(),
                  optim        = lavoptim,
                  implied      = implied,          # list *
                  vcov         = lavvcov,             #*
                  test         = TEST,                #* blav_model_test
                  h1           = h1,
                  loglik       = loglik_lav,
                  baseline     = list(),
                  internal     = list(),
                  external     = list()              # can add extra info from latent
    )

    tt <- lavaan:::lav_object_independence(object = result)

    baseline <- list(partable = tt@ParTable,
                     test = tt@test)

    result@baseline <- baseline

  }

  return(result)

}


