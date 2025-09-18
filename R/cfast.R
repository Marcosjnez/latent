# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 05/09/2025
#'
#' @title
#' Fit a Confirmatory Factor Analysis (CFA) model with lavaan syntax.
#'
#' @usage
#'
#' cfast(data, model = NULL, cor = "pearson", estimator = "ml",
#' group = NULL, sample.cov = NULL, nobs = NULL,
#' missing = "pairwise.complete.obs", W = NULL, std.lv = FALSE,
#' positive = FALSE, do.fit = TRUE, control = NULL)
#'
#' @examples
#'
#' \dontrun{
#' # The famous Holzinger and Swineford (1939) example
#' HS.model <- ' visual  =~ x1 + x2 + x3
#'               textual =~ x4 + x5 + x6
#'               speed   =~ x7 + x8 + x9 '
#'
#' fit <- cfast(model = HS.model, data = HolzingerSwineford1939)
#' summary(fit, fit.measures = TRUE)
#'}
#'
#' @export
cfast <- function(data, model = NULL, cor = "pearson", estimator = "uls",
                  group = NULL, sample.cov = NULL, nobs = NULL,
                  missing = "pairwise.complete.obs", W = NULL, std.lv = FALSE,
                  positive = FALSE, do.fit = TRUE, control = NULL) {

  # Check the arguments to control_optimizer and create defaults:
  control <- cfast_control(control)

  if(is.null(nobs)) {
    if(nrow(data) == ncol(data)) {
      stop("Please, provide the number of observations in nobs.")
    }
    nobs <- nrow(data)
  }

  # Extract the lavaan model:
  model_syntax <- model
  extract_fit <- lavaan::cfa(model = model_syntax, data = data,
                             sample.cov = sample.cov, std.lv = std.lv,
                             do.fit = FALSE, group = group)
  item_names <- unique(extract_fit@ParTable$rhs[extract_fit@ParTable$op == "=~"])
  # Model for the parameters:
  model <- getmodel_fromlavaan(extract_fit)

  if(!is.list(model[[1]])) {
    model <- list(model)
  }
  ngroups <- length(model)

  # Get the correlation and weight matrices:
  correl <- vector("list", length = ngroups)
  if(ngroups > 1) {
    data_split <- split(data, data[[group]])
    for(i in 1:ngroups) {
      correl[[i]] <- correlation(data = data_split[[i]], item_names = item_names,
                                 cor = cor, estimator = estimator, missing = missing)
    }
  } else {
    correl[[1]] <- correlation(data = data, item_names = item_names, cor = cor,
                               estimator = estimator, missing = missing)
  }

  # Data and structure information:
  nitems <- length(item_names)
  npatterns <- 0.5*nitems*(nitems+1)
  nfactors <- ncol(model[[1]]$lambda)
  data_list <- vector("list")
  data_list$data <- data
  data_list$nobs <- nobs
  data_list$ngroups <- ngroups
  data_list$nitems <- nitems
  data_list$npatterns <- npatterns
  data_list$nfactors <- nfactors
  data_list$correl <- correl
  data_list$positive <- positive

  ## store original call
  mc  <- match.call()

  #### Create the model ####

  # Get the model specification:
  full_model <- get_full_cfa_model(data_list = data_list,
                                   model = model,
                                   control = control)
  list2env(full_model, envir = environment())

  #### Create the structures ####

  # Generate the structures for optimization:
  structures <- get_cfa_structures(data_list = data_list,
                                   full_model = full_model,
                                   control = control)
  list2env(structures, envir = environment())

  #### Fit the model ####

  # Model information:
  modelInfo <- list(nobs = nobs,
                    nparam = nparam,
                    npatterns = npatterns,
                    df = npatterns - nparam,
                    ntrans = ntrans,
                    parameters_labels = parameters_labels,
                    transparameters_labels = transparameters_labels,
                    cfa_param = cfa_param,
                    cfa_trans = cfa_trans)

  # Data for the optimization algorithms:
  Optim <- list(data = data,
                data_list = data_list,
                control_manifold = control_manifold,
                control_transform = control_transform,
                control_estimator = control_estimator,
                control = control)

  if(!do.fit) {

    lcfa_list <- new("lcfa",
                     version            = as.character( packageVersion('latent') ),
                     call               = mc, # matched call
                     timing             = numeric(), # timing information
                     modelInfo          = modelInfo, # modelInfo
                     Optim              = Optim, # Optim
                     parameters         = list(),
                     transformed_pars   = list(),
                     loglik             = numeric(), # loglik values
                     penalized_loglik   = numeric(),
                     loss               = numeric(),
                     penalized_loss     = numeric()
    )

    return(lcfa_list)

  }

  control$cores <- min(control$rstarts, control$cores)
  # Perform the optimization (fit the model):
  x <- optimizer(control_manifold = control_manifold,
                 control_transform = control_transform,
                 control_estimator = control_estimator,
                 control_optimizer = control)

  # Collect all the information about the optimization:

  Optim$opt <- x
  elapsed <- x$elapsed

  if(estimator == "uls" || estimator == "dwls") {
    loss <- x$f
    penalized_loss <- x$f
    loglik <- numeric()
    penalized_loglik <- numeric()
  } else if(estimator == "ml") {
    loss <- x$f
    penalized_loss <- x$f
    loglik <- -loss
    penalized_loglik <- -penalized_loss
  }

  #### Process the outputs ####

  matrices <- outputs <- vector("list", length = ngroups)

  for(i in 1:ngroups) {

    p <- nitems
    q <- nfactors

    # Arrange lambda parameter estimates:
    matrices[[i]]$lambda <- matrix(x$outputs$estimators$matrices[[i]][[1]], p, q)

    # Arrange psi parameter estimates:
    matrices[[i]]$psi <- matrix(x$outputs$estimators$matrices[[i]][[2]], q, q)

    # Arrange theta parameter estimates:
    matrices[[i]]$theta <- matrix(x$outputs$estimators$matrices[[i]][[3]], p, p)
    # uniquenesses_hat[[i]] <- diag(psi_hat[[i]])

    # # Model matrix:
    # outputs[[i]]$model <- matrix(x$outputs$estimators$matrices[[i]][[4]], p, p)
    #
    # # Residual matrix:
    # outputs[[i]]$residuals <- matrix(x$outputs$estimators$matrices[[i]][[5]], p, p)
    #
    # # Weight matrix:
    # if(control_estimator[[i]]$estimator == "uls" ||
    #    control_estimator[[i]]$estimator == "dwls") {
    #   outputs[[i]]$W <- matrix(x$outputs$estimators$matrices[[i]][[6]], p, p)
    # }
    #
    # # Uniquenesses:
    # outputs[[i]]$uniquenesses <- c(x$outputs$estimators$vectors[[i]][[1]])

  }

  #### Return ####

  lcfa_list <- new("lcfa",
                   version            = as.character( packageVersion('latent') ),
                   call               = mc, # matched call
                   timing             = elapsed, # timing information
                   modelInfo          = modelInfo, # modelInfo
                   Optim              = Optim, # Optim
                   parameters         = matrices,
                   transformed_pars   = matrices,
                   loglik             = loglik, # loglik values
                   penalized_loglik   = penalized_loglik,
                   loss               = loss,
                   penalized_loss     = penalized_loss
  )

  # class(lcfa_list) <- "lcfa.list"

  return(lcfa_list)

}

