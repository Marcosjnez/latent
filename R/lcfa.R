# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 27/10/2025
#'
#' @title
#' Fit a Confirmatory Factor Analysis (CFA) model with lavaan syntax.
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
lcfa <- function(data, model = NULL, cor = "pearson",
                 estimator = "ml", group = NULL,
                 sample.cov = NULL, nobs = NULL, W = NULL,
                 positive = FALSE, penalties = TRUE,
                 missing = "pairwise.complete.obs",
                 std.lv = TRUE, do.fit = TRUE, control = NULL) {

  # Check the arguments to control_optimizer and create defaults:
  control$penalties <- penalties
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
                             sample.cov = sample.cov,
                             sample.nobs = nobs,
                             std.lv = std.lv,
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
  } else {
    data_split <- list()
    data_split[[1]] <- data
  }

  # if(is.null(sample.cov)) {
  # }

  for(ng in 1:ngroups) {

    if(is.null(sample.cov)) {

      correl[[ng]] <- correlation(data = data_split[[ng]], item_names = item_names,
                                  cor = cor, estimator = estimator, missing = missing)
    } else {

      correl[[ng]]$R <- sample.cov
      p <- nrow(sample.cov)
      correl[[ng]]$W <- matrix(1, nrow = p, ncol = p)

    }

  }

  # Data and structure information:
  nitems <- length(item_names)
  npatterns <- 0.5*nitems*(nitems+1) * ngroups
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
  data_list$estimator <- estimator

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
                    nparam = nparam - rest,
                    npatterns = npatterns,
                    dof = npatterns - nparam + rest,
                    ntrans = ntrans,
                    parameters_labels = parameters_labels,
                    transparameters_labels = transparameters_labels,
                    cfa_param = cfa_param,
                    cfa_trans = cfa_trans,
                    control_manifold = control_manifold,
                    control_transform = control_transform,
                    control_estimator = control_estimator,
                    control = control)

  # Data for the optimization algorithms:
  Optim <- list(data = data,
                data_list = data_list)

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
    loss <- x$f / ngroups
    penalized_loss <- x$f / ngroups
    loglik <- numeric()
    penalized_loglik <- numeric()
  } else if(estimator == "ml") {
    loss <- x$f / ngroups
    penalized_loss <- x$f / ngroups
    loglik <- -loss / ngroups
    penalized_loglik <- -penalized_loss / ngroups
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

    # Model matrix:
    outputs[[i]]$model <- matrix(x$outputs$estimators$matrices[[i]][[4]], p, p)

    # Residual matrix:
    outputs[[i]]$residuals <- matrix(x$outputs$estimators$matrices[[i]][[5]], p, p)
    #
    # # Weight matrix:
    # if(control_estimator[[i]]$estimator == "uls" ||
    #    control_estimator[[i]]$estimator == "dwls") {
    #   outputs[[i]]$W <- matrix(x$outputs$estimators$matrices[[i]][[6]], p, p)
    # }
    #
    # Uniquenesses:
    outputs[[i]]$uniquenesses <- c(x$outputs$estimators$vectors[[i]][[1]])

  }

  indices_trans <- match(modelInfo$transparameters_labels,
                         unlist(modelInfo$cfa_trans))
  vv <- rep(0, times = length(unlist(modelInfo$cfa_trans)))
  vv[indices_trans] <- Optim$opt$transparameters
  transformed_pars <- fill_list_with_vector(modelInfo$cfa_trans,
                                            vv)
  transformed_pars <- allnumeric(transformed_pars)

  #### Return ####

  lcfa_list <- new("lcfa",
                   version            = as.character( packageVersion('latent') ),
                   call               = mc, # matched call
                   timing             = elapsed, # timing information
                   modelInfo          = modelInfo, # modelInfo
                   Optim              = Optim, # Optim
                   parameters         = matrices,
                   transformed_pars   = transformed_pars,
                   loglik             = loglik, # loglik values
                   penalized_loglik   = penalized_loglik,
                   loss               = loss,
                   penalized_loss     = penalized_loss
  )

  # class(lcfa_list) <- "lcfa.list"

  return(lcfa_list)

}

