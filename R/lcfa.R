# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 31/10/2025
#'
#' @title
#' Fit a Confirmatory Factor Analysis (CFA) model with lavaan syntax.
#'
#' @usage
#'
#' lcfa(data, model = NULL, cor = "pearson",
#' estimator = "ml", group = NULL,
#' sample.cov = NULL, nobs = NULL,
#' positive = FALSE, penalties = TRUE,
#' missing = "pairwise.complete.obs",
#' std.lv = FALSE, do.fit = TRUE, control = NULL, ...)
#'
#' @param data data frame or matrix.
#' @param model lavaan's model syntax.
#' @param cor Correlation types: "pearson" and "poly". Defaults to "pearson".
#' @param estimator Available estimators: "ml", "uls", and "dwls". Defaults to "ml".
#' @param group .
#' @param sample.cov Covariance matrix between the items. Defaults to NULL.
#' @param nobs Number of observations. Defaults to NULL.
#' @param positive Force a positive-definite solution. Defaults to FALSE.
#' @param penalties list of penalty terms for the parameters.
#' @param missing Method to handle missing data.
#' @param std.lv Provide the parameters of the standardized model.
#' @param do.fit TRUE to fit the model and FALSE to return only the model setup. Defaults to TRUE.
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
lcfa <- function(data, model = NULL, cor = "pearson",
                 estimator = "ml", group = NULL,
                 sample.cov = NULL, nobs = NULL,
                 positive = FALSE, penalties = TRUE,
                 missing = "pairwise.complete.obs",
                 std.lv = TRUE, do.fit = TRUE, control = NULL,
                 ...) {

  # Check the arguments to control_optimizer and create defaults:
  control$penalties <- penalties
  control <- lcfa_control(control)

  # Extract the lavaan model:
  model_syntax <- model
  LAV <- lavaan::cfa(model = model_syntax, data = data,
                     sample.cov = sample.cov,
                     sample.nobs = nobs,
                     std.lv = std.lv,
                     do.fit = FALSE, group = group,
                     ...)

  LAV@Options$positive <- positive

  # extract slots from dummy lavaan object
  lavpartable    <- LAV@ParTable
  lavmodel       <- LAV@Model
  lavdata        <- LAV@Data
  lavoptions     <- LAV@Options
  lavsamplestats <- LAV@SampleStats
  lavcache       <- LAV@Cache
  timing         <- LAV@timing
  ngroups        <- LAV@Model@ngroups
  item_names     <- LAV@Data@ov.names
  nobs           <- LAV@Data@nobs
  group_label    <- LAV@Data@group.label
  item_label     <- LAV@Data@ov.names
  factor_label   <- replicate(ngroups, list(LAV@Model@dimNames[[1]][[2]]))
  X              <- LAV@Data@X

  # Rename columns:
  for(i in 1:ngroups) {
    colnames(X[[i]]) <- item_names[[i]]
  }

  # Get the correlation and weight matrices:
  correl <- vector("list", length = ngroups)
  names(X) <- names(correl) <- group_label

  # Model for the parameters:
  model <- getmodel_fromlavaan(LAV)

  if(ngroups == 1) model <- list(model)

  # Estimate the correlation matrix for each group:
  for(i in 1:ngroups) {

    if(is.null(sample.cov)) {

      correl[[i]] <- correlation(data = X[[i]],
                                 item_names = item_names[[i]],
                                 cor = cor,
                                 estimator = estimator,
                                 missing = missing)

    } else {

      correl[[i]]$R <- sample.cov
      p <- nrow(sample.cov)
      correl[[i]]$W <- matrix(1, nrow = p, ncol = p)

    }

  }

  # Data and structure information:
  nitems <- as.list(lavmodel@nvar)
  npatterns <- lapply(nitems, FUN = \(p) 0.5*p*(p+1))
  nfactors <- lapply(model, FUN = \(x) ncol(x$lambda))

  data_list <- vector("list")
  data_list$ngroups <- ngroups
  data_list$data <- X
  data_list$nobs <- nobs
  data_list$nitems <- nitems
  data_list$npatterns <- npatterns
  data_list$nfactors <- nfactors
  data_list$correl <- correl
  data_list$positive <- positive
  data_list$estimator <- estimator
  data_list$group_label <- group_label
  data_list$item_label <- item_label
  data_list$factor_label <- factor_label

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

  #### Collect all the model information ####

  # Model information:
  modelInfo <- list(nobs = nobs,
                    nparam = nparam - rest,
                    npatterns = npatterns,
                    dof = sum(unlist(npatterns)) - nparam + rest,
                    ntrans = ntrans,
                    parameters_labels = parameters_labels,
                    transparameters_labels = transparameters_labels,
                    cfa_param = cfa_param,
                    cfa_trans = cfa_trans,
                    control_manifold = control_manifold,
                    control_transform = control_transform,
                    control_estimator = control_estimator,
                    control = control)

  #### Fit the model ####

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
  # Fit the model:
  x <- optimizer(control_manifold = control_manifold,
                 control_transform = control_transform,
                 control_estimator = control_estimator,
                 control_optimizer = control)

  # Collect all the information about the optimization:

  Optim$opt <- x
  elapsed <- x$elapsed

  #### Estimated model structures ####

  # Create the structures of untransformed parameters:
  indices_pars <- match(modelInfo$parameters_labels,
                        unlist(modelInfo$cfa_param))
  vv <- rep(0, times = length(unlist(modelInfo$cfa_param)))
  vv[indices_pars] <- Optim$opt$parameters
  parameters <- fill_list_with_vector(modelInfo$cfa_param, vv)
  parameters <- allnumeric(parameters)
  # FIXED PARAMETERS?

  # Create the structures of transformed parameters:
  indices_trans <- match(modelInfo$transparameters_labels,
                         unlist(modelInfo$cfa_trans))
  vv <- rep(0, times = length(unlist(modelInfo$cfa_trans)))
  vv[indices_trans] <- Optim$opt$transparameters
  transformed_pars <- fill_list_with_vector(modelInfo$cfa_trans, vv)
  transformed_pars <- allnumeric(transformed_pars)

  #### Process the fit information ####

  # Get the indices of the estimator structures "cfa_dwls" and "cfa_ml":
  all_estimators <- unlist(lapply(modelInfo$control_estimator, FUN = \(x) x$estimator))
  indices_cfa <- which(all_estimators == "cfa_dwls" | all_estimators == "cfa_ml")

  # Get the indices of the estimator structures "logdetmat" (penalties):
  indices_logdetmat <- which(all_estimators == "logdetmat")

  # Initialize the objects to be returned:
  loss <- penalized_loss <- loglik <- penalized_loglik <- penalty <-
    vector("list", length = ngroups)

  # For each group, extract the loss, penalized loss, loglik and penalized loglik
  for(i in 1:ngroups) {

    k <- indices_cfa[i]

    loss[[i]] <- c(x$outputs$estimators$doubles[[k]][[1]])
    loglik[[i]] <- c(x$outputs$estimators$doubles[[k]][[2]])

    # If there are penalties, add the penalties to the loss or loglik:
    if(length(indices_logdetmat) > 0) {

      l <- indices_logdetmat[i]
      penalty[[i]] <- c(x$outputs$estimators$doubles[[l]][[1]])
      penalized_loss[[i]] <- loss[[i]] + penalty[[i]]
      penalized_loglik[[i]] <- loglik[[i]] + penalty[[i]]

    } else {

      penalized_loss[[i]] <- loss[[i]]
      penalized_loglik[[i]] <- loglik[[i]]

    }

  }

  loss <- sum(unlist(loss))
  penalized_loss <- sum(unlist(penalized_loss))
  loglik <- sum(unlist(loglik))
  penalized_loglik <- sum(unlist(penalized_loglik))

  #### Return ####

  lcfa_list <- new("lcfa",
                   version            = as.character( packageVersion('latent') ),
                   call               = mc, # matched call
                   timing             = elapsed, # timing information
                   modelInfo          = modelInfo, # modelInfo
                   Optim              = Optim, # Optim
                   parameters         = parameters,
                   transformed_pars   = transformed_pars,
                   loglik             = loglik, # loglik values
                   penalized_loglik   = penalized_loglik,
                   loss               = loss,
                   penalized_loss     = penalized_loss
  )

  # class(lcfa_list) <- "lcfa.list"

  return(lcfa_list)

  ## for lavaan like object
  ## need to fill in these objects into the lav dummy objects
  ## will take me longer to find how to match these slots
  # lcfa <- new("lcfa",
  #             version      = as.character( packageVersion('latent') ),
  #             call         = mc,                  # match.call
  #             timing       = timing,              # list
  #             Options      = lavoptions,          # list *
  #             ParTable     = lavpartable,         # list *
  #             pta          = LAV@pta,             # list
  #             Data         = lavdata,             # S4 class
  #             SampleStats  = lavsamplestats,      # S4 class
  #             Model        = lavmodel,            # S4 class *
  #             Cache        = lavcache,            # list
  #             Fit          = lavfit,              # S4 class *
  #             boot         = list(),
  #             optim        = lavoptim,
  #             implied      = lavimplied,          # list *
  #             vcov         = lavvcov,                    *
  #             test         = TEST,                       *
  #             h1           = h1,
  #             internal     = internal,
  #             external     = extslot              # can add extra info from latent
  # )


}

