#' @title Maximum likelihood estimation of positive-definite polychoric correlation matrices
#'
#' @description
#' \code{lpoly} estimates a polychoric correlation matrix from ordinal data using
#' maximum likelihood. The function can use either a one-step or a two-step
#' estimation approach and optionally enforce positive-semidefiniteness or
#' positive-definiteness through constrained estimation or penalties.
#'
#' @usage
#' lpoly(data,
#'       method = "two-step",
#'       model = NULL,
#'       positive = FALSE,
#'       penalties = FALSE,
#'       do.fit = TRUE,
#'       message = FALSE,
#'       control = NULL,
#'       ...)
#'
#' @param data A data frame or matrix containing the raw ordinal data.
#' @param method Character string indicating the estimation method. Possible
#'   values are \code{"one-step"} and \code{"two-step"}. Default is
#'   \code{"two-step"}.
#' @param model Optional model object used internally for initialization or
#'   custom model setup. Default is \code{NULL}.
#' @param positive Logical. If \code{TRUE}, the estimated polychoric correlation
#'   matrix is forced to be positive semidefinite. Default is \code{FALSE}.
#' @param penalties Logical. If \code{TRUE}, penalties are added to the
#'   objective function to encourage a positive-definite solution. Default is
#'   \code{FALSE}.
#' @param do.fit Logical. If \code{TRUE}, the model is fitted. If \code{FALSE},
#'   only the model setup is returned. Default is \code{TRUE}.
#' @param message Logical. If \code{TRUE}, progress messages are printed during
#'   estimation. Default is \code{FALSE}.
#' @param control A list of control parameters for the optimization algorithm.
#'   This may include starting values, convergence tolerances, maximum number of
#'   iterations, and other optimizer-specific options.
#' @param ... Additional arguments passed to internal optimization and model
#'   setup routines.
#'
#' @details
#' \code{lpoly} estimates positive-definite or positive-semidefinite
#' polychoric correlation matrices from ordinal data. The function is designed
#' for situations in which the unrestricted polychoric matrix is not
#' guaranteed to be admissible, for example because of sampling variability or
#' sparse response patterns.
#'
#' Two estimation strategies are available:
#' \itemize{
#'   \item \code{"two-step"}: thresholds are estimated first and the correlation
#'   matrix is estimated in a second step.
#'   \item \code{"one-step"}: thresholds and correlations are estimated jointly.
#' }
#'
#' If \code{positive = TRUE}, the estimated matrix is constrained to be
#' positive semidefinite. If \code{penalties = TRUE}, penalty terms are added
#' to the objective function to encourage positive-definiteness.
#'
#' If \code{do.fit = FALSE}, the function returns the model setup without
#' running the optimizer.
#'
#' @return
#' A list containing the fitted model and related information. Typical elements
#' include:
#' \itemize{
#'   \item \code{version}: Version number of \pkg{latent} used when the model
#'   was estimated.
#'   \item \code{call}: Matched call used to estimate the model.
#'   \item \code{ModelInfo}: Information about the model specification and data.
#'   \item \code{Optim}: Output of the optimization routine.
#'   \item \code{parameters}: Structure containing the model parameters.
#'   \item \code{transparameters}: Structure containing transformed model
#'   parameters.
#'   \item \code{loglik}: Log-likelihood of the fitted model.
#'   \item \code{penalized_loglik}: Penalized log-likelihood of the fitted model.
#' }
#'
#' @examples
#' \dontrun{
#' fit <- lpoly(data = values)
#'
#' fit_psd <- lpoly(data = values, positive = TRUE)
#'
#' setup_only <- lpoly(data = values, do.fit = FALSE)
#' }
#'
#' @export
lpoly <- function(data,
                  method = "two-step",
                  model = NULL,
                  positive = FALSE,
                  penalties = FALSE,
                  do.fit = TRUE,
                  message = FALSE,
                  control = NULL,
                  ...) {

  dots <- list(...)

  ## store original call
  mc  <- match.call()

  if(method == "two-step") {
    positive <- FALSE
    penalties <- FALSE
  }

  # Check the arguments to control_optimizer and create defaults:
  control$penalties <- penalties
  control$positive <- positive
  control <- lpoly_control(control)

  #### Create the data_list ####

  data_list <- create_lpoly_datalist(data, control)
  list2env(data_list, envir = environment())

  #### Create the model ####

  full_model <- create_lpoly_model(data_list = data_list,
                                   model = model,
                                   control = control)
  list2env(full_model, envir = environment())

  #### Create the modelInfo ####

  modelInfo <- create_lpoly_modelInfo(data_list = data_list,
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

  if(method == "one-step") {

    modelInfo$control_optimizer$cores <- min(modelInfo$control_optimizer$rstarts,
                                             modelInfo$control_optimizer$cores)

    # Fit the model:
    Optim <- optimizer(control_manifold = modelInfo$control_manifold,
                       control_transform = modelInfo$control_transform,
                       control_estimator = modelInfo$control_estimator,
                       control_optimizer = modelInfo$control_optimizer)

  } else if(method == "two-step") {
    Optim <- polyfast(as.matrix(data), cores = parallel::detectCores())
    threslds <- lapply(Optim$thresholds, FUN = \(x) matrix(x[!is.infinite(x)],
                                                           ncol = 1L))
    Optim$parameters <- c(unlist(threslds),
                          Optim$correlation[lower.tri(Optim$correlation,
                                                      diag = FALSE)])
    Optim$transparameters <- c(unlist(threslds),
                               Optim$correlation[lower.tri(Optim$correlation,
                                                           diag = TRUE)])
    Optim$f <- 0
  } else {
    stop("Unknown method")
  }

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

  result <- new("latent",
                version            = as.character( packageVersion('latent') ),
                call               = mc,
                timing             = elapsed,
                dataList           = data_list,
                modelInfo          = modelInfo,
                Optim              = Optim,
                parameters         = parameters,
                transformed_pars   = transformed_pars,
                loglik             = loglik,
                penalized_loglik   = penalized_loglik,
                loss               = loss,
                penalized_loss     = penalized_loss
  )

  #### Return ####

  return(result)

}

create_lpoly_datalist <- function(data, control) {

  # Estimate the polychoric correlations without positive-definite constraints:
  polychorics <- polyfast(as.matrix(data))
  S <- polychorics$correlation # Polychoric correlation matrix
  taus <- polychorics$thresholds # Thresholds
  n <- polychorics$contingency_tables # Contingency tables
  nobs <- nrow(data)
  nitems <- ncol(data)
  item_label <- colnames(data)
  npatterns <- length(unlist(n))

  K <- vector(length = nitems)
  item_cat <- vector("list", length = nitems)
  for(i in 1:nitems) {
    K[i] <- length(taus[[i]])-2L
    item_cat[[i]] <- unique(data[, i])
  }

  data_list <- vector("list")
  data_list$data <- data
  data_list$nobs <- nobs
  data_list$nitems <- nitems
  data_list$npatterns <- npatterns
  data_list$item_label <- item_label
  data_list$n <- n
  data_list$S <- S
  data_list$taus <- taus
  data_list$K <- K
  data_list$item_cat <- item_cat

  return(data_list)

}

create_lpoly_model <- function(data_list, model, control) {

  # Generate the model syntax and initial parameter values

  list2env(data_list, envir = environment())

  # Initialize the objects to store the initial parameters:
  param <- trans <- vector("list")
  fixed <- nonfixed <- fixed_values_list <- vector("list")

  #### Model for the transformed parameters ####

  taus_item <- paste("taus.", item_label, control$subfix, sep = "")

  # Transformed parameters:
  list_struct <- vector("list")
  k <- 1L
  for(i in 1:nitems) {

    # Taus:
    list_struct[[k]] <- list(name = taus_item[i],
                             type = "matrix",
                             dim = c(K[i], 1),
                             rownames = 1:K[i],
                             colnames = item_label[i])
    k <- k+1L

  }

  if(control$positive) {

    # X:
    list_struct[[k]] <- list(name = "X",
                             type = "matrix",
                             dim = c(nitems, nitems),
                             rownames = item_label,
                             colnames = item_label)
    k <- k+1L

  }

  # Covariance matrix:
  list_struct[[k]] <- list(name = "S",
                           type = "matrix",
                           dim = c(nitems, nitems),
                           rownames = item_label,
                           colnames = item_label,
                           symmetric = TRUE)
  k <- k+1L

  trans <- create_parameters(list_struct)

  #### Model for the parameters ####

  param <- trans
  diag(param$S) <- "1"

  #### Create the initial values for the parameters ####

  init_param <- vector("list", length = control$rstarts)
  # Thresholds without the infite values:
  threslds <- lapply(taus, FUN = \(x) x[!is.infinite(x)])

  for(rs in 1:control$rstarts) {

    init_param[[rs]] <- vector("list")

    for(i in 1:nitems) {
      init_param[[rs]][[taus_item[i]]] <- matrix(threslds[[i]], ncol = 1L)
    }

    if(control$positive) {
      # Initial values for the square root of correlations:
      init_param[[rs]]$X <- real_sqrtmat(data_list$S)
      # Initial values for the polychoric correlations:
      init_param[[rs]]$S <- crossprod(init_param[[rs]]$X)
    } else {
      init_param[[rs]]$S <- data_list$S
    }

  }


  #### Custom initial values ####

  # # Replace initial starting values by custom starting values:
  #
  # if(!is.null(control$start)) {
  #
  #   nm <- names(control$start)
  #   nm <- nm[!vapply(control$start, is.null, logical(1))]
  #
  #   for (i in seq_len(control$rstarts)) {
  #     common_nm <- intersect(nm, names(init_param[[i]]))
  #     for (j in common_nm) {
  #       init_param[[i]][[j]] <- insert_object(init_param[[i]][[j]],
  #                                             control$start[[j]])
  #     }
  #   }
  #
  # }

  #### Return ####

  result <- list(param = param,
                 trans = trans,
                 init_param = init_param)

  return(result)

}

create_lpoly_modelInfo <- function(data_list, full_model, control) {

  # Generate control_manifold, control_transform, and control_estimator

  list2env(data_list, envir = environment())
  list2env(full_model, envir = environment())

  taus_item <- paste("taus.", item_label, control$subfix, sep = "")

  #### Manifolds ####

  if(control$positive) {

    manifolds <- list(
      list(manifold = "euclidean",
           parameters = list(param[taus_item])),
      list(manifold = "oblq", parameters = "X",
           extra = list(p = nitems, q = nitems))
    )

  } else {

    manifolds <- list(
      list(manifold = "euclidean",
           parameters = list(param[taus_item])),
      list(manifold = "euclidean", parameters = "S",
           extra = list(p = nitems, q = nitems))
    )

  }

  control_manifold <- create_manifolds(manifolds = manifolds,
                                       structures = param)

  #### Transformations ####

  transform <- list()

  if(control$positive) {

    transforms <- list(
      list(transform = "crossprod",
           parameters_in = "X",
           parameters_out = "S",
           extra = list(p = nitems))
    )

  } else {

    transforms <- list()

  }

  control_transform <- create_transforms(transforms = transforms,
                                         structures = trans)

  #### Estimators ####

  estimators <- list()

  estimators[[1]] <- list(estimator = "polycor",
                          parameters = c(list(S = trans$S),
                                         trans[taus_item]),
                          extra = list(n = data_list$n,
                                       p = nitems,
                                       N = data_list$nobs))

  if(control$reg) {

    lower_indices <- which(lower.tri(trans$S, diag = TRUE))
    estimators[[2]] <- list(estimator = "logdetmat",
                            parameters = "S",
                            extra = list(lower_indices = lower_indices-1L,
                                         p = nitems,
                                         logdetw = control$penalties$logdet$w))

  }

  control_estimator <- create_estimators(estimators = estimators,
                                         structures = trans)

  #### Pass the initial values to vectors ####

  idx_transformed <- unlist(lapply(control_transform,
                                   FUN = \(x) unlist(x$indices_out)+1L))
  inits <- create_init(trans, param, init_param,
                       idx_transformed = idx_transformed, control)

  parameters <- inits$parameters
  parameters_labels <- names(parameters[[1]])
  nparam <- length(parameters_labels)

  transparameters <- inits$transparameters
  transparameters_labels <- names(transparameters[[1]])
  ntrans <- length(transparameters_labels)

  trans2param <- match(parameters_labels, transparameters_labels)

  #### Set up the optimizer ####

  # Create defaults for the control of the optimizer:
  control_optimizer <- control
  control_optimizer$parameters <- parameters
  control_optimizer$transparameters <- transparameters
  control_optimizer$init_param <- init_param
  control_optimizer$transparam2param <- trans2param-1L

  #### Collect all the model information ####

  modelInfo <- list(param = param,
                    trans = trans,
                    nparam = nparam,
                    ntrans = ntrans,
                    parameters_labels = parameters_labels,
                    transparameters_labels = transparameters_labels,
                    dof = npatterns - nparam,
                    control_manifold = control_manifold,
                    control_transform = control_transform,
                    control_estimator = control_estimator,
                    control_optimizer = control_optimizer)

  #### Return ####

  return(modelInfo)

}

asymptotic_poly <- function(fit, model = NULL) {

  control_manifold <- fit@modelInfo$control_manifold
  control_transform <- fit@modelInfo$control_transform
  control_estimator <- fit@modelInfo$control_estimator
  control_optimizer <- fit@modelInfo$control_optimizer

  control_optimizer$parameters[[1]] <- fit@Optim$parameters
  control_optimizer$transparameters[[1]] <- fit@Optim$transparameters

  x <- get_hess(control_manifold, control_transform,
                control_estimator, control_optimizer,
                cores = parallel::detectCores())
  ACOV <- solve(x$h)

  return(ACOV)

}
