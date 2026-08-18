# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 17/08/2026
#'
#' Polychoric Correlation Matrix
#'
#' Estimate a polychoric correlation matrix from ordinal data using one-step or
#' two-step estimation.
#'
#' @usage
#' lpoly(data, method = "two-step", model = NULL,
#'       positive = FALSE, penalties = FALSE,
#'       start = NULL, do.fit = TRUE, message = FALSE,
#'       control = NULL, ...)
#'
#' @param data A data frame or matrix containing ordinal variables coded
#'   numerically.
#' @param method Character string indicating the estimation method. Available
#'   options are \code{"one-step"} and \code{"two-step"}.
#' @param model Optional named list used with \code{method = "one-step"} to
#'   fix parameters or impose equality constraints. Names should correspond to
#'   parameter blocks such as \code{S} or the threshold blocks. Numeric values
#'   fix parameters, repeated character labels impose equality constraints, and
#'   \code{NA} leaves the corresponding parameter unchanged.
#' @param positive Logical. If \code{TRUE}, the correlation matrix is
#'   parameterized through an oblique manifold so that the resulting matrix is
#'   positive semidefinite. This option is available for one-step estimation.
#' @param penalties Logical value or named list controlling regularization. If
#'   \code{TRUE}, the default log-determinant penalty is used. Two-step
#'   estimation does not use penalties.
#' @param start Optional named list of starting values used with
#'   \code{method = "one-step"}. Names should correspond to parameter blocks in
#'   the model. Partial matrices/vectors and \code{NA} values are handled in the
#'   same way as in \code{lca()}.
#' @param do.fit Logical. If \code{FALSE}, return the prepared but unfitted
#'   \code{"latent"} object.
#' @param message Logical. Print progress messages during estimation.
#' @param control Optional list of optimization controls.
#' @param ... Additional arguments reserved for future extensions.
#'
#' @details
#' With \code{method = "two-step"}, thresholds are obtained from the marginal
#' distributions and the polychoric correlations are then estimated pairwise.
#' With \code{method = "one-step"}, thresholds and correlations are optimized
#' jointly using the general optimization infrastructure of \pkg{latent}.
#'
#' The two-step method always uses the unconstrained pairwise estimator and
#' therefore does not support \code{model}, \code{start}, \code{positive}, or
#' \code{penalties}. Model constraints and custom starts are available for
#' \code{method = "one-step"}.
#'
#' @return An S4 object of class \code{"latent"}. The object contains the
#'   processed data in \code{dataList}, the model and optimization structures in
#'   \code{modelInfo}, optimizer output in \code{Optim}, and the estimated
#'   parameter structures in \code{parameters} and \code{transformed_pars}.
#'
#' @examples
#' \dontrun{
#' fit <- lpoly(data = values)
#' fit_one_step <- lpoly(data = values, method = "one-step")
#' fit_positive <- lpoly(data = values, method = "one-step", positive = TRUE)
#' }
#'
#' @export
lpoly <- function(data,
                  method = "two-step",
                  model = NULL,
                  positive = FALSE,
                  penalties = FALSE,
                  start = NULL,
                  do.fit = TRUE,
                  message = FALSE,
                  control = NULL,
                  ...) {

  #### Check input arguments ####

  if(!is.data.frame(data) && !is.matrix(data)) {
    stop("data must be a data.frame or matrix")
  }

  if(is.matrix(data)) {
    data <- as.data.frame(data)
  }

  if(nrow(data) == 0L || ncol(data) == 0L) {
    stop("data must contain at least one observation and one variable")
  }

  if(is.null(colnames(data)) ||
     any(colnames(data) == "") ||
     anyDuplicated(colnames(data))) {
    stop("data must have unique, non-empty column names")
  }

  numeric_data <- vapply(data, is.numeric, logical(1L))
  if(!all(numeric_data)) {
    stop("All variables in data must be numeric ordinal variables")
  }

  nlevels <- vapply(data, FUN = function(x) {
    length(unique(x[!is.na(x)]))
  }, FUN.VALUE = integer(1L))

  if(any(nlevels < 2L)) {
    stop("Every variable must contain at least two observed categories")
  }

  method <- match.arg(tolower(method), c("one-step", "two-step"))

  if(!is.logical(positive) || length(positive) != 1L || is.na(positive)) {
    stop("positive must be TRUE or FALSE")
  }

  if(!is.logical(do.fit) || length(do.fit) != 1L || is.na(do.fit)) {
    stop("do.fit must be TRUE or FALSE")
  }

  if(!is.logical(message) || length(message) != 1L || is.na(message)) {
    stop("message must be TRUE or FALSE")
  }

  if(!is.null(model)) {

    if(!is.list(model)) {
      stop("model must be NULL or a named list")
    }

    if(is.null(names(model)) ||
       any(names(model) == "") ||
       anyDuplicated(names(model))) {
      stop("model must have unique, non-empty names")
    }

  }

  if(!is.null(start) && !is.list(start)) {
    stop("start must be NULL or a named list")
  }

  if(!is.null(control) && !is.list(control)) {
    stop("control must be NULL or a list")
  }

  if(is.null(control)) {
    control <- list()
  }

  # Two-step estimation uses the ordinary unconstrained pairwise estimator.
  if(method == "two-step") {

    if(!is.null(model)) {
      stop("model constraints are only available with method = 'one-step'")
    }

    if(!is.null(start)) {
      stop("start is only available with method = 'one-step'")
    }

    positive <- FALSE
    penalties <- FALSE

  }

  #### Store original call ####

  mc <- match.call(expand.dots = TRUE)
  args <- lapply(as.list(mc)[-1L], eval, envir = parent.frame())

  #### Check control parameters ####

  control$method <- method
  control$positive <- positive
  control$penalties <- penalties
  control$start <- start
  control <- lpoly_control(control)

  #### Create the dataList ####

  dataList <- create_lpoly_dataList(data = data,
                                    control = control)
  dataList$args <- args

  #### Create the model ####

  full_model <- create_lpoly_model(dataList = dataList,
                                   model = model,
                                   control = control)

  #### Create the manifold, transformation, and estimator structures ####

  modelInfo <- create_lpoly_modelInfo(dataList = dataList,
                                      full_model = full_model,
                                      control = control)

  #### Fit the model ####

  if(!do.fit) {

    result <- new("latent",
                  version          = as.character(packageVersion("latent")),
                  call             = mc,
                  timing           = numeric(),
                  dataList         = dataList,
                  modelInfo        = modelInfo,
                  Optim            = list(),
                  parameters       = list(),
                  transformed_pars = list(),
                  extra            = list())

    #### Result ####

    return(result)

  }

  if(message) {
    print_lpoly_message("Fitting the model")
  }

  Optim <- fit_lpoly(dataList = dataList,
                     modelInfo = modelInfo,
                     method = method)

  #### Process the outputs ####

  transformed_pars <- fill_in(modelInfo$trans,
                              Optim$transparameters)

  parameters <- transformed_pars[names(modelInfo$param)]

  #### Standard errors ####

  Optim$SE <- compute_se_lpoly(dataList = dataList,
                               modelInfo = modelInfo,
                               Optim = Optim)

  #### latent object ####

  result <- new("latent",
                version          = as.character(packageVersion("latent")),
                call             = mc,
                timing           = Optim$elapsed,
                dataList         = dataList,
                modelInfo        = modelInfo,
                Optim            = Optim,
                parameters       = parameters,
                transformed_pars = transformed_pars,
                extra            = list())

  #### Result ####

  return(result)

}

#### Function to create the dataList ####

create_lpoly_dataList <- function(data, control) {

  #### Polychoric sample statistics ####

  polychorics <- polyfast(as.matrix(data), cores = 1L)

  S <- polychorics$correlation
  taus <- polychorics$thresholds
  n <- polychorics$contingency_tables

  item_label <- colnames(data)
  rownames(S) <- colnames(S) <- item_label
  names(taus) <- item_label

  nobs <- nrow(data)
  nitems <- ncol(data)
  npatterns <- length(unlist(n))

  K <- vapply(taus, FUN = function(x) {
    length(x)-2L
  }, FUN.VALUE = integer(1L))

  item_cat <- lapply(data, FUN = function(x) {
    sort(unique(x[!is.na(x)]))
  })

  #### Result ####

  result <- list(data = data,
                 nobs = nobs,
                 nitems = nitems,
                 npatterns = npatterns,
                 item_label = item_label,
                 n = n,
                 polychorics = polychorics,
                 S = S,
                 taus = taus,
                 K = K,
                 item_cat = item_cat)

  return(result)

}

#### Function to create the model ####

create_lpoly_model <- function(dataList, model = NULL, control) {

  # Generate the model syntax and initial parameter values.

  #### Parameter-block names ####

  data_param <- create_lpoly_data_param(dataList = dataList,
                                        control = control)

  #### Model for the transformed parameters ####

  trans <- model_lpoly(dataList = dataList,
                       data_param = data_param,
                       control = control)

  #### Model for the parameters ####

  param <- constraints_lpoly(dataList = dataList,
                             data_param = data_param,
                             trans = trans,
                             control = control)

  #### Fixed parameters and equality constraints ####

  model_constraints <- apply_model_lpoly(model = model,
                                         param = param,
                                         trans = trans,
                                         data_param = data_param,
                                         control = control)

  param <- model_constraints$param
  trans <- model_constraints$trans

  #### Create the initial values for the parameters ####

  init_param <- start_lpoly(dataList = dataList,
                            data_param = data_param,
                            param = param,
                            trans = trans,
                            control = control)

  #### Custom initial values ####

  init_param <- custom_init_param(control$start, init_param)

  #### Result ####

  result <- list(param = param,
                 trans = trans,
                 init_param = init_param,
                 data_param = data_param,
                 control = control)

  return(result)

}

#### Function to create the modelInfo ####

create_lpoly_modelInfo <- function(dataList, full_model, control) {

  # Generate the manifold, transformation, estimator, and optimizer structures.

  list2env(full_model, envir = environment())

  #### Manifolds ####

  manifolds <- manifolds_lpoly(dataList = dataList,
                               data_param = data_param,
                               param = param,
                               control = control)

  control_manifold <- create_manifolds(manifolds = manifolds,
                                       structures = param)

  #### Transformations ####

  transforms <- transformations_lpoly(dataList = dataList,
                                       data_param = data_param,
                                       trans = trans,
                                       control = control)

  control_transform <- create_transforms(transforms = transforms,
                                         structures = trans)

  #### Estimators ####

  estimators <- estimators_lpoly(dataList = dataList,
                                 data_param = data_param,
                                 trans = trans,
                                 control = control)

  control_estimator <- create_estimators(estimators = estimators,
                                         structures = trans)

  #### Pass the initial values to vectors ####

  inits <- create_init(trans, param, init_param,
                       control_transform = control_transform,
                       control = control)
  list2env(inits, envir = environment())

  #### Set up the optimizer ####

  control_optimizer <- control
  control_optimizer$parameters <- parameters
  control_optimizer$transparameters <- transparameters
  control_optimizer$init_param <- init_param
  control_optimizer$transparam2param <- trans2param-1L

  #### Result ####

  modelInfo <- list(param = param,
                    trans = trans,
                    nparam = nparam,
                    ntrans = ntrans,
                    parameters_labels = parameters_labels,
                    transparameters_labels = transparameters_labels,
                    dof = dataList$npatterns-nparam,
                    control_manifold = control_manifold,
                    control_transform = control_transform,
                    control_estimator = control_estimator,
                    control_optimizer = control_optimizer)

  return(modelInfo)

}

#### Function to create the control list of optimization parameters ####

lpoly_control <- function(control) {

  # Auxiliary function for lpoly.R.

  penalty_defaults <- list(
    logdet = list(w = 1e-03)
  )

  if(isFALSE(control$penalties)) {

    control$reg <- FALSE

  } else if(isTRUE(control$penalties)) {

    control$reg <- TRUE
    control$penalties <- penalty_defaults

  } else if(is.list(control$penalties)) {

    unknown <- setdiff(names(control$penalties), names(penalty_defaults))

    if(length(unknown) > 0L) {
      stop("Unknown penalty name(s): ", paste(unknown, collapse = ", "))
    }

    control$penalties <- utils::modifyList(penalty_defaults,
                                           control$penalties)

    if(!is.numeric(control$penalties$logdet$w) ||
       length(control$penalties$logdet$w) != 1L ||
       !is.finite(control$penalties$logdet$w) ||
       control$penalties$logdet$w <= 0) {
      stop("The logdet penalty weight must be a positive number")
    }

    control$reg <- TRUE

  } else {

    stop("penalties should be TRUE, FALSE, or a named list")

  }

  if(is.null(control$start)) {
    control$start <- NULL
  }

  if(is.null(control$subfix)) {
    control$subfix <- ""
  }

  if(is.null(control$opt)) {
    control$opt <- "lbfgs"
  }

  if(is.null(control$step_maxit)) {
    control$step_maxit <- 30L
  } else if(control$step_maxit < 1L) {
    stop("step_maxit must be an integer greater than 0")
  }

  if(is.null(control$c1)) {
    control$c1 <- 0.5
  } else if(control$c1 < 0) {
    stop("c1 must be a positive number")
  }

  if(is.null(control$c2)) {
    control$c2 <- 0.5
  } else if(control$c2 < 0) {
    stop("c2 must be a positive number")
  }

  if(is.null(control$step_eps)) {
    control$step_eps <- 1e-09
  } else if(control$step_eps < 0) {
    stop("step_eps must be a positive number")
  }

  if(is.null(control$df_eps)) {
    control$df_eps <- 1e-09
  } else if(control$df_eps < 0) {
    stop("df_eps must be a positive number")
  }

  if(is.null(control$M)) {
    control$M <- 100L
  } else if(control$M < 1L) {
    stop("M must be a positive integer")
  }

  if(is.null(control$eps)) {
    control$eps <- 1e-04
  } else if(control$eps < 0) {
    stop("eps must be a positive number")
  }

  if(is.null(control$ss_fac)) {
    control$ss_fac <- 2
  } else if(control$ss_fac <= 1) {
    stop("ss_fac must be larger than 1")
  }

  if(is.null(control$maxit)) {
    control$maxit <- 1000L
  } else if(control$maxit < 1L) {
    stop("maxit must be a positive integer")
  }

  if(is.null(control$tcg_maxit)) {
    control$tcg_maxit <- 10L
  } else if(control$tcg_maxit < 1L) {
    stop("tcg_maxit must be a positive integer")
  }

  if(is.null(control$ss)) {
    # Step sizes should be small so thresholds remain close to sensible bounds.
    control$ss <- 0.001
  } else if(control$ss <= 0) {
    stop("ss must be a positive number")
  }

  # lpoly currently uses a single start and a single optimizer core.
  control$rstarts <- 1L
  control$cores <- 1L

  #### Result ####

  return(control)

}

#### Auxiliary functions for create_lpoly_model ####

create_lpoly_data_param <- function(dataList, control) {

  taus_item <- paste0("taus", dataList$item_label, control$subfix)
  X_matrix <- paste0("X", control$subfix)
  S_matrix <- paste0("S", control$subfix)

  #### Result ####

  result <- list(taus_item = taus_item,
                 X_matrix = X_matrix,
                 S_matrix = S_matrix)

  return(result)

}

model_lpoly <- function(dataList, data_param, control) {

  list2env(data_param, envir = environment())

  list_struct <- list()
  k <- 1L

  #### Thresholds ####

  for(i in seq_len(dataList$nitems)) {

    list_struct[[k]] <- list(name = taus_item[i],
                             type = "matrix",
                             dim = c(dataList$K[i], 1L),
                             rownames = seq_len(dataList$K[i]),
                             colnames = dataList$item_label[i])
    k <- k+1L

  }

  #### Positive-semidefinite parameterization ####

  if(control$positive) {

    list_struct[[k]] <- list(name = X_matrix,
                             type = "matrix",
                             dim = c(dataList$nitems, dataList$nitems),
                             rownames = dataList$item_label,
                             colnames = dataList$item_label)
    k <- k+1L

  }

  #### Correlation matrix ####

  list_struct[[k]] <- list(name = S_matrix,
                           type = "matrix",
                           dim = c(dataList$nitems, dataList$nitems),
                           rownames = dataList$item_label,
                           colnames = dataList$item_label,
                           symmetric = TRUE)

  #### Result ####

  result <- create_parameters(list_struct)

  return(result)

}

constraints_lpoly <- function(dataList, data_param, trans, control) {

  list2env(data_param, envir = environment())

  param <- trans

  # A correlation matrix has unit diagonal.
  diag(param[[S_matrix]]) <- "1"

  #### Result ####

  return(param)

}

apply_model_lpoly <- function(model, param, trans, data_param, control) {

  list2env(data_param, envir = environment())

  #### User-specified parameter blocks ####

  if(!is.null(model)) {

    unknown <- setdiff(names(model), names(param))

    if(length(unknown) > 0L) {
      stop("Unknown model parameter block(s): ",
           paste(unknown, collapse = ", "))
    }

    if(control$positive) {

      constrained_matrix <- intersect(names(model), c(X_matrix, S_matrix))
      constrained_matrix <-
        constrained_matrix[!vapply(model[constrained_matrix],
                                   is.null, logical(1L))]

      if(length(constrained_matrix) > 0L) {
        stop("Constraints on X or S are not currently supported when positive = TRUE")
      }

    }

    nm <- intersect(names(model), names(param))
    nm <- nm[!vapply(model[nm], is.null, logical(1L))]

    for(j in nm) {

      if(j == S_matrix) {

        param[[j]] <- insert_symmetric_model_lpoly(param[[j]],
                                                   model[[j]],
                                                   object_name = j)

      } else {

        param[[j]] <- insert_partial_object(param[[j]], model[[j]],
                                            object_name = j)

      }

    }

  }

  #### Identification constraints ####

  # The diagonal of a correlation matrix must remain fixed to one.
  diag(param[[S_matrix]]) <- "1"

  #### Equality constraints ####

  # Match lca(): equality labels must be copied into trans so create_init()
  # sees repeated labels as one free parameter. Numeric values remain only in
  # param and are therefore treated as fixed values.
  for(nm in names(param)) {

    x <- param[[nm]]

    keep <- !is.na(x) & grepl("[A-Za-z]", x)

    if(any(keep)) {
      trans[[nm]][keep] <- x[keep]
    }

  }

  #### Result ####

  result <- list(param = param,
                 trans = trans)

  return(result)

}

insert_symmetric_model_lpoly <- function(X, Y, object_name = NULL) {

  if(is.null(Y)) {

    #### Result ####

    return(X)

  }

  if(is.data.frame(Y)) {
    Y <- as.matrix(Y)
  }

  if(is.null(dim(Y)) || length(dim(Y)) != 2L) {
    stop("Object '", object_name, "' must be supplied as a matrix")
  }

  #### Align the supplied values ####

  supplied <- matrix(NA_character_,
                     nrow = nrow(X),
                     ncol = ncol(X),
                     dimnames = dimnames(X))

  supplied <- insert_partial_object(supplied, Y,
                                    object_name = object_name)

  #### Enforce symmetric constraints ####

  for(j in seq_len(ncol(X))) {

    for(i in seq_len(j-1L)) {

      upper <- supplied[i, j]
      lower <- supplied[j, i]

      upper_set <- !is.na(upper)
      lower_set <- !is.na(lower)

      if(upper_set && lower_set &&
         !identical(as.character(upper), as.character(lower))) {
        stop("Conflicting symmetric constraints in '", object_name,
             "' for [", rownames(X)[i], ", ", colnames(X)[j],
             "] and [", rownames(X)[j], ", ", colnames(X)[i], "]")
      }

      if(upper_set || lower_set) {

        value <- if(lower_set) lower else upper

        X[i, j] <- value
        X[j, i] <- value

      }

    }

  }

  #### Diagonal ####

  diagonal <- diag(supplied)
  replace <- !is.na(diagonal)

  if(any(replace)) {
    idx <- which(replace)
    X[cbind(idx, idx)] <- diagonal[replace]
  }

  #### Result ####

  return(X)

}

start_lpoly <- function(dataList, data_param, param, trans, control) {

  list2env(data_param, envir = environment())

  threslds <- lapply(dataList$taus, FUN = function(x) {
    x[!is.infinite(x)]
  })

  init_param <- vector("list", length = control$rstarts)

  for(rs in seq_len(control$rstarts)) {

    init_param[[rs]] <- list()

    #### Thresholds ####

    for(i in seq_len(dataList$nitems)) {

      init_param[[rs]][[taus_item[i]]] <-
        matrix(threslds[[i]], ncol = 1L,
               dimnames = dimnames(trans[[taus_item[i]]]))

    }

    #### Correlation matrix ####

    if(control$positive) {

      init_param[[rs]][[X_matrix]] <- real_sqrtmat(dataList$S)
      dimnames(init_param[[rs]][[X_matrix]]) <- dimnames(trans[[X_matrix]])

      init_param[[rs]][[S_matrix]] <-
        crossprod(init_param[[rs]][[X_matrix]])
      dimnames(init_param[[rs]][[S_matrix]]) <- dimnames(trans[[S_matrix]])

    } else {

      init_param[[rs]][[S_matrix]] <- dataList$S
      dimnames(init_param[[rs]][[S_matrix]]) <- dimnames(trans[[S_matrix]])

    }

  }

  #### Result ####

  return(init_param)

}

#### Auxiliary functions for create_lpoly_modelInfo ####

manifolds_lpoly <- function(dataList, data_param, param, control) {

  list2env(data_param, envir = environment())

  if(control$positive) {

    manifolds <- list(
      list(manifold = "euclidean",
           parameters = list(param[taus_item])),
      list(manifold = "oblq",
           parameters = X_matrix,
           extra = list(p = dataList$nitems,
                        q = dataList$nitems))
    )

  } else {

    manifolds <- list(
      list(manifold = "euclidean",
           parameters = list(param[taus_item])),
      list(manifold = "euclidean",
           parameters = S_matrix)
    )

  }

  #### Result ####

  return(manifolds)

}

transformations_lpoly <- function(dataList, data_param, trans, control) {

  list2env(data_param, envir = environment())

  transforms <- list()

  if(control$positive) {

    transforms[[1L]] <- list(transform = "crossprod",
                             parameters_in = X_matrix,
                             parameters_out = S_matrix,
                             extra = list(p = dataList$nitems))

  }

  #### Result ####

  return(transforms)

}

estimators_lpoly <- function(dataList, data_param, trans, control) {

  list2env(data_param, envir = environment())

  estimators <- list()
  k <- 1L

  estimators[[k]] <- list(estimator = "polycor",
                          parameters = c(list(trans[[S_matrix]]),
                                         trans[taus_item]),
                          extra = list(n = dataList$n,
                                       p = dataList$nitems,
                                       N = dataList$nobs))
  k <- k+1L

  if(control$reg) {

    lower_indices <- which(lower.tri(trans[[S_matrix]], diag = TRUE))

    estimators[[k]] <- list(estimator = "logdetR",
                            parameters = S_matrix,
                            extra = list(lower_indices = lower_indices-1L,
                                         p = dataList$nitems,
                                         logdetw = control$penalties$logdet$w))

  }

  #### Result ####

  return(estimators)

}

#### Auxiliary functions for fitting and post-estimation ####

fit_lpoly <- function(dataList, modelInfo, method) {

  if(method == "one-step") {

    control_optimizer <- modelInfo$control_optimizer
    control_optimizer$cores <- min(control_optimizer$rstarts,
                                   control_optimizer$cores)

    Optim <- optimizer(control_manifold = modelInfo$control_manifold,
                       control_transform = modelInfo$control_transform,
                       control_estimator = modelInfo$control_estimator,
                       control_optimizer = control_optimizer)

  } else {

    cores <- parallel::detectCores()
    if(is.na(cores) || cores < 1L) cores <- 1L

    Optim <- polyfast(as.matrix(dataList$data), cores = cores)

    rownames(Optim$correlation) <- colnames(Optim$correlation) <-
      dataList$item_label
    names(Optim$thresholds) <- dataList$item_label

    threslds <- lapply(Optim$thresholds, FUN = function(x) {
      matrix(x[!is.infinite(x)], ncol = 1L)
    })

    Optim$parameters <- c(unlist(threslds, use.names = FALSE),
                          Optim$correlation[lower.tri(Optim$correlation,
                                                      diag = FALSE)])

    Optim$transparameters <- c(unlist(threslds, use.names = FALSE),
                               Optim$correlation[lower.tri(Optim$correlation,
                                                           diag = TRUE)])

    Optim$f <- 0

  }

  names(Optim$parameters) <- modelInfo$parameters_labels
  names(Optim$transparameters) <- modelInfo$transparameters_labels

  if(!is.null(Optim$g) && length(Optim$g) == length(modelInfo$parameters_labels)) {
    names(Optim$g) <- modelInfo$parameters_labels
  }

  if(!is.null(Optim$rg) && length(Optim$rg) == length(modelInfo$parameters_labels)) {
    names(Optim$rg) <- modelInfo$parameters_labels
  }

  if(!is.null(Optim$dir) && length(Optim$dir) == length(modelInfo$parameters_labels)) {
    names(Optim$dir) <- modelInfo$parameters_labels
  }

  #### Result ####

  return(Optim)

}

compute_se_lpoly <- function(dataList, modelInfo, Optim) {

  #### Hessian ####

  modelInfo$control_optimizer$parameters[[1L]] <- Optim$parameters
  modelInfo$control_optimizer$transparameters[[1L]] <-
    Optim$transparameters

  H <- get_hess(control_manifold = modelInfo$control_manifold,
                control_transform = modelInfo$control_transform,
                control_estimator = modelInfo$control_estimator,
                control_optimizer = modelInfo$control_optimizer,
                cores = 1L)$h

  rownames(H) <- colnames(H) <- modelInfo$parameters_labels

  #### Variance-covariance matrix ####

  # The polychoric objective is averaged over observations. Therefore, the
  # inverse Hessian is the N-scaled asymptotic covariance and must be divided by
  # the sample size to obtain the sampling variance-covariance matrix.
  VCOV <- approx_Hinv(H)/dataList$nobs
  VCOV <- (VCOV+t(VCOV))/2
  rownames(VCOV) <- colnames(VCOV) <- modelInfo$parameters_labels

  se <- sqrt(diag(VCOV))
  names(se) <- modelInfo$parameters_labels

  #### Result ####

  result <- list(H = H,
                 VCOV = VCOV,
                 se = se)

  return(result)

}

print_lpoly_message <- function(msg) {

  w <- nchar(msg)+4L

  cat("\n", "+", strrep("-", w), "+\n",
      "|  ", msg, "  |\n",
      "+", strrep("-", w), "+\n\n", sep = "")

  #### Result ####

  return(invisible(NULL))

}
