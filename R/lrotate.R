# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 13/08/2026
#'
#' @title
#' Rotate a factor loading matrix
#'
#' @description
#' \code{lrotate} rotates one or more factor loading matrices using an
#' orthogonal or oblique projection and a selected rotation criterion.
#'
#' @usage
#' lrotate(lambda, projection = "oblq", rotation = "oblimin",
#'         group = NULL, positive = FALSE, penalties = TRUE,
#'         do.fit = TRUE, control = NULL, ...)
#'
#' @param lambda A matrix or a list of loading matrices, one for each group.
#' @param projection Character string. Available projections are
#'   \code{"orth"}, \code{"oblq"}, and \code{"poblq"}.
#' @param rotation Character string identifying the rotation criterion.
#' @param group Optional grouping information retained in the fitted object.
#' @param positive Logical. Retained for compatibility with the factor-analysis
#'   interface.
#' @param penalties Logical or list of penalty settings retained in the
#'   optimization control.
#' @param do.fit Logical. If \code{TRUE}, fit the rotation. If \code{FALSE},
#'   return only the model specification.
#' @param control List of optimization-control arguments.
#' @param ... Additional arguments required by the selected projection or
#'   rotation criterion.
#'
#' @return An object of class \code{"latent"}.
#'
#' @examples
#' \dontrun{
#' fit <- lrotate(lambda = list(lambda),
#'                projection = "oblq",
#'                rotation = "oblimin")
#' }
#'
#' @export
lrotate <- function(lambda, projection = "oblq", rotation = "oblimin",
                    group = NULL, positive = FALSE, penalties = TRUE,
                    do.fit = TRUE, control = NULL, ...) {

  #### Check input arguments ####

  lambda <- check_lambda_lrotate(lambda)

  projection <- tolower(projection)
  rotation <- tolower(rotation)

  supported_projection <- c("orth", "oblq", "poblq")

  if(!(projection %in% supported_projection)) {
    stop("Unknown projection: ", projection)
  }

  supported_rotation <- c("cf", "geomin", "lclf", "oblimin",
                          "target", "varimax", "varimin", "xtarget")

  if(!(rotation %in% supported_rotation)) {
    stop("Unknown rotation criterion: ", rotation)
  }

  if(length(positive) != 1L || !is.logical(positive) || is.na(positive)) {
    stop("positive must be TRUE or FALSE")
  }

  if(length(do.fit) != 1L || !is.logical(do.fit) || is.na(do.fit)) {
    stop("do.fit must be TRUE or FALSE")
  }

  #### Store original call ####

  mc <- match.call()

  #### Rotation-specific arguments ####

  dots <- list(...)
  dots <- rotation_defaults_lrotate(rotation = rotation,
                                    dots = dots)

  #### Check control parameters ####

  control$penalties <- penalties
  control$positive <- positive
  control$estimator <- rotation
  control$projection <- projection
  control <- lrotate_control(control)

  #### Create the dataList ####

  dataList <- create_lrotate_dataList(lambda = lambda,
                                      projection = projection,
                                      rotation = rotation,
                                      group = group,
                                      positive = positive)

  #### Create the model ####

  full_model <- create_lrotate_model(dataList = dataList,
                                     control = control)

  #### Create the manifold, transformation, and estimator structures ####

  modelInfo <- create_lrotate_modelInfo(dataList = dataList,
                                        full_model = full_model,
                                        control = control,
                                        dots = dots)

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

  Optim <- fit_lrotate(modelInfo = modelInfo)

  #### Process the outputs ####

  transformed_pars <- fill_in(modelInfo$trans,
                              Optim$transparameters)

  parameters <- transformed_pars[names(modelInfo$param)]

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

#### Function to check the loading matrices ####

check_lambda_lrotate <- function(lambda) {

  if(is.matrix(lambda)) {
    lambda <- list(lambda)
  }

  if(!is.list(lambda) || length(lambda) == 0L) {
    stop("lambda must be a matrix or a non-empty list of matrices")
  }

  ngroups <- length(lambda)

  group_label <- names(lambda)

  if(is.null(group_label)) {

    group_label <- paste("group", seq_len(ngroups), sep = "")

  } else {

    empty_names <- is.na(group_label) | group_label == ""

    if(any(empty_names)) {
      group_label[empty_names] <-
        paste("group", which(empty_names), sep = "")
    }

    if(anyDuplicated(group_label)) {
      stop("lambda must have unique group names")
    }

  }

  for(i in seq_len(ngroups)) {

    if(is.data.frame(lambda[[i]])) {
      lambda[[i]] <- as.matrix(lambda[[i]])
    }

    if(!is.matrix(lambda[[i]]) || !is.numeric(lambda[[i]])) {
      stop("Every element of lambda must be a numeric matrix")
    }

    if(nrow(lambda[[i]]) < 1L || ncol(lambda[[i]]) < 1L) {
      stop("Every loading matrix must contain at least one row and one column")
    }

    if(anyNA(lambda[[i]]) || any(!is.finite(lambda[[i]]))) {
      stop("Loading matrices cannot contain missing or non-finite values")
    }

    item_label <- rownames(lambda[[i]])
    factor_label <- colnames(lambda[[i]])

    if(is.null(item_label)) {
      item_label <- paste("item", seq_len(nrow(lambda[[i]])), sep = "")
    }

    if(is.null(factor_label)) {
      factor_label <- paste("factor", seq_len(ncol(lambda[[i]])), sep = "")
    }

    if(any(item_label == "") || anyDuplicated(item_label)) {
      stop("Loading-matrix row names must be unique and non-empty")
    }

    if(any(factor_label == "") || anyDuplicated(factor_label)) {
      stop("Loading-matrix column names must be unique and non-empty")
    }

    rownames(lambda[[i]]) <- item_label
    colnames(lambda[[i]]) <- factor_label

  }

  names(lambda) <- group_label

  #### Result ####

  return(lambda)

}

#### Function to create defaults for rotation-specific arguments ####

rotation_defaults_lrotate <- function(rotation, dots) {

  # These defaults make the two commonly used criteria directly usable from
  # lrotate()/lefa() without requiring an otherwise undocumented argument.

  if(rotation == "oblimin" && is.null(dots$gamma)) {
    dots$gamma <- 0
  }

  if(rotation == "geomin" && is.null(dots$epsilon)) {
    dots$epsilon <- 0.01
  }

  #### Result ####

  return(dots)

}

#### Function to create the dataList ####

create_lrotate_dataList <- function(lambda, projection, rotation,
                                    group, positive) {

  ngroups <- length(lambda)

  item_label <- lapply(lambda, rownames)
  factor_label <- lapply(lambda, colnames)
  nitems <- lapply(lambda, nrow)
  nfactors <- lapply(lambda, ncol)

  orthogonal <- projection == "orth"

  dataList <- list(
    data = lambda,
    ngroups = ngroups,
    group = group,
    group_label = names(lambda),
    item_label = item_label,
    factor_label = factor_label,
    nitems = nitems,
    nfactors = nfactors,
    projection = projection,
    rotation = rotation,
    estimator = rotation,
    positive = positive,
    orthogonal = orthogonal
  )

  #### Result ####

  return(dataList)

}

#### Function to create the model ####

create_lrotate_model <- function(dataList, control) {

  data_param <- create_lrotate_data_param(dataList = dataList)

  #### Model for the transformed parameters ####

  trans <- model_lrotate(dataList = dataList,
                         data_param = data_param)

  #### Model for the parameters ####

  param <- constraints_lrotate(trans = trans,
                               dataList = dataList,
                               data_param = data_param)

  #### Create the initial values for the parameters ####

  init_param <- start_lrotate(trans = trans,
                              dataList = dataList,
                              data_param = data_param,
                              control = control)

  #### Custom initial values ####

  init_param <- custom_init_param(control$start, init_param)

  # Recompute every transformed quantity from X after custom starts have been
  # inserted. The unrotated loadings remain fixed at the supplied values.
  init_param <- refresh_start_lrotate(init_param = init_param,
                                      trans = trans,
                                      dataList = dataList,
                                      data_param = data_param)

  #### Result ####

  result <- list(param = param,
                 trans = trans,
                 init_param = init_param,
                 data_param = data_param)

  return(result)

}

#### Function to create the parameter-block names ####

create_lrotate_data_param <- function(dataList) {

  groups <- seq_len(dataList$ngroups)

  result <- list(
    ulambda_group = paste("ulambda.group", groups, sep = ""),
    X_group = paste("X.group", groups, sep = ""),
    Xinv_group = paste("Xinv.group", groups, sep = ""),
    lambda_group = paste("lambda.group", groups, sep = ""),
    psi_group = paste("psi.group", groups, sep = "")
  )

  #### Result ####

  return(result)

}

#### Function to create matrix parameter labels ####

labels_lrotate <- function(group, object, nrow, ncol) {

  labels <- paste("g", group, ".", object, "[",
                  rep(seq_len(nrow), times = ncol), ",",
                  rep(seq_len(ncol), each = nrow), "]",
                  sep = "")

  result <- matrix(labels, nrow = nrow, ncol = ncol)

  #### Result ####

  return(result)

}

#### Function to create the transformed-parameter model ####

model_lrotate <- function(dataList, data_param) {

  list2env(data_param, envir = environment())

  list_struct <- list()
  k <- 1L

  for(i in seq_len(dataList$ngroups)) {

    p <- dataList$nitems[[i]]
    q <- dataList$nfactors[[i]]
    item_names <- dataList$item_label[[i]]
    factor_names <- dataList$factor_label[[i]]

    #### Unrotated loadings ####

    list_struct[[k]] <- list(
      name = ulambda_group[i],
      type = "matrix",
      dim = c(p, q),
      rownames = item_names,
      colnames = factor_names,
      labels = labels_lrotate(i, "ulambda", p, q)
    )
    k <- k+1L

    #### Rotation matrix ####

    list_struct[[k]] <- list(
      name = X_group[i],
      type = "matrix",
      dim = c(q, q),
      rownames = factor_names,
      colnames = factor_names,
      labels = labels_lrotate(i, "X", q, q)
    )
    k <- k+1L

    #### Inverse rotation matrix ####

    if(!dataList$orthogonal) {

      list_struct[[k]] <- list(
        name = Xinv_group[i],
        type = "matrix",
        dim = c(q, q),
        rownames = factor_names,
        colnames = factor_names,
        labels = labels_lrotate(i, "Xinv", q, q)
      )
      k <- k+1L

    }

    #### Rotated loadings ####

    list_struct[[k]] <- list(
      name = lambda_group[i],
      type = "matrix",
      dim = c(p, q),
      rownames = item_names,
      colnames = factor_names,
      labels = labels_lrotate(i, "lambda", p, q)
    )
    k <- k+1L

    #### Latent correlations ####

    list_struct[[k]] <- list(
      name = psi_group[i],
      type = "matrix",
      dim = c(q, q),
      rownames = factor_names,
      colnames = factor_names,
      symmetric = TRUE,
      labels = labels_lrotate(i, "psi", q, q)
    )
    k <- k+1L

  }

  trans <- create_parameters(list_struct)

  #### Result ####

  return(trans)

}

#### Function to create the parameter constraints ####

constraints_lrotate <- function(trans, dataList, data_param) {

  list2env(data_param, envir = environment())

  param <- list()

  for(i in seq_len(dataList$ngroups)) {

    # The unrotated loading matrix is fixed at the supplied solution.
    param[[ulambda_group[i]]] <- dataList$data[[i]]
    dimnames(param[[ulambda_group[i]]]) <-
      dimnames(trans[[ulambda_group[i]]])

    # X is the only free parameter block for the rotation.
    param[[X_group[i]]] <- trans[[X_group[i]]]

  }

  #### Result ####

  return(param)

}

#### Function to create starting values ####

start_lrotate <- function(trans, dataList, data_param, control) {

  list2env(data_param, envir = environment())

  init_param <- vector("list", length = control$rstarts)

  for(rs in seq_len(control$rstarts)) {

    init_param[[rs]] <- list()

    for(i in seq_len(dataList$ngroups)) {

      q <- dataList$nfactors[[i]]

      #### Fixed unrotated loadings ####

      init_param[[rs]][[ulambda_group[i]]] <- dataList$data[[i]]
      dimnames(init_param[[rs]][[ulambda_group[i]]]) <-
        dimnames(trans[[ulambda_group[i]]])

      #### Rotation matrix ####

      X <- rorth(q, q)
      dimnames(X) <- dimnames(trans[[X_group[i]]])
      init_param[[rs]][[X_group[i]]] <- X

      #### Derived transformed quantities ####

      init_param[[rs]] <- update_rotation_lrotate(
        x = init_param[[rs]],
        group_index = i,
        trans = trans,
        dataList = dataList,
        data_param = data_param
      )

    }

  }

  #### Result ####

  return(init_param)

}

#### Function to refresh starting transformed quantities ####

refresh_start_lrotate <- function(init_param, trans,
                                  dataList, data_param) {

  list2env(data_param, envir = environment())

  for(rs in seq_along(init_param)) {

    for(i in seq_len(dataList$ngroups)) {

      # The unrotated loadings are data, not free starting values.
      init_param[[rs]][[ulambda_group[i]]] <- dataList$data[[i]]
      dimnames(init_param[[rs]][[ulambda_group[i]]]) <-
        dimnames(trans[[ulambda_group[i]]])

      init_param[[rs]] <- update_rotation_lrotate(
        x = init_param[[rs]],
        group_index = i,
        trans = trans,
        dataList = dataList,
        data_param = data_param
      )

    }

  }

  #### Result ####

  return(init_param)

}

#### Function to update quantities implied by a rotation matrix ####

update_rotation_lrotate <- function(x, group_index, trans,
                                    dataList, data_param) {

  list2env(data_param, envir = environment())

  i <- group_index
  U <- dataList$data[[i]]
  X <- x[[X_group[i]]]

  if(is.null(X) ||
     !identical(dim(X), c(dataList$nfactors[[i]],
                          dataList$nfactors[[i]]))) {
    stop("Invalid starting value for '", X_group[i], "'")
  }

  dimnames(X) <- dimnames(trans[[X_group[i]]])
  x[[X_group[i]]] <- X

  if(dataList$orthogonal) {

    rotated_lambda <- U %*% X

  } else {

    Xinv <- solve(X)
    dimnames(Xinv) <- dimnames(trans[[Xinv_group[i]]])
    x[[Xinv_group[i]]] <- Xinv

    rotated_lambda <- U %*% t(Xinv)

  }

  dimnames(rotated_lambda) <- dimnames(trans[[lambda_group[i]]])
  x[[lambda_group[i]]] <- rotated_lambda

  psi <- crossprod(X)
  dimnames(psi) <- dimnames(trans[[psi_group[i]]])
  x[[psi_group[i]]] <- psi

  #### Result ####

  return(x)

}

#### Function to create the modelInfo ####

create_lrotate_modelInfo <- function(dataList, full_model,
                                     control, dots) {

  list2env(full_model, envir = environment())

  #### Manifolds ####

  manifolds <- manifolds_lrotate(dataList = dataList,
                                 data_param = data_param,
                                 dots = dots)

  control_manifold <- create_manifolds(manifolds = manifolds,
                                       structures = param)

  #### Transformations ####

  transforms <- transformations_lrotate(dataList = dataList,
                                        data_param = data_param)

  control_transform <- create_transforms(transforms = transforms,
                                         structures = trans)

  #### Estimators ####

  estimators <- estimators_lrotate(dataList = dataList,
                                   data_param = data_param,
                                   dots = dots)

  control_estimator <- create_estimators(estimators = estimators,
                                         structures = trans)

  #### Pass the initial values to vectors ####

  inits <- create_init(trans, param, init_param,
                       control_transform = control_transform, control)

  list2env(inits, envir = environment())

  #### Set up the optimizer ####

  control_optimizer <- control
  control_optimizer$parameters <- parameters
  control_optimizer$transparameters <- transparameters
  control_optimizer$init_param <- init_param
  control_optimizer$transparam2param <- trans2param-1L

  #### Result ####

  modelInfo <- list(
    param = param,
    trans = trans,
    nparam = nparam,
    ntrans = ntrans,
    parameters_labels = parameters_labels,
    transparameters_labels = transparameters_labels,
    dof = NA_integer_,
    rotation = dataList$rotation,
    projection = dataList$projection,
    data_param = data_param,
    control_manifold = control_manifold,
    control_transform = control_transform,
    control_estimator = control_estimator,
    control_optimizer = control_optimizer
  )

  return(modelInfo)

}

#### Function to create the manifolds ####

manifolds_lrotate <- function(dataList, data_param, dots) {

  X_group <- data_param$X_group

  manifolds <- vector("list", length = dataList$ngroups)

  for(i in seq_len(dataList$ngroups)) {

    q <- dataList$nfactors[[i]]
    extra <- group_dots_lrotate(dots = dots,
                                group_index = i,
                                ngroups = dataList$ngroups)

    extra$p <- q
    extra$q <- q

    manifolds[[i]] <- list(
      manifold = dataList$projection,
      parameters = X_group[i],
      extra = extra
    )

  }

  #### Result ####

  return(manifolds)

}

#### Function to create the transformations ####

transformations_lrotate <- function(dataList, data_param) {

  list2env(data_param, envir = environment())

  transforms <- list()
  k <- 1L

  for(i in seq_len(dataList$ngroups)) {

    p <- dataList$nitems[[i]]
    q <- dataList$nfactors[[i]]

    if(dataList$orthogonal) {

      #### Rotated loadings ####

      transforms[[k]] <- list(
        transform = "XY",
        parameters_in = c(ulambda_group[i], X_group[i]),
        parameters_out = lambda_group[i],
        extra = list(p = p, q = q)
      )
      k <- k+1L

    } else {

      #### Inverse rotation matrix ####

      transforms[[k]] <- list(
        transform = "matrix_inverse",
        parameters_in = X_group[i],
        parameters_out = Xinv_group[i],
        extra = list(p = q)
      )
      k <- k+1L

      #### Rotated loadings ####

      transforms[[k]] <- list(
        transform = "XYt",
        parameters_in = c(ulambda_group[i], Xinv_group[i]),
        parameters_out = lambda_group[i],
        extra = list(p = p, q = q)
      )
      k <- k+1L

    }

    #### Latent correlations ####

    transforms[[k]] <- list(
      transform = "crossprod",
      parameters_in = X_group[i],
      parameters_out = psi_group[i],
      extra = list(p = q)
    )
    k <- k+1L

  }

  #### Result ####

  return(transforms)

}

#### Function to create the estimators ####

estimators_lrotate <- function(dataList, data_param, dots) {

  lambda_group <- data_param$lambda_group
  psi_group <- data_param$psi_group

  estimators <- vector("list", length = dataList$ngroups)

  for(i in seq_len(dataList$ngroups)) {

    extra <- group_dots_lrotate(dots = dots,
                                group_index = i,
                                ngroups = dataList$ngroups)

    extra$p <- dataList$nitems[[i]]
    extra$q <- dataList$nfactors[[i]]

    estimators[[i]] <- list(
      estimator = dataList$rotation,
      parameters = c(lambda_group[i], psi_group[i]),
      extra = extra
    )

  }

  #### Result ####

  return(estimators)

}

#### Function to select group-specific extra arguments ####

group_dots_lrotate <- function(dots, group_index, ngroups) {

  extra <- dots

  # Matrix-valued projection/criterion arguments may be supplied either once
  # for every group or as a list with one object per group.
  group_objects <- c("constraints", "target", "weight",
                     "psitarget", "psiweight")

  for(nm in intersect(names(extra), group_objects)) {

    object <- extra[[nm]]

    if(is.list(object) &&
       !is.data.frame(object) &&
       length(object) == ngroups) {
      extra[[nm]] <- object[[group_index]]
    }

  }

  #### Result ####

  return(extra)

}

#### Function to fit the rotation ####

fit_lrotate <- function(modelInfo) {

  control_optimizer <- modelInfo$control_optimizer

  control_optimizer$cores <-
    min(control_optimizer$rstarts,
        control_optimizer$cores)

  Optim <- optimizer(
    control_manifold = modelInfo$control_manifold,
    control_transform = modelInfo$control_transform,
    control_estimator = modelInfo$control_estimator,
    control_optimizer = control_optimizer
  )

  names(Optim$parameters) <- modelInfo$parameters_labels
  names(Optim$transparameters) <- modelInfo$transparameters_labels

  #### Result ####

  return(Optim)

}

#### Function to create the control list ####

lrotate_control <- function(control) {

  # Keep the optimizer defaults used by the CFA/EFA machinery. Rotation is
  # normally more stable with the Newton optimizer, matching the previous
  # implementation.

  control <- lcfa_control(control)

  if(control$opt == "lbfgs") {
    control$opt <- "newton"
  }

  if(is.null(control$start)) {
    control$start <- NULL
  }

  #### Result ####

  return(control)

}
