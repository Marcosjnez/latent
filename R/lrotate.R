# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 23/08/2026
#'
#' @title
#' Rotate a fitted CFA model
#'
#' @description
#' \code{lrotate} rotates the factor loading matrices of a fitted \code{lcfa}
#' object using an orthogonal or oblique projection and a selected rotation
#' criterion.
#'
#' @usage
#' lrotate(fit, projection = "oblq", rotation = "oblimin",
#'         do.fit = TRUE, control = NULL, ...)
#'
#' @param fit A fitted object inheriting from class \code{"lcfa"}.
#' @param projection Character string. Available projections are
#'   \code{"orth"}, \code{"oblq"}, and \code{"poblq"}.
#' @param rotation Character string identifying the rotation criterion.
#' @param do.fit Logical. If \code{TRUE}, fit the rotation. If \code{FALSE},
#'   return the unrestricted rotation specification used for derivative
#'   calculations.
#' @param control List of optimization-control arguments.
#' @param ... Additional arguments required by the selected projection or
#'   rotation criterion.
#'
#' @details
#' Let \eqn{X} be the rotation matrix and let \eqn{\Lambda_0}, \eqn{\Psi_0},
#' and \eqn{\alpha_0} denote the unrotated factor loadings, factor covariance
#' matrix, and factor means. The rotated quantities are
#' \deqn{\Lambda_r=\Lambda_0X^{-T},}
#' \deqn{\Psi_r=X^T\Psi_0X,}
#' and
#' \deqn{\alpha_r=X^T\alpha_0.}
#' For an orthogonal projection, \eqn{X^{-T}=X}. If \eqn{\Psi_0} is a fixed
#' identity matrix, \eqn{\Psi_r} is computed as \eqn{X^TX}.
#'
#' @return An object of class \code{"multistep"}. The fitted \code{lcfa}
#' object is stored in \code{extra} as the preceding estimation step.
#'
#' @examples
#' \dontrun{
#' fit_cfa <- lcfa(data = HolzingerSwineford1939,
#'                 model = model,
#'                 std.lv = TRUE)
#' fit_rotation <- lrotate(fit = fit_cfa,
#'                         projection = "oblq",
#'                         rotation = "oblimin")
#' }
#'
#' @export
lrotate <- function(fit, projection = "oblq", rotation = "oblimin",
                    do.fit = TRUE, control = NULL, ...) {

  #### Check input arguments ####

  check_fit_lrotate(fit)

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

  if(length(do.fit) != 1L || !is.logical(do.fit) || is.na(do.fit)) {
    stop("do.fit must be TRUE or FALSE")
  }

  if(!is.null(control) && !is.list(control)) {
    stop("control must be NULL or a list")
  }

  if(is.null(control)) {
    control <- list()
  }

  #### Store original call ####

  mc <- match.call()

  #### Rotation-specific arguments ####

  dots <- list(...)
  dots <- rotation_defaults_lrotate(rotation = rotation,
                                    dots = dots)

  #### Check control parameters ####

  control$penalties <- FALSE
  control$positive <- FALSE
  control$estimator <- rotation
  control$projection <- projection
  control <- lrotate_control(control)
  control$free_previous <- !do.fit

  #### Create the dataList ####

  dataList <- create_lrotate_dataList(fit = fit,
                                      projection = projection,
                                      rotation = rotation)

  dataList$args <- c(list(fit = fit,
                          projection = projection,
                          rotation = rotation,
                          do.fit = do.fit,
                          control = control),
                     dots)

  #### Create the model ####

  full_model <- create_lrotate_model(dataList = dataList,
                                     control = control)

  #### Create the manifold, transformation, and estimator structures ####

  modelInfo <- create_lrotate_modelInfo(dataList = dataList,
                                        full_model = full_model,
                                        control = control,
                                        dots = dots)

  modelInfo$propagate_uncertainty <- TRUE
  modelInfo$step_labels <- modelInfo$parameters_labels

  #### Fit the model ####

  if(!do.fit) {

    result <- new("multistep",
                  version          = as.character(packageVersion("latent")),
                  call             = mc,
                  timing           = numeric(),
                  dataList         = dataList,
                  modelInfo        = modelInfo,
                  Optim            = list(),
                  parameters       = list(),
                  transformed_pars = list(),
                  extra            = list(fit))

    #### Result ####

    return(result)

  }

  Optim <- fit_lrotate(modelInfo = modelInfo)

  #### Process the outputs ####

  transformed_pars <- fill_in(modelInfo$trans,
                              Optim$transparameters)

  parameters <- transformed_pars[names(modelInfo$param)]

  #### latent object ####

  result <- new("multistep",
                version          = as.character(packageVersion("latent")),
                call             = mc,
                timing           = Optim$elapsed,
                dataList         = dataList,
                modelInfo        = modelInfo,
                Optim            = Optim,
                parameters       = parameters,
                transformed_pars = transformed_pars,
                extra            = list(fit))

  #### Result ####

  return(result)

}

#### Function to check the fitted CFA object ####

check_fit_lrotate <- function(fit) {

  if(!inherits(fit, "lcfa")) {
    stop("fit must be a fitted object inheriting from class 'lcfa'")
  }

  if(length(fit@Optim$transparameters) == 0L ||
     length(fit@transformed_pars) == 0L) {
    stop("fit must contain a fitted lcfa model")
  }

  data_param <- fit@dataList$data_param
  required <- c("lambda_group", "psi_group", "alpha_group")

  if(is.null(data_param) ||
     !all(required %in% names(data_param))) {
    stop("fit does not contain the CFA parameter-block information required for rotation")
  }

  #### Result ####

  return(invisible(NULL))

}

#### Function to identify a fixed identity factor covariance matrix ####

identity_psi_lrotate <- function(fit, psi_name, psi) {

  if(!is.matrix(psi) || nrow(psi) != ncol(psi)) {

    #### Result ####

    return(FALSE)

  }

  q <- nrow(psi)
  identity <- diag(q)
  dimnames(identity) <- dimnames(psi)

  if(!isTRUE(all.equal(psi, identity,
                       tolerance = sqrt(.Machine$double.eps)))) {

    #### Result ####

    return(FALSE)

  }

  if(!(psi_name %in% names(fit@modelInfo$param))) {

    #### Result ####

    return(FALSE)

  }

  model_psi <- fit@modelInfo$param[[psi_name]]
  numeric_psi <- suppressWarnings(as.numeric(model_psi))

  if(anyNA(numeric_psi)) {

    #### Result ####

    return(FALSE)

  }

  numeric_psi <- matrix(numeric_psi,
                        nrow = q,
                        ncol = q,
                        dimnames = dimnames(psi))

  result <- isTRUE(all.equal(numeric_psi, identity,
                             tolerance = sqrt(.Machine$double.eps)))

  #### Result ####

  return(result)

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

create_lrotate_dataList <- function(fit, projection, rotation) {

  source_data_param <- fit@dataList$data_param
  ngroups <- fit@dataList$ngroups

  ulambda_group <- source_data_param$lambda_group
  upsi_group <- source_data_param$psi_group
  ualpha_group <- source_data_param$alpha_group

  if(length(ulambda_group) != ngroups ||
     length(upsi_group) != ngroups ||
     length(ualpha_group) != ngroups) {
    stop("The fitted lcfa object contains incompatible group-specific parameter blocks")
  }

  lambda <- fit@transformed_pars[ulambda_group]
  psi <- fit@transformed_pars[upsi_group]
  alpha <- fit@transformed_pars[ualpha_group]

  valid_lambda <- vapply(lambda,
                         FUN = \(x) is.matrix(x) && is.numeric(x),
                         FUN.VALUE = logical(1L))
  valid_psi <- vapply(psi,
                      FUN = \(x) is.matrix(x) && is.numeric(x),
                      FUN.VALUE = logical(1L))
  valid_alpha <- vapply(alpha,
                        FUN = \(x) is.matrix(x) && is.numeric(x),
                        FUN.VALUE = logical(1L))

  if(!all(valid_lambda) || !all(valid_psi) || !all(valid_alpha)) {
    stop("The fitted lcfa object contains invalid lambda, psi, or alpha matrices")
  }

  item_label <- lapply(lambda, rownames)
  factor_label <- lapply(lambda, colnames)
  nitems <- lapply(lambda, nrow)
  nfactors <- lapply(lambda, ncol)

  for(i in seq_len(ngroups)) {

    q <- nfactors[[i]]

    if(!identical(dim(psi[[i]]), c(q, q)) ||
       !identical(dim(alpha[[i]]), c(q, 1L))) {
      stop("The fitted lcfa object contains incompatible psi or alpha dimensions")
    }

  }

  identity_psi <- vapply(seq_len(ngroups), FUN = \(i) {
    identity_psi_lrotate(fit = fit,
                         psi_name = upsi_group[i],
                         psi = psi[[i]])
  }, FUN.VALUE = logical(1L))

  orthogonal <- projection == "orth"

  dataList <- list(
    data = lambda,
    ngroups = ngroups,
    group = fit@dataList$group,
    group_label = fit@dataList$group_label,
    item_label = item_label,
    factor_label = factor_label,
    nitems = nitems,
    nfactors = nfactors,
    nobs = fit@dataList$nobs,
    projection = projection,
    rotation = rotation,
    estimator = rotation,
    positive = isTRUE(fit@dataList$positive),
    orthogonal = orthogonal,
    identity_psi = identity_psi,
    meanstructure = fit@dataList$meanstructure,
    se_type = fit@dataList$se_type,
    source_dof = fit@modelInfo$dof,
    source_nparam = fit@modelInfo$nparam,
    source_data_param = source_data_param,
    source_param = fit@modelInfo$param,
    source_trans = fit@modelInfo$trans,
    source_parameters = fit@parameters,
    source_transformed_pars = fit@transformed_pars,
    source_control_manifold = fit@modelInfo$control_manifold,
    source_control_transform = fit@modelInfo$control_transform
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
                               data_param = data_param,
                               control = control)

  #### Create the initial values for the parameters ####

  init_param <- start_lrotate(trans = trans,
                              dataList = dataList,
                              data_param = data_param,
                              control = control)

  #### Custom initial values ####

  init_param <- custom_init_param(control$start, init_param)

  # Recompute every transformed quantity from X after custom starts have been
  # inserted. The unrotated CFA parameters remain fixed at the supplied values.
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
  source_data_param <- dataList$source_data_param

  result <- list(
    ulambda_group = source_data_param$lambda_group,
    upsi_group = source_data_param$psi_group,
    ualpha_group = source_data_param$alpha_group,
    X_group = paste("X.group", groups, sep = ""),
    Xinv_group = paste("Xinv.group", groups, sep = ""),
    lambda_group = paste("lambda_rotated.group", groups, sep = ""),
    psi_group = paste("psi_rotated.group", groups, sep = ""),
    alpha_group = paste("alpha_rotated.group", groups, sep = "")
  )

  #### Result ####

  return(result)

}

#### Function to create matrix parameter labels ####

labels_lrotate <- function(group, object, nrow, ncol) {

  labels <- paste("rotation.g", group, ".", object, "[",
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

    #### Rotated factor covariance matrix ####

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

    #### Rotated factor means ####

    list_struct[[k]] <- list(
      name = alpha_group[i],
      type = "matrix",
      dim = c(q, 1L),
      rownames = factor_names,
      colnames = "intrcp",
      labels = labels_lrotate(i, "alpha", q, 1L)
    )
    k <- k+1L

  }

  rotation_trans <- create_parameters(list_struct)

  duplicated_names <- intersect(names(dataList$source_trans),
                                names(rotation_trans))

  if(length(duplicated_names) > 0L) {
    stop("The rotation parameter-block names conflict with the fitted lcfa model: ",
         paste(duplicated_names, collapse = ", "))
  }

  source_labels <- unique(c(unlist(dataList$source_trans,
                                   use.names = FALSE)))
  rotation_labels <- unique(c(unlist(rotation_trans,
                                     use.names = FALSE)))
  duplicated_labels <- intersect(source_labels, rotation_labels)

  if(length(duplicated_labels) > 0L) {
    stop("The rotation parameter labels conflict with the fitted lcfa model: ",
         paste(duplicated_labels, collapse = ", "))
  }

  trans <- c(dataList$source_trans, rotation_trans)

  #### Result ####

  return(trans)

}

#### Function to create the parameter constraints ####

constraints_lrotate <- function(trans, dataList, data_param, control) {

  list2env(data_param, envir = environment())

  source_names <- names(dataList$source_param)

  if(isTRUE(control$free_previous)) {

    param <- dataList$source_param

  } else {

    missing_source <- setdiff(source_names,
                              names(dataList$source_parameters))

    if(length(missing_source) > 0L) {
      stop("The fitted lcfa object is missing parameter estimate(s): ",
           paste(missing_source, collapse = ", "))
    }

    param <- dataList$source_parameters[source_names]

  }

  for(i in seq_len(dataList$ngroups)) {

    # X is the only new free parameter block for the rotation.
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

    init_param[[rs]] <- dataList$source_transformed_pars

    for(i in seq_len(dataList$ngroups)) {

      q <- dataList$nfactors[[i]]

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

    init_param[[rs]][names(dataList$source_transformed_pars)] <-
      dataList$source_transformed_pars

    for(i in seq_len(dataList$ngroups)) {

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
  U <- x[[ulambda_group[i]]]
  psi_0 <- x[[upsi_group[i]]]
  alpha_0 <- x[[ualpha_group[i]]]
  X <- x[[X_group[i]]]

  if(is.null(X) ||
     !identical(dim(X), c(dataList$nfactors[[i]],
                          dataList$nfactors[[i]]))) {
    stop("Invalid starting value for '", X_group[i], "'")
  }

  dimnames(X) <- dimnames(trans[[X_group[i]]])
  x[[X_group[i]]] <- X

  if(dataList$orthogonal) {

    rotated_lambda <- U%*%X

  } else {

    Xinv <- solve(X)
    dimnames(Xinv) <- dimnames(trans[[Xinv_group[i]]])
    x[[Xinv_group[i]]] <- Xinv

    rotated_lambda <- U%*%t(Xinv)

  }

  dimnames(rotated_lambda) <- dimnames(trans[[lambda_group[i]]])
  x[[lambda_group[i]]] <- rotated_lambda

  if(dataList$identity_psi[i]) {
    rotated_psi <- crossprod(X)
  } else {
    rotated_psi <- t(X)%*%psi_0%*%X
  }

  rotated_psi <- 0.5*(rotated_psi+t(rotated_psi))
  dimnames(rotated_psi) <- dimnames(trans[[psi_group[i]]])
  x[[psi_group[i]]] <- rotated_psi

  rotated_alpha <- t(X)%*%alpha_0
  dimnames(rotated_alpha) <- dimnames(trans[[alpha_group[i]]])
  x[[alpha_group[i]]] <- rotated_alpha

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

  control_manifold_rotation <-
    create_manifolds(manifolds = manifolds,
                     structures = param)

  if(isTRUE(control$free_previous)) {
    control_manifold <- c(dataList$source_control_manifold,
                          control_manifold_rotation)
  } else {
    control_manifold <- control_manifold_rotation
  }

  #### Transformations ####

  transforms <- transformations_lrotate(dataList = dataList,
                                        data_param = data_param)

  control_transform_rotation <-
    create_transforms(transforms = transforms,
                      structures = trans)

  control_transform <- c(dataList$source_control_transform,
                         control_transform_rotation)

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
    dof = dataList$source_dof,
    rotation = dataList$rotation,
    projection = dataList$projection,
    data_param = data_param,
    source_nparam = dataList$source_nparam,
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

    #### Rotated factor covariance matrix ####

    if(dataList$identity_psi[i]) {

      transforms[[k]] <- list(
        transform = "crossprod",
        parameters_in = X_group[i],
        parameters_out = psi_group[i],
        extra = list(p = q)
      )

    } else {

      transforms[[k]] <- list(
        transform = "XtYX",
        parameters_in = c(X_group[i], upsi_group[i]),
        parameters_out = psi_group[i],
        extra = list(p = q, q = q)
      )

    }

    k <- k+1L

    #### Rotated factor means ####

    transforms[[k]] <- list(
      transform = "XtY",
      parameters_in = c(X_group[i], ualpha_group[i]),
      parameters_out = alpha_group[i],
      extra = list(p = q, q = q, r = 1L)
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
