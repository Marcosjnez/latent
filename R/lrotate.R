# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 23/08/2026
#'
#' @title
#' Rotate factor loading and covariance matrices
#'
#' @description
#' \code{lrotate} rotates the factor loading and factor covariance matrices
#' supplied directly or extracted from a fitted \code{lcfa} object using an
#' orthogonal or oblique projection and a selected rotation criterion.
#'
#' @usage
#' lrotate(fit = NULL, lambda = NULL, psi = NULL,
#'         projection = "oblq", rotation = "oblimin",
#'         se = TRUE, do.fit = TRUE, control = NULL, ...)
#'
#' @param fit Optional fitted object inheriting from class \code{"lcfa"}.
#' @param lambda Optional loading matrix or list of loading matrices. This is an
#'   alternative to supplying \code{fit}.
#' @param psi Optional factor covariance matrix or list of factor covariance
#'   matrices corresponding to \code{lambda}. If omitted, identity matrices are
#'   used. This argument cannot be used together with \code{fit}.
#' @param projection Character string. Available projections are
#'   \code{"orth"}, \code{"oblq"}, and \code{"poblq"}.
#' @param rotation Character string identifying the rotation criterion.
#' @param se Logical. If \code{TRUE} and \code{fit} is supplied, propagate
#'   standard errors from the fitted \code{lcfa} model to the rotated
#'   parameters. Standard errors are not available when matrices are supplied
#'   directly.
#' @param do.fit Logical. If \code{TRUE}, fit the rotation. If \code{FALSE},
#'   return the model specification. With a fitted \code{lcfa} input, the
#'   unrestricted specification used for derivative calculations is returned.
#' @param control List of optimization-control arguments.
#' @param ... Additional arguments required by the selected projection or
#'   rotation criterion.
#'
#' @details
#' Exactly one of \code{fit} and \code{lambda} must be supplied. Let \eqn{X}
#' be the rotation matrix and let \eqn{\Lambda_0}, \eqn{\Psi_0}, and
#' \eqn{\alpha_0} denote the unrotated factor loadings, factor covariance
#' matrix, and factor means. The rotated quantities are
#' \deqn{\Lambda_r=\Lambda_0X^{-T},}
#' \deqn{\Psi_r=X^T\Psi_0X,}
#' and
#' \deqn{\alpha_r=X^T\alpha_0.}
#' For an orthogonal projection, \eqn{X^{-T}=X}. If \eqn{\Psi_0} is a fixed
#' identity matrix, \eqn{\Psi_r} is computed as \eqn{X^TX}.
#'
#' When \code{fit} is supplied, the returned object inherits from
#' \code{"multistep"} and the fitted \code{lcfa} object is stored in
#' \code{extra}. When matrices are supplied directly, the returned object
#' inherits only from \code{"latent"}; sampling uncertainty is not propagated
#' because no fitted source model is available.
#'
#' @return An object inheriting from class \code{"latent"}.
#'
#' @examples
#' \dontrun{
#' fit_cfa <- lcfa(data = HolzingerSwineford1939,
#'                 model = model,
#'                 std.lv = TRUE)
#' fit_rotation <- lrotate(fit = fit_cfa,
#'                         projection = "oblq",
#'                         rotation = "oblimin")
#'
#' direct_rotation <- lrotate(lambda = lambda,
#'                            psi = psi,
#'                            projection = "oblq",
#'                            rotation = "oblimin")
#' }
#'
#' @export
lrotate <- function(fit = NULL, lambda = NULL, psi = NULL,
                    projection = "oblq", rotation = "oblimin",
                    se = TRUE, do.fit = TRUE, control = NULL, ...) {

  #### Check input arguments ####

  fit_input <- !is.null(fit)
  matrix_input <- !is.null(lambda)

  if(fit_input == matrix_input) {
    stop("Supply exactly one of fit or lambda")
  }

  if(fit_input) {

    check_fit_lrotate(fit)

    if(!is.null(psi)) {
      stop("psi cannot be supplied together with fit")
    }

  } else {

    lambda <- check_lambda_lrotate(lambda)
    psi <- check_psi_lrotate(psi = psi,
                             lambda = lambda)

  }

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

  if(length(se) != 1L || !is.logical(se) || is.na(se)) {
    stop("se must be TRUE or FALSE")
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
  control$free_previous <- fit_input && !do.fit

  #### Create the dataList ####

  dataList <- create_lrotate_dataList(fit = fit,
                                      lambda = lambda,
                                      psi = psi,
                                      projection = projection,
                                      rotation = rotation)

  input_args <- if(fit_input) {
    list(fit = fit)
  } else {
    list(lambda = lambda,
         psi = psi)
  }

  dataList$args <- c(input_args,
                     list(projection = projection,
                          rotation = rotation,
                          se = se,
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

  modelInfo$propagate_uncertainty <- fit_input
  modelInfo$step_labels <- modelInfo$parameters_labels

  object_class <- if(fit_input) "multistep" else "latent"
  extra <- if(fit_input) list(efa = fit) else list()

  #### Fit the model ####

  if(!do.fit) {

    result <- new(object_class,
                  version          = as.character(packageVersion("latent")),
                  call             = mc,
                  timing           = numeric(),
                  dataList         = dataList,
                  modelInfo        = modelInfo,
                  Optim            = list(),
                  parameters       = list(),
                  transformed_pars = list(),
                  extra            = extra)

    #### Result ####

    return(result)

  }

  Optim <- fit_lrotate(modelInfo = modelInfo)

  #### Process the outputs ####

  transformed_pars <- fill_in(modelInfo$trans,
                              Optim$transparameters)

  parameters <- transformed_pars[names(modelInfo$param)]

  #### latent object ####

  result <- new(object_class,
                version          = as.character(packageVersion("latent")),
                call             = mc,
                timing           = Optim$elapsed,
                dataList         = dataList,
                modelInfo        = modelInfo,
                Optim            = Optim,
                parameters       = parameters,
                transformed_pars = transformed_pars,
                extra            = extra)

  #### Standard errors ####

  if(fit_input && se) {

    data_param <- modelInfo$data_param
    rotated_parameters <- modelInfo$trans[
      unique(c(data_param$lambda_group,
               data_param$psi_group,
               data_param$alpha_group))
    ]

    result@Optim$SE <- se_lrotate_multistep(
      fit = result,
      parameters = rotated_parameters
    )

  }

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

#### Function to check loading matrices ####

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

    if(ngroups == 1L) {
      group_label <- ""
    } else {
      group_label <- paste("group", seq_len(ngroups), sep = "")
    }

  } else {

    if(ngroups == 1L) {

      if(is.na(group_label) || group_label == "") {
        group_label <- ""
      }

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

#### Function to check factor covariance matrices ####

check_psi_lrotate <- function(psi, lambda) {

  ngroups <- length(lambda)
  group_label <- names(lambda)

  if(is.null(psi)) {

    psi <- lapply(lambda, FUN = \(x) {
      factor_label <- colnames(x)
      result <- diag(ncol(x))
      dimnames(result) <- list(factor_label, factor_label)
      return(result)
    })

  } else if(is.matrix(psi)) {

    if(ngroups != 1L) {
      stop("psi must be a list with one matrix per group")
    }

    psi <- list(psi)

  } else if(!is.list(psi) || length(psi) == 0L) {

    stop("psi must be NULL, a matrix, or a non-empty list of matrices")

  }

  if(length(psi) != ngroups) {
    stop("lambda and psi must contain the same number of groups")
  }

  psi_names <- names(psi)

  if(!is.null(psi_names) &&
     length(psi_names) == ngroups &&
     all(!is.na(psi_names)) &&
     all(psi_names != "")) {

    if(anyDuplicated(psi_names)) {
      stop("psi must have unique group names")
    }

    if(ngroups > 1L) {

      if(!setequal(psi_names, group_label)) {
        stop("The group names of psi must match the group names of lambda")
      }

      psi <- psi[group_label]

    }

  }

  names(psi) <- group_label

  for(i in seq_len(ngroups)) {

    if(is.data.frame(psi[[i]])) {
      psi[[i]] <- as.matrix(psi[[i]])
    }

    q <- ncol(lambda[[i]])
    factor_label <- colnames(lambda[[i]])

    if(!is.matrix(psi[[i]]) ||
       !is.numeric(psi[[i]]) ||
       !identical(dim(psi[[i]]), c(q, q))) {
      stop("Every factor covariance matrix must be a numeric q by q matrix")
    }

    if(anyNA(psi[[i]]) || any(!is.finite(psi[[i]]))) {
      stop("Factor covariance matrices cannot contain missing or non-finite values")
    }

    if(!isSymmetric(psi[[i]], tol = sqrt(.Machine$double.eps))) {
      stop("Every factor covariance matrix must be symmetric")
    }

    rn <- rownames(psi[[i]])
    cn <- colnames(psi[[i]])

    if(is.null(rn) && is.null(cn)) {

      dimnames(psi[[i]]) <- list(factor_label, factor_label)

    } else {

      if(is.null(rn) || is.null(cn) ||
         any(rn == "") || any(cn == "") ||
         anyDuplicated(rn) || anyDuplicated(cn) ||
         !setequal(rn, factor_label) ||
         !setequal(cn, factor_label)) {
        stop("Factor covariance matrix names must match the loading-matrix factor names")
      }

      psi[[i]] <- psi[[i]][factor_label, factor_label, drop = FALSE]

    }

  }

  #### Result ####

  return(psi)

}

#### Function to identify an identity matrix ####

identity_matrix_lrotate <- function(psi) {

  if(!is.matrix(psi) || nrow(psi) != ncol(psi)) {

    #### Result ####

    return(FALSE)

  }

  identity <- diag(nrow(psi))
  dimnames(identity) <- dimnames(psi)

  result <- isTRUE(all.equal(psi, identity,
                             tolerance = sqrt(.Machine$double.eps)))

  #### Result ####

  return(result)

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

#### Function to create group-specific parameter-block names ####

group_names_lrotate <- function(name, group_label) {

  if(length(group_label) < 2L) {
    result <- rep(name, length(group_label))
  } else {
    result <- paste(name, group_label, sep = ".")
  }

  #### Result ####

  return(result)

}

#### Function to create the dataList ####

create_lrotate_dataList <- function(fit = NULL, lambda = NULL, psi = NULL,
                                    projection, rotation) {

  if(!is.null(fit)) {

    result <- create_lrotate_dataList_fit(
      fit = fit,
      projection = projection,
      rotation = rotation
    )

  } else {

    result <- create_lrotate_dataList_matrices(
      lambda = lambda,
      psi = psi,
      projection = projection,
      rotation = rotation
    )

  }

  #### Result ####

  return(result)

}

#### Function to create the dataList from a fitted CFA model ####

create_lrotate_dataList_fit <- function(fit, projection, rotation) {

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
    source = "lcfa",
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

#### Function to create the dataList from matrices ####

create_lrotate_dataList_matrices <- function(lambda, psi,
                                             projection, rotation) {

  ngroups <- length(lambda)
  group_label <- names(lambda)
  item_label <- lapply(lambda, rownames)
  factor_label <- lapply(lambda, colnames)
  nitems <- lapply(lambda, nrow)
  nfactors <- lapply(lambda, ncol)

  alpha <- lapply(seq_len(ngroups), FUN = \(i) {
    matrix(0,
           nrow = nfactors[[i]],
           ncol = 1L,
           dimnames = list(factor_label[[i]], "intrcp"))
  })
  names(alpha) <- group_label

  source_data_param <- list(
    lambda_group = group_names_lrotate("lambda", group_label),
    psi_group = group_names_lrotate("psi", group_label),
    alpha_group = group_names_lrotate("alpha", group_label)
  )

  source_model <- create_lrotate_source_model(
    lambda = lambda,
    psi = psi,
    alpha = alpha,
    source_data_param = source_data_param
  )

  identity_psi <- vapply(psi,
                         FUN = identity_matrix_lrotate,
                         FUN.VALUE = logical(1L))

  orthogonal <- projection == "orth"

  dataList <- list(
    data = lambda,
    source = "matrices",
    ngroups = ngroups,
    group = NULL,
    group_label = group_label,
    item_label = item_label,
    factor_label = factor_label,
    nitems = nitems,
    nfactors = nfactors,
    nobs = vector("list", ngroups),
    projection = projection,
    rotation = rotation,
    estimator = rotation,
    positive = FALSE,
    orthogonal = orthogonal,
    identity_psi = identity_psi,
    meanstructure = FALSE,
    se_type = NULL,
    source_dof = NA_integer_,
    source_nparam = 0L,
    source_data_param = source_data_param,
    source_param = source_model$param,
    source_trans = source_model$trans,
    source_parameters = source_model$param,
    source_transformed_pars = source_model$param,
    source_control_manifold = list(),
    source_control_transform = list()
  )

  #### Result ####

  return(dataList)

}

#### Function to create fixed source matrix structures ####

create_lrotate_source_model <- function(lambda, psi, alpha,
                                        source_data_param) {

  list2env(source_data_param, envir = environment())

  list_struct <- list()
  k <- 1L

  for(i in seq_along(lambda)) {

    p <- nrow(lambda[[i]])
    q <- ncol(lambda[[i]])
    item_label <- rownames(lambda[[i]])
    factor_label <- colnames(lambda[[i]])

    list_struct[[k]] <- list(
      name = lambda_group[i],
      type = "matrix",
      dim = c(p, q),
      rownames = item_label,
      colnames = factor_label,
      labels = labels_lrotate(i, "lambda_unrotated", p, q)
    )
    k <- k+1L

    list_struct[[k]] <- list(
      name = psi_group[i],
      type = "matrix",
      dim = c(q, q),
      rownames = factor_label,
      colnames = factor_label,
      symmetric = TRUE,
      labels = labels_lrotate(i, "psi_unrotated", q, q)
    )
    k <- k+1L

    list_struct[[k]] <- list(
      name = alpha_group[i],
      type = "matrix",
      dim = c(q, 1L),
      rownames = factor_label,
      colnames = "intrcp",
      labels = labels_lrotate(i, "alpha_unrotated", q, 1L)
    )
    k <- k+1L

  }

  trans <- create_parameters(list_struct)
  param <- list()

  for(i in seq_along(lambda)) {

    param[[lambda_group[i]]] <- lambda[[i]]
    dimnames(param[[lambda_group[i]]]) <-
      dimnames(trans[[lambda_group[i]]])

    param[[psi_group[i]]] <- psi[[i]]
    dimnames(param[[psi_group[i]]]) <-
      dimnames(trans[[psi_group[i]]])

    param[[alpha_group[i]]] <- alpha[[i]]
    dimnames(param[[alpha_group[i]]]) <-
      dimnames(trans[[alpha_group[i]]])

  }

  #### Result ####

  result <- list(param = param,
                 trans = trans)

  return(result)

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

  source_data_param <- dataList$source_data_param
  group_label <- dataList$group_label

  result <- list(
    ulambda_group = source_data_param$lambda_group,
    upsi_group = source_data_param$psi_group,
    ualpha_group = source_data_param$alpha_group,
    X_group = group_names_lrotate("X", group_label),
    Xinv_group = group_names_lrotate("Xinv", group_label),
    lambda_group = group_names_lrotate("lambda_rotated", group_label),
    psi_group = group_names_lrotate("psi_rotated", group_label),
    alpha_group = group_names_lrotate("alpha_rotated", group_label)
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
