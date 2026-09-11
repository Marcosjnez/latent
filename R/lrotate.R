# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 11/09/2026
#'
#' @title
#' Rotate factor loading and covariance matrices
#'
#' @description
#' \code{lrotate} rotates the factor loading and factor covariance matrices
#' supplied directly or extracted from a fitted \code{lcfa} object using an
#' orthogonal, oblique, or orthoblique projection and one or more rotation criteria.
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
#' @param rotation A criterion name, a character vector of criterion names, or
#'   a named list of parameter lists. Vector entries apply to the full rotated
#'   loading matrix by default and their losses are summed. For a named list,
#'   names identify the criteria (duplicate names are allowed), and each list
#'   contains that component's arguments, optionally including \code{items}
#'   and \code{factors} to select rows and columns. Missing selectors use all
#'   rows or columns. See Details for defaults and target subsetting.
#' @param se Logical. If \code{TRUE} and \code{fit} is supplied, propagate
#'   standard errors from the fitted \code{lcfa} model to the rotated
#'   parameters. Standard errors are not available when matrices are supplied
#'   directly.
#' @param do.fit Logical. If \code{TRUE}, fit the rotation. If \code{FALSE},
#'   return the model specification. With a fitted \code{lcfa} input, the
#'   unrestricted specification used for derivative calculations is returned.
#' @param control List of optimization-control arguments.
#' @param ... Additional projection or rotation arguments shared by the
#'   components to which they apply. Component-specific arguments take
#'   precedence. Missing (or NULL) loading weights default to \code{1-target},
#'   and covariance weights to \code{1-psitarget}. Explicit weights, including
#'   zero matrices, are preserved. Group-specific lists are supported.
#'
#' @details
#' Parameter labels are generated from their block names, just as in the CFA
#' parameter constructor. Single-group labels have no group suffix, for example
#' \code{X[1,1]} and \code{lambda_rotated[1,1]}. Multiple groups use the
#' corresponding group names to keep parameter labels distinct.
#'
#' Fitted objects retain the estimated factor order and signs. Sorting is
#' available only through \code{latInspect(fit, sort = TRUE)} and does not
#' modify the fitted object or its standard-error calculations.
#'
#' All components are optimized simultaneously over the same rotation matrix.
#'   With row sets \eqn{I_s} and factor sets \eqn{J_s}, the total criterion is
#'   \deqn{Q(\Lambda,\Psi)=\sum_s Q_s(\Lambda_{I_s,J_s},\Psi_{J_s,J_s}).}
#'   Only \code{xtarget} uses the selected factor covariance matrix; other
#'   criteria use the loading submatrix only. Overlapping selections are
#'   allowed and their objective, gradient, and Hessian contributions are added.
#'   There is no automatic rescaling or averaging of the component criteria.
#'
#' \code{items} and \code{factors} accept positive integer positions, matrix
#'   row/column names, or logical vectors of the full corresponding length.
#'   Duplicated or empty selections are rejected. Selector order is retained.
#'   A target or weight matrix may have full loading-matrix dimensions, in which
#'   case it is subset automatically, or the selected submatrix's dimensions,
#'   in which case it is used in the supplied selector order. When both sizes
#'   coincide, it is treated as a full matrix. For \code{xtarget}, the same rule
#'   applies to the principal factor submatrices of \code{psitarget} and
#'   \code{psiweight}; \code{items} affects only the loading part.
#'
#' Defaults are resolved independently for each component and group.
#'   \code{geomin} defaults to \code{epsilon = 0.01}; \code{oblimin}
#'   defaults to \code{gamma = 0}. \code{alpha} is an alias for oblimin's
#'   \code{gamma}, not a component weight. Conflicting values supplied at the
#'   same level are rejected. The existing required arguments of other
#'   criteria remain required, including \code{k} for \code{cf},
#'   \code{epsilon} for \code{lclf}, and \code{w} for \code{xtarget}.
#'   A local NULL requests the criterion default rather than inheriting a
#'   shared value. Group-specific selectors or parameters may be supplied as
#'   lists with one entry per group, in the input group order.
#'
#' For target criteria the loss uses squared weighted residuals. To multiply a
#'   target loss by \eqn{c}, multiply its weight matrix by \eqn{\sqrt{c}}.
#'   For \code{xtarget}, scale both weight matrices to scale the entire term;
#'   its \code{w} remains the relative covariance-target weight. Sparse
#'   selections must jointly provide a sufficiently identified rotation for
#'   standard errors; selecting a submatrix alone does not guarantee this.
#'
#' \code{dataList$rotation} is a printable criterion label, while
#'   \code{dataList$rotation_spec} retains the vector or named-list
#'   specification. \code{modelInfo$rotation_components} records the resolved
#'   row/factor positions and parameters for each component in each group.
#'
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
#' With \code{projection = "poblq"}, either \code{constraints} or \code{oblique}
#' must be supplied through \code{...}. The former uses arbitrary structural
#' constraints. A numeric \code{oblique} vector gives the sizes of consecutive
#' oblique blocks. Alternatively, a list gives their explicit one-based factor
#' positions, for example \code{list(c(1, 3, 5, 9), c(2, 4, 6, 7, 8))}.
#' Blocks must be non-empty and disjoint; singleton blocks are allowed.
#' Unlisted factors are mutually orthogonal and orthogonal to all listed blocks.
#' The factor order, parameter labels, and target/weight indexing are unchanged.
#' These constraints apply to \eqn{X^TX}; the formula for \eqn{\Psi_r} above
#' remains unchanged when the input factor covariance is not the identity.
#' \code{constraints} and \code{oblique} cannot be used together.
#'
#' A flat list of \code{oblique} position vectors specifies blocks shared by
#' all groups, even when its length equals the number of groups. Group-specific
#' blocks require an outer list with one block list per group, in input group
#' order: for example \code{list(list(c(1, 3), c(2, 4)), list(1:2, 3:4))}.
#' Previously used flat lists of group-specific block-size vectors must instead
#' be expressed as nested lists of explicit positions. A numeric size vector
#' shared by all groups retains its existing meaning.
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
#'
#' mixed_rotation <- lrotate(lambda = lambda,
#'                            rotation = c("oblimin", "target", "geomin"),
#'                            target = target, weight = 1-target)
#'
#' subset_rotation <- lrotate(lambda = lambda,
#'                             rotation = list(
#'                               oblimin = list(alpha = 0, items = 1:5, factors = 1:3),
#'                               oblimin = list(alpha = 0.5, items = 6:15, factors = 4:5),
#'                               target = list(target = target, weight = 1-target),
#'                               geomin = list(items = 16:20)))
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
  rotation <- check_rotation_lrotate(rotation)
  rotation_names <- rotation_names_lrotate(rotation)
  rotation_label <- paste(rotation_names, collapse = " + ")

  supported_projection <- c("orth", "oblq", "poblq")

  if(!(projection %in% supported_projection)) {
    stop("Unknown projection: ", projection)
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

  if("sort" %in% names(dots)) {
    stop("sort is only available in latInspect().")
  }
  check_rotation_dots_lrotate(dots)

  # Defaults are applied separately after resolving each component and group.
  # In particular, a component's target must not inherit weights calculated
  # from a different target supplied through ... .

  check_poblq_arguments_lrotate(projection = projection,
                                dots = dots)

  #### Check control parameters ####

  control$penalties <- FALSE
  control$positive <- FALSE
  control$estimator <- rotation_label
  control$projection <- projection
  control <- lrotate_control(control)
  control$free_previous <- fit_input && !do.fit

  #### Create the dataList ####

  dataList <- create_lrotate_dataList(fit = fit,
                                      lambda = lambda,
                                      psi = psi,
                                      projection = projection,
                                      rotation = rotation_label)
  dataList$rotation_spec <- rotation

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

    result@Optim$SE <- se.multistep(
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

#### Functions to describe and validate rotation components ####

rotation_parameters_lrotate <- function(rotation) {

  parameters <- switch(rotation,
                        cf = "k",
                        geomin = "epsilon",
                        lclf = "epsilon",
                        oblimin = c("gamma", "alpha"),
                        target = c("target", "weight"),
                        varimax = character(0L),
                        varimin = character(0L),
                        xtarget = c("target", "weight", "w", "psitarget", "psiweight"),
                        stop("Unknown rotation criterion: ", rotation))
  result <- c(parameters, "items", "factors")

  #### Result ####

  return(result)

}

check_rotation_dots_lrotate <- function(dots) {

  if(length(dots) > 0L &&
     (is.null(names(dots)) || anyNA(names(dots)) || any(names(dots) == "") ||
      anyDuplicated(names(dots)))) {
    stop("Rotation parameters must have unique, non-empty names")
  }

  #### Result ####

  return(invisible(NULL))

}

check_rotation_lrotate <- function(rotation) {

  if(is.character(rotation) && is.null(dim(rotation))) {

    if(length(rotation) == 0L || anyNA(rotation) || any(rotation == "")) {
      stop("rotation must contain at least one non-missing criterion name")
    }
    rotation <- tolower(rotation)
    criteria <- unname(rotation)

  } else if(is.list(rotation) && !is.data.frame(rotation) && is.null(dim(rotation))) {

    criteria <- names(rotation)
    if(length(rotation) == 0L || is.null(criteria) || anyNA(criteria) ||
       any(criteria == "")) {
      stop("rotation must be a non-empty named list of criterion parameter lists")
    }
    criteria <- tolower(criteria)
    names(rotation) <- criteria

    for(i in seq_along(rotation)) {
      component <- rotation[[i]]
      if(!is.list(component) || is.data.frame(component) || !is.null(dim(component))) {
        stop("Every rotation component must be a parameter list; use list() for defaults")
      }
      check_rotation_dots_lrotate(component)
      unknown <- setdiff(names(component), rotation_parameters_lrotate(criteria[i]))
      if(length(unknown) > 0L) {
        stop("Unknown parameter(s) for rotation component ", i, " ('", criteria[i],
             "'): ", paste(unknown, collapse = ", "))
      }
    }

  } else {

    stop("rotation must be a character vector or a named list of parameter lists")

  }

  supported <- c("cf", "geomin", "lclf", "oblimin",
                 "target", "varimax", "varimin", "xtarget")
  unknown <- setdiff(criteria, supported)
  if(length(unknown) > 0L) {
    stop("Unknown rotation criterion: ", paste(unknown, collapse = ", "))
  }

  #### Result ####

  return(rotation)

}

rotation_names_lrotate <- function(rotation) {

  result <- if(is.character(rotation)) unname(rotation) else names(rotation)

  #### Result ####

  return(result)

}

rotation_components_lrotate <- function(rotation) {

  if(is.character(rotation)) {
    result <- rep(list(list()), length(rotation))
    names(result) <- unname(rotation)
  } else {
    result <- rotation
  }

  #### Result ####

  return(result)

}

rotation_component_dots_lrotate <- function(criterion, component, dots) {

  allowed <- rotation_parameters_lrotate(criterion)
  extra <- dots[intersect(names(dots), allowed)]

  # Treat the oblimin aliases at the same precedence level. A local alpha
  # overrides a shared gamma (and vice versa), rather than silently losing it.
  if(criterion == "oblimin" && any(c("alpha", "gamma") %in% names(component))) {
    extra$alpha <- extra$gamma <- NULL
  }
  for(nm in names(component)) extra[nm] <- component[nm]

  #### Result ####

  return(extra)

}

rotation_indices_lrotate <- function(index, size, labels, argument) {

  if(is.null(index)) {
    result <- seq_len(size)
  } else if(is.character(index) && is.null(dim(index))) {
    if(anyNA(index) || any(index == "") || is.null(labels)) {
      stop(argument, " must contain valid matrix names")
    }
    result <- match(index, labels)
    if(anyNA(result)) stop("Unknown ", argument, ": ", paste(index[is.na(result)], collapse = ", "))
  } else if(is.logical(index) && is.null(dim(index))) {
    if(length(index) != size || anyNA(index)) {
      stop("Logical ", argument, " must have length ", size, " and no missing values")
    }
    result <- which(index)
  } else if(is.numeric(index) && !is.complex(index) && is.null(dim(index))) {
    if(any(!is.finite(index)) || any(index < 1 | index > size | index != trunc(index))) {
      stop(argument, " must contain integer positions between 1 and ", size)
    }
    result <- as.integer(index)
  } else {
    stop(argument, " must be NULL, positive integer positions, names, or a logical vector")
  }

  if(length(result) == 0L || anyDuplicated(result)) {
    stop(argument, " must select at least one entry without duplicates")
  }

  #### Result ####

  return(result)

}

rotation_matrix_lrotate <- function(x, rows, columns, p, q, argument) {

  if(!is.matrix(x) || !(is.numeric(x) || is.logical(x)) || is.complex(x) ||
     any(!is.finite(x))) {
    stop(argument, " must be a finite numeric or logical matrix")
  }

  if(identical(dim(x), c(p, q))) {
    result <- x[rows, columns, drop = FALSE]
  } else if(identical(dim(x), c(length(rows), length(columns)))) {
    result <- x
  } else {
    stop(argument, " must have full dimensions ", p, " by ", q,
         " or selected dimensions ", length(rows), " by ", length(columns))
  }
  storage.mode(result) <- "double"

  #### Result ####

  return(result)

}

rotation_extra_lrotate <- function(criterion, extra, items, factors, p, q) {

  scalar <- switch(criterion, cf = "k", geomin = "epsilon", lclf = "epsilon",
                    oblimin = "gamma", xtarget = "w", character(0L))
  for(nm in scalar) {
    value <- extra[[nm]]
    if(length(value) != 1L || !is.numeric(value) || is.complex(value) ||
       !is.finite(value)) {
      stop("Rotation '", criterion, "' requires one finite numeric '", nm, "'")
    }
    if(nm == "epsilon" && value <= 0) stop("epsilon must be strictly positive")
  }

  if(criterion %in% c("target", "xtarget")) {
    for(nm in c("target", "weight")) {
      extra[[nm]] <- rotation_matrix_lrotate(extra[[nm]], items, factors, p, q, nm)
    }
  }
  if(criterion == "xtarget") {
    for(nm in c("psitarget", "psiweight")) {
      extra[[nm]] <- rotation_matrix_lrotate(extra[[nm]], factors, factors, q, q, nm)
    }
  }

  #### Result ####

  return(extra)

}

#### Function to create defaults for rotation-specific arguments ####

rotation_defaults_lrotate <- function(rotation, dots) {

  # These defaults make the two commonly used criteria directly usable from
  # lrotate()/lefa() without requiring an otherwise undocumented argument.

  if(rotation == "oblimin") {

    # alpha is accepted as an alias; the native criterion uses gamma.
    if(!is.null(dots$alpha)) {
      if(!is.null(dots$gamma) &&
         !isTRUE(all.equal(dots$alpha, dots$gamma, check.attributes = FALSE))) {
        stop("Supply only one of alpha and gamma for an oblimin component")
      }
      dots$gamma <- dots$alpha
    }
    dots$alpha <- NULL
    if(is.null(dots$gamma)) dots$gamma <- 0

  }

  if(rotation == "geomin" && is.null(dots$epsilon)) {
    dots$epsilon <- 0.01
  }

  if(rotation %in% c("target", "xtarget") && is.null(dots$weight) &&
     !is.null(dots$target)) {
    dots$weight <- complement_target_weights_lrotate(dots$target)
  }

  if(rotation == "xtarget" && is.null(dots$psiweight) &&
     !is.null(dots$psitarget)) {
    dots$psiweight <- complement_target_weights_lrotate(dots$psitarget)
  }

  #### Result ####

  return(dots)

}

complement_target_weights_lrotate <- function(target) {

  if(is.null(target)) {
    result <- NULL
  } else if(is.list(target) && !is.data.frame(target)) {
    result <- lapply(target, complement_target_weights_lrotate)
  } else {
    if(!is.numeric(target) && !is.logical(target)) {
      stop("Targets must be numeric/logical matrices or group-specific lists.")
    }
    result <- 1-target
  }

  #### Result ####

  return(result)

}

#### Function to check partially oblique projection arguments ####

check_poblq_arguments_lrotate <- function(projection, dots) {

  constraints <- !is.null(dots$constraints)
  oblique <- !is.null(dots$oblique)

  if(constraints && oblique) {
    stop("constraints and oblique cannot be supplied together")
  }

  if(projection == "poblq" && !constraints && !oblique) {
    stop("projection = 'poblq' requires either constraints or oblique")
  }

  if(projection != "poblq" && oblique) {
    stop("oblique can only be used with projection = 'poblq'")
  }

  #### Result ####

  return(invisible(NULL))

}

#### Function to check one group's oblique block specification ####

check_oblique_lrotate <- function(oblique, q) {

  if(is.numeric(oblique) && !is.complex(oblique) && is.null(dim(oblique))) {

    if(length(oblique) == 0L || any(!is.finite(oblique)) ||
       any(oblique < 1 | oblique != trunc(oblique)) || sum(oblique) > q) {
      stop("oblique block sizes must be positive integers with sum <= ", q)
    }

  } else if(is.list(oblique) && !is.data.frame(oblique) && is.null(dim(oblique))) {

    valid <- vapply(oblique, FUN = \(block) {
      result <- is.numeric(block) && !is.complex(block) && is.null(dim(block)) &&
        length(block) > 0L && all(is.finite(block)) &&
        all(block >= 1 & block <= q & block == trunc(block))

      #### Result ####

      return(result)
    }, FUN.VALUE = logical(1L))

    if(length(oblique) == 0L || !all(valid)) {
      stop("Each oblique block must contain integer factor positions between 1 and ", q)
    }
    if(anyDuplicated(unlist(oblique, use.names = FALSE))) {
      stop("A factor cannot appear more than once within or across oblique blocks")
    }

  } else {

    stop("oblique must be a numeric vector of block sizes or a list of factor-position vectors")

  }

  #### Result ####

  return(oblique)

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
      colnames = factor_label
    )
    k <- k+1L

    list_struct[[k]] <- list(
      name = psi_group[i],
      type = "matrix",
      dim = c(q, q),
      rownames = factor_label,
      colnames = factor_label,
      symmetric = TRUE
    )
    k <- k+1L

    list_struct[[k]] <- list(
      name = alpha_group[i],
      type = "matrix",
      dim = c(q, 1L),
      rownames = factor_label,
      colnames = "intrcp"
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
      colnames = factor_names
    )
    k <- k+1L

    #### Inverse rotation matrix ####

    if(!dataList$orthogonal) {

      list_struct[[k]] <- list(
        name = Xinv_group[i],
        type = "matrix",
        dim = c(q, q),
        rownames = factor_names,
        colnames = factor_names
      )
      k <- k+1L

    }

    #### Rotated loadings ####

    list_struct[[k]] <- list(
      name = lambda_group[i],
      type = "matrix",
      dim = c(p, q),
      rownames = item_names,
      colnames = factor_names
    )
    k <- k+1L

    #### Rotated factor covariance matrix ####

    list_struct[[k]] <- list(
      name = psi_group[i],
      type = "matrix",
      dim = c(q, q),
      rownames = factor_names,
      colnames = factor_names,
      symmetric = TRUE
    )
    k <- k+1L

    #### Rotated factor means ####

    list_struct[[k]] <- list(
      name = alpha_group[i],
      type = "matrix",
      dim = c(q, 1L),
      rownames = factor_names,
      colnames = "intrcp"
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
                                   dots = dots, trans = trans)

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
    rotation_spec = dataList$rotation_spec,
    rotation_components = lapply(estimators, FUN = \(x) x$component),
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

    manifold <- dataList$projection

    if(manifold == "poblq" && !is.null(extra$oblique)) {
      extra$oblique <- check_oblique_lrotate(extra$oblique, q)
      manifold <- "poblq_blocks"
    }

    extra$p <- q
    extra$q <- q

    manifolds[[i]] <- list(
      manifold = manifold,
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

estimators_lrotate <- function(dataList, data_param, dots, trans) {

  lambda_group <- data_param$lambda_group
  psi_group <- data_param$psi_group
  rotation <- dataList$rotation_spec
  if(is.null(rotation)) rotation <- dataList$rotation
  components <- rotation_components_lrotate(rotation)

  estimators <- vector("list", length = dataList$ngroups*length(components))
  k <- 1L

  for(i in seq_len(dataList$ngroups)) {

    p <- dataList$nitems[[i]]
    q <- dataList$nfactors[[i]]

    for(j in seq_along(components)) {

      criterion <- names(components)[j]
      extra <- rotation_component_dots_lrotate(criterion, components[[j]], dots)
      extra <- group_dots_lrotate(dots = extra, group_index = i,
                                  ngroups = dataList$ngroups)
      extra <- rotation_defaults_lrotate(rotation = criterion, dots = extra)

      items <- rotation_indices_lrotate(extra$items, p, dataList$item_label[[i]],
                                         "items")
      factors <- rotation_indices_lrotate(extra$factors, q,
                                           dataList$factor_label[[i]], "factors")
      extra$items <- extra$factors <- NULL
      extra$p <- length(items)
      extra$q <- length(factors)
      extra <- rotation_extra_lrotate(criterion, extra, items, factors, p, q)

      # Supply label submatrices to the existing estimator-index machinery.
      # No additional parameters or transformations are introduced. Overlapping
      # components contribute additively to the same gradient/Hessian entries.
      parameters <- list(trans[[lambda_group[i]]][items, factors, drop = FALSE],
                          trans[[psi_group[i]]][factors, factors, drop = FALSE])

      estimators[[k]] <- list(
        estimator = criterion,
        parameters = parameters,
        extra = extra,
        component = list(group_index = i, group_label = dataList$group_label[i],
                          component_index = j, criterion = criterion,
                          items = items, factors = factors, extra = extra)
      )
      k <- k+1L

    }

  }

  #### Result ####

  return(estimators)

}

#### Function to select group-specific extra arguments ####

group_dots_lrotate <- function(dots, group_index, ngroups) {

  extra <- dots

  # Group-specific projection/criterion arguments may be supplied either once
  # for every group or as a list with one object per group.
  group_objects <- c("constraints", "oblique", "target", "weight",
                     "psitarget", "psiweight", "items", "factors",
                     "gamma", "alpha", "epsilon", "k", "w")

  for(nm in intersect(names(extra), group_objects)) {

    object <- extra[[nm]]

    if(is.list(object) && !is.data.frame(object)) {

      # A flat oblique list contains factor positions, not group entries.
      # Group-specific blocks use an additional (outer) list level.
      if(nm == "oblique" && !any(vapply(object, is.list, logical(1L)))) {
        next
      }

      if(length(object) != ngroups) {
        stop("Group-specific '", nm, "' must have one entry per group")
      }
      extra[nm] <- list(object[[group_index]])
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
