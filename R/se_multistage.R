# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 25/08/2026
#'
#' Standard Errors for Multistage Latent Models
#'
#' Compute standard errors for a multistage latent model by propagating the
#' joint variance-covariance matrix of parameters estimated in previous stages.
#'
#' @param fit A fitted object inheriting from class \code{"multistage"}.
#' @param type Character string specifying the covariance estimator used at each
#'   estimation stage. Available options are \code{"information"} and
#'   \code{"robust"}.
#' @param parameters Optional parameter specification identifying the parameters
#'   or transformed parameters for which standard errors should be returned. If
#'   \code{NULL}, the freely estimated parameters of the top-level model are
#'   used.
#' @param digits Non-negative integer specifying the number of decimal places
#'   used in the formatted parameter tables. If \code{NULL}, values are not
#'   rounded.
#' @param ... Additional arguments reserved for future extensions.
#'
#' @details
#' Previous fitted models are obtained from the \code{extra} slot. The sequence
#' of models is ordered from the ultimate ancestor to the top-level model and
#' processed iteratively.
#'
#' Let \eqn{A} be the joint variance-covariance matrix of the parameters fixed
#' from previous stages, \eqn{P_2} the inverse-Hessian operator of the parameters
#' estimated in the current stage, \eqn{C} the cross-Hessian between
#' previous-stage and current-stage parameters obtained from the corresponding
#' unrestricted model, and \eqn{V_2} the conditional variance-covariance matrix
#' of the current-stage estimator. At each stage the joint covariance matrix is
#' enlarged as
#' \deqn{
#' V =
#' \begin{pmatrix}
#' A & -A C P_2 \\
#' -P_2^\top C^\top A &
#' V_2 + P_2 C^\top A C P_2^\top
#' \end{pmatrix}.
#' }
#'
#' By default, \eqn{P_2=H_2^{-1}}. If the current-stage control contains
#' \code{se_method = "KKT"} and active equality constraints are available,
#' \eqn{P_2} is the parameter block of the inverse KKT matrix.
#'
#' If \code{type = "robust"}, robust covariance matrices are used for
#' \eqn{V_2} whenever a class-specific robust method is available. The ordinary
#' objective Hessian and, when requested, its KKT-constrained inverse operator
#' are still used for propagation of uncertainty between estimation stages.
#'
#' The top-level unrestricted model must contain all parameter labels carried
#' forward from previous stages. This is required to obtain the complete
#' cross-Hessian between the previously estimated parameters and the parameters
#' estimated in the new stage.
#'
#' @return A list containing the formatted estimates and standard errors,
#'   variance-covariance matrix, Hessian of the top-level estimated parameters,
#'   final propagation matrix \code{B}, transformation Jacobian, and the final
#'   joint variance-covariance matrix before propagation to transformed
#'   parameters.
#'
#' @method se multistage
#' @export
se.multistage <- function(fit, type = "information", parameters = NULL,
                          digits = 4L, ...) {

  #### Check inputs ####

  if(!inherits(fit, "multistage")) {
    stop("fit must inherit from class 'multistage'.")
  }

  if(length(fit@Optim$parameters) == 0L) {
    stop("The multistage object has not been fitted.")
  }

  type <- match.arg(tolower(type), c("information", "robust"))

  if(!is.null(digits) &&
     (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits < 0L || digits != as.integer(digits))) {
    stop("digits must be NULL or a non-negative integer.")
  }

  if(is.null(parameters)) {

    # parameters <- fit@modelInfo$parameters_labels
    parameters <- fit@modelInfo$trans[names(fit@modelInfo$param)]

  } else {

    selected_parameters <- unique(unlist(parameters))

    if(!all(selected_parameters %in%
            fit@modelInfo$transparameters_labels)) {
      stop("Unknown parameters.")
    }

  }

  #### Order the estimation stages ####

  stages <- multistage_models(fit)

  if(length(stages) < 2L) {
    stop("No previous estimation stage was found in fit@extra.")
  }

  #### Ultimate ancestor ####

  ancestor <- stages[[1L]]

  v <- conditional_vcov_multistage(ancestor, type = type)
  joint_labels <- ancestor@modelInfo$parameters_labels

  rownames(v) <- colnames(v) <- joint_labels

  H2 <- NULL
  B <- NULL
  C <- NULL
  A <- NULL
  fit_full <- NULL

  #### Increase the joint covariance matrix stage by stage ####

  for(i in 2:length(stages)) {

    stage <- stages[[i]]
    current_labels <- stage@modelInfo$parameters_labels

    #### Unrestricted model for this stage ####

    fit_full <- unrestricted_model_multistage(stage)
    full_labels <- fit_full@modelInfo$parameters_labels

    # Every parameter carried from the preceding iteration must still be
    # represented in the new unrestricted model. A parameter may be estimated
    # again in the current stage, in which case its old estimate is replaced by
    # the current-stage estimate below.
    missing_ancestors <- setdiff(joint_labels, full_labels)

    if(length(missing_ancestors) > 0L) {
      stop("The unrestricted model does not contain previous-stage parameter(s): ",
           paste(missing_ancestors, collapse = ", "))
    }

    missing_current <- setdiff(current_labels, full_labels)

    if(length(missing_current) > 0L) {
      stop("The unrestricted model does not contain current-stage parameter(s): ",
           paste(missing_current, collapse = ", "))
    }

    # Parameters in H_full that are not estimated in the current stage are the
    # parameters fixed from previous stages.
    previous_labels <- full_labels[
      !(full_labels %in% current_labels)
    ]

    missing_previous <- setdiff(previous_labels, joint_labels)

    if(length(missing_previous) > 0L) {
      stop("The covariance matrix from previous stages does not contain parameter(s): ",
           paste(missing_previous, collapse = ", "))
    }

    if(length(previous_labels) == 0L) {
      stop("No parameters estimated in previous stages were identified in the unrestricted model.")
    }

    #### Previous-stage joint covariance ####

    A <- v[previous_labels, previous_labels, drop = FALSE]
    A <- (A+t(A))/2

    #### Conditional covariance of the current stage ####

    v2 <- conditional_vcov_multistage(stage, type = type)
    v2 <- v2[current_labels, current_labels, drop = FALSE]
    v2 <- (v2+t(v2))/2

    #### Current-stage Hessian ####

    H2 <- hessian(stage)
    H2 <- H2[current_labels, current_labels, drop = FALSE]
    H2 <- validate_covariance_matrix(
      H2,
      labels = current_labels,
      object_name = paste0("Hessian of multistage step ", i)
    )
    H2_inv <- invert_hessian_latent(
      fit = stage,
      H = H2,
      labels = current_labels,
      object_name = paste0("Hessian of multistage step ", i)
    )

    #### Full Hessian and cross derivatives ####

    H_full <- hessian(fit_full)
    H_full <- H_full[full_labels, full_labels, drop = FALSE]
    H_full <- validate_covariance_matrix(
      H_full,
      labels = full_labels,
      object_name = paste0("unrestricted Hessian of multistage step ", i)
    )

    C <- H_full[previous_labels, current_labels, drop = FALSE]

    #### Additional variance from previous stages ####

    B <- t(C) %*% A %*% C
    B <- (B+t(B))/2

    propagated <- H2_inv %*% B %*% H2_inv
    propagated <- (propagated+t(propagated))/2

    v22 <- v2 + propagated
    v22 <- (v22+t(v22))/2

    #### Covariance between previous and current parameters ####

    v12 <- -A %*% C %*% H2_inv

    #### Increase the joint covariance matrix ####

    new_labels <- c(previous_labels, current_labels)

    v_new <- matrix(0,
                    nrow = length(new_labels),
                    ncol = length(new_labels),
                    dimnames = list(new_labels, new_labels))

    v_new[previous_labels, previous_labels] <- A
    v_new[previous_labels, current_labels] <- v12
    v_new[current_labels, previous_labels] <- t(v12)
    v_new[current_labels, current_labels] <- v22

    v <- (v_new+t(v_new))/2
    joint_labels <- new_labels

  }

  #### Top-level variance-covariance matrix ####

  # The final unrestricted model contains the complete parameter vector used by
  # the iterative covariance calculation.
  final_labels <- fit_full@modelInfo$parameters_labels

  if(!setequal(final_labels, joint_labels)) {
    stop("The final joint covariance matrix does not match the parameters of the unrestricted top-level model.")
  }

  v <- v[final_labels, final_labels, drop = FALSE]

  VCOV <- vcov(fit_full,
               v = v,
               parameters = parameters)

  #### Parameter table ####

  est <- fill_in(parameters,
                 fit_full@Optim$transparameters,
                 miss = NA)

  table_se <- fill_in(parameters,
                      VCOV$se,
                      miss = NA)

  table <- combine_est_se(est, table_se, digits = digits)

  #### Result ####

  result <- list(table = table,
                 table_se = table_se,
                 se = c(VCOV$se),
                 vcov = VCOV$vcov,
                 VCOV = VCOV$vcov,
                 H = H2,
                 B = B,
                 A = A,
                 C = C,
                 jacob = VCOV$jacob,
                 joint_vcov = v,
                 type = type)

  return(result)

}

#### Conditional covariance matrix of one estimation stage ####

conditional_vcov_multistage <- function(fit, type = "information") {

  if(type == "information") {
    stage_vcov <- information(fit)$VCOV
  } else {
    stage_vcov <- robust(fit)$VCOV
  }

  labels <- fit@modelInfo$parameters_labels

  stage_vcov <- validate_covariance_matrix(
    stage_vcov,
    labels = labels,
    object_name = "conditional multistage variance-covariance matrix"
  )

  #### Result ####

  return(stage_vcov)

}

#### Order models from the ultimate ancestor to the top-level model ####

multistage_models <- function(fit) {

  stages <- list(fit)
  current <- fit

  while(inherits(current, "multistage")) {

    previous <- current@extra
    previous <- previous[
      vapply(previous,
             FUN = \(x) inherits(x, "latent"),
             FUN.VALUE = logical(1L))
    ]

    previous <- unique_latent_models_multistage(previous)

    if(length(previous) == 0L) {
      stop("A multistage object does not contain a previous latent model in its extra slot.")
    }

    previous <- immediate_previous_model_multistage(previous)

    duplicated <- any(vapply(stages,
                             FUN = \(x) identical(x, previous),
                             FUN.VALUE = logical(1L)))

    if(duplicated) {
      stop("A cycle was detected among the models stored in the extra slots.")
    }

    stages <- c(list(previous), stages)
    current <- previous

  }

  #### Result ####

  return(stages)

}

#### Identify the immediately preceding model ####

immediate_previous_model_multistage <- function(models) {

  if(length(models) == 1L) {

    #### Result ####

    return(models[[1L]])

  }

  candidate <- logical(length(models))

  for(i in seq_along(models)) {

    candidate[i] <- TRUE

    for(j in seq_along(models)) {

      if(i == j) {
        next
      }

      if(!is_ancestor_multistage(models[[i]], models[[j]])) {
        candidate[i] <- FALSE
        break
      }

    }

  }

  if(sum(candidate) != 1L) {
    stop("The extra slot contains more than one independent previous model. ",
         "A unique sequential estimation path is required.")
  }

  result <- models[[which(candidate)]]

  #### Result ####

  return(result)

}

#### Check whether one model is an ancestor of another ####

is_ancestor_multistage <- function(fit, ancestor) {

  queue <- list(fit)
  visited <- list()
  result <- FALSE

  while(length(queue) > 0L) {

    current <- queue[[1L]]
    queue <- queue[-1L]

    already_visited <- length(visited) > 0L &&
      any(vapply(visited,
                 FUN = \(x) identical(x, current),
                 FUN.VALUE = logical(1L)))

    if(already_visited) {
      next
    }

    visited[[length(visited)+1L]] <- current

    if(!inherits(current, "multistage")) {
      next
    }

    previous <- current@extra
    previous <- previous[
      vapply(previous,
             FUN = \(x) inherits(x, "latent"),
             FUN.VALUE = logical(1L))
    ]

    previous <- unique_latent_models_multistage(previous)

    if(any(vapply(previous,
                  FUN = \(x) identical(x, ancestor),
                  FUN.VALUE = logical(1L)))) {
      result <- TRUE
      break
    }

    queue <- c(queue, previous)

  }

  #### Result ####

  return(result)

}

#### Remove duplicated latent models ####

unique_latent_models_multistage <- function(models) {

  result <- list()

  for(model in models) {

    duplicated <- length(result) > 0L &&
      any(vapply(result,
                 FUN = \(x) identical(x, model),
                 FUN.VALUE = logical(1L)))

    if(!duplicated) {
      result[[length(result)+1L]] <- model
    }

  }

  #### Result ####

  return(result)

}

#### Reconstruct an unrestricted version of one estimation stage ####

unrestricted_model_multistage <- function(fit) {

  args <- fit@dataList$args

  if(is.null(args) || !is.list(args)) {
    stop("The original fitting arguments are missing from dataList$args.")
  }

  if("model" %in% names(args)) {
    args$model <- remove_previous_models_multistage(args$model)
  }

  #### Fitting function ####

  if(inherits(fit, "llca")) {

    fit_function <- lca

  } else if(inherits(fit, "lcfa")) {

    fit_function <- lcfa

  } else {

    fit_function <- tryCatch(
      eval(fit@call[[1L]], envir = parent.frame()),
      error = function(e) NULL
    )

    if(!is.function(fit_function)) {
      stop("Could not identify the fitting function used to create this latent object.")
    }

  }

  function_args <- names(formals(fit_function))

  if(!("do.fit" %in% function_args)) {
    stop("The fitting function does not provide a do.fit argument required to reconstruct the unrestricted model.")
  }

  args$do.fit <- FALSE

  if("adjustment" %in% function_args) {
    args$adjustment <- "none"
  }

  if("se" %in% function_args) {
    args$se <- FALSE
  }

  if("verbose" %in% function_args) {
    args$verbose <- FALSE
  }

  if("message" %in% function_args) {
    args$message <- FALSE
  }

  result <- do.call(fit_function, args)

  if(!inherits(result, "latent")) {
    stop("The unrestricted model is not a latent object.")
  }

  #### Post-estimation estimator ####

  # An unfitted EM LCA object still contains the EM estimator. Hessians used for
  # multistage inference must be evaluated with the observed-data likelihood.
  if(inherits(result, "llca") &&
     length(result@modelInfo$control_estimator) > 0L) {
    result@modelInfo$control_estimator[[1L]]$estimator <- "lca"
  }

  #### Insert fitted parameter values ####

  parameters <- fit@Optim$transparameters[
    result@modelInfo$parameters_labels
  ]

  transparameters <- fit@Optim$transparameters[
    result@modelInfo$transparameters_labels
  ]

  if(anyNA(parameters)) {

    missing <- result@modelInfo$parameters_labels[is.na(parameters)]

    stop("The unrestricted model contains free parameter(s) that are not present ",
         "in the fitted multistage model: ",
         paste(missing, collapse = ", "))

  }

  if(anyNA(transparameters)) {

    missing <- result@modelInfo$transparameters_labels[is.na(transparameters)]

    stop("The unrestricted model contains transformed parameter(s) that are not ",
         "present in the fitted multistage model: ",
         paste(missing, collapse = ", "))

  }

  result@Optim$parameters <- parameters
  result@Optim$transparameters <- transparameters

  #### Result ####

  return(result)

}

#### Remove previous fitted models from the model argument ####

remove_previous_models_multistage <- function(model) {

  if(is.null(model)) {

    result <- NULL

  } else if(inherits(model, "latent")) {

    result <- NULL

  } else if(is.list(model)) {

    keep <- !vapply(model,
                    FUN = \(x) inherits(x, "latent"),
                    FUN.VALUE = logical(1L))

    result <- model[keep]

    if(length(result) == 0L) {
      result <- NULL
    }

  } else {

    result <- model

  }

  #### Result ####

  return(result)

}
