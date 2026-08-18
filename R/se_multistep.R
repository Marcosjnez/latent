# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 18/08/2026
#'
#' Standard Errors for Multistep Latent Models
#'
#' Compute standard errors for deterministic multistep estimators by propagating
#' the joint variance-covariance matrix of parameters estimated in previous
#' steps.
#'
#' @param fit A fitted object inheriting from class \code{"multistep"}.
#' @param type Optional compatibility argument. If supplied, it is validated but
#'   does not replace the sample-statistic covariance method selected when the
#'   multistep object was fitted.
#' @param parameters Optional parameter specification identifying the parameters
#'   or transformed parameters for which standard errors should be returned. If
#'   \code{NULL}, the parameter structures of the top-level model are used.
#' @param digits Non-negative integer specifying the number of decimal places
#'   used in formatted parameter tables. If \code{NULL}, values are not rounded.
#' @param ... Additional arguments reserved for future multistep methods.
#'
#' @details
#' A multistep estimator is deterministic conditional on the parameter estimates
#' passed from preceding steps. Consequently, no inverse-Hessian variance is
#' added for the new step.
#'
#' Let \eqn{A} be the joint variance-covariance matrix of parameters estimated
#' in all preceding steps, \eqn{H_2} the Hessian of parameters estimated in the
#' current step, and \eqn{C} the cross-Hessian between preceding and current
#' parameters in the corresponding unrestricted model. The covariance matrix is
#' enlarged at every step as
#' \deqn{
#' V =
#' \begin{pmatrix}
#' A & -A C H_2^{-1} \\
#' -H_2^{-1} C^\top A &
#' H_2^{-1} C^\top A C H_2^{-1}
#' \end{pmatrix}.
#' }
#'
#' The final joint covariance matrix is passed to \code{vcov.latent()}, which
#' propagates it through all requested parameter transformations.
#'
#' @return A list containing parameter tables, standard errors, the selected
#'   variance-covariance matrix, the final joint variance-covariance matrix, and
#'   the Hessian and cross-derivative matrices used in the final step.
#'
#' @method se multistep
#' @export
se.multistep <- function(fit, type = NULL, parameters = NULL,
                           digits = 4L, ...) {

  #### Check inputs ####

  if(!inherits(fit, "multistep")) {
    stop("fit must inherit from class 'multistep'.")
  }

  if(length(fit@Optim$transparameters) == 0L) {
    stop("The multistep object has not been fitted.")
  }

  if(!is.null(type)) {
    type <- match.arg(tolower(type),
                      c("standard", "information", "robust"))
  }

  if(is.null(parameters)) {

    parameters <- fit@modelInfo$trans[names(fit@modelInfo$param)]

  } else {

    selected_parameters <- unique(unlist(parameters))

    if(!all(selected_parameters %in%
            fit@modelInfo$transparameters_labels)) {
      stop("Unknown parameters.")
    }

  }

  if(!is.null(digits) &&
     (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits < 0L || digits != as.integer(digits))) {
    stop("digits must be NULL or a non-negative integer.")
  }

  #### Order model-estimation steps ####

  stages <- multistep_models(fit)

  v <- NULL
  joint_labels <- character(0L)
  step_results <- vector("list", length(stages))

  H2 <- matrix(numeric(0L), nrow = 0L, ncol = 0L)
  A <- matrix(numeric(0L), nrow = 0L, ncol = 0L)
  B <- matrix(numeric(0L), nrow = 0L, ncol = 0L)
  C <- matrix(numeric(0L), nrow = 0L, ncol = 0L)
  fit_full <- NULL

  #### Increase the joint covariance matrix step by step ####

  for(i in seq_along(stages)) {

    stage <- stages[[i]]
    current_labels <- stage@modelInfo$parameters_labels

    fit_full <- unrestricted_model_multistep(stage)
    full_labels <- fit_full@modelInfo$parameters_labels

    missing_current <- setdiff(current_labels, full_labels)

    if(length(missing_current) > 0L) {
      stop("The unrestricted model does not contain current-step parameter(s): ",
           paste(missing_current, collapse = ", "))
    }

    previous_labels <- full_labels[
      !(full_labels %in% current_labels)
    ]

    #### Joint covariance from preceding steps ####

    if(i == 1L) {

      A <- initial_vcov_multistep(stage = stage,
                                  parameters = previous_labels)

    } else {

      missing_previous <- setdiff(previous_labels, joint_labels)

      if(length(missing_previous) > 0L) {
        stop("The accumulated covariance matrix does not contain preceding-step ",
             "parameter(s): ",
             paste(missing_previous, collapse = ", "))
      }

      A <- v[previous_labels, previous_labels, drop = FALSE]

    }

    A <- (A+t(A))/2

    #### Hessian of the current step ####

    H2 <- hessian(stage)
    H2 <- H2[current_labels, current_labels, drop = FALSE]
    H2 <- (H2+t(H2))/2

    H2_inv <- invert_hessian_multistep(H2,
                                       object_name = paste0("step ", i))

    #### Full Hessian and cross derivatives ####

    H_full <- hessian(fit_full)
    H_full <- H_full[full_labels, full_labels, drop = FALSE]
    H_full <- (H_full+t(H_full))/2

    if(length(previous_labels) > 0L) {

      C <- H_full[previous_labels, current_labels, drop = FALSE]

      B <- t(C) %*% A %*% C
      B <- (B+t(B))/2

      v12 <- -A %*% C %*% H2_inv
      v22 <- H2_inv %*% B %*% H2_inv
      v22 <- (v22+t(v22))/2

    } else {

      C <- matrix(0,
                  nrow = 0L,
                  ncol = length(current_labels),
                  dimnames = list(character(0L), current_labels))

      B <- matrix(0,
                  nrow = length(current_labels),
                  ncol = length(current_labels),
                  dimnames = list(current_labels, current_labels))

      v12 <- matrix(0,
                    nrow = 0L,
                    ncol = length(current_labels),
                    dimnames = list(character(0L), current_labels))

      v22 <- B

    }

    #### Enlarged joint covariance ####

    new_labels <- c(previous_labels, current_labels)

    v_new <- matrix(0,
                    nrow = length(new_labels),
                    ncol = length(new_labels),
                    dimnames = list(new_labels, new_labels))

    if(length(previous_labels) > 0L) {
      v_new[previous_labels, previous_labels] <- A
      v_new[previous_labels, current_labels] <- v12
      v_new[current_labels, previous_labels] <- t(v12)
    }

    v_new[current_labels, current_labels] <- v22

    v <- (v_new+t(v_new))/2
    joint_labels <- new_labels

    step_results[[i]] <- list(
      class = class(stage)[1L],
      previous_parameters = previous_labels,
      current_parameters = current_labels,
      H = H2,
      C = C,
      B = B,
      joint_vcov = v
    )

  }

  #### Top-level variance-covariance matrix ####

  final_labels <- fit_full@modelInfo$parameters_labels

  if(!setequal(final_labels, joint_labels)) {
    stop("The final joint covariance matrix does not match the unrestricted ",
         "top-level model.")
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

  free_labels <- fit@modelInfo$parameters_labels
  free_VCOV <- v[free_labels, free_labels, drop = FALSE]

  sample_se <- fit@dataList$se_type

  if(is.null(sample_se)) {
    sample_se <- NA_character_
  }

  #### Result ####

  result <- list(table = table,
                 table_se = table_se,
                 se = c(VCOV$se),
                 vcov = VCOV$vcov,
                 VCOV = VCOV$vcov,
                 free_VCOV = free_VCOV,
                 jacob = VCOV$jacob,
                 H = H2,
                 B = B,
                 A = A,
                 C = C,
                 joint_vcov = v,
                 steps = step_results,
                 type = "multistep",
                 sample_se = sample_se)

  return(result)

}

# Explicit concrete-class method to make S3 dispatch deterministic for S4
# multiple inheritance.
#'
#' @rdname se.multistep
#' @method se multistep_lcfa
#' @export
se.multistep_lcfa <- function(fit, ...) {

  result <- se.multistep(fit, ...)

  #### Result ####

  return(result)

}

#### Initial covariance of preceding fitted models ####

initial_vcov_multistep <- function(stage, parameters) {

  result <- matrix(0,
                   nrow = length(parameters),
                   ncol = length(parameters),
                   dimnames = list(parameters, parameters))

  if(length(parameters) == 0L) {

    #### Result ####

    return(result)

  }

  assigned <- setNames(rep(FALSE, length(parameters)), parameters)

  #### Covariance blocks stored in previous latent objects ####

  models <- collect_latent_models_multistep(stage@extra)

  for(model in models) {

    model_labels <- intersect(parameters,
                              model@modelInfo$parameters_labels)

    if(length(model_labels) == 0L) {
      next
    }

    model_VCOV <- stored_vcov_multistep(model)
    model_VCOV <- model_VCOV[model_labels, model_labels, drop = FALSE]

    overlap <- model_labels[assigned[model_labels]]

    if(length(overlap) > 0L) {
      existing <- result[overlap, overlap, drop = FALSE]
      proposed <- model_VCOV[overlap, overlap, drop = FALSE]

      if(!isTRUE(all.equal(existing, proposed,
                           tolerance = sqrt(.Machine$double.eps)))) {
        stop("More than one preceding model supplies incompatible covariance ",
             "information for parameter(s): ",
             paste(overlap, collapse = ", "))
      }
    }

    result[model_labels, model_labels] <- model_VCOV
    assigned[model_labels] <- TRUE

  }

  #### Numeric covariance blocks stored in the current model ####

  blocks <- sample_vcov_blocks_multistep(stage)

  for(block in blocks) {

    if(!is.matrix(block)) {
      block <- as.matrix(block)
    }

    if(is.null(rownames(block)) || is.null(colnames(block))) {
      next
    }

    block_labels <- intersect(parameters, rownames(block))
    block_labels <- intersect(block_labels, colnames(block))

    if(length(block_labels) == 0L) {
      next
    }

    block <- block[block_labels, block_labels, drop = FALSE]
    overlap <- block_labels[assigned[block_labels]]

    if(length(overlap) > 0L) {
      existing <- result[overlap, overlap, drop = FALSE]
      proposed <- block[overlap, overlap, drop = FALSE]

      if(!isTRUE(all.equal(existing, proposed,
                           tolerance = sqrt(.Machine$double.eps)))) {
        stop("Stored model objects and numeric sample-statistic covariance ",
             "matrices disagree for parameter(s): ",
             paste(overlap, collapse = ", "))
      }
    }

    # Insert the complete block, including covariances between labels that
    # were already represented and labels supplied for the first time.
    result[block_labels, block_labels] <- block
    assigned[block_labels] <- TRUE

  }

  missing <- names(assigned)[!assigned]

  if(length(missing) > 0L) {
    stop("No preceding-step covariance matrix was found for parameter(s): ",
         paste(missing, collapse = ", "))
  }

  result <- (result+t(result))/2

  #### Result ####

  return(result)

}

#### Retrieve a stored covariance matrix ####

stored_vcov_multistep <- function(fit) {

  labels <- fit@modelInfo$parameters_labels

  VCOV <- tryCatch(fit@Optim$SE$VCOV,
                   error = function(e) NULL)

  if(is.null(VCOV)) {
    VCOV <- information(fit)$VCOV
  }

  if(!is.matrix(VCOV)) {
    VCOV <- as.matrix(VCOV)
  }

  if(nrow(VCOV) != length(labels) ||
     ncol(VCOV) != length(labels) ||
     !isSymmetric(VCOV) ||
     any(!is.finite(VCOV))) {
    stop("A preceding latent object contains an invalid VCOV matrix.")
  }

  if(!is.null(rownames(VCOV)) || !is.null(colnames(VCOV))) {

    if(is.null(rownames(VCOV)) || is.null(colnames(VCOV)) ||
       !all(labels %in% rownames(VCOV)) ||
       !all(labels %in% colnames(VCOV))) {
      stop("The names of a preceding VCOV matrix do not match its ",
           "parameter labels.")
    }

    VCOV <- VCOV[labels, labels, drop = FALSE]

  }

  rownames(VCOV) <- colnames(VCOV) <- labels
  VCOV <- (VCOV+t(VCOV))/2

  #### Result ####

  return(VCOV)

}

#### Sample-statistic covariance blocks stored in an lcfa object ####

sample_vcov_blocks_multistep <- function(fit) {

  if(!inherits(fit, "lcfa") ||
     is.null(fit@dataList$data_param)) {

    #### Result ####

    return(list())

  }

  data_param <- fit@dataList$data_param
  blocks <- list()

  if(isTRUE(fit@modelInfo$control_optimizer$meanstructure)) {
    blocks <- c(blocks,
                flatten_matrix_list_multistep(data_param$VCOV_means))
  }

  blocks <- c(blocks,
              flatten_matrix_list_multistep(data_param$VCOV_cov))

  #### Result ####

  return(blocks)

}

flatten_matrix_list_multistep <- function(x) {

  if(is.matrix(x)) {

    result <- list(x)

  } else if(is.list(x)) {

    result <- unlist(lapply(x, flatten_matrix_list_multistep),
                     recursive = FALSE)

  } else {

    result <- list()

  }

  #### Result ####

  return(result)

}

#### Collect leaf latent objects ####

collect_latent_models_multistep <- function(objects) {

  result <- list()

  visit <- function(object) {

    if(!inherits(object, "latent")) {
      return(invisible(NULL))
    }

    children <- object@extra[
      vapply(object@extra,
             FUN = \(x) inherits(x, "latent"),
             FUN.VALUE = logical(1L))
    ]

    if(length(children) > 0L) {

      for(child in children) {
        visit(child)
      }

    } else {

      result[[length(result)+1L]] <<- object

    }

    return(invisible(NULL))

  }

  for(object in objects) {
    visit(object)
  }

  result <- unique_latent_models_multistep(result)

  #### Result ####

  return(result)

}

#### Order multistep model stages ####

multistep_models <- function(fit) {

  stages <- list(fit)
  current <- fit

  repeat {

    previous <- current@extra[
      vapply(current@extra,
             FUN = \(x) inherits(x, "multistep"),
             FUN.VALUE = logical(1L))
    ]

    previous <- unique_latent_models_multistep(previous)

    if(length(previous) == 0L) {
      break
    }

    previous <- immediate_previous_model_multistep(previous)

    duplicated <- any(vapply(stages,
                             FUN = \(x) identical(x, previous),
                             FUN.VALUE = logical(1L)))

    if(duplicated) {
      stop("A cycle was detected among the models stored in extra.")
    }

    stages <- c(list(previous), stages)
    current <- previous

  }

  #### Result ####

  return(stages)

}

immediate_previous_model_multistep <- function(models) {

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

      if(!is_ancestor_multistep(models[[i]], models[[j]])) {
        candidate[i] <- FALSE
        break
      }

    }

  }

  if(sum(candidate) != 1L) {
    stop("The extra slot contains more than one independent multistep model. ",
         "A unique sequential path is required.")
  }

  result <- models[[which(candidate)]]

  #### Result ####

  return(result)

}

is_ancestor_multistep <- function(fit, ancestor) {

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

    previous <- current@extra[
      vapply(current@extra,
             FUN = \(x) inherits(x, "multistep"),
             FUN.VALUE = logical(1L))
    ]

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

unique_latent_models_multistep <- function(models) {

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

#### Reconstruct an unrestricted model ####

unrestricted_model_multistep <- function(fit) {

  if(inherits(fit, "lcfa")) {

    result <- unrestricted_lcfa_multistep(fit)

    #### Result ####

    return(result)

  }

  args <- fit@dataList$args

  if(is.null(args) || !is.list(args)) {

    call_args <- as.list(fit@call)[-1L]

    args <- tryCatch(
      lapply(call_args, eval, envir = parent.frame()),
      error = function(e) NULL
    )

  }

  if(is.null(args) || !is.list(args)) {
    stop("The original fitting arguments are missing. Store them in ",
         "dataList$args when creating a multistep object.")
  }

  if("model" %in% names(args)) {
    args$model <- remove_previous_models_multistep(args$model)
  }

  fit_name <- tail(as.character(fit@call[[1L]]), 1L)
  fit_function <- get0(fit_name, mode = "function", inherits = TRUE)

  if(is.null(fit_function)) {
    stop("Could not identify the fitting function used to create this object.")
  }

  function_args <- names(formals(fit_function))

  if(!("do.fit" %in% function_args)) {
    stop("The fitting function does not provide the do.fit argument required ",
         "to reconstruct an unrestricted model.")
  }

  args$do.fit <- FALSE

  if("se" %in% function_args) {
    args$se <- FALSE
  }

  if("adjustment" %in% function_args) {
    args$adjustment <- "none"
  }

  if("message" %in% function_args) {
    args$message <- FALSE
  }

  if("verbose" %in% function_args) {
    args$verbose <- FALSE
  }

  result <- do.call(fit_function, args)

  if(!inherits(result, "latent")) {
    stop("The unrestricted model is not a latent object.")
  }

  result <- insert_multistep_estimates(result, fit)

  #### Result ####

  return(result)

}

unrestricted_lcfa_multistep <- function(fit) {

  dataList <- fit@dataList
  control_free <- fit@modelInfo$control_optimizer

  control_free$free_taus <- TRUE
  control_free$free_S <- TRUE
  control_free$free_M <- TRUE
  control_free$rstarts <- 1L
  control_free$cores <- 1L
  control_free$start <- NULL

  full_model <- create_lcfa_model(dataList = dataList,
                                  model = dataList$model,
                                  control = control_free)

  modelInfo <- create_lcfa_modelInfo(dataList = dataList,
                                     full_model = full_model,
                                     control = control_free)

  Optim <- fit@Optim
  Optim$parameters <-
    fit@Optim$transparameters[modelInfo$parameters_labels]
  Optim$transparameters <-
    fit@Optim$transparameters[modelInfo$transparameters_labels]

  if(anyNA(Optim$parameters) || anyNA(Optim$transparameters)) {
    stop("The unrestricted lcfa model contains parameters that could not be ",
         "matched to the fitted multistep object.")
  }

  transformed_pars <- fill_in(modelInfo$trans,
                              Optim$transparameters)
  parameters <- transformed_pars[names(modelInfo$param)]

  result <- new("lcfa",
                version          = fit@version,
                call             = fit@call,
                timing           = fit@timing,
                dataList         = dataList,
                modelInfo        = modelInfo,
                Optim            = Optim,
                parameters       = parameters,
                transformed_pars = transformed_pars,
                extra            = list())

  #### Result ####

  return(result)

}

insert_multistep_estimates <- function(result, fit) {

  parameters <- fit@Optim$transparameters[
    result@modelInfo$parameters_labels
  ]

  transparameters <- fit@Optim$transparameters[
    result@modelInfo$transparameters_labels
  ]

  if(anyNA(parameters)) {

    missing <- result@modelInfo$parameters_labels[is.na(parameters)]

    stop("The unrestricted model contains free parameter(s) that are not ",
         "present in the fitted multistep model: ",
         paste(missing, collapse = ", "))

  }

  if(anyNA(transparameters)) {

    missing <- result@modelInfo$transparameters_labels[is.na(transparameters)]

    stop("The unrestricted model contains transformed parameter(s) that are ",
         "not present in the fitted multistep model: ",
         paste(missing, collapse = ", "))

  }

  result@Optim$parameters <- parameters
  result@Optim$transparameters <- transparameters
  result@transformed_pars <- fill_in(result@modelInfo$trans,
                                     transparameters)
  result@parameters <-
    result@transformed_pars[names(result@modelInfo$param)]

  #### Result ####

  return(result)

}

remove_previous_models_multistep <- function(model) {

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

invert_hessian_multistep <- function(H, object_name = "model") {

  if(length(H) == 0L) {

    #### Result ####

    return(H)

  }

  eigenvalues <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  tolerance <- sqrt(.Machine$double.eps)*
    max(1, max(abs(eigenvalues)))

  if(any(eigenvalues <= tolerance)) {
    warning("The Hessian of ", object_name,
            " is not positive definite; an approximate inverse was used.")
  }

  result <- approx_Hinv(H)

  if(any(!is.finite(result))) {
    stop("The Hessian of ", object_name, " could not be inverted.")
  }

  result <- (result+t(result))/2
  dimnames(result) <- dimnames(H)

  #### Result ####

  return(result)

}
