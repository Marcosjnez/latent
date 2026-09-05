# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 05/09/2026
#'
#' Saturated Multivariate-Normal Moments with Incomplete Data
#'
#' Estimate an unrestricted vector of means and variance-covariance matrix by
#' maximizing the multivariate-normal observed-data likelihood. Missing-data
#' patterns are represented internally by their means, covariance matrices, and
#' frequencies, but only one global mean vector and covariance matrix are
#' estimated per substantive group.
#'
#' @usage
#' lmvnorm(data, group = NULL, variables = NULL,
#'         se = TRUE, do.fit = TRUE, message = FALSE,
#'         control = NULL, ...)
#'
#' @param data A data frame or numeric matrix containing the observed variables.
#' @param group Optional character string identifying a grouping variable.
#' @param variables Optional character vector of observed-variable names, or a
#'   list containing one character vector per group. If \code{NULL}, all numeric
#'   columns other than the grouping variable are used.
#' @param se Logical or character. \code{TRUE}, \code{"standard"}, and
#'   \code{"information"} compute the information covariance matrix.
#'   \code{FALSE} skips its computation. Robust covariance estimation is not
#'   yet implemented for this estimator.
#' @param do.fit Logical. If \code{FALSE}, return the prepared but unfitted
#'   \code{"latent"} object.
#' @param message Logical. Print progress messages.
#' @param control Optional list of optimizer controls. Custom starting values
#'   may be supplied through \code{control$start}.
#' @param ... Additional arguments reserved for future extensions.
#'
#' @details
#' For every missingness pattern \eqn{r}, the likelihood contribution is
#' evaluated from the pattern frequency, mean vector, and covariance matrix.
#' The covariance matrix uses the maximum-likelihood divisor \eqn{n_r}; for a
#' singleton pattern it is the zero matrix.
#'
#' The optimized objective is the total negative observed-data log-likelihood.
#' Pattern means and covariance matrices are fixed transformed parameters, so
#' the information covariance is the inverse Hessian on its natural
#' finite-sample scale.
#'
#' @return An S4 object of class \code{"latent"}. Its public parameter blocks
#'   are named \code{means} and \code{S} for a single group, with group suffixes
#'   for multiple groups. The joint finite-sample covariance of all estimated
#'   moments is stored in \code{Optim$SE$VCOV}.
#'
#' @examples
#' \dontrun{
#' X <- HolzingerSwineford1939[, paste0("x", 1:9)]
#' X$x1[1:20] <- NA
#' fit <- lmvnorm(X)
#' fit@transformed_pars$means
#' fit@transformed_pars$S
#' }
#'
#' @export
lmvnorm <- function(data, group = NULL, variables = NULL,
                    se = TRUE, do.fit = TRUE, message = FALSE,
                    control = NULL, ...) {

  #### Check input arguments ####

  if(!is.data.frame(data) && !is.matrix(data)) {
    stop("data must be a data.frame or matrix")
  }

  if(is.matrix(data)) {
    data <- as.data.frame(data)
  }

  if(nrow(data) < 1L || ncol(data) < 1L) {
    stop("data must contain at least one observation and one variable")
  }

  if(is.null(colnames(data)) ||
     any(colnames(data) == "") ||
     anyDuplicated(colnames(data))) {
    stop("data must have unique, non-empty column names")
  }

  if(!is.null(group) &&
     (!is.character(group) || length(group) != 1L || is.na(group))) {
    stop("group must be NULL or the name of one grouping variable")
  }

  if(!is.null(group) && !(group %in% colnames(data))) {
    stop("The grouping variable is not present in data")
  }

  if(!is.null(control) && !is.list(control)) {
    stop("control must be NULL or a list")
  }

  if(is.null(control)) {
    control <- list()
  }

  if(length(do.fit) != 1L || !is.logical(do.fit) || is.na(do.fit)) {
    stop("do.fit must be TRUE or FALSE")
  }

  if(length(message) != 1L || !is.logical(message) || is.na(message)) {
    stop("message must be TRUE or FALSE")
  }

  se_control <- normalize_lmvnorm_se(se)

  #### Store original call ####

  mc <- match.call(expand.dots = TRUE)
  args <- lapply(as.list(mc)[-1L], eval, envir = parent.frame())

  #### Check control parameters ####

  control <- lmvnorm_control(control)

  #### Create the dataList ####

  dataList <- create_lmvnorm_dataList(data = data,
                                      group = group,
                                      variables = variables,
                                      args = args)

  #### Create the model ####

  full_model <- create_lmvnorm_model(dataList = dataList,
                                     control = control)

  #### Create the manifold, transformation, and estimator structures ####

  modelInfo <- create_lmvnorm_modelInfo(dataList = dataList,
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
    print_lmvnorm_message("Fitting saturated multivariate-normal moments")
  }

  Optim <- fit_lmvnorm(modelInfo = modelInfo)

  #### Standard errors ####

  if(se_control$compute) {

    if(message) {
      print_lmvnorm_message("Computing standard errors")
    }

    Optim$SE <- compute_se_lmvnorm(dataList = dataList,
                                   modelInfo = modelInfo,
                                   Optim = Optim)

  }

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

#### Function to create the dataList ####

create_lmvnorm_dataList <- function(data, group, variables, args) {

  #### Groups ####

  if(is.null(group)) {

    group_label <- ""
    group_index <- rep(1L, nrow(data))

  } else {

    if(anyNA(data[[group]])) {
      stop("The grouping variable cannot contain missing values")
    }

    group_label <- unique(as.character(data[[group]]))

    if(any(group_label == "") || anyDuplicated(group_label)) {
      stop("Group labels must be unique and non-empty")
    }

    group_index <- match(as.character(data[[group]]), group_label)

  }

  ngroups <- length(group_label)

  variables <- normalize_lmvnorm_variables(data = data,
                                            group = group,
                                            group_label = group_label,
                                            variables = variables)

  #### Group data and missing-data patterns ####

  data_per_group <- vector("list", length = ngroups)
  patterns <- vector("list", length = ngroups)
  nobs_group <- vector("list", length = ngroups)
  initial_means <- vector("list", length = ngroups)
  initial_covariance <- vector("list", length = ngroups)

  for(i in seq_len(ngroups)) {

    X <- data[group_index == i, variables[[i]], drop = FALSE]
    X <- X[!apply(is.na(X), MARGIN = 1L, FUN = all), , drop = FALSE]

    if(nrow(X) < 2L) {
      stop("Each group must contain at least two observations with observed data")
    }

    observed_counts <- colSums(!is.na(X))

    if(any(observed_counts < 2L)) {
      variables_i <- names(observed_counts)[observed_counts < 2L]
      stop("Every modeled variable must contain at least two observed values: ",
           paste(variables_i, collapse = ", "))
    }

    observed_variances <- vapply(X,
                                 FUN = \(x) stats::var(x, na.rm = TRUE),
                                 FUN.VALUE = numeric(1L))

    if(any(!is.finite(observed_variances) |
           observed_variances <= 0)) {
      variables_i <- names(observed_variances)[
        !is.finite(observed_variances) | observed_variances <= 0
      ]
      stop("Every modeled variable must have a positive observed variance: ",
           paste(variables_i, collapse = ", "))
    }

    coverage <- crossprod(!is.na(as.matrix(X)))

    if(any(coverage[lower.tri(coverage)] == 0L)) {

      missing_pairs <- which(coverage == 0L & lower.tri(coverage),
                             arr.ind = TRUE)

      pair_names <- apply(missing_pairs, MARGIN = 1L,
                          FUN = \(z) paste(colnames(X)[z], collapse = " / "))

      stop("Every pair of modeled variables must be observed together at ",
           "least once. Zero-coverage pair(s): ",
           paste(pair_names, collapse = ", "))

    }

    data_per_group[[i]] <- X
    patterns[[i]] <- create_lmvnorm_patterns(X)
    nobs_group[[i]] <- nrow(X)
    initial_means[[i]] <- colMeans(X, na.rm = TRUE)
    initial_covariance[[i]] <- initial_lmvnorm_covariance(X)

  }

  #### Result ####

  result <- list(data = data,
                 data_per_group = data_per_group,
                 group = group,
                 group_label = group_label,
                 ngroups = ngroups,
                 variables = variables,
                 item_names = variables,
                 nitems = lapply(variables, length),
                 nobs = sum(unlist(nobs_group)),
                 nobs_group = nobs_group,
                 npatterns = lapply(patterns, length),
                 patterns = patterns,
                 initial_means = initial_means,
                 initial_covariance = initial_covariance,
                 args = args)

  return(result)

}

normalize_lmvnorm_variables <- function(data, group, group_label,
                                        variables) {

  ngroups <- length(group_label)

  if(is.null(variables)) {

    candidates <- setdiff(colnames(data), group)
    candidates <- candidates[vapply(data[candidates], is.numeric,
                                    FUN.VALUE = logical(1L))]

    if(length(candidates) == 0L) {
      stop("No numeric observed variables were found in data")
    }

    variables <- replicate(ngroups, candidates, simplify = FALSE)

  } else if(is.character(variables)) {

    variables <- replicate(ngroups, variables, simplify = FALSE)

  } else if(is.list(variables)) {

    if(!is.null(names(variables)) && ngroups > 1L &&
       all(group_label %in% names(variables))) {
      variables <- variables[group_label]
    }

    if(length(variables) != ngroups) {
      stop("variables must contain one character vector per group")
    }

  } else {

    stop("variables must be NULL, a character vector, or a list")

  }

  for(i in seq_len(ngroups)) {

    vars <- variables[[i]]

    if(!is.character(vars) || length(vars) == 0L ||
       anyNA(vars) || any(vars == "") || anyDuplicated(vars)) {
      stop("Each variables entry must contain unique, non-empty names")
    }

    if(!is.null(group) && group %in% vars) {
      stop("The grouping variable cannot also be a modeled variable")
    }

    unknown <- setdiff(vars, colnames(data))

    if(length(unknown) > 0L) {
      stop("Unknown observed variable(s): ",
           paste(unknown, collapse = ", "))
    }

    nonnumeric <- vars[!vapply(data[vars], is.numeric,
                               FUN.VALUE = logical(1L))]

    if(length(nonnumeric) > 0L) {
      stop("All modeled variables must be numeric: ",
           paste(nonnumeric, collapse = ", "))
    }

    variables[[i]] <- vars

  }

  names(variables) <- if(ngroups > 1L) group_label else NULL

  #### Result ####

  return(variables)

}

create_lmvnorm_patterns <- function(data) {

  missing_matrix <- is.na(data)
  pattern_key <- apply(missing_matrix, MARGIN = 1L,
                       FUN = \(x) paste(as.integer(x), collapse = ""))
  pattern_levels <- unique(pattern_key)
  frequencies <- tabulate(match(pattern_key, pattern_levels),
                          nbins = length(pattern_levels))
  ord <- order(-frequencies, seq_along(pattern_levels))
  pattern_levels <- pattern_levels[ord]

  result <- vector("list", length = length(pattern_levels))

  for(i in seq_along(pattern_levels)) {

    rows <- which(pattern_key == pattern_levels[i])
    observed <- !missing_matrix[rows[1L], ]
    observed_names <- colnames(data)[observed]
    X <- as.matrix(data[rows, observed, drop = FALSE])
    n <- nrow(X)
    means <- colMeans(X)

    if(n > 1L) {
      centered <- sweep(X, MARGIN = 2L, STATS = means, FUN = "-")
      S <- crossprod(centered)/n
    } else {
      S <- matrix(0, nrow = ncol(X), ncol = ncol(X))
    }

    rownames(S) <- colnames(S) <- observed_names
    names(means) <- observed_names

    result[[i]] <- list(observed = observed,
                        observed_names = observed_names,
                        nobs = n,
                        means = means,
                        covariance = S)

  }

  #### Result ####

  return(result)

}

initial_lmvnorm_covariance <- function(data) {

  S <- suppressWarnings(stats::cov(data,
                                   use = "pairwise.complete.obs"))
  S <- as.matrix(S)
  dimnames(S) <- list(colnames(data), colnames(data))

  S[!is.finite(S)] <- 0
  S <- (S+t(S))/2

  for(i in seq_len(ncol(data))) {

    variance_i <- suppressWarnings(stats::var(data[[i]], na.rm = TRUE))

    if(!is.finite(variance_i) || variance_i <= 0) {
      variance_i <- 1
    }

    if(!is.finite(S[i, i]) || S[i, i] <= 0) {
      S[i, i] <- variance_i
    }

  }

  result <- make_positive_definite_lmvnorm(S)

  #### Result ####

  return(result)

}

make_positive_definite_lmvnorm <- function(S, eps = 1e-05) {

  S <- (S+t(S))/2
  eig <- eigen(S, symmetric = TRUE)
  eig$values <- pmax(eig$values, eps)

  result <- eig$vectors %*%
    diag(eig$values, nrow = length(eig$values)) %*%
    t(eig$vectors)
  result <- (result+t(result))/2
  dimnames(result) <- dimnames(S)

  #### Result ####

  return(result)

}

#### Function to create the model ####

create_lmvnorm_model <- function(dataList, control) {

  data_param <- create_lmvnorm_data_param(dataList = dataList)

  #### Model for the transformed parameters ####

  trans <- model_lmvnorm(dataList = dataList,
                         data_param = data_param)

  #### Model for the parameters ####

  param <- constraints_lmvnorm(dataList = dataList,
                               data_param = data_param,
                               trans = trans)

  #### Create the initial values for the parameters ####

  init_param <- start_lmvnorm(dataList = dataList,
                              data_param = data_param,
                              trans = trans,
                              control = control)

  #### Custom initial values ####

  init_param <- custom_init_param(control$start, init_param)

  #### Result ####

  result <- list(param = param,
                 trans = trans,
                 init_param = init_param,
                 data_param = data_param)

  return(result)

}

create_lmvnorm_data_param <- function(dataList) {

  if(dataList$ngroups < 2L) {
    suffix <- ""
  } else {
    suffix <- paste0(".", dataList$group_label)
  }

  pattern_means_group <- vector("list", length = dataList$ngroups)
  pattern_S_group <- vector("list", length = dataList$ngroups)

  for(i in seq_len(dataList$ngroups)) {

    pattern_id <- seq_along(dataList$patterns[[i]])
    pattern_means_group[[i]] <-
      paste0("means", suffix[i], ".pattern", pattern_id)
    pattern_S_group[[i]] <-
      paste0("S", suffix[i], ".pattern", pattern_id)

  }

  #### Result ####

  result <- list(means_group = paste0("means", suffix),
                 S_group = paste0("S", suffix),
                 pattern_means_group = pattern_means_group,
                 pattern_S_group = pattern_S_group)

  return(result)

}

model_lmvnorm <- function(dataList, data_param) {

  list2env(data_param, envir = environment())

  list_struct <- list()
  k <- 1L

  for(i in seq_len(dataList$ngroups)) {

    items <- dataList$item_names[[i]]
    p <- length(items)

    list_struct[[k]] <- list(name = means_group[i],
                             type = "matrix",
                             dim = c(p, 1L),
                             rownames = items,
                             colnames = "intrcpt")
    k <- k+1L

    list_struct[[k]] <- list(name = S_group[i],
                             type = "matrix",
                             dim = c(p, p),
                             rownames = items,
                             colnames = items,
                             symmetric = TRUE)
    k <- k+1L

    for(j in seq_along(dataList$patterns[[i]])) {

      items_pattern <- dataList$patterns[[i]][[j]]$observed_names
      p_pattern <- length(items_pattern)

      list_struct[[k]] <- list(name = pattern_means_group[[i]][j],
                               type = "matrix",
                               dim = c(p_pattern, 1L),
                               rownames = items_pattern,
                               colnames = "intrcpt")
      k <- k+1L

      list_struct[[k]] <- list(name = pattern_S_group[[i]][j],
                               type = "matrix",
                               dim = c(p_pattern, p_pattern),
                               rownames = items_pattern,
                               colnames = items_pattern,
                               symmetric = TRUE)
      k <- k+1L

    }

  }

  #### Result ####

  result <- create_parameters(list_struct)

  return(result)

}

constraints_lmvnorm <- function(dataList, data_param, trans) {

  list2env(data_param, envir = environment())

  param <- trans

  for(i in seq_len(dataList$ngroups)) {

    for(j in seq_along(dataList$patterns[[i]])) {

      pattern <- dataList$patterns[[i]][[j]]

      means <- matrix(pattern$means,
                      ncol = 1L,
                      dimnames = dimnames(
                        trans[[pattern_means_group[[i]][j]]]
                      ))

      S <- pattern$covariance
      dimnames(S) <- dimnames(trans[[pattern_S_group[[i]][j]]])

      param[[pattern_means_group[[i]][j]]] <- means
      param[[pattern_S_group[[i]][j]]] <- S

    }

  }

  #### Result ####

  return(param)

}

start_lmvnorm <- function(dataList, data_param, trans, control) {

  list2env(data_param, envir = environment())

  init_param <- vector("list", length = control$rstarts)

  for(rs in seq_len(control$rstarts)) {

    init_param[[rs]] <- list()

    for(i in seq_len(dataList$ngroups)) {

      means <- matrix(dataList$initial_means[[i]],
                      ncol = 1L,
                      dimnames = dimnames(trans[[means_group[i]]]))

      S <- dataList$initial_covariance[[i]]
      dimnames(S) <- dimnames(trans[[S_group[i]]])

      init_param[[rs]][[means_group[i]]] <- means
      init_param[[rs]][[S_group[i]]] <- S

    }

  }

  #### Result ####

  return(init_param)

}

#### Function to create the modelInfo ####

create_lmvnorm_modelInfo <- function(dataList, full_model, control) {

  list2env(full_model, envir = environment())

  #### Manifolds ####

  manifolds <- manifolds_lmvnorm(dataList = dataList,
                                 data_param = data_param)

  control_manifold <- create_manifolds(manifolds = manifolds,
                                       structures = param)

  #### Transformations ####

  transforms <- transformations_lmvnorm()

  control_transform <- create_transforms(transforms = transforms,
                                         structures = trans)

  #### Estimators ####

  estimators <- estimators_lmvnorm(dataList = dataList,
                                   data_param = data_param,
                                   trans = trans)

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
                    step_labels = parameters_labels,
                    propagate_uncertainty = TRUE,
                    dof = 0L,
                    data_param = data_param,
                    control_manifold = control_manifold,
                    control_transform = control_transform,
                    control_estimator = control_estimator,
                    control_optimizer = control_optimizer)

  return(modelInfo)

}

manifolds_lmvnorm <- function(dataList, data_param) {

  list2env(data_param, envir = environment())

  manifolds <- list()
  k <- 1L

  for(i in seq_len(dataList$ngroups)) {

    manifolds[[k]] <- list(manifold = "euclidean",
                           parameters = means_group[i])
    k <- k+1L

    manifolds[[k]] <- list(manifold = "euclidean",
                           parameters = S_group[i])
    k <- k+1L

  }

  #### Result ####

  return(manifolds)

}

transformations_lmvnorm <- function() {

  #### Result ####

  return(list())

}

estimators_lmvnorm <- function(dataList, data_param, trans) {

  list2env(data_param, envir = environment())

  estimators <- list()
  k <- 1L

  for(i in seq_len(dataList$ngroups)) {

    for(j in seq_along(dataList$patterns[[i]])) {

      pattern <- dataList$patterns[[i]][[j]]
      pick <- pattern$observed_names

      model_covariance <-
        trans[[S_group[i]]][pick, pick, drop = FALSE]
      sample_covariance <- trans[[pattern_S_group[[i]][j]]]
      model_means <-
        trans[[means_group[i]]][pick, , drop = FALSE]
      sample_means <- trans[[pattern_means_group[[i]][j]]]

      estimators[[k]] <- list(
        estimator = "cfa_means_ml",
        parameters = list(model_covariance,
                          sample_covariance,
                          model_means,
                          sample_means),
        extra = list(
          p = length(pick),
          n = pattern$nobs,
          w = pattern$nobs/dataList$nobs,
          double_names = c("loss", "loglik"),
          matrix_names = c("covariance_residuals", "mean_residuals")
        )
      )
      k <- k+1L

    }

  }

  #### Result ####

  return(estimators)

}

#### Function to fit the model ####

fit_lmvnorm <- function(modelInfo) {

  control_optimizer <- modelInfo$control_optimizer
  control_optimizer$cores <-
    min(control_optimizer$rstarts,
        control_optimizer$cores)

  Optim <- optimizer(control_manifold = modelInfo$control_manifold,
                     control_transform = modelInfo$control_transform,
                     control_estimator = modelInfo$control_estimator,
                     control_optimizer = control_optimizer)

  names(Optim$parameters) <- modelInfo$parameters_labels
  names(Optim$transparameters) <- modelInfo$transparameters_labels

  if(!is.null(Optim$g)) names(Optim$g) <- modelInfo$parameters_labels
  if(!is.null(Optim$rg)) names(Optim$rg) <- modelInfo$parameters_labels
  if(!is.null(Optim$dir)) names(Optim$dir) <- modelInfo$parameters_labels

  #### Result ####

  return(Optim)

}

#### Function to compute standard errors ####

compute_se_lmvnorm <- function(dataList, modelInfo, Optim) {

  # control_optimizer <- modelInfo$control_optimizer
  # control_optimizer$parameters[[1L]] <- Optim$parameters
  # control_optimizer$transparameters[[1L]] <- Optim$transparameters

  labels <- modelInfo$parameters_labels

  # H <- get_hess(
  #   modelInfo$control_manifold,
  #   modelInfo$control_transform,
  #   modelInfo$control_estimator,
  #   control_optimizer
  # )$h
  H <- get_hess(fit)$h
  H <- validate_covariance_matrix(
    H,
    labels = labels,
    object_name = "saturated-moment Hessian"
  )

  VCOV <- invert_information_matrix(
    H,
    labels = labels,
    object_name = "saturated-moment Hessian"
  )

  se <- standard_errors_from_vcov(
    VCOV,
    object_name = "saturated-moment variance-covariance matrix"
  )

  #### Result ####

  result <- list(
    H = H,
    VCOV = VCOV,
    se = se,
    type = "information"
  )

  return(result)

}

normalize_lmvnorm_se <- function(se) {

  if(isTRUE(se)) {

    result <- list(type = "information",
                   compute = TRUE)

  } else if(isFALSE(se)) {

    result <- list(type = "information",
                   compute = FALSE)

  } else if(is.character(se) &&
            length(se) == 1L && !is.na(se)) {

    se <- match.arg(tolower(se),
                    c("standard", "information", "robust"))

    if(se == "robust") {
      stop("Robust covariance estimation is not yet implemented in lmvnorm")
    }

    result <- list(type = "information",
                   compute = TRUE)

  } else {

    stop("se must be TRUE, FALSE, 'standard', or 'information'")

  }

  #### Result ####

  return(result)

}

#### Function to create the control list ####

lmvnorm_control <- function(control) {

  if(is.null(control$opt)) control$opt <- "lbfgs"
  if(is.null(control$step_maxit)) control$step_maxit <- 30L
  if(is.null(control$c1)) control$c1 <- 0.5
  if(is.null(control$c2)) control$c2 <- 0.5
  if(is.null(control$step_eps)) control$step_eps <- 1e-09
  if(is.null(control$df_eps)) control$df_eps <- 1e-09
  if(is.null(control$M)) control$M <- 100L
  if(is.null(control$eps)) control$eps <- 1e-06
  if(is.null(control$ss_fac)) control$ss_fac <- 2
  if(is.null(control$maxit)) control$maxit <- 1000L
  if(is.null(control$rstarts)) control$rstarts <- 1L
  if(is.null(control$cores)) control$cores <- 1L
  if(is.null(control$tcg_maxit)) control$tcg_maxit <- 10L
  if(is.null(control$ss)) control$ss <- 0.001
  if(is.null(control$start)) control$start <- NULL

  if(control$step_maxit < 1L) stop("step_maxit must be positive")
  if(control$c1 < 0) stop("c1 must be non-negative")
  if(control$c2 < 0) stop("c2 must be non-negative")
  if(control$step_eps < 0) stop("step_eps must be non-negative")
  if(control$df_eps < 0) stop("df_eps must be non-negative")
  if(control$M < 1L) stop("M must be positive")
  if(control$eps < 0) stop("eps must be non-negative")
  if(control$ss_fac <= 1) stop("ss_fac must be larger than 1")
  if(control$maxit < 1L) stop("maxit must be positive")
  if(control$rstarts < 1L) stop("rstarts must be positive")
  if(control$cores < 1L) stop("cores must be positive")
  if(control$tcg_maxit < 1L) stop("tcg_maxit must be positive")
  if(control$ss <= 0) stop("ss must be positive")

  control$rstarts <- as.integer(control$rstarts)
  control$cores <- as.integer(control$cores)

  #### Result ####

  return(control)

}

#### Function to print progress messages ####

print_lmvnorm_message <- function(msg) {

  w <- nchar(msg)+4L

  cat("\n", "+", strrep("-", w), "+\n",
      "|  ", msg, "  |\n",
      "+", strrep("-", w), "+\n\n", sep = "")

  #### Result ####

  return(invisible(NULL))

}
