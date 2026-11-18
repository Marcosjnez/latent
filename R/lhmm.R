# Hidden Markov extension of the existing LCA model construction.
# A separate class prevents LCA-specific post-estimation methods from being
# applied before their temporal versions are implemented.
setClass("llhmm", contains = "latent")

#' Hidden Markov Models
#'
#' Extend \code{lca()} with a common transition matrix and covariate-dependent
#' initial state probabilities. Indicator and outcome distributions, parameter
#' constraints, starts, penalties, and numerical optimizers use the LCA helpers.
#'
#' @inheritParams lca
#' @param time Name of the time variable. Observations are sorted within each
#'   trajectory. With \code{NULL}, call \code{lca()} with the same arguments.
#' @param id Name of the trajectory variable. With \code{NULL}, all observations
#'   form one trajectory.
#'
#' @details
#' In temporal models, covariates at the first observation of each trajectory
#' predict its initial state probabilities. Later covariate values are unused.
#' Transition probabilities are common across trajectories and time; columns
#' correspond to the previous state and rows to the next state.
#'
#' Weights contain one positive value per trajectory, in order of first
#' appearance in the input data. Outcomes contribute at their observed time
#' points, conditional on the current state. An outcome observed once should
#' appear only once, with missing values at the other times.
#'
#' Each adjacent pair of rows represents one transition. Time orders rows;
#' elapsed-time gaps are not modeled. Repeated response patterns are retained.
#'
#' Temporal fitting currently uses one-step estimation. When omitted,
#' \code{adjustment} is \code{"none"} for temporal data and retains the LCA
#' default otherwise. Temporal EM and multistep adjustments are not implemented.
#' If every trajectory has one observation, fitting delegates to \code{lca()}.
#'
#' Use \code{start$transition} for a positive column-stochastic starting matrix,
#' or \code{start$logtransition} for reference logits. The first row of
#' \code{logtransition} is fixed to zero by default. Parameter constraints use
#' \code{model$logtransition}; other blocks have the same names as in LCA.
#'
#' @return An object of class \code{"llhmm"}, inheriting from \code{"latent"},
#'   or an \code{"llhmmlist"} for several state counts. Non-temporal calls return
#'   the existing LCA objects. Temporal post-estimation methods are separate.
#' @export
lhmm <- function(data,
                 nclasses = 1L,
                 gaussian = NULL,
                 multinomial = NULL,
                 covariates = NULL,
                 outcomes = NULL,
                 model = NULL,
                 weights = NULL,
                 adjustment = "bk",
                 classification = "modal",
                 penalties = TRUE,
                 start = NULL,
                 control = NULL,
                 do.fit = TRUE,
                 verbose = TRUE,
                 time = NULL,
                 id = NULL) {

  #### Temporal structure and the ordinary LCA case ####

  mc <- match.call(expand.dots = TRUE)

  if(!is.null(time)) {
    if(!is.data.frame(data) || nrow(data) == 0L) {
      stop("data must be a non-empty data.frame.")
    }
    if(!is.character(time) || length(time) != 1L || !time %in% names(data) ||
       (!is.null(id) && (!is.character(id) || length(id) != 1L ||
                        !id %in% names(data)))) {
      stop("time and id must name columns of data; id may be NULL.")
    }
    ids <- if(is.null(id)) rep.int(1L, nrow(data)) else data[[id]]
    if(anyNA(ids) || anyNA(data[[time]])) stop("id and time cannot be missing.")
    trajectory <- match(ids, unique(ids))
    lengths <- tabulate(trajectory)
    if(all(lengths == 1L)) time <- NULL
  }

  if(is.null(time)) {
    return(lca(data = data, nclasses = nclasses,
               gaussian = gaussian, multinomial = multinomial,
               covariates = covariates, outcomes = outcomes, model = model,
               weights = weights, adjustment = adjustment,
               classification = classification, penalties = penalties,
               start = start, control = control, do.fit = do.fit,
               verbose = verbose))
  }

  #### Check the temporal inputs ####

  if(missing(adjustment)) adjustment <- "none"
  if(tolower(adjustment) != "none") {
    stop("Temporal models currently require adjustment = 'none'.")
  }
  if(is.null(control)) control <- list()
  if(!is.list(control)) stop("control must be NULL or a list.")
  if(identical(control$opt, "em")) stop("EM is not yet implemented for lhmm.")
  if(!is.numeric(nclasses) || !length(nclasses) ||
     any(!is.finite(nclasses) | nclasses < 1L | nclasses != floor(nclasses))) {
    stop("nclasses must contain positive integers.")
  }
  responses <- c(gaussian, multinomial, unlist(outcomes, use.names = FALSE))
  variables <- c(responses, covariates)
  if(!length(responses)) stop("Supply at least one indicator or outcome.")
  if(any(!variables %in% names(data))) stop("Some model variables are not in data.")
  if(anyDuplicated(variables)) stop("Variables cannot have multiple model roles.")
  if(!is.null(model) && !is.list(model)) {
    stop("For temporal models, model must be a parameter specification list.")
  }
  if(is.null(weights)) weights <- rep(1, length(lengths))
  if(!is.numeric(weights) || length(weights) != length(lengths) ||
     any(!is.finite(weights) | weights <= 0)) {
    stop("weights must contain one positive finite value per trajectory.")
  }

  #### Class enumeration ####

  if(length(nclasses) > 1L) {
    fits <- vector("list", length(nclasses))
    for(i in seq_along(nclasses)) {
      if(verbose) message("Model nclasses=", nclasses[i])
      fits[[i]] <- lhmm(data = data, nclasses = nclasses[i],
                        gaussian = gaussian, multinomial = multinomial,
                        covariates = covariates, outcomes = outcomes,
                        model = model, weights = weights, adjustment = adjustment,
                        classification = classification, penalties = penalties,
                        start = start, control = control, do.fit = do.fit,
                        verbose = verbose, time = time, id = id)
    }
    names(fits) <- paste0("nclasses=", nclasses)
    class(fits) <- "llhmmlist"
    return(fits)
  }

  #### Ordered observations and initial-state covariates ####

  ordering <- order(trajectory, data[[time]])
  data <- data[ordering, , drop = FALSE]
  trajectory <- trajectory[ordering]
  if(anyDuplicated(data.frame(id = trajectory, time = data[[time]]))) {
    stop("Time values must be unique within each trajectory.")
  }
  first <- c(1L, head(cumsum(lengths), -1L) + 1L)
  initial_data <- data[first, , drop = FALSE]
  if(anyNA(initial_data[, covariates, drop = FALSE])) {
    stop("Initial-state covariates cannot be missing at trajectory starts.")
  }

  # Reuse all LCA emission preparation without collapsing rows or deleting
  # interior observations because of covariates that the HMM does not use.
  dataList <- create_lca_dataList(data = data, nclasses = nclasses,
                                 gaussian = gaussian, multinomial = multinomial,
                                 outcomes = outcomes, covariates = NULL,
                                 model = model, collapse = FALSE)
  design <- make_design_matrix(covariates = covariates, data = initial_data)
  rownames(design) <- as.character(ids[ordering][first])
  dataList$design <- dataList$cov_patterns <- design
  dataList$pdesign <- ncol(design)
  dataList$covariates_names <- covariates
  dataList$variables <- unique(c(dataList$variables, covariates))
  dataList$lengths <- lengths
  dataList$trajectory_weights <- weights
  dataList$initial_rows <- first
  dataList$trajectory_ids <- ids[ordering][first]
  dataList$ordering <- ordering
  dataList$inverse_ordering <- order(ordering)
  dataList$time <- time
  dataList$id <- id
  dataList$npossible_patterns <- NA_real_

  #### Control settings and transition starts ####

  if(!is.null(start$transition)) {
    A <- start$transition
    if(!is.null(start$logtransition)) {
      stop("Supply start$transition or start$logtransition, not both.")
    }
    if(!is.matrix(A) || !is.numeric(A) ||
       !all(dim(A) == c(nclasses, nclasses)) ||
       any(!is.finite(A) | A <= 0) || any(abs(colSums(A) - 1) > 1e-8)) {
      stop("start$transition must be a positive column-stochastic matrix.")
    }
    start$logtransition <- sweep(log(A), 2L, log(A[1L, ]), "-")
    start$transition <- NULL
  }
  control$penalties <- penalties
  control$start <- start
  control <- lca_control(control)
  control$estimator <- "lhmm"

  #### Reuse the LCA model and add the transition blocks ####

  full_model <- create_lhmm_model(dataList, nclasses, model, control)
  modelInfo <- create_lca_modelInfo(dataList, full_model, control)

  #### Fit and collect the results, as in lca() ####

  Optim <- list()
  parameters <- transformed_pars <- list()
  elapsed <- numeric()

  if(do.fit) {
    modelInfo$control_optimizer$cores <- min(control$rstarts, control$cores)
    Optim <- optimizer(control_manifold = modelInfo$control_manifold,
                       control_transform = modelInfo$control_transform,
                       control_estimator = modelInfo$control_estimator,
                       control_optimizer = modelInfo$control_optimizer)

    names(Optim$parameters) <- modelInfo$parameters_labels
    names(Optim$transparameters) <- modelInfo$transparameters_labels
    names(Optim$g) <- names(Optim$rg) <- names(Optim$dir) <- modelInfo$parameters_labels
    transformed_pars <- fill_in(modelInfo$trans, Optim$transparameters)
    parameters <- transformed_pars[names(modelInfo$param)]
    elapsed <- Optim$elapsed
  }

  new("llhmm",
       version = as.character(packageVersion("latent")),
       call = mc, timing = elapsed, dataList = dataList, modelInfo = modelInfo,
       Optim = Optim, parameters = parameters,
       transformed_pars = transformed_pars, extra = list())

}

#### HMM additions to the shared LCA parameter construction ####

create_lhmm_model <- function(dataList, nclasses, model, control) {

  # Apply custom starts once, after the transition block has been added.
  lca_settings <- control
  lca_settings$start <- NULL
  result <- create_lca_model(dataList, nclasses, model, lca_settings)
  if(length(result$previous_models)) {
    stop("Multistage fitting is not yet implemented for temporal models.")
  }

  transition_structures <- lapply(c("logtransition", "transition"), function(nm) {
    list(name = nm, type = "matrix", dim = c(nclasses, nclasses),
          rownames = dataList$class_names, colnames = dataList$class_names)
  })
  trans <- create_parameters(transition_structures)
  param <- list(logtransition = trans$logtransition)
  param$logtransition[1L, ] <- "0"
  transition_model <- apply_model_lca(model, list(), param, trans, control)
  result$param <- c(result$param, transition_model$param)
  result$trans <- c(result$trans, transition_model$trans)

  # A modest diagonal preference gives persistent states an initial advantage.
  initial <- diag(2, nclasses)
  initial <- sweep(initial, 2L, initial[1L, ], "-")
  for(i in seq_along(result$init_param)) {
    result$init_param[[i]]$logtransition <- initial
  }
  result$init_param <- custom_init_param(control$start, result$init_param)
  result$control$start <- control$start

  result

}
