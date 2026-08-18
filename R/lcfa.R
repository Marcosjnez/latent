# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 18/08/2026
#'
#' Confirmatory Factor Analysis
#'
#' Fit confirmatory factor analysis models using lavaan model syntax and the
#' optimization infrastructure of \pkg{latent}.
#'
#' @usage
#' lcfa(data = NULL, model = NULL, estimator = "ml",
#'      ordered = FALSE, group = NULL,
#'      sample.cov = NULL, sample.mean = NULL, sample.nobs = NULL,
#'      positive = FALSE, penalties = FALSE,
#'      missing = "pairwise.complete.obs",
#'      std.lv = FALSE, std.ov = FALSE,
#'      acov = "standard", meanstructure = TRUE,
#'      parameterization = NULL,
#'      likelihood = NULL, se = TRUE,
#'      control = NULL, message = FALSE,
#'      do.fit = TRUE, ...)
#'
#' @param data Optional data frame or matrix containing the observed variables.
#'   If NULL, sample.cov and sample.nobs must be supplied.
#' @param model Confirmatory factor model specified using lavaan syntax.
#' @param estimator Estimation method. Available options include \code{"ml"},
#'   \code{"uls"}, and \code{"dwls"}.
#' @param ordered Logical value indicating whether indicators are ordinal. The
#'   character value \code{"yule"} requests Yule correlations.
#' @param group Optional character string identifying the grouping variable.
#' @param sample.cov Optional sample covariance matrix or list of covariance
#'   matrices. Used when data is NULL.
#' @param sample.mean Optional sample mean vector or list of vectors. Required
#'   when data is NULL and meanstructure = TRUE.
#' @param sample.nobs Optional number of observations, or one value per group,
#'   used when data is NULL.
#' @param positive Logical. If \code{TRUE}, positive-definite covariance
#'   structures are imposed through the corresponding manifold parameterization.
#' @param penalties Logical value or list controlling regularization.
#' @param missing Missing-data method.
#' @param std.lv Logical. Standardize latent variables.
#' @param std.ov Logical. Standardize observed variables.
#' @param acov Character string selecting the variance-covariance estimator
#'   for the sample statistics. Available options are \code{"standard"} and
#'   \code{"robust"}.
#' @param meanstructure Logical. Estimate the observed-variable mean structure.
#' @param parameterization Optional parameterization specification.
#' @param likelihood Character string controlling the normal/Wishart likelihood
#'   convention.
#' @param se Logical. Compute standard errors before returning the fitted object.
#' @param control Optional list of optimization controls.
#' @param message Logical. Print progress messages.
#' @param do.fit Logical. If \code{FALSE}, return the prepared but unfitted
#'   \code{"lcfa"} object.
#' @param ... Additional arguments passed to lavaan and the sample-statistic
#'   estimators where applicable.
#'
#' @return An S4 object of class \code{"lcfa"}, which inherits from
#'   \code{"latent"}.
#'
#' @examples
#' \dontrun{
#' HS.model <- '
#'   visual  =~ x1 + x2 + x3
#'   textual =~ x4 + x5 + x6
#'   speed   =~ x7 + x8 + x9
#' '
#'
#' fit <- lcfa(model = HS.model, data = HolzingerSwineford1939)
#' summary(fit, digits = 3L)
#'
#' S <- cov(HolzingerSwineford1939[, paste0("x", 1:9)])
#' M <- colMeans(HolzingerSwineford1939[, paste0("x", 1:9)])
#' fit_cov <- lcfa(model = HS.model, sample.cov = S, sample.mean = M,
#'                 sample.nobs = nrow(HolzingerSwineford1939))
#' }
#'
#' @export
lcfa <- function(data = NULL, model = NULL, estimator = "ml",
                 ordered = FALSE, group = NULL,
                 sample.cov = NULL, sample.mean = NULL, sample.nobs = NULL,
                 positive = FALSE, penalties = FALSE,
                 missing = "pairwise.complete.obs",
                 std.lv = FALSE, std.ov = FALSE,
                 acov = "standard", meanstructure = TRUE,
                 parameterization = NULL,
                 likelihood = NULL, se = TRUE,
                 control = NULL, message = FALSE,
                 do.fit = TRUE,
                 ...) {

  #### Check input arguments ####

  if(is.null(data) && is.null(sample.cov)) {
    stop("Either data or sample.cov must be provided")
  }

  if(!is.null(data) && !is.data.frame(data) && !is.matrix(data)) {
    stop("data must be NULL, a data.frame, or a matrix")
  }

  if(is.matrix(data)) {
    data <- as.data.frame(data)
  }

  if(is.null(model)) {
    stop("model must contain a confirmatory factor model specified with lavaan syntax")
  }

  if(!is.null(control) && !is.list(control)) {
    stop("control must be NULL or a list")
  }

  if(is.null(control)) {
    control <- list()
  }

  estimator <- tolower(estimator)
  missing <- tolower(missing)
  VCOV_type <- match.arg(tolower(acov), c("standard", "robust"))

  if(isTRUE(ordered)) {

    cor <- "poly"
    std.ov <- TRUE

    if(positive) {
      control$deltaparam <- FALSE
      std.lv <- TRUE
    } else {
      control$deltaparam <- TRUE
    }

  } else if(is.character(ordered) &&
            length(ordered) == 1L &&
            tolower(ordered) == "yule") {

    cor <- "yule"
    std.ov <- TRUE
    control$deltaparam <- TRUE

  } else if(isFALSE(ordered)) {

    cor <- "pearson"

  } else {

    stop("ordered must be FALSE, TRUE, or 'yule'")

  }

  # Preserve the likelihood convention currently used by lcfa.
  if(is.null(likelihood)) {
    likelihood <- "normal"
  }

  if(missing == "fiml") {

    if(is.null(data)) {
      stop("missing = 'fiml' requires raw data")
    }

    meanstructure <- TRUE
    std.ov <- FALSE

  }

  if(is.null(data)) {

    if(meanstructure && is.null(sample.mean)) {
      stop("sample.mean= argument is missing, but model contains mean/intercept parameters.")
    }

    if(is.null(sample.nobs)) {
      stop("sample.nobs must be supplied when data is NULL")
    }

    if(!isFALSE(ordered)) {
      stop("ordered models require raw data")
    }

    if(VCOV_type == "robust") {
      stop("acov = 'robust' requires raw data")
    }

  }

  if(meanstructure) {
    if(estimator %in% c("ml", "fml")) estimator <- "means_fml"
    if(estimator == "uls") estimator <- "means_uls"
    if(estimator == "dwls") estimator <- "means_dwls"
  }

  #### Store original call ####

  mc <- match.call(expand.dots = TRUE)
  args <- lapply(as.list(mc)[-1L], eval, envir = parent.frame())

  #### Check control parameters ####

  control$ordered <- ordered
  control$std.lv <- std.lv
  control$std.ov <- std.ov
  control$positive <- positive
  control$penalties <- penalties
  control$estimator <- estimator
  control$meanstructure <- meanstructure
  control$missing <- missing
  control$VCOV <- VCOV_type
  control <- lcfa_control(control)

  #### Create the dataList ####

  dataList <- create_lcfa_dataList(data = data,
                                   model = model,
                                   cor = cor,
                                   estimator = estimator,
                                   ordered = ordered,
                                   group = group,
                                   sample.cov = sample.cov,
                                   sample.mean = sample.mean,
                                   sample.nobs = sample.nobs,
                                   positive = positive,
                                   penalties = penalties,
                                   missing = missing,
                                   std.lv = std.lv,
                                   std.ov = std.ov,
                                   VCOV = VCOV_type,
                                   message = message,
                                   likelihood = likelihood,
                                   meanstructure = meanstructure,
                                   args = args,
                                   control = control,
                                   ...)

  #### Create the model ####

  full_model <- create_lcfa_model(dataList = dataList,
                                  model = model,
                                  control = control)
  list2env(full_model, envir = environment())
  dataList$data_param <- data_param

  #### Create the manifold, transformation, and estimator structures ####

  modelInfo <- create_lcfa_modelInfo(dataList = dataList,
                                     full_model = full_model,
                                     control = control)

  #### Fit the model ####

  if(!do.fit) {

    previous_models <- c(dataList$fit_means,
                         dataList$fit_cov)
    previous_models <- previous_models[
      vapply(previous_models,
             FUN = \(x) inherits(x, "latent"),
             FUN.VALUE = logical(1L))
    ]

    result <- new("lcfa",
                  version          = as.character(packageVersion("latent")),
                  call             = mc,
                  timing           = numeric(),
                  dataList         = dataList,
                  modelInfo        = modelInfo,
                  Optim            = list(),
                  parameters       = list(),
                  transformed_pars = list(),
                  extra            = previous_models)

    #### Result ####

    return(result)

  }

  if(message) {
    print_lcfa_message("Fitting the model")
  }

  modelInfo$control_optimizer$cores <-
    min(modelInfo$control_optimizer$rstarts,
        modelInfo$control_optimizer$cores)

  Optim <- optimizer(control_manifold = modelInfo$control_manifold,
                     control_transform = modelInfo$control_transform,
                     control_estimator = modelInfo$control_estimator,
                     control_optimizer = modelInfo$control_optimizer)

  names(Optim$parameters) <- modelInfo$parameters_labels
  names(Optim$transparameters) <- modelInfo$transparameters_labels

  if(!is.null(Optim$g)) names(Optim$g) <- modelInfo$parameters_labels
  if(!is.null(Optim$rg)) names(Optim$rg) <- modelInfo$parameters_labels
  if(!is.null(Optim$dir)) names(Optim$dir) <- modelInfo$parameters_labels

  #### Process the outputs ####

  transformed_pars <- fill_in(modelInfo$trans,
                              Optim$transparameters)

  parameters <- transformed_pars[names(modelInfo$param)]

  #### Previous sample-statistic models ####

  previous_models <- c(dataList$fit_means,
                       dataList$fit_cov)
  previous_models <- previous_models[
    vapply(previous_models,
           FUN = \(x) inherits(x, "latent"),
           FUN.VALUE = logical(1L))
  ]

  #### latent object ####

  result <- new("lcfa",
                version          = as.character(packageVersion("latent")),
                call             = mc,
                timing           = Optim$elapsed,
                dataList         = dataList,
                modelInfo        = modelInfo,
                Optim            = Optim,
                parameters       = parameters,
                transformed_pars = transformed_pars,
                extra            = previous_models)

  #### Standard errors ####

  if(isTRUE(se)) {

    if(message) {
      print_lcfa_message("Computing standard errors")
    }

    # parameters_se <- modelInfo$trans[names(modelInfo$param)]
    parameters <- modelInfo$parameters_labels
    type_se <- if(VCOV_type == "robust") "robust" else "information"

    result@Optim$SE <- se(result,
                          type = type_se,
                          parameters = parameters)

  }

  #### Result ####

  return(result)

}

#### Function to create the dataList ####

create_lcfa_dataList <- function(data = NULL, model = NULL, cor = "pearson",
                                 estimator = "ml", ordered = FALSE,
                                 group = NULL, sample.cov = NULL,
                                 sample.mean = NULL, sample.nobs = NULL,
                                 positive = FALSE, penalties = TRUE,
                                 missing = "pairwise.complete.obs",
                                 std.lv = TRUE, std.ov = FALSE,
                                 VCOV = "standard", message = FALSE,
                                 likelihood = NULL, meanstructure = TRUE,
                                 args = NULL, control = NULL,
                                 ...) {

  cor <- tolower(cor)
  estimator <- tolower(estimator)
  VCOV_type <- match.arg(tolower(VCOV), c("standard", "robust"))
  missing <- tolower(missing)

  sample_stats_only <- is.null(data)

  #### Groups ####

  if(sample_stats_only) {

    if(!is.null(group)) {
      stop("group cannot be used when data is NULL; provide sample.cov as a list for multiple groups")
    }

    if(is.data.frame(sample.cov)) {
      sample.cov <- as.matrix(sample.cov)
    }

    if(!is.list(sample.cov)) {
      sample.cov <- list(sample.cov)
    }

    ngroups <- length(sample.cov)

    if(ngroups < 1L) {
      stop("sample.cov must contain at least one covariance matrix")
    }

    group_label <- names(sample.cov)

    if(is.null(group_label)) {
      group_label <- rep("", ngroups)
    }

    if(ngroups == 1L) {

      group_label <- ""

    } else {

      empty <- is.na(group_label) | group_label == ""
      group_label[empty] <- as.character(which(empty))

      if(anyDuplicated(group_label)) {
        stop("sample.cov must have unique group names when supplied as a named list")
      }

      names(sample.cov) <- group_label

    }

  } else if(is.null(group)) {

    ngroups <- 1L
    group <- "group"
    group_label <- ""
    data$group <- group_label

  } else {

    if(!group %in% colnames(data)) {
      stop("The grouping variable is not present in data")
    }

    group_label <- unique(data[[group]])
    ngroups <- length(group_label)

  }

  item_names <- extract_item_names_lavaan(model, ngroups = ngroups)

  for(i in seq_len(ngroups)) {
    if(length(item_names[[i]]) == 0L) {
      stop("No observed variables were identified for group ", i)
    }
  }

  #### Raw data ####

  if(!sample_stats_only) {

    keep <- rep(FALSE, nrow(data))

    for(i in seq_len(ngroups)) {

      group_i <- data[[group]] == group_label[i]
      items_i <- item_names[[i]]

      not_all_na_i <- !apply(is.na(data[group_i, items_i, drop = FALSE]),
                             MARGIN = 1L, FUN = all)
      keep[group_i] <- not_all_na_i

    }

    data <- data[keep, , drop = FALSE]

    if(nrow(data) == 0L) {
      stop("No observations remain after removing cases with all model variables missing")
    }

  }

  #### Sample statistics ####

  X <- vector("list", length = ngroups)
  nobs_list <- vector("list", length = ngroups)
  sample.mean_list <- vector("list", length = ngroups)
  sample.cov_input <- vector("list", length = ngroups)
  NVCOV <- vector("list", length = ngroups)
  VCOV_cov <- vector("list", length = ngroups)
  NVCOV_means <- vector("list", length = ngroups)
  VCOV_means <- vector("list", length = ngroups)
  WLS.V <- vector("list", length = ngroups)
  thresholds <- vector("list", length = ngroups)
  fit_cov <- vector("list", length = ngroups)
  fit_means <- vector("list", length = ngroups)

  if(sample_stats_only) {

    nobs_list <- normalize_lcfa_sample_nobs(sample.nobs = sample.nobs,
                                            ngroups = ngroups,
                                            group_label = group_label)

    sample.mean_list <- normalize_lcfa_sample_mean(sample.mean = sample.mean,
                                                   item_names = item_names,
                                                   ngroups = ngroups,
                                                   group_label = group_label,
                                                   meanstructure = meanstructure,
                                                   std.ov = std.ov)

    for(i in seq_len(ngroups)) {

      S_input <- validate_lcfa_sample_covariance(sample.cov[[i]],
                                                 item_names = item_names[[i]])

      sample.cov_input[[i]] <- S_input

      S <- prepare_lcfa_sample_covariance(S = S_input,
                                          sample.nobs = nobs_list[[i]],
                                          std.ov = std.ov,
                                          likelihood = likelihood)

      sample.cov[[i]] <- S

      # asymptotic_normal() returns the N-scaled covariance matrix.
      NVCOV[[i]] <- asymptotic_normal(S,
                                      cov = !std.ov,
                                      diag = FALSE)
      VCOV_cov[[i]] <- NVCOV[[i]]/nobs_list[[i]]

      # lavaan and the DWLS discrepancy use the N-scaled covariance matrix.
      WLS.V[[i]] <- diag(1/diag(NVCOV[[i]]))

      if(std.ov) {

        NVCOV_means[[i]] <- matrix(0,
                                    nrow = nrow(S_input),
                                    ncol = ncol(S_input),
                                    dimnames = dimnames(S_input))

      } else {

        # Under normal sampling, Var(sqrt(N)*sample.mean) is the observed
        # covariance matrix. Retain its off-diagonal elements.
        NVCOV_means[[i]] <- S_input

      }

      VCOV_means[[i]] <- NVCOV_means[[i]]/nobs_list[[i]]

      thresholds[[i]] <- list()

    }

  } else {

    sample.cov <- normalize_lcfa_group_input(sample.cov, ngroups)

    for(i in seq_len(ngroups)) {

      X[[i]] <- extract_lcfa_group_data(data = data,
                                        group = group,
                                        group_label = group_label,
                                        item_names = item_names,
                                        i = i,
                                        ngroups = ngroups)

      control_i <- control

      if(ngroups < 2L) {
        control_i$subfix <- ""
      } else {
        control_i$subfix <- paste0(".", group_label[i])
      }

      sample_stats <- estimate_lcfa_sample_statistics(data = X[[i]],
                                                      model = model,
                                                      cor = cor,
                                                      std.ov = std.ov,
                                                      VCOV = VCOV_type,
                                                      likelihood = likelihood,
                                                      missing = missing,
                                                      control = control_i,
                                                      ...)

      fit_means[[i]] <- sample_stats$fit_means
      fit_cov[[i]] <- sample_stats$fit_cov

      nobs_list[[i]] <- fit_cov[[i]]@dataList$nobs
      sample.cov[[i]] <- fit_cov[[i]]@transformed_pars$S
      sample.cov_input[[i]] <- sample.cov[[i]]
      VCOV_cov[[i]] <- fit_cov[[i]]@Optim$SE$VCOV
      NVCOV[[i]] <- VCOV_cov[[i]]*fit_cov[[i]]@dataList$nobs

      VCOV_means[[i]] <- fit_means[[i]]@Optim$SE$VCOV
      NVCOV_means[[i]] <-
        VCOV_means[[i]]*fit_means[[i]]@dataList$nobs

      WLS.V[[i]] <- diag(1/diag(NVCOV[[i]]))

      idx_taus <- startsWith(names(fit_cov[[i]]@transformed_pars), "taus")
      thresholds[[i]] <- fit_cov[[i]]@transformed_pars[idx_taus]

      idx_means <- startsWith(names(fit_means[[i]]@transformed_pars), "means")

      if(any(idx_means)) {

        M <- fit_means[[i]]@transformed_pars[[which(idx_means)[1L]]]
        sample.mean_list[[i]] <- c(M)
        names(sample.mean_list[[i]]) <- rownames(M)

      } else {

        sample.mean_list[[i]] <- setNames(rep(0, length(item_names[[i]])),
                                          item_names[[i]])

      }

    }

  }

  #### Lavaan model structure ####

  sample.cov_lav <- if(sample_stats_only) sample.cov_input else sample.cov

  if(ngroups > 1L && sample_stats_only) {
    names(sample.cov_lav) <- group_label
    names(sample.mean_list) <- group_label
    names(nobs_list) <- group_label
  }

  LAV <- lavaan::cfa(model = model,
                     sample.cov = unwrap_lcfa_group_input(sample.cov_lav, ngroups),
                     sample.mean = if(meanstructure) {
                       unwrap_lcfa_group_input(sample.mean_list, ngroups)
                     } else {
                       NULL
                     },
                     sample.nobs = unwrap_lcfa_group_input(nobs_list, ngroups),
                     group = group,
                     NACOV = if(sample_stats_only) {
                       NULL
                     } else {
                       unwrap_lcfa_group_input(NVCOV, ngroups)
                     },
                     WLS.V = if(sample_stats_only) {
                       NULL
                     } else {
                       unwrap_lcfa_group_input(WLS.V, ngroups)
                     },
                     ordered = ordered,
                     std.lv = std.lv,
                     std.ov = std.ov,
                     meanstructure = meanstructure,
                     do.fit = FALSE,
                     warn = FALSE,
                     ...)

  LAV@Options$positive <- positive

  lavmodel <- LAV@Model
  item_label <- LAV@Data@ov.names
  nobs_list <- LAV@Data@nobs
  factor_label <- replicate(ngroups,
                            list(LAV@Model@dimNames[[1L]][[2L]]))

  if(sample_stats_only && ngroups > 1L) {

    lav_group_label <- tryCatch(lavaan::lavInspect(LAV, "group.label"),
                                error = function(e) NULL)

    if(length(lav_group_label) == ngroups) {
      group_label <- lav_group_label
    }

  }

  model_out <- getmodel_fromlavaan(LAV)

  if(ngroups == 1L) {
    model_out <- list(model_out)
  }

  nitems <- as.list(lavmodel@nvar)

  if(meanstructure) {
    npatterns <- lapply(nitems, FUN = \(p) 0.5*p*(p+1)+p)
  } else {
    npatterns <- lapply(nitems, FUN = \(p) 0.5*p*(p+1))
  }

  nfactors <- lapply(model_out, FUN = \(x) ncol(x$lambda))

  if(is.null(args)) {
    args <- as.list(match.call(expand.dots = TRUE))[-1L]
  }

  #### Store objects in dataList ####

  dataList <- list(ngroups = ngroups,
                   data = data,
                   data_per_group = X,
                   sample_stats_only = sample_stats_only,
                   nobs = nobs_list,
                   sample.nobs = nobs_list,
                   nitems = nitems,
                   npatterns = npatterns,
                   nfactors = nfactors,
                   positive = positive,
                   estimator = estimator,
                   cor = cor,
                   VCOV_type = VCOV_type,
                   group = group,
                   group_label = group_label,
                   item_label = item_label,
                   factor_label = factor_label,
                   LAV = LAV,
                   args = args,
                   model = model_out,
                   sample.cov = sample.cov,
                   sample.cov.input = sample.cov_input,
                   sample.mean = sample.mean_list,
                   NVCOV = NVCOV,
                   VCOV = VCOV_cov,
                   NVCOV_means = NVCOV_means,
                   VCOV_means = VCOV_means,
                   WLS.V = WLS.V,
                   thresholds = thresholds,
                   fit_means = fit_means,
                   fit_cov = fit_cov)

  #### Result ####

  return(dataList)

}

#### Function to create the model ####

create_lcfa_model <- function(dataList, model = NULL, control) {

  #### Parameters from the sample statistics ####

  data_param <- create_lcfa_data_param(dataList = dataList,
                                       control = control)

  #### Model for the transformed parameters ####

  trans <- model_lcfa(dataList = dataList,
                      data_param = data_param,
                      control = control)

  #### Model for the parameters ####

  constraints <- constraints_lcfa(dataList = dataList,
                                  data_param = data_param,
                                  trans = trans,
                                  control = control)

  param <- constraints$param
  trans <- constraints$trans

  #### Create the initial values for the parameters ####

  init_param <- start_lcfa(dataList = dataList,
                           data_param = data_param,
                           param = param,
                           trans = trans,
                           fixed = constraints$fixed,
                           fixed_values_list = constraints$fixed_values_list,
                           control = control)

  #### Custom initial values ####

  init_param <- custom_init_param(control$start, init_param)

  #### Result ####

  result <- list(param = param,
                 trans = trans,
                 init_param = init_param,
                 target_psi = constraints$target_psi,
                 target_theta = constraints$target_theta,
                 data_param = data_param,
                 control = control)

  return(result)

}

#### Function to create the modelInfo ####

create_lcfa_modelInfo <- function(dataList, full_model, control) {

  list2env(full_model, envir = environment())

  #### Manifolds ####

  control_manifold <- manifolds_lcfa(dataList = dataList,
                                     data_param = data_param,
                                     param = param,
                                     target_psi = target_psi,
                                     target_theta = target_theta,
                                     control = control)

  #### Transformations ####

  control_transform <- transformations_lcfa(dataList = dataList,
                                             data_param = data_param,
                                             trans = trans,
                                             control = control)

  #### Estimators ####

  control_estimator <- estimators_lcfa(dataList = dataList,
                                       data_param = data_param,
                                       trans = trans,
                                       control = control)

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

  #### Collect all the model information ####

  modelInfo <- list(param = param,
                    trans = trans,
                    nparam = nparam,
                    ntrans = ntrans,
                    parameters_labels = parameters_labels,
                    transparameters_labels = transparameters_labels,
                    dof = sum(unlist(dataList$npatterns))-nparam,
                    control_manifold = control_manifold,
                    control_transform = control_transform,
                    control_estimator = control_estimator,
                    control_optimizer = control_optimizer)

  #### Result ####

  return(modelInfo)

}

#### Function to create the control list of optimization parameters ####

lcfa_control <- function(control) {

  # Auxiliary function for lcfa.R

  #### Penalties ####

  penalty_defaults <- list(
    logdet = list(w = 1e-03)
  )

  if(!control$positive) {
    control$penalties <- FALSE
  }

  if(isFALSE(control$penalties)) {

    control$reg <- FALSE

  } else if(isTRUE(control$penalties)) {

    control$reg <- TRUE
    control$penalties <- penalty_defaults

  } else if(is.list(control$penalties)) {

    unknown_penalties <- setdiff(names(control$penalties),
                                 names(penalty_defaults))

    if(length(unknown_penalties) > 0L) {
      stop("Unknown penalty name(s): ",
           paste(unknown_penalties, collapse = ", "))
    }

    control$penalties <- utils::modifyList(penalty_defaults,
                                           control$penalties)

    if(!is.numeric(control$penalties$logdet$w) ||
       length(control$penalties$logdet$w) != 1L ||
       !is.finite(control$penalties$logdet$w) ||
       control$penalties$logdet$w <= 0) {
      stop("The logdet penalty w must be a positive number")
    }

    control$reg <- TRUE

  } else {

    stop("penalties should be TRUE, FALSE, or a named list")

  }

  #### Model controls ####

  if(is.null(control$free_S)) control$free_S <- FALSE
  if(is.null(control$free_taus)) control$free_taus <- FALSE
  if(is.null(control$free_M)) control$free_M <- FALSE
  if(is.null(control$deltaparam)) control$deltaparam <- FALSE
  if(is.null(control$start)) control$start <- NULL

  #### Optimizer ####

  if(is.null(control$opt)) {

    if(control$positive) {
      control$opt <- "grad"
    } else {
      control$opt <- "lbfgs"
    }

  }

  if(is.null(control$rstarts)) {

    if(control$positive) {
      control$rstarts <- 10L
    } else {
      control$rstarts <- 1L
    }

  } else if(control$rstarts < 1L ||
            !all(control$rstarts == as.integer(control$rstarts))) {
    stop("rstarts must be a positive integer")
  }

  if(is.null(control$step_maxit)) {
    control$step_maxit <- 30L
  } else if(control$step_maxit < 1L) {
    stop("step_maxit must be a positive integer")
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
  } else if(control$M < 0L) {
    stop("M must be a positive integer")
  }

  if(is.null(control$eps)) {
    control$eps <- 1e-06
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
  } else if(control$maxit < 0L) {
    stop("maxit must be a positive integer")
  }

  if(is.null(control$cores)) {
    control$cores <- 1L
  } else if(control$cores < 1L) {
    stop("cores must be a positive integer")
  }

  if(is.null(control$tcg_maxit)) {
    control$tcg_maxit <- 10L
  } else if(control$tcg_maxit < 1L) {
    stop("tcg_maxit must be a positive integer")
  }

  #### Result ####

  return(control)

}

#### Auxiliary functions for create_lcfa_dataList ####

print_lcfa_message <- function(msg) {

  w <- nchar(msg)+4L
  cat("\n", "+", strrep("-", w), "+\n",
      "|  ", msg, "  |\n",
      "+", strrep("-", w), "+\n\n", sep = "")

  #### Result ####

  return(invisible(NULL))

}

normalize_lcfa_group_input <- function(x, ngroups) {

  if(is.null(x)) {

    result <- vector("list", ngroups)

  } else if(ngroups == 1L && !is.list(x)) {

    result <- list(x)

  } else {

    result <- x

  }

  #### Result ####

  return(result)

}

unwrap_lcfa_group_input <- function(x, ngroups) {

  if(length(x) == 0L || all(vapply(x, is.null, logical(1L)))) {

    result <- NULL

  } else if(ngroups == 1L) {

    result <- x[[1L]]

  } else {

    result <- x

  }

  #### Result ####

  return(result)

}

normalize_lcfa_sample_nobs <- function(sample.nobs, ngroups,
                                        group_label = NULL) {

  input_names <- names(sample.nobs)

  if(is.list(sample.nobs)) {
    sample.nobs <- unlist(sample.nobs, use.names = TRUE)
  }

  if(!is.null(input_names) && is.null(names(sample.nobs))) {
    names(sample.nobs) <- input_names
  }

  if(ngroups > 1L &&
     !is.null(group_label) &&
     !is.null(names(sample.nobs))) {

    if(!setequal(names(sample.nobs), group_label)) {
      stop("Names of sample.nobs must match the groups in sample.cov")
    }

    sample.nobs <- sample.nobs[group_label]

  }

  if(!is.numeric(sample.nobs) ||
     length(sample.nobs) != ngroups ||
     any(!is.finite(sample.nobs)) ||
     any(sample.nobs < 2L) ||
     any(sample.nobs != as.integer(sample.nobs))) {
    stop("sample.nobs must contain one integer greater than 1 for each group")
  }

  result <- as.list(as.integer(sample.nobs))

  if(ngroups > 1L && !is.null(group_label)) {
    names(result) <- group_label
  }

  #### Result ####

  return(result)

}

normalize_lcfa_sample_mean <- function(sample.mean, item_names, ngroups,
                                       group_label = NULL,
                                       meanstructure, std.ov) {

  if(is.null(sample.mean)) {

    if(meanstructure) {
      stop("sample.mean= argument is missing, but model contains mean/intercept parameters.")
    }

    sample.mean <- lapply(item_names, FUN = \(x) {
      setNames(rep(0, length(x)), x)
    })

  } else if(ngroups == 1L && !is.list(sample.mean)) {

    sample.mean <- list(sample.mean)

  } else if(!is.list(sample.mean) || length(sample.mean) != ngroups) {

    stop("sample.mean must be a vector or a list with one vector per group")

  }

  if(ngroups > 1L &&
     !is.null(group_label) &&
     !is.null(names(sample.mean))) {

    if(!setequal(names(sample.mean), group_label)) {
      stop("Names of sample.mean must match the groups in sample.cov")
    }

    sample.mean <- sample.mean[group_label]

  }

  result <- vector("list", ngroups)

  for(i in seq_len(ngroups)) {

    x <- sample.mean[[i]]

    if(is.matrix(x)) {

      if(!any(dim(x) == 1L)) {
        stop("Each sample.mean element must be a vector")
      }

      if(ncol(x) == 1L) {
        x_names <- rownames(x)
      } else {
        x_names <- colnames(x)
      }

      x <- c(x)

      if(!is.null(x_names)) {
        names(x) <- x_names
      }

    }

    if(!is.numeric(x) || length(x) != length(item_names[[i]]) ||
       any(!is.finite(x))) {
      stop("Each sample.mean element must contain one finite numeric value for each observed variable")
    }

    if(is.null(names(x))) {

      names(x) <- item_names[[i]]

    } else {

      if(any(names(x) == "") || anyDuplicated(names(x))) {
        stop("sample.mean must have unique, non-empty names")
      }

      unknown <- setdiff(names(x), item_names[[i]])
      missing <- setdiff(item_names[[i]], names(x))

      if(length(unknown) > 0L || length(missing) > 0L) {
        stop("Names of sample.mean must match the observed variables in the model")
      }

      x <- x[item_names[[i]]]

    }

    if(std.ov) {
      x[] <- 0
    }

    result[[i]] <- x

  }

  if(ngroups > 1L && !is.null(group_label)) {
    names(result) <- group_label
  }

  #### Result ####

  return(result)

}

validate_lcfa_sample_covariance <- function(S, item_names) {

  if(is.data.frame(S)) {
    S <- as.matrix(S)
  }

  if(!is.matrix(S) || !is.numeric(S) || nrow(S) != ncol(S)) {
    stop("sample.cov must be a numeric square matrix or a list of numeric square matrices")
  }

  if(any(!is.finite(S))) {
    stop("sample.cov cannot contain missing or non-finite values")
  }

  rn <- rownames(S)
  cn <- colnames(S)

  if(is.null(rn) || is.null(cn) ||
     any(rn == "") || any(cn == "") ||
     anyDuplicated(rn) || anyDuplicated(cn)) {
    stop("sample.cov must have unique, non-empty rownames and colnames")
  }

  unknown_rows <- setdiff(rn, item_names)
  unknown_cols <- setdiff(cn, item_names)
  missing_rows <- setdiff(item_names, rn)
  missing_cols <- setdiff(item_names, cn)

  if(length(unknown_rows) > 0L || length(unknown_cols) > 0L ||
     length(missing_rows) > 0L || length(missing_cols) > 0L) {
    stop("rownames and colnames of sample.cov must match the observed variables in the model")
  }

  S <- S[item_names, item_names, drop = FALSE]

  if(!isSymmetric(S, tol = sqrt(.Machine$double.eps))) {
    stop("sample.cov must be symmetric")
  }

  #### Result ####

  return(S)

}

prepare_lcfa_sample_covariance <- function(S, sample.nobs,
                                           std.ov, likelihood) {

  if(std.ov) {

    if(any(diag(S) <= 0)) {
      stop("Observed variables cannot be standardized because at least one variance is non-positive")
    }

    item_names <- rownames(S)
    inv_sqrtdiagS <- diag(1/sqrt(diag(S)), nrow = nrow(S))
    S <- inv_sqrtdiagS %*% S %*% inv_sqrtdiagS
    rownames(S) <- colnames(S) <- item_names

  }

  if(likelihood == "normal") {
    S <- S*(sample.nobs-1L)/sample.nobs
  }

  #### Result ####

  return(S)

}

extract_lcfa_group_data <- function(data, group, group_label,
                                    item_names, i, ngroups) {

  if(ngroups > 1L) {
    result <- data[data[[group]] == group_label[i],
                   item_names[[i]], drop = FALSE]
  } else {
    result <- data[, item_names[[i]], drop = FALSE]
  }

  #### Result ####

  return(result)

}

estimate_lcfa_sample_statistics <- function(data, model, cor,
                                            std.ov, VCOV, likelihood,
                                            missing, control, ...) {

  #### Means ####

  fit_means <- lmean(data = data,
                     std.ov = std.ov,
                     control = control,
                     do.fit = TRUE)

  #### Covariances ####

  if(cor == "pearson") {

    fit_cov <- lpearson(data = data,
                        std.ov = std.ov,
                        VCOV = VCOV,
                        likelihood = likelihood,
                        missing = missing,
                        control = control,
                        do.fit = TRUE)

    if(missing == "fiml") {

      patterns <- split_by_missing_pattern(data)
      npatterns <- length(patterns)

      fit_means@extra <- vector("list", length = npatterns)
      fit_cov@extra <- vector("list", length = npatterns)

      subfix <- control$subfix

      for(j in seq_len(npatterns)) {

        control_j <- control
        control_j$subfix <- paste0(subfix, ".pattern", j)

        fit_means@extra[[j]] <- lmean(data = patterns[[j]]$data,
                                      std.ov = std.ov,
                                      do.fit = TRUE,
                                      control = control_j,
                                      ...)

        fit_cov@extra[[j]] <- lpearson(data = patterns[[j]]$data,
                                       model = model,
                                       std.ov = std.ov,
                                       VCOV = VCOV,
                                       likelihood = likelihood,
                                       missing = "pairwise.complete.obs",
                                       do.fit = TRUE,
                                       control = control_j,
                                       ...)

      }

    }

  } else if(cor == "poly") {

    fit_cov <- lpoly(data = data,
                     method = "two-step",
                     control = control,
                     do.fit = TRUE)

  } else if(cor == "yule") {

    fit_cov <- lyule(data = data,
                     control = control,
                     do.fit = TRUE)

  } else {

    stop("Unknown correlation type")

  }

  #### Result ####

  result <- list(fit_means = fit_means,
                 fit_cov = fit_cov)

  return(result)

}

split_by_missing_pattern <- function(data) {

  if(!is.data.frame(data) && !is.matrix(data)) {
    stop("data must be a data.frame or matrix.")
  }

  if(nrow(data) == 0L) {

    #### Result ####

    return(list())

  }

  miss <- is.na(data)
  pattern_key <- apply(miss, MARGIN = 1L,
                       FUN = \(x) paste(as.integer(x), collapse = ""))
  pattern_levels <- unique(pattern_key)
  counts <- tabulate(match(pattern_key, pattern_levels),
                     nbins = length(pattern_levels))
  ord <- order(-counts, seq_along(pattern_levels))
  pattern_levels <- pattern_levels[ord]

  result <- vector("list", length(pattern_levels))

  for(k in seq_along(pattern_levels)) {

    idx_rows <- which(pattern_key == pattern_levels[k])
    idx_vars <- !miss[idx_rows[1L], ]

    result[[k]] <- list(data = data[idx_rows, idx_vars, drop = FALSE],
                        vars = idx_vars,
                        nobs = length(idx_rows))

  }

  names(result) <- NULL

  #### Result ####

  return(result)

}

#### Auxiliary functions for create_lcfa_model ####

create_lcfa_data_param <- function(dataList, control) {

  ngroups <- dataList$ngroups

  #### Parameter-block names ####

  if(ngroups < 2L) {
    sep <- ""
  } else {
    sep <- "."
  }

  lambda_group <- paste("lambda", dataList$group_label, sep = sep)
  psi_group <- paste("psi", dataList$group_label, sep = sep)
  theta_group <- paste("theta", dataList$group_label, sep = sep)
  xpsi_group <- paste("xpsi", dataList$group_label, sep = sep)
  xtheta_group <- paste("xtheta", dataList$group_label, sep = sep)
  model_group <- paste("model", dataList$group_label, sep = sep)
  nu_group <- paste("nu", dataList$group_label, sep = sep)
  delta_group <- paste("delta", dataList$group_label, sep = sep)
  tau_group <- paste("tau", dataList$group_label, sep = sep)

  #### Sample-statistic parameters ####

  S_group <- vector("list", length = ngroups)
  taus_group <- vector("list", length = ngroups)
  M_group <- vector("list", length = ngroups)
  means_params <- vector("list", length = ngroups)
  means_params_labels <- vector("list", length = ngroups)
  cov_params <- vector("list", length = ngroups)
  cov_params_labels <- vector("list", length = ngroups)
  VCOV_means <- vector("list", length = ngroups)
  VCOV_cov <- vector("list", length = ngroups)
  NVCOV_means <- vector("list", length = ngroups)
  NVCOV_cov <- vector("list", length = ngroups)
  nobs_ij <- vector("list", length = ngroups)

  if(isTRUE(dataList$sample_stats_only)) {

    for(i in seq_len(ngroups)) {

      if(ngroups < 2L) {
        subfix <- ""
      } else {
        subfix <- paste0(".", dataList$group_label[i])
      }

      items <- dataList$item_label[[i]]
      p <- length(items)

      #### Means ####

      means_name <- paste0("means", subfix)

      means_struct <- create_parameters(list(
        list(name = means_name,
             type = "matrix",
             dim = c(p, 1L),
             rownames = items,
             colnames = "intrcpt")
      ))

      means_params_labels[[i]] <- means_struct
      means_params[[i]] <- list()
      means_params[[i]][[means_name]] <-
        matrix(dataList$sample.mean[[i]][items],
               nrow = p,
               ncol = 1L,
               dimnames = list(items, "intrcpt"))

      mean_labels <- c(means_struct[[means_name]])
      VCOV_means[[i]] <- list(dataList$VCOV_means[[i]])
      NVCOV_means[[i]] <- list(dataList$NVCOV_means[[i]])

      rownames(VCOV_means[[i]][[1L]]) <-
        colnames(VCOV_means[[i]][[1L]]) <- mean_labels
      rownames(NVCOV_means[[i]][[1L]]) <-
        colnames(NVCOV_means[[i]][[1L]]) <- mean_labels

      #### Covariance matrix ####

      S_name <- paste0("S", subfix)

      cov_struct <- create_parameters(list(
        list(name = S_name,
             type = "matrix",
             dim = c(p, p),
             rownames = items,
             colnames = items,
             symmetric = TRUE)
      ))

      cov_params_labels[[i]] <- cov_struct
      cov_params[[i]] <- list()
      cov_params[[i]][[S_name]] <- dataList$sample.cov[[i]]
      dimnames(cov_params[[i]][[S_name]]) <- dimnames(cov_struct[[S_name]])

      S_labels <- c(cov_struct[[S_name]][lower.tri(cov_struct[[S_name]],
                                                   diag = !control$std.ov)])

      VCOV_cov[[i]] <- list(dataList$VCOV[[i]])
      NVCOV_cov[[i]] <- list(dataList$NVCOV[[i]])

      rownames(VCOV_cov[[i]][[1L]]) <-
        colnames(VCOV_cov[[i]][[1L]]) <- S_labels
      rownames(NVCOV_cov[[i]][[1L]]) <-
        colnames(NVCOV_cov[[i]][[1L]]) <- S_labels

      nobs_ij[[i]] <- dataList$nobs[[i]]

      M_group[[i]] <- means_name
      S_group[[i]] <- S_name
      taus_group[[i]] <- character(0L)

    }

  } else {

    for(i in seq_len(ngroups)) {

      if(control$missing == "fiml") {

        means_params[[i]] <- unlist(lapply(dataList$fit_means[[i]]@extra,
                                           FUN = \(x) x@parameters),
                                    recursive = FALSE)
        means_params_labels[[i]] <-
          unlist(lapply(dataList$fit_means[[i]]@extra,
                        FUN = \(x) x@modelInfo$trans),
                 recursive = FALSE)
        VCOV_means[[i]] <- lapply(dataList$fit_means[[i]]@extra,
                                  FUN = \(x) x@Optim$SE$VCOV)
        NVCOV_means[[i]] <- Map(
          FUN = \(V, n) V*n,
          VCOV_means[[i]],
          lapply(dataList$fit_means[[i]]@extra,
                 FUN = \(x) x@dataList$nobs)
        )

        cov_params[[i]] <- unlist(lapply(dataList$fit_cov[[i]]@extra,
                                         FUN = \(x) x@parameters),
                                  recursive = FALSE)
        cov_params_labels[[i]] <-
          unlist(lapply(dataList$fit_cov[[i]]@extra,
                        FUN = \(x) x@modelInfo$trans),
                 recursive = FALSE)
        VCOV_cov[[i]] <- lapply(dataList$fit_cov[[i]]@extra,
                                FUN = \(x) x@Optim$SE$VCOV)
        NVCOV_cov[[i]] <- Map(
          FUN = \(V, n) V*n,
          VCOV_cov[[i]],
          lapply(dataList$fit_cov[[i]]@extra,
                 FUN = \(x) x@dataList$nobs)
        )
        nobs_ij[[i]] <- lapply(dataList$fit_cov[[i]]@extra,
                               FUN = \(x) x@dataList$nobs)

      } else {

        means_params[[i]] <- dataList$fit_means[[i]]@parameters
        means_params_labels[[i]] <- dataList$fit_means[[i]]@modelInfo$trans
        VCOV_means[[i]] <- list(dataList$fit_means[[i]]@Optim$SE$VCOV)
        NVCOV_means[[i]] <- list(
          VCOV_means[[i]][[1L]]*dataList$fit_means[[i]]@dataList$nobs
        )

        cov_params[[i]] <- dataList$fit_cov[[i]]@parameters
        cov_params_labels[[i]] <- dataList$fit_cov[[i]]@modelInfo$trans
        VCOV_cov[[i]] <- list(dataList$fit_cov[[i]]@Optim$SE$VCOV)
        NVCOV_cov[[i]] <- list(
          VCOV_cov[[i]][[1L]]*dataList$fit_cov[[i]]@dataList$nobs
        )
        nobs_ij[[i]] <- dataList$fit_cov[[i]]@dataList$nobs

      }

      means_names <- names(means_params[[i]])
      cov_names <- names(cov_params[[i]])

      M_group[[i]] <- means_names[startsWith(means_names, "means")]
      S_group[[i]] <- cov_names[startsWith(cov_names, "S")]
      taus_group[[i]] <- cov_names[startsWith(cov_names, "taus")]

    }

  }

  #### Result ####

  result <- list(lambda_group = lambda_group,
                 theta_group = theta_group,
                 psi_group = psi_group,
                 xtheta_group = xtheta_group,
                 xpsi_group = xpsi_group,
                 model_group = model_group,
                 nu_group = nu_group,
                 delta_group = delta_group,
                 tau_group = tau_group,
                 M_group = M_group,
                 S_group = S_group,
                 taus_group = taus_group,
                 means_params = means_params,
                 means_params_labels = means_params_labels,
                 cov_params = cov_params,
                 cov_params_labels = cov_params_labels,
                 VCOV_means = VCOV_means,
                 VCOV_cov = VCOV_cov,
                 NVCOV_means = NVCOV_means,
                 NVCOV_cov = NVCOV_cov,
                 nobs_ij = nobs_ij)

  return(result)

}

model_lcfa <- function(dataList, data_param, control) {

  list2env(data_param, envir = environment())

  ngroups <- dataList$ngroups
  list_struct <- list()
  k <- 1L

  #### Group-specific CFA parameters ####

  for(i in seq_len(ngroups)) {

    list_struct[[k]] <- list(name = lambda_group[i],
                             type = "matrix",
                             dim = c(dataList$nitems[[i]],
                                     dataList$nfactors[[i]]),
                             rownames = dataList$item_label[[i]],
                             colnames = dataList$factor_label[[i]])
    k <- k+1L

    if(dataList$positive) {

      list_struct[[k]] <- list(name = xtheta_group[i],
                               type = "matrix",
                               dim = c(dataList$nitems[[i]],
                                       dataList$nitems[[i]]),
                               rownames = dataList$item_label[[i]],
                               colnames = dataList$item_label[[i]])
      k <- k+1L

      list_struct[[k]] <- list(name = xpsi_group[i],
                               type = "matrix",
                               dim = c(dataList$nfactors[[i]],
                                       dataList$nfactors[[i]]),
                               rownames = dataList$factor_label[[i]],
                               colnames = dataList$factor_label[[i]])
      k <- k+1L

    }

    list_struct[[k]] <- list(name = theta_group[i],
                             type = "matrix",
                             dim = c(dataList$nitems[[i]],
                                     dataList$nitems[[i]]),
                             rownames = dataList$item_label[[i]],
                             colnames = dataList$item_label[[i]],
                             symmetric = TRUE)
    k <- k+1L

    list_struct[[k]] <- list(name = psi_group[i],
                             type = "matrix",
                             dim = c(dataList$nfactors[[i]],
                                     dataList$nfactors[[i]]),
                             rownames = dataList$factor_label[[i]],
                             colnames = dataList$factor_label[[i]],
                             symmetric = TRUE)
    k <- k+1L

    list_struct[[k]] <- list(name = model_group[i],
                             type = "matrix",
                             dim = c(dataList$nitems[[i]],
                                     dataList$nitems[[i]]),
                             rownames = dataList$item_label[[i]],
                             colnames = dataList$item_label[[i]],
                             symmetric = TRUE)
    k <- k+1L

    if(control$meanstructure) {

      list_struct[[k]] <- list(name = nu_group[i],
                               type = "matrix",
                               dim = c(dataList$nitems[[i]], 1L),
                               rownames = dataList$item_label[[i]],
                               colnames = "intrcp")
      k <- k+1L

    }

    if(control$deltaparam) {

      list_struct[[k]] <- list(name = delta_group[i],
                               type = "matrix",
                               dim = c(dataList$nitems[[i]], 1L),
                               rownames = dataList$item_label[[i]],
                               colnames = "latent.var")
      k <- k+1L

    }

  }

  trans <- create_parameters(list_struct)

  if(control$meanstructure) {
    trans <- c(trans, unlist(means_params_labels, recursive = FALSE))
  }

  trans <- c(trans, unlist(cov_params_labels, recursive = FALSE))

  #### Result ####

  return(trans)

}

constraints_lcfa <- function(dataList, data_param, trans, control) {

  list2env(data_param, envir = environment())

  ngroups <- dataList$ngroups
  lavaan_model <- dataList$model

  fixed <- vector("list")
  nonfixed <- vector("list")
  fixed_values_list <- vector("list")
  target_psi <- vector("list", length = ngroups)
  target_theta <- vector("list", length = ngroups)

  #### Replace latent labels by lavaan labels ####

  for(i in seq_len(ngroups)) {

    group_names <- c(lambda_group[i], theta_group[i],
                     psi_group[i], nu_group[i])

    model_blocks <- list(lavaan_model[[i]]$lambda,
                         lavaan_model[[i]]$theta,
                         lavaan_model[[i]]$psi,
                         lavaan_model[[i]]$nu)
    names(model_blocks) <- group_names

    nonfixed[group_names] <- lapply(model_blocks, FUN = \(x) {
      which(is.na(suppressWarnings(as.numeric(x))))
    })

    fixed[group_names] <- lapply(model_blocks, FUN = \(x) {
      which(!is.na(suppressWarnings(as.numeric(x))))
    })

    fixed_values_list[group_names] <- lapply(model_blocks, FUN = \(x) {
      numerals <- suppressWarnings(as.numeric(x))
      idx <- which(!is.na(numerals))
      return(numerals[idx])
    })

    trans[[lambda_group[i]]][nonfixed[[lambda_group[i]]]] <-
      lavaan_model[[i]]$lambda[nonfixed[[lambda_group[i]]]]
    trans[[theta_group[i]]][nonfixed[[theta_group[i]]]] <-
      lavaan_model[[i]]$theta[nonfixed[[theta_group[i]]]]
    trans[[psi_group[i]]][nonfixed[[psi_group[i]]]] <-
      lavaan_model[[i]]$psi[nonfixed[[psi_group[i]]]]

    if(control$meanstructure) {
      trans[[nu_group[i]]][nonfixed[[nu_group[i]]]] <-
        lavaan_model[[i]]$nu[nonfixed[[nu_group[i]]]]
    }

  }

  #### Model for the parameters ####

  param <- list()

  for(i in seq_len(ngroups)) {

    # Factor loadings:
    param[[lambda_group[i]]] <- trans[[lambda_group[i]]]
    param[[lambda_group[i]]][fixed[[lambda_group[i]]]] <-
      fixed_values_list[[lambda_group[i]]]

    # Residual and latent covariance matrices:
    if(dataList$positive) {

      param[[xtheta_group[i]]] <- trans[[xtheta_group[i]]]
      param[[xpsi_group[i]]] <- trans[[xpsi_group[i]]]

    } else {

      param[[theta_group[i]]] <- trans[[theta_group[i]]]
      param[[theta_group[i]]][fixed[[theta_group[i]]]] <-
        fixed_values_list[[theta_group[i]]]

      if(control$deltaparam) {
        diag(param[[theta_group[i]]]) <- "1"
      }

      param[[psi_group[i]]] <- trans[[psi_group[i]]]
      param[[psi_group[i]]][fixed[[psi_group[i]]]] <-
        fixed_values_list[[psi_group[i]]]

    }

    # Sample covariance matrices:
    if(control$free_S) {
      param[S_group[[i]]] <- trans[S_group[[i]]]
    } else {
      param[S_group[[i]]] <- cov_params[[i]][S_group[[i]]]
    }

    # Sample thresholds:
    if(control$free_taus) {
      param[taus_group[[i]]] <- trans[taus_group[[i]]]
    } else {
      param[taus_group[[i]]] <- cov_params[[i]][taus_group[[i]]]
    }

    # Polychoric correlations have unit diagonal:
    if(dataList$cor %in% c("poly", "polys", "polychoric", "polychorics")) {
      param[S_group[[i]]] <- lapply(param[S_group[[i]]], FUN = \(x) {
        diag(x) <- 1
        return(x)
      })
    }

    # Mean structure:
    if(control$meanstructure) {

      if(control$free_M) {
        param[M_group[[i]]] <- trans[M_group[[i]]]
      } else {
        param[M_group[[i]]] <- means_params[[i]][M_group[[i]]]
      }

      if(control$std.ov) {

        param[[nu_group[i]]] <- matrix(0,
                                       nrow = dataList$nitems[[i]],
                                       ncol = 1L,
                                       dimnames = list(dataList$item_label[[i]],
                                                       "intrcp"))

      } else {

        param[[nu_group[i]]] <- trans[[nu_group[i]]]
        param[[nu_group[i]]][fixed[[nu_group[i]]]] <-
          fixed_values_list[[nu_group[i]]]

      }

    }

    # Delta parameterization:
    if(control$deltaparam) {

      if(control$std.lv) {

        param[[delta_group[i]]] <- trans[[delta_group[i]]]

      } else {

        param[[delta_group[i]]] <- matrix(1,
                                          nrow = dataList$nitems[[i]],
                                          ncol = 1L,
                                          dimnames = list(dataList$item_label[[i]],
                                                          "latent.var"))

      }

    }

    #### Positive-definite constraint targets ####

    if(dataList$positive) {

      target_theta[[i]] <- matrix(0,
                                  nrow = dataList$nitems[[i]],
                                  ncol = dataList$nitems[[i]],
                                  dimnames = dimnames(trans[[theta_group[i]]]))
      target_theta[[i]][nonfixed[[theta_group[i]]]] <- 1

      target_psi[[i]] <- matrix(0,
                                nrow = dataList$nfactors[[i]],
                                ncol = dataList$nfactors[[i]],
                                dimnames = dimnames(trans[[psi_group[i]]]))
      target_psi[[i]][nonfixed[[psi_group[i]]]] <- 1

    }

  }

  #### Result ####

  result <- list(param = param,
                 trans = trans,
                 fixed = fixed,
                 nonfixed = nonfixed,
                 fixed_values_list = fixed_values_list,
                 target_psi = target_psi,
                 target_theta = target_theta)

  return(result)

}

start_lcfa <- function(dataList, data_param, param, trans,
                       fixed, fixed_values_list, control) {

  list2env(data_param, envir = environment())

  ngroups <- dataList$ngroups
  init_param <- vector("list", length = control$rstarts)

  for(rs in seq_len(control$rstarts)) {

    init_param[[rs]] <- list()

    for(i in seq_len(ngroups)) {

      #### Factor loadings ####

      init_param[[rs]][[lambda_group[i]]] <-
        rorth(dataList$nitems[[i]], dataList$nfactors[[i]])
      dimnames(init_param[[rs]][[lambda_group[i]]]) <-
        dimnames(trans[[lambda_group[i]]])
      init_param[[rs]][[lambda_group[i]]][fixed[[lambda_group[i]]]] <-
        fixed_values_list[[lambda_group[i]]]

      #### Covariance matrices ####

      if(dataList$positive) {

        init_param[[rs]][[xtheta_group[i]]] <-
          rorth(dataList$nitems[[i]], dataList$nitems[[i]])
        dimnames(init_param[[rs]][[xtheta_group[i]]]) <-
          dimnames(trans[[xtheta_group[i]]])

        init_param[[rs]][[xpsi_group[i]]] <-
          rorth(dataList$nfactors[[i]], dataList$nfactors[[i]])
        dimnames(init_param[[rs]][[xpsi_group[i]]]) <-
          dimnames(trans[[xpsi_group[i]]])

        init_param[[rs]][[theta_group[i]]] <-
          crossprod(init_param[[rs]][[xtheta_group[i]]])
        dimnames(init_param[[rs]][[theta_group[i]]]) <-
          dimnames(trans[[theta_group[i]]])

        init_param[[rs]][[psi_group[i]]] <-
          crossprod(init_param[[rs]][[xpsi_group[i]]])
        dimnames(init_param[[rs]][[psi_group[i]]]) <-
          dimnames(trans[[psi_group[i]]])

      } else {

        init_param[[rs]][[theta_group[i]]] <-
          diag(runif(dataList$nitems[[i]]))
        dimnames(init_param[[rs]][[theta_group[i]]]) <-
          dimnames(trans[[theta_group[i]]])
        init_param[[rs]][[theta_group[i]]][fixed[[theta_group[i]]]] <-
          fixed_values_list[[theta_group[i]]]

        init_param[[rs]][[psi_group[i]]] <-
          diag(dataList$nfactors[[i]])
        dimnames(init_param[[rs]][[psi_group[i]]]) <-
          dimnames(trans[[psi_group[i]]])
        init_param[[rs]][[psi_group[i]]][fixed[[psi_group[i]]]] <-
          fixed_values_list[[psi_group[i]]]

      }

      #### Mean structure ####

      if(control$meanstructure) {

        if(isTRUE(dataList$sample_stats_only)) {
          init_means <- dataList$sample.mean[[i]][dataList$item_label[[i]]]
        } else {
          init_means <- colMeans(dataList$data_per_group[[i]], na.rm = TRUE)
        }

        init_param[[rs]][[nu_group[i]]] <-
          matrix(init_means,
                 ncol = 1L,
                 dimnames = dimnames(trans[[nu_group[i]]]))

      }

      #### Delta parameterization ####

      if(control$deltaparam) {

        init_param[[rs]][[delta_group[i]]] <-
          matrix(1,
                 nrow = dataList$nitems[[i]],
                 ncol = 1L,
                 dimnames = dimnames(trans[[delta_group[i]]]))

      }

    }

  }

  #### Result ####

  return(init_param)

}

#### Auxiliary functions for create_lcfa_modelInfo ####

manifolds_lcfa <- function(dataList, data_param, param,
                           target_psi, target_theta, control) {

  list2env(data_param, envir = environment())

  manifolds <- list()
  k <- 1L

  add_euclidean <- function(parameters) {

    parameters <- intersect(parameters, names(param))

    if(length(parameters) > 0L) {
      manifolds[[k]] <<- list(manifold = "euclidean",
                              parameters = parameters)
      k <<- k+1L
    }

    #### Result ####

    return(invisible(NULL))

  }

  for(i in seq_len(dataList$ngroups)) {

    add_euclidean(lambda_group[i])
    add_euclidean(nu_group[i])
    add_euclidean(delta_group[i])
    add_euclidean(M_group[[i]])
    add_euclidean(S_group[[i]])

    if(dataList$positive) {

      if(xpsi_group[i] %in% names(param)) {
        manifolds[[k]] <- list(manifold = "poblq",
                               parameters = xpsi_group[i],
                               extra = list(p = dataList$nfactors[[i]],
                                            q = dataList$nfactors[[i]],
                                            constraints = target_psi[[i]]))
        k <- k+1L
      }

      if(xtheta_group[i] %in% names(param)) {
        manifolds[[k]] <- list(manifold = "poblq",
                               parameters = xtheta_group[i],
                               extra = list(p = dataList$nitems[[i]],
                                            q = dataList$nitems[[i]],
                                            constraints = target_theta[[i]]))
        k <- k+1L
      }

    } else {

      add_euclidean(psi_group[i])
      add_euclidean(theta_group[i])

    }

  }

  control_manifold <- create_manifolds(manifolds = manifolds,
                                       structures = param)

  #### Result ####

  return(control_manifold)

}

transformations_lcfa <- function(dataList, data_param, trans, control) {

  list2env(data_param, envir = environment())

  transforms <- list()
  k <- 1L

  for(i in seq_len(dataList$ngroups)) {

    #### Positive-definite covariance matrices ####

    if(dataList$positive) {

      transforms[[k]] <- list(transform = "crossprod",
                              parameters_in = xpsi_group[i],
                              parameters_out = psi_group[i],
                              extra = list(p = nrow(trans[[psi_group[i]]])))
      k <- k+1L

      transforms[[k]] <- list(transform = "crossprod",
                              parameters_in = xtheta_group[i],
                              parameters_out = theta_group[i],
                              extra = list(p = nrow(trans[[theta_group[i]]])))
      k <- k+1L

    }

    #### Delta parameterization ####

    if(control$deltaparam) {

      transforms[[k]] <- list(transform = "deltaparam",
                              parameters_in = c(delta_group[i],
                                                lambda_group[i],
                                                psi_group[i]),
                              parameters_out = list(diag(trans[[theta_group[i]]])),
                              extra = list(p = nrow(trans[[theta_group[i]]]),
                                           q = nrow(trans[[psi_group[i]]])))
      k <- k+1L

    }

    #### Model-implied covariance matrix ####

    transforms[[k]] <- list(transform = "factor_cor",
                            parameters_in = c(lambda_group[i],
                                              psi_group[i],
                                              theta_group[i]),
                            parameters_out = model_group[i],
                            extra = list(p = dataList$nitems[[i]],
                                         q = dataList$nfactors[[i]]))
    k <- k+1L

  }

  control_transform <- create_transforms(transforms = transforms,
                                         structures = trans)

  #### Result ####

  return(control_transform)

}

estimators_lcfa <- function(dataList, data_param, trans, control) {

  list2env(data_param, envir = environment())

  estimators <- list()
  k <- 1L

  #### CFA discrepancy functions ####

  for(i in seq_len(dataList$ngroups)) {

    cfa_estimator <- switch(tolower(dataList$estimator),
                            uls = "cfa_dwls",
                            dwls = "cfa_dwls",
                            ml = "cfa_fml",
                            fml = "cfa_fml",
                            means_fml = "cfa_means_fml",
                            means_dwls = "cfa_means_dwls",
                            means_uls = "cfa_means_dwls",
                            stop("Unknown estimator: ", dataList$estimator))

    cov_params_i <- cov_params[[i]]
    cov_params_i <- cov_params_i[startsWith(names(cov_params_i), "S")]

    for(j in seq_along(cov_params_i)) {

      pick <- rownames(cov_params_i[[j]])
      S_group_ij <- S_group[[i]][[j]]
      M_group_ij <- M_group[[i]][[j]]
      p <- nrow(cov_params_i[[j]])

      if(control$estimator %in%
         c("uls", "means_uls", "ml", "fml", "means_fml")) {

        W_cov <- matrix(1, nrow = p, ncol = p)

      } else {

        idx <- startsWith(rownames(NVCOV_cov[[i]][[j]]), "S")
        W_cov <- matrix(NA_real_, nrow = p, ncol = p)
        W_cov[lower.tri(W_cov, diag = !control$std.ov)] <-
          diag(NVCOV_cov[[i]][[j]][idx, idx, drop = FALSE])
        W_cov[upper.tri(W_cov)] <- t(W_cov)[upper.tri(W_cov)]
        W_cov <- 1/W_cov

        if(control$std.ov) diag(W_cov) <- 0

      }

      w_means <- diag(NVCOV_means[[i]][[j]])

      model_parameters <- c(trans[[model_group[i]]][pick, pick])
      sample_covariance <- c(trans[[S_group_ij]])

      if(control$meanstructure) {
        model_means <- c(trans[[nu_group[i]]][pick, ])
        sample_means <- c(trans[[M_group_ij]])
      } else {
        model_means <- numeric(0L)
        sample_means <- numeric(0L)
      }

      estimators[[k]] <- list(estimator = cfa_estimator,
                              parameters = list(model_parameters,
                                                sample_covariance,
                                                model_means,
                                                sample_means),
                              extra = list(W = W_cov,
                                           w = nobs_ij[[i]][[j]]/
                                             sum(unlist(dataList$nobs)),
                                           w_means = w_means,
                                           q = nrow(trans[[psi_group[i]]]),
                                           p = p,
                                           n = nobs_ij[[i]][[j]],
                                           double_names = S_group_ij,
                                           matrix_names = c(S_group_ij,
                                                            M_group_ij)))
      k <- k+1L

    }

  }

  #### Penalties ####

  if(control$reg) {

    for(i in seq_len(dataList$ngroups)) {

      lower_indices <- which(lower.tri(trans[[psi_group[i]]], diag = TRUE))
      estimators[[k]] <- list(estimator = "logdetR",
                              parameters = psi_group[i],
                              extra = list(lower_indices = lower_indices-1L,
                                           p = nrow(trans[[psi_group[i]]]),
                                           logdetw = control$penalties$logdet$w,
                                           double_names = "logdetR psi"))
      k <- k+1L

      lower_indices <- which(lower.tri(trans[[theta_group[i]]], diag = TRUE))
      estimators[[k]] <- list(estimator = "logdetR",
                              parameters = theta_group[i],
                              extra = list(lower_indices = lower_indices-1L,
                                           p = nrow(trans[[theta_group[i]]]),
                                           logdetw = control$penalties$logdet$w,
                                           double_names = "logdetR theta"))
      k <- k+1L

    }

  }

  control_estimator <- create_estimators(estimators = estimators,
                                         structures = trans)

  #### Result ####

  return(control_estimator)

}
