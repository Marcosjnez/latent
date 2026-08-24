# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 24/08/2026
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
#'      meanstructure = TRUE,
#'      parameterization = NULL,
#'      likelihood = NULL, se = TRUE,
#'      control = NULL, message = FALSE,
#'      do.fit = TRUE, ...)
#'
#' @param data Optional data frame or matrix containing the observed variables.
#'   If NULL, sample.cov and sample.nobs must be supplied.
#' @param model Confirmatory factor model specified using lavaan syntax.
#' @param estimator Estimation method. Available options include \code{"ml"},
#'   \code{"uls"}, and \code{"dwls"}. The value \code{"fiml"} requests
#'   direct pattern-likelihood FIML unless \code{missing = "fiml"} is
#'   supplied explicitly.
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
#' @param missing Missing-data method. \code{"ml"} uses direct
#'   pattern-likelihood FIML. \code{"fiml"} first estimates saturated
#'   incomplete-data moments with \code{lmvnorm()} and then fits CFA as a
#'   deterministic multistep estimator.
#' @param std.lv Logical. Standardize latent variables.
#' @param std.ov Logical. Standardize observed variables. With direct or
#'   saturated-moment FIML, raw variables are standardized within
#'   substantive groups before missingness patterns are constructed;
#'   observed-variable means remain freely estimated.
#' @param meanstructure Logical. Estimate the observed-variable mean structure.
#' @param parameterization Optional parameterization specification.
#' @param likelihood Character string controlling the normal/Wishart likelihood
#'   convention.
#' @param se Logical or character. \code{TRUE}, \code{"standard"}, and
#'   \code{"information"} use standard sampling covariance matrices for the
#'   sample statistics. \code{"robust"} requests robust sampling covariance
#'   matrices where implemented. \code{FALSE} skips computation of the final
#'   CFA standard errors. Ordinary ML and direct FIML use information from the
#'   likelihood; multistep analyses propagate the covariance of their parent
#'   statistics with \code{se.multistep()}.
#' @param control Optional list of optimization controls.
#' @param message Logical. Print progress messages.
#' @param do.fit Logical. If \code{FALSE}, return the prepared but unfitted
#'   \code{"lcfa"} or \code{"multistep_lcfa"} object.
#' @param ... Additional arguments passed to lavaan and the sample-statistic
#'   estimators where applicable.
#'
#' @details
#' The model-implied observed means are computed as
#' \deqn{\widehat{\mu}=\nu+\Lambda\alpha,}
#' where \eqn{\nu} contains observed-variable intercepts and \eqn{\alpha}
#' contains latent-factor means. For ordinal models, standardized model
#' thresholds are computed from the unstandardized thresholds, model-implied
#' means, and model-implied variances.
#'
#' Direct FIML creates one likelihood contribution for every missingness
#' pattern and substantive group. Saturated-moment FIML instead stores one
#' \code{lmvnorm} source object in \code{extra}; its uncertainty is propagated
#' automatically.
#'
#' @return An S4 object of class \code{"lcfa"} for ordinary likelihood
#'   analyses, or \code{"multistep_lcfa"} when at least one parent object is
#'   marked for uncertainty propagation.
#'
#' @examples
#' \dontrun{
#' HS.model <- '
#'   visual  =~ x1 + x2 + x3
#'   textual =~ x4 + x5 + x6
#'   speed   =~ x7 + x8 + x9
#' '
#' fit <- lcfa(model = HS.model, data = HolzingerSwineford1939)
#' }
#'
#' @export
lcfa <- function(data = NULL, model = NULL, estimator = "ml",
                 ordered = FALSE, group = NULL,
                 sample.cov = NULL, sample.mean = NULL, sample.nobs = NULL,
                 positive = FALSE, penalties = FALSE,
                 missing = "pairwise.complete.obs",
                 std.lv = FALSE, std.ov = FALSE,
                 meanstructure = TRUE,
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

  missing_was_supplied <- !missing(missing)
  estimator <- tolower(estimator)
  missing <- tolower(missing)

  if(estimator == "fiml") {
    estimator <- "ml"
    if(!missing_was_supplied) missing <- "ml"
  }

  se_control <- normalize_lcfa_se(se)
  sample_se <- se_control$sample
  compute_se <- se_control$compute

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

  if(is.null(likelihood)) {
    likelihood <- "normal"
  }

  if(!is.character(likelihood) ||
     length(likelihood) != 1L ||
     is.na(likelihood)) {
    stop("likelihood must be a single character string")
  }

  likelihood <- tolower(likelihood)

  if(missing %in% c("ml", "fiml")) {

    if(is.null(data)) {
      stop("missing = '", missing, "' requires raw data")
    }

    if(!isFALSE(ordered)) {
      stop("missing = '", missing, "' currently requires continuous variables")
    }

    if(likelihood != "normal") {
      stop("missing = '", missing, "' requires likelihood = 'normal'")
    }

    meanstructure <- TRUE

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

    if(sample_se == "robust") {
      stop("se = 'robust' requires raw data")
    }

  }

  if(sample_se == "robust" &&
     !(missing %in% c("ml", "fiml")) &&
     (meanstructure || cor != "pearson")) {
    warning("Robust sampling covariance is currently implemented only for ",
            "Pearson covariance/correlation statistics. Sample means, ",
            "polychoric correlations, and Yule correlations use their ",
            "standard covariance matrices.")
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
  control$sample_se <- sample_se
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
                                   se = sample_se,
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

  #### Previous sample-statistic models ####

  previous_models <- previous_models_lcfa(dataList)

  propagate <- vapply(previous_models,
                      FUN = \(x) isTRUE(x@modelInfo$propagate_uncertainty),
                      FUN.VALUE = logical(1L))

  numeric_multistep <- isTRUE(dataList$sample_stats_only) &&
    dataList$estimator %in% c("uls", "dwls")

  object_class <- if(any(propagate) || numeric_multistep) {
    "multistep_lcfa"
  } else {
    "lcfa"
  }

  if(sample_se == "robust" && identical(object_class, "lcfa")) {
    stop("Robust standard errors are not yet implemented for ordinary lcfa ",
         "likelihood models. Use information standard errors, or a multistep ",
         "analysis whose parent statistic supplies a robust VCOV matrix.")
  }

  modelInfo$propagate_uncertainty <- TRUE
  modelInfo$step_labels <- modelInfo$parameters_labels

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
                extra            = previous_models)

  #### Standard errors ####

  if(compute_se) {

    if(message) {
      print_lcfa_message("Computing standard errors")
    }

    if(inherits(result, "multistep")) {

      result@Optim$SE <- se.multistep(
        fit = result,
        parameters = modelInfo$parameters_labels
      )

    } else {

      result@Optim$SE <- se.latent(
        fit = result,
        type = if(sample_se == "robust") "robust" else "information",
        parameters = modelInfo$parameters_labels
      )

    }

    result@Optim$SE$sample_se <- sample_se

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
                                 se = "standard", message = FALSE,
                                 likelihood = NULL, meanstructure = TRUE,
                                 args = NULL, control = NULL,
                                 ...) {

  cor <- tolower(cor)
  estimator <- tolower(estimator)
  se_type <- match.arg(tolower(se), c("standard", "robust"))
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
    group_label <- ""

  } else {

    if(!group %in% colnames(data)) {
      stop("The grouping variable is not present in data")
    }

    if(anyNA(data[[group]])) {
      stop("The grouping variable cannot contain missing values")
    }

    group_label <- unique(as.character(data[[group]]))

    if(any(group_label == "")) {
      stop("Group labels must be non-empty")
    }

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

      group_i <- lcfa_group_rows(data = data,
                                 group = group,
                                 group_label = group_label,
                                 i = i,
                                 ngroups = ngroups)
      items_i <- item_names[[i]]

      not_all_na_i <- !apply(is.na(data[group_i, items_i, drop = FALSE]),
                             MARGIN = 1L, FUN = all)
      keep[group_i] <- not_all_na_i

    }

    data <- data[keep, , drop = FALSE]

    if(nrow(data) == 0L) {
      stop("No observations remain after removing cases with all model variables missing")
    }

    if(std.ov && missing %in% c("ml", "fiml")) {
      data <- standardize_lcfa_raw_data(data = data,
                                        group = group,
                                        group_label = group_label,
                                        item_names = item_names,
                                        ngroups = ngroups)
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
  fit_moments <- NULL
  fiml_moments <- FALSE
  direct_fiml <- FALSE
  direct_patterns <- vector("list", length = ngroups)

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

      NVCOV[[i]] <- asymptotic_normal(S,
                                      cov = !std.ov,
                                      diag = FALSE)
      VCOV_cov[[i]] <- NVCOV[[i]]/nobs_list[[i]]
      WLS.V[[i]] <- diag(1/diag(NVCOV[[i]]))

      if(std.ov) {

        NVCOV_means[[i]] <- matrix(0,
                                    nrow = nrow(S_input),
                                    ncol = ncol(S_input),
                                    dimnames = dimnames(S_input))

      } else {

        NVCOV_means[[i]] <- S_input

      }

      VCOV_means[[i]] <- NVCOV_means[[i]]/nobs_list[[i]]
      thresholds[[i]] <- list()

    }

  } else {

    sample.cov <- normalize_lcfa_group_input(sample.cov, ngroups)

    #### Data per group ####

    for(i in seq_len(ngroups)) {

      X[[i]] <- extract_lcfa_group_data(data = data,
                                        group = group,
                                        group_label = group_label,
                                        item_names = item_names,
                                        i = i,
                                        ngroups = ngroups)

    }

    has_missing <- any(vapply(X, FUN = anyNA, FUN.VALUE = logical(1L)))
    direct_fiml <- missing == "ml" && has_missing
    fiml_moments <- missing == "fiml" && has_missing

    if(fiml_moments) {

      if(cor != "pearson") {
        stop("The saturated multivariate-normal moment estimator requires continuous variables")
      }

      if(se_type == "robust") {
        stop("Robust covariance estimation is not yet implemented for ",
             "incomplete multivariate-normal moments")
      }

      #### Saturated incomplete-data moments ####

      control_moments <- control
      control_moments$start <- NULL
      control_moments$rstarts <- 1L
      control_moments$cores <- 1L

      if(ngroups < 2L) {
        data_moments <- X[[1L]]
        group_moments <- NULL
        variables_moments <- item_names[[1L]]
      } else {
        data_moments <- data
        group_moments <- group
        variables_moments <- item_names
      }

      fit_moments <- lmvnorm(data = data_moments,
                             group = group_moments,
                             variables = variables_moments,
                             se = "information",
                             do.fit = TRUE,
                             message = message,
                             control = control_moments)

      moment_param <- fit_moments@modelInfo$data_param
      moment_VCOV <- fit_moments@Optim$SE$VCOV

      for(i in seq_len(ngroups)) {

        M_name <- moment_param$means_group[i]
        S_name <- moment_param$S_group[i]

        M <- fit_moments@transformed_pars[[M_name]]
        S <- fit_moments@transformed_pars[[S_name]]

        mean_labels <- fit_moments@modelInfo$parameters_labels[
          fit_moments@modelInfo$parameters_labels %in%
            c(fit_moments@modelInfo$trans[[M_name]])
        ]
        S_labels <- fit_moments@modelInfo$parameters_labels[
          fit_moments@modelInfo$parameters_labels %in%
            c(fit_moments@modelInfo$trans[[S_name]])
        ]

        nobs_list[[i]] <- fit_moments@dataList$nobs_group[[i]]
        sample.mean_list[[i]] <- c(M)
        names(sample.mean_list[[i]]) <- rownames(M)
        sample.cov[[i]] <- S
        sample.cov_input[[i]] <- S

        VCOV_means[[i]] <- moment_VCOV[mean_labels, mean_labels,
                                       drop = FALSE]
        VCOV_cov[[i]] <- moment_VCOV[S_labels, S_labels,
                                     drop = FALSE]
        NVCOV_means[[i]] <- VCOV_means[[i]]*nobs_list[[i]]
        NVCOV[[i]] <- VCOV_cov[[i]]*nobs_list[[i]]
        WLS.V[[i]] <- diag(1/diag(NVCOV[[i]]))
        thresholds[[i]] <- list()

      }

    } else if(direct_fiml) {

      #### Direct pattern-likelihood FIML ####

      for(i in seq_len(ngroups)) {

        direct_patterns[[i]] <- create_lmvnorm_patterns(X[[i]])
        nobs_list[[i]] <- nrow(X[[i]])
        sample.mean_list[[i]] <- colMeans(X[[i]], na.rm = TRUE)

        S_pairwise <- initial_lmvnorm_covariance(X[[i]])
        S_pairwise <- S_pairwise*(nobs_list[[i]]-1L)/nobs_list[[i]]
        dimnames(S_pairwise) <- list(item_names[[i]], item_names[[i]])

        sample.cov[[i]] <- S_pairwise
        sample.cov_input[[i]] <- S_pairwise
        NVCOV[[i]] <- asymptotic_normal(S_pairwise,
                                        cov = TRUE,
                                        diag = FALSE)
        VCOV_cov[[i]] <- NVCOV[[i]]/nobs_list[[i]]
        NVCOV_means[[i]] <- S_pairwise
        VCOV_means[[i]] <- S_pairwise/nobs_list[[i]]
        WLS.V[[i]] <- diag(1/diag(NVCOV[[i]]))
        thresholds[[i]] <- list()

      }

    } else {

      #### Complete-data or ordinary sample statistics ####

      for(i in seq_len(ngroups)) {

        control_i <- control

        if(ngroups < 2L) {
          control_i$subfix <- ""
        } else {
          control_i$subfix <- paste0(".", group_label[i])
        }

        missing_i <- if(missing %in% c("ml", "fiml")) {
          "pairwise.complete.obs"
        } else {
          missing
        }

        direct_complete_ml <- cor == "pearson" &&
          control$estimator %in% c("ml", "fml") &&
          !anyNA(X[[i]])

        source_se <- !direct_complete_ml
        control_i$propagate_uncertainty <- source_se

        sample_stats <- estimate_lcfa_sample_statistics(
          data = X[[i]],
          model = model,
          cor = cor,
          std.ov = std.ov,
          se = se_type,
          compute_se = source_se,
          likelihood = likelihood,
          missing = missing_i,
          control = control_i,
          ...
        )

        fit_means[[i]] <- sample_stats$fit_means
        fit_cov[[i]] <- sample_stats$fit_cov

        nobs_list[[i]] <- fit_cov[[i]]@dataList$nobs
        S_names <- names(fit_cov[[i]]@transformed_pars)
        S_names <- S_names[startsWith(S_names, "S")]

        if(length(S_names) != 1L) {
          stop("Exactly one sample covariance/correlation matrix was expected per group.")
        }

        sample.cov[[i]] <- fit_cov[[i]]@transformed_pars[[S_names]]
        sample.cov_input[[i]] <- sample.cov[[i]]

        if(source_se) {

          VCOV_cov[[i]] <- stored_sample_vcov_lcfa(fit_cov[[i]])
          NVCOV[[i]] <- VCOV_cov[[i]]*fit_cov[[i]]@dataList$nobs
          VCOV_means[[i]] <- stored_sample_vcov_lcfa(fit_means[[i]])
          NVCOV_means[[i]] <-
            VCOV_means[[i]]*fit_means[[i]]@dataList$nobs

        } else {

          NVCOV[[i]] <- asymptotic_normal(sample.cov[[i]],
                                          cov = !std.ov,
                                          diag = FALSE)
          VCOV_cov[[i]] <- NVCOV[[i]]/fit_cov[[i]]@dataList$nobs

          if(std.ov) {
            NVCOV_means[[i]] <- matrix(
              0, nrow = nrow(sample.cov[[i]]),
              ncol = ncol(sample.cov[[i]]),
              dimnames = dimnames(sample.cov[[i]])
            )
          } else {
            NVCOV_means[[i]] <- sample.cov[[i]]
          }

          VCOV_means[[i]] <-
            NVCOV_means[[i]]/fit_means[[i]]@dataList$nobs

        }

        WLS.V[[i]] <- diag(1/diag(NVCOV[[i]]))

        idx_taus <- startsWith(names(fit_cov[[i]]@transformed_pars),
                               "taus")
        thresholds[[i]] <- fit_cov[[i]]@transformed_pars[idx_taus]

        idx_means <- startsWith(names(fit_means[[i]]@transformed_pars),
                                "means")

        if(any(idx_means)) {

          M <- fit_means[[i]]@transformed_pars[[which(idx_means)[1L]]]
          sample.mean_list[[i]] <- c(M)
          names(sample.mean_list[[i]]) <- rownames(M)

        } else {

          sample.mean_list[[i]] <-
            setNames(rep(0, length(item_names[[i]])), item_names[[i]])

        }

      }

    }

  }

  #### Lavaan model structure ####

  sample.cov_lav <- if(sample_stats_only) sample.cov_input else sample.cov
  sample.th_lav <- NULL

  if(cor == "poly") {

    sample.th_lav <- lapply(seq_len(ngroups), FUN = \(i) {
      lcfa_lavaan_sample_thresholds(
        threshold_blocks = thresholds[[i]],
        item_names = item_names[[i]]
      )
    })

    threshold_indices <- lapply(sample.th_lav, FUN = \(x) {
      attr(x, "th.idx")
    })

    if(ngroups == 1L) {

      sample.th_lav <- sample.th_lav[[1L]]

    } else {

      sample.th_lav <- lapply(sample.th_lav, FUN = \(x) {
        attr(x, "th.idx") <- NULL
        return(x)
      })
      names(sample.th_lav) <- group_label
      names(threshold_indices) <- group_label
      attr(sample.th_lav, "th.idx") <- threshold_indices

    }

  }

  if(ngroups > 1L && (sample_stats_only || fiml_moments || direct_fiml)) {
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
                     sample.th = sample.th_lav,
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
  nfactors <- lapply(model_out, FUN = \(x) ncol(x$lambda))

  #### Normalize mean and threshold structures ####

  for(i in seq_len(ngroups)) {

    p <- nitems[[i]]
    q <- nfactors[[i]]
    items <- item_label[[i]]
    factors <- factor_label[[i]]

    if(is.null(model_out[[i]]$nu) || length(model_out[[i]]$nu) == 0L) {
      model_out[[i]]$nu <- matrix(0, nrow = p, ncol = 1L,
                                  dimnames = list(items, "intrcp"))
    }

    if(is.null(model_out[[i]]$alpha) || length(model_out[[i]]$alpha) == 0L) {
      model_out[[i]]$alpha <- matrix(0, nrow = q, ncol = 1L,
                                     dimnames = list(factors, "intrcp"))
    }

    if(cor == "poly") {

      if(is.null(model_out[[i]]$tau) || length(model_out[[i]]$tau) == 0L) {
        stop("The ordinal CFA model did not produce threshold parameters.")
      }

      tau_names <- if(ngroups == 1L) {
        names(sample.th_lav)
      } else {
        names(sample.th_lav[[i]])
      }

      if(is.null(dim(model_out[[i]]$tau))) {
        model_out[[i]]$tau <- matrix(model_out[[i]]$tau, ncol = 1L)
      }

      if(length(tau_names) != nrow(model_out[[i]]$tau)) {
        stop("The lavaan threshold matrix and sample thresholds have different lengths.")
      }

      rownames(model_out[[i]]$tau) <- tau_names

      if(is.null(colnames(model_out[[i]]$tau))) {
        colnames(model_out[[i]]$tau) <- "threshold"
      }

    }

  }

  #### Number of modeled sample statistics ####

  if(cor == "poly") {

    npatterns <- lapply(seq_len(ngroups), FUN = \(i) {
      p <- nitems[[i]]
      p*(p-1L)/2L+length(model_out[[i]]$tau)
    })

  } else if(meanstructure) {

    npatterns <- lapply(nitems, FUN = \(p) 0.5*p*(p+1)+p)

  } else {

    npatterns <- lapply(nitems, FUN = \(p) 0.5*p*(p+1))

  }

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
                   likelihood = likelihood,
                   cor = cor,
                   se_type = se_type,
                   meanstructure = meanstructure,
                   missing = missing,
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
                   fit_cov = fit_cov,
                   fit_moments = fit_moments,
                   fiml_moments = fiml_moments,
                   direct_fiml = direct_fiml,
                   direct_patterns = direct_patterns)

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

lcfa_group_rows <- function(data, group, group_label,
                            i, ngroups) {

  if(ngroups == 1L || is.null(group)) {
    result <- seq_len(nrow(data))
  } else {
    result <- which(as.character(data[[group]]) ==
                      as.character(group_label[i]))
  }

  if(length(result) == 0L) {
    stop("Group '", group_label[i], "' contains no observations.")
  }

  #### Result ####

  return(result)

}

standardize_lcfa_raw_data <- function(data, group, group_label,
                                      item_names, ngroups) {

  for(i in seq_len(ngroups)) {

    group_i <- lcfa_group_rows(data = data,
                               group = group,
                               group_label = group_label,
                               i = i,
                               ngroups = ngroups)
    items_i <- item_names[[i]]
    X_i <- data[group_i, items_i, drop = FALSE]

    means_i <- colMeans(X_i, na.rm = TRUE)
    sds_i <- vapply(X_i, FUN = stats::sd, FUN.VALUE = numeric(1L),
                    na.rm = TRUE)

    invalid <- !is.finite(means_i) |
      !is.finite(sds_i) |
      sds_i <= 0

    if(any(invalid)) {
      stop("Observed variables cannot be standardized because at least one ",
           "variable has a non-finite mean or a non-positive/non-finite ",
           "standard deviation in group '", group_label[i], "': ",
           paste(items_i[invalid], collapse = ", "))
    }

    X_i <- sweep(X_i, MARGIN = 2L, STATS = means_i, FUN = "-")
    X_i <- sweep(X_i, MARGIN = 2L, STATS = sds_i, FUN = "/")
    data[group_i, items_i] <- X_i

  }

  #### Result ####

  return(data)

}

extract_lcfa_group_data <- function(data, group, group_label,
                                    item_names, i, ngroups) {

  rows <- lcfa_group_rows(data = data,
                          group = group,
                          group_label = group_label,
                          i = i,
                          ngroups = ngroups)

  result <- data[rows, item_names[[i]], drop = FALSE]

  #### Result ####

  return(result)

}

estimate_lcfa_sample_statistics <- function(data, model, cor,
                                            std.ov, se, compute_se = TRUE,
                                            likelihood, missing, control, ...) {

  #### Means ####

  control$propagate_uncertainty <- isTRUE(compute_se)

  fit_means <- lmean(data = data,
                     std.ov = std.ov,
                     se = compute_se,
                     control = control,
                     do.fit = TRUE)

  #### Covariances ####

  if(cor == "pearson") {

    fit_cov <- lpearson(data = data,
                        std.ov = std.ov,
                        VCOV = if(se == "robust") "robust" else "standard",
                        likelihood = likelihood,
                        missing = missing,
                        se = compute_se,
                        control = control,
                        do.fit = TRUE)

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

  #### Sampling covariance methods ####

  fit_means@modelInfo$propagate_uncertainty <- isTRUE(compute_se)
  fit_cov@modelInfo$propagate_uncertainty <- isTRUE(compute_se)

  if(length(fit_means@Optim$SE) > 0L) {
    fit_means@Optim$SE$type <- "standard"
  }

  if(length(fit_cov@Optim$SE) > 0L) {
    if(cor == "pearson") {
      fit_cov@Optim$SE$type <- se
    } else {
      fit_cov@Optim$SE$type <- "standard"
    }
  }

  #### Result ####

  result <- list(fit_means = fit_means,
                 fit_cov = fit_cov)

  return(result)

}

stored_sample_vcov_lcfa <- function(fit) {

  labels <- fit@modelInfo$parameters_labels
  VCOV <- tryCatch(fit@Optim$SE$VCOV,
                   error = function(e) NULL)

  if(is.null(VCOV)) {
    VCOV <- matrix(0,
                   nrow = length(labels),
                   ncol = length(labels),
                   dimnames = list(labels, labels))
  }

  if(!is.matrix(VCOV)) {
    VCOV <- as.matrix(VCOV)
  }

  if(!is.null(rownames(VCOV)) &&
     !is.null(colnames(VCOV)) &&
     all(labels %in% rownames(VCOV)) &&
     all(labels %in% colnames(VCOV))) {
    VCOV <- VCOV[labels, labels, drop = FALSE]
  }

  rownames(VCOV) <- colnames(VCOV) <- labels

  #### Result ####

  return(VCOV)

}

normalize_lcfa_se <- function(se) {

  if(isTRUE(se)) {

    result <- list(sample = "standard", compute = TRUE)

  } else if(isFALSE(se)) {

    result <- list(sample = "standard", compute = FALSE)

  } else if(is.character(se) &&
            length(se) == 1L &&
            !is.na(se)) {

    se <- match.arg(tolower(se),
                    c("standard", "information", "robust"))

    result <- list(sample = if(se == "robust") "robust" else "standard",
                   compute = TRUE)

  } else {

    stop("se must be TRUE, FALSE, 'standard', 'information', or 'robust'")

  }

  #### Result ####

  return(result)

}

previous_models_lcfa <- function(dataList) {

  if(inherits(dataList$fit_moments, "latent")) {

    #### Result ####

    return(list(dataList$fit_moments))

  }

  models <- list()

  append_model <- function(object) {

    if(!inherits(object, "latent")) {
      return(invisible(NULL))
    }

    if(inherits(object, "multistep")) {
      models[[length(models)+1L]] <<- object
      return(invisible(NULL))
    }

    children <- object@extra[
      vapply(object@extra,
             FUN = \(x) inherits(x, "latent"),
             FUN.VALUE = logical(1L))
    ]

    if(length(children) > 0L) {
      for(child in children) append_model(child)
    } else {
      models[[length(models)+1L]] <<- object
    }

    return(invisible(NULL))

  }

  if(isTRUE(dataList$meanstructure) && dataList$cor == "pearson") {
    for(object in dataList$fit_means) append_model(object)
  }

  for(object in dataList$fit_cov) append_model(object)

  unique_models <- list()

  for(model in models) {

    duplicated <- length(unique_models) > 0L &&
      any(vapply(unique_models,
                 FUN = \(x) identical(x, model),
                 FUN.VALUE = logical(1L)))

    if(!duplicated) {
      unique_models[[length(unique_models)+1L]] <- model
    }

  }

  #### Result ####

  return(unique_models)

}

#### Threshold helpers ####

lcfa_lavaan_sample_thresholds <- function(threshold_blocks, item_names) {

  if(!is.list(threshold_blocks) || length(threshold_blocks) == 0L) {
    stop("Sample threshold blocks are unavailable for the ordinal CFA model.")
  }

  block_variables <- vapply(seq_along(threshold_blocks), FUN = \(b) {

    block <- threshold_blocks[[b]]
    cn <- colnames(block)

    if(length(cn) == 1L && !is.na(cn) && nzchar(cn)) {
      return(cn)
    }

    block_names <- names(threshold_blocks)
    nm <- if(is.null(block_names) || length(block_names) < b) {
      NA_character_
    } else {
      block_names[b]
    }

    if(is.na(nm) || !nzchar(nm)) {
      stop("Every sample threshold block must identify its observed variable.")
    }

    nm <- sub("^taus", "", nm)
    nm <- sub("\\..*$", "", nm)
    return(nm)

  }, FUN.VALUE = character(1L))

  if(anyDuplicated(block_variables)) {
    stop("More than one sample threshold block was found for an observed variable.")
  }

  values <- vector("list", length(item_names))
  indices <- vector("list", length(item_names))

  for(j in seq_along(item_names)) {

    variable <- item_names[j]
    b <- match(variable, block_variables)

    if(is.na(b)) {
      stop("No sample threshold block was found for variable '", variable, "'.")
    }

    values[[j]] <- as.numeric(c(threshold_blocks[[b]]))

    if(length(values[[j]]) == 0L || any(!is.finite(values[[j]]))) {
      stop("Every ordinal variable must have at least one finite sample threshold.")
    }

    threshold_names <- paste0(variable, "|t", seq_along(values[[j]]))
    indices[[j]] <- setNames(rep.int(j, length(values[[j]])),
                             threshold_names)

  }

  result <- unlist(values, use.names = FALSE)
  threshold_indices <- unlist(indices, use.names = TRUE)
  names(result) <- names(threshold_indices)
  attr(result, "th.idx") <- threshold_indices

  #### Result ####

  return(result)

}

lcfa_threshold_items <- function(model_tau, item_names) {

  if(is.null(model_tau) || length(model_tau) == 0L) {

    #### Result ####

    return(integer(0L))

  }

  threshold_names <- rownames(model_tau)

  if(is.null(threshold_names) ||
     any(is.na(threshold_names)) ||
     any(threshold_names == "")) {
    stop("The lavaan threshold matrix must have non-empty row names.")
  }

  threshold_variables <- sub("\\|.*$", "", threshold_names)
  result <- match(threshold_variables, item_names)

  if(anyNA(result)) {
    stop("At least one model threshold could not be matched to an observed variable.")
  }

  #### Result ####

  return(result)

}

lcfa_order_thresholds <- function(threshold_blocks, model_tau) {

  if(is.null(model_tau) || length(model_tau) == 0L) {

    #### Result ####

    return(numeric(0L))

  }

  if(!is.list(threshold_blocks) || length(threshold_blocks) == 0L) {
    stop("Sample threshold blocks are unavailable for the ordinal CFA model.")
  }

  threshold_variables <- sub("\\|.*$", "", rownames(model_tau))

  block_variables <- vapply(seq_along(threshold_blocks), FUN = \(b) {

    block <- threshold_blocks[[b]]
    cn <- colnames(block)

    if(length(cn) == 1L && !is.na(cn) && nzchar(cn)) {
      return(cn)
    }

    block_names <- names(threshold_blocks)
    nm <- if(is.null(block_names) || length(block_names) < b) {
      NA_character_
    } else {
      block_names[b]
    }

    if(is.na(nm) || !nzchar(nm)) {
      stop("Every sample threshold block must identify its observed variable.")
    }

    nm <- sub("^taus", "", nm)
    nm <- sub("\\..*$", "", nm)
    return(nm)

  }, FUN.VALUE = character(1L))

  if(anyDuplicated(block_variables)) {
    stop("More than one sample threshold block was found for an observed variable.")
  }

  seen <- setNames(integer(length(unique(threshold_variables))),
                   unique(threshold_variables))
  result <- vector("list", length(threshold_variables))

  for(h in seq_along(threshold_variables)) {

    variable <- threshold_variables[h]
    b <- match(variable, block_variables)

    if(is.na(b)) {
      stop("No sample threshold block was found for variable '", variable, "'.")
    }

    seen[variable] <- seen[variable]+1L
    values <- c(threshold_blocks[[b]])
    position <- seen[variable]

    if(position > length(values)) {
      stop("The model contains more thresholds than the sample block for variable '",
           variable, "'.")
    }

    result[[h]] <- values[position]

  }

  result <- unlist(result, use.names = FALSE)

  #### Result ####

  return(result)

}

lcfa_diagonal_weights <- function(NVCOV, labels, zero_allowed = FALSE,
                                  object_name = "sample statistic") {

  if(length(labels) == 0L) {

    #### Result ####

    return(numeric(0L))

  }

  if(is.null(rownames(NVCOV)) || is.null(colnames(NVCOV))) {
    stop("The N-scaled covariance matrix for ", object_name,
         " must have row and column names.")
  }

  idx <- match(labels, rownames(NVCOV))

  if(anyNA(idx)) {
    stop("At least one ", object_name,
         " label could not be matched to its N-scaled covariance matrix.")
  }

  variances <- diag(NVCOV)[idx]

  if(any(!is.finite(variances)) || any(variances < 0)) {
    stop("The N-scaled variances of ", object_name,
         " must be finite and non-negative.")
  }

  if(!zero_allowed && any(variances <= 0)) {
    stop("The N-scaled variances of ", object_name,
         " must be strictly positive.")
  }

  result <- numeric(length(variances))
  positive <- variances > 0
  result[positive] <- 1/variances[positive]

  #### Result ####

  return(result)

}

#### Auxiliary functions for create_lcfa_model ####

create_lcfa_data_param <- function(dataList, control) {

  ngroups <- dataList$ngroups
  sep <- if(ngroups < 2L) "" else "."

  #### Model parameter-block names ####

  lambda_group <- paste("lambda", dataList$group_label, sep = sep)
  alpha_group <- paste("alpha", dataList$group_label, sep = sep)
  psi_group <- paste("psi", dataList$group_label, sep = sep)
  theta_group <- paste("theta", dataList$group_label, sep = sep)
  logvars_group <- paste("logvars", dataList$group_label, sep = sep)
  xpsi_group <- paste("xpsi", dataList$group_label, sep = sep)
  xtheta_group <- paste("xtheta", dataList$group_label, sep = sep)
  model_group <- paste("model", dataList$group_label, sep = sep)
  nu_group <- paste("nu", dataList$group_label, sep = sep)
  meanshat_group <- paste("meanshat", dataList$group_label, sep = sep)
  delta_group <- paste("delta", dataList$group_label, sep = sep)
  kappa_group <- paste("kappa", dataList$group_label, sep = sep)
  tauhat_group <- paste("tauhat", dataList$group_label, sep = sep)

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

      subfix <- if(ngroups < 2L) "" else paste0(".", dataList$group_label[i])
      items <- dataList$item_label[[i]]
      p <- length(items)

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
               nrow = p, ncol = 1L,
               dimnames = list(items, "intrcpt"))

      mean_labels <- c(means_struct[[means_name]])
      VCOV_means[[i]] <- list(dataList$VCOV_means[[i]])
      NVCOV_means[[i]] <- list(dataList$NVCOV_means[[i]])
      rownames(VCOV_means[[i]][[1L]]) <-
        colnames(VCOV_means[[i]][[1L]]) <- mean_labels
      rownames(NVCOV_means[[i]][[1L]]) <-
        colnames(NVCOV_means[[i]][[1L]]) <- mean_labels

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

  } else if(isTRUE(dataList$direct_fiml)) {

    for(i in seq_len(ngroups)) {

      subfix <- if(ngroups < 2L) "" else paste0(".", dataList$group_label[i])

      means_params[[i]] <- list()
      means_params_labels[[i]] <- list()
      cov_params[[i]] <- list()
      cov_params_labels[[i]] <- list()
      VCOV_means[[i]] <- list()
      VCOV_cov[[i]] <- list()
      NVCOV_means[[i]] <- list()
      NVCOV_cov[[i]] <- list()
      M_group[[i]] <- character(0L)
      S_group[[i]] <- character(0L)
      taus_group[[i]] <- character(0L)
      nobs_ij[[i]] <- list()

      for(j in seq_along(dataList$direct_patterns[[i]])) {

        pattern <- dataList$direct_patterns[[i]][[j]]
        items <- pattern$observed_names
        p <- length(items)
        means_name <- paste0("means", subfix, ".pattern", j)
        S_name <- paste0("S", subfix, ".pattern", j)

        means_struct <- create_parameters(list(
          list(name = means_name,
               type = "matrix",
               dim = c(p, 1L),
               rownames = items,
               colnames = "intrcpt")
        ))
        cov_struct <- create_parameters(list(
          list(name = S_name,
               type = "matrix",
               dim = c(p, p),
               rownames = items,
               colnames = items,
               symmetric = TRUE)
        ))

        means_params_labels[[i]] <- c(means_params_labels[[i]], means_struct)
        cov_params_labels[[i]] <- c(cov_params_labels[[i]], cov_struct)
        means_params[[i]][[means_name]] <-
          matrix(pattern$means, ncol = 1L,
                 dimnames = list(items, "intrcpt"))
        cov_params[[i]][[S_name]] <- pattern$covariance
        dimnames(cov_params[[i]][[S_name]]) <- list(items, items)

        mean_labels <- c(means_struct[[means_name]])
        S_labels <- c(cov_struct[[S_name]][lower.tri(cov_struct[[S_name]],
                                                     diag = TRUE)])
        VCOV_means[[i]][[j]] <- matrix(
          0, nrow = length(mean_labels), ncol = length(mean_labels),
          dimnames = list(mean_labels, mean_labels)
        )
        NVCOV_means[[i]][[j]] <- VCOV_means[[i]][[j]]
        VCOV_cov[[i]][[j]] <- matrix(
          0, nrow = length(S_labels), ncol = length(S_labels),
          dimnames = list(S_labels, S_labels)
        )
        NVCOV_cov[[i]][[j]] <- VCOV_cov[[i]][[j]]

        M_group[[i]] <- c(M_group[[i]], means_name)
        S_group[[i]] <- c(S_group[[i]], S_name)
        nobs_ij[[i]][[j]] <- pattern$nobs

      }

    }

  } else if(isTRUE(dataList$fiml_moments)) {

    fit_moments <- dataList$fit_moments
    moment_param <- fit_moments@modelInfo$data_param
    moment_VCOV <- fit_moments@Optim$SE$VCOV

    for(i in seq_len(ngroups)) {

      M_name <- moment_param$means_group[i]
      S_name <- moment_param$S_group[i]
      means_params[[i]] <- fit_moments@parameters[M_name]
      means_params_labels[[i]] <- fit_moments@modelInfo$trans[M_name]
      cov_params[[i]] <- fit_moments@parameters[S_name]
      cov_params_labels[[i]] <- fit_moments@modelInfo$trans[S_name]

      mean_labels <- fit_moments@modelInfo$parameters_labels[
        fit_moments@modelInfo$parameters_labels %in%
          c(fit_moments@modelInfo$trans[[M_name]])
      ]
      S_labels <- fit_moments@modelInfo$parameters_labels[
        fit_moments@modelInfo$parameters_labels %in%
          c(fit_moments@modelInfo$trans[[S_name]])
      ]

      VCOV_means[[i]] <- list(moment_VCOV[mean_labels, mean_labels,
                                                    drop = FALSE])
      VCOV_cov[[i]] <- list(moment_VCOV[S_labels, S_labels, drop = FALSE])
      NVCOV_means[[i]] <- list(VCOV_means[[i]][[1L]]*dataList$nobs[[i]])
      NVCOV_cov[[i]] <- list(VCOV_cov[[i]][[1L]]*dataList$nobs[[i]])
      nobs_ij[[i]] <- dataList$nobs[[i]]
      M_group[[i]] <- M_name
      S_group[[i]] <- S_name
      taus_group[[i]] <- character(0L)

    }

  } else {

    for(i in seq_len(ngroups)) {

      means_params[[i]] <- dataList$fit_means[[i]]@parameters
      means_params_labels[[i]] <- dataList$fit_means[[i]]@modelInfo$trans
      VCOV_means[[i]] <- list(dataList$VCOV_means[[i]])
      NVCOV_means[[i]] <- list(
        VCOV_means[[i]][[1L]]*dataList$fit_means[[i]]@dataList$nobs
      )
      cov_params[[i]] <- dataList$fit_cov[[i]]@parameters
      cov_params_labels[[i]] <- dataList$fit_cov[[i]]@modelInfo$trans
      VCOV_cov[[i]] <- list(dataList$VCOV[[i]])
      NVCOV_cov[[i]] <- list(
        VCOV_cov[[i]][[1L]]*dataList$fit_cov[[i]]@dataList$nobs
      )
      nobs_ij[[i]] <- dataList$fit_cov[[i]]@dataList$nobs

      means_names <- names(means_params[[i]])
      cov_names <- names(cov_params[[i]])
      M_group[[i]] <- means_names[startsWith(means_names, "means")]
      S_group[[i]] <- cov_names[startsWith(cov_names, "S")]
      taus_group[[i]] <- cov_names[startsWith(cov_names, "taus")]

    }

  }

  #### Result ####

  result <- list(lambda_group = lambda_group,
                 alpha_group = alpha_group,
                 theta_group = theta_group,
                 logvars_group = logvars_group,
                 psi_group = psi_group,
                 xtheta_group = xtheta_group,
                 xpsi_group = xpsi_group,
                 model_group = model_group,
                 nu_group = nu_group,
                 meanshat_group = meanshat_group,
                 delta_group = delta_group,
                 kappa_group = kappa_group,
                 tauhat_group = tauhat_group,
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

  list_struct <- list()
  k <- 1L

  for(i in seq_len(dataList$ngroups)) {

    p <- dataList$nitems[[i]]
    q <- dataList$nfactors[[i]]
    items <- dataList$item_label[[i]]
    factors <- dataList$factor_label[[i]]

    list_struct[[k]] <- list(name = lambda_group[i],
                             type = "matrix",
                             dim = c(p, q),
                             rownames = items,
                             colnames = factors)
    k <- k+1L

    list_struct[[k]] <- list(name = alpha_group[i],
                             type = "matrix",
                             dim = c(q, 1L),
                             rownames = factors,
                             colnames = "intrcp")
    k <- k+1L

    if(dataList$positive) {

      list_struct[[k]] <- list(name = xtheta_group[i],
                               type = "matrix",
                               dim = c(p, p),
                               rownames = items,
                               colnames = items)
      k <- k+1L

      list_struct[[k]] <- list(name = xpsi_group[i],
                               type = "matrix",
                               dim = c(q, q),
                               rownames = factors,
                               colnames = factors)
      k <- k+1L

    }

    list_struct[[k]] <- list(name = theta_group[i],
                             type = "matrix",
                             dim = c(p, p),
                             rownames = items,
                             colnames = items,
                             symmetric = TRUE)
    k <- k+1L

    if(!dataList$positive && !control$deltaparam) {

      list_struct[[k]] <- list(name = logvars_group[i],
                               type = "vector",
                               dim = p,
                               rownames = items)
      k <- k+1L

    }

    list_struct[[k]] <- list(name = psi_group[i],
                             type = "matrix",
                             dim = c(q, q),
                             rownames = factors,
                             colnames = factors,
                             symmetric = TRUE)
    k <- k+1L

    list_struct[[k]] <- list(name = model_group[i],
                             type = "matrix",
                             dim = c(p, p),
                             rownames = items,
                             colnames = items,
                             symmetric = TRUE)
    k <- k+1L

    list_struct[[k]] <- list(name = nu_group[i],
                             type = "matrix",
                             dim = c(p, 1L),
                             rownames = items,
                             colnames = "intrcp")
    k <- k+1L

    list_struct[[k]] <- list(name = meanshat_group[i],
                             type = "matrix",
                             dim = c(p, 1L),
                             rownames = items,
                             colnames = "mean")
    k <- k+1L

    if(dataList$cor == "poly") {

      tau <- dataList$model[[i]]$tau

      list_struct[[k]] <- list(name = kappa_group[i],
                               type = "matrix",
                               dim = dim(tau),
                               rownames = rownames(tau),
                               colnames = colnames(tau))
      k <- k+1L

      list_struct[[k]] <- list(name = tauhat_group[i],
                               type = "matrix",
                               dim = dim(tau),
                               rownames = rownames(tau),
                               colnames = colnames(tau))
      k <- k+1L

    }

    if(control$deltaparam) {

      list_struct[[k]] <- list(name = delta_group[i],
                               type = "matrix",
                               dim = c(p, 1L),
                               rownames = items,
                               colnames = "latent.var")
      k <- k+1L

    }

  }

  trans <- create_parameters(list_struct)
  trans <- c(trans, unlist(means_params_labels, recursive = FALSE))
  trans <- c(trans, unlist(cov_params_labels, recursive = FALSE))

  #### Result ####

  return(trans)

}

constraints_lcfa <- function(dataList, data_param, trans, control) {

  list2env(data_param, envir = environment())

  fixed <- vector("list")
  nonfixed <- vector("list")
  fixed_values_list <- vector("list")
  target_psi <- vector("list", length = dataList$ngroups)
  target_theta <- vector("list", length = dataList$ngroups)

  #### Replace latent labels by lavaan labels ####

  for(i in seq_len(dataList$ngroups)) {

    group_names <- c(lambda_group[i], theta_group[i], psi_group[i],
                     nu_group[i], alpha_group[i])
    model_blocks <- list(dataList$model[[i]]$lambda,
                         dataList$model[[i]]$theta,
                         dataList$model[[i]]$psi,
                         dataList$model[[i]]$nu,
                         dataList$model[[i]]$alpha)

    if(dataList$cor == "poly") {
      group_names <- c(group_names, kappa_group[i])
      model_blocks <- c(model_blocks, list(dataList$model[[i]]$tau))
    }

    names(model_blocks) <- group_names

    nonfixed[group_names] <- lapply(model_blocks, FUN = \(x) {
      which(is.na(suppressWarnings(as.numeric(x))))
    })
    fixed[group_names] <- lapply(model_blocks, FUN = \(x) {
      which(!is.na(suppressWarnings(as.numeric(x))))
    })
    fixed_values_list[group_names] <- lapply(model_blocks, FUN = \(x) {
      numerals <- suppressWarnings(as.numeric(x))
      numerals[!is.na(numerals)]
    })

    for(nm in group_names) {
      trans[[nm]][nonfixed[[nm]]] <-
        model_blocks[[nm]][nonfixed[[nm]]]
    }

    if(!dataList$positive && !control$deltaparam) {

      p <- dataList$nitems[[i]]
      diag_indices <- seq.int(1L, p*p, by = p+1L)
      free_diag <- diag_indices %in% nonfixed[[theta_group[i]]]
      theta_diag_labels <- diag(trans[[theta_group[i]]])

      trans[[logvars_group[i]]][free_diag] <-
        paste0("log(", theta_diag_labels[free_diag], ")")

    }

  }

  #### Model for the parameters ####

  param <- list()

  for(i in seq_len(dataList$ngroups)) {

    p <- dataList$nitems[[i]]
    q <- dataList$nfactors[[i]]
    fiml_missing <- isTRUE(dataList$direct_fiml) ||
      isTRUE(dataList$fiml_moments)

    param[[lambda_group[i]]] <- trans[[lambda_group[i]]]
    param[[lambda_group[i]]][fixed[[lambda_group[i]]]] <-
      fixed_values_list[[lambda_group[i]]]

    if(dataList$positive) {

      param[[xtheta_group[i]]] <- trans[[xtheta_group[i]]]
      param[[xpsi_group[i]]] <- trans[[xpsi_group[i]]]

    } else {

      param[[theta_group[i]]] <- trans[[theta_group[i]]]
      param[[theta_group[i]]][fixed[[theta_group[i]]]] <-
        fixed_values_list[[theta_group[i]]]

      if(control$deltaparam) {

        diag(param[[theta_group[i]]]) <- "1"

      } else {

        diag_indices <- seq.int(1L, p*p, by = p+1L)
        free_diag <- diag_indices %in% nonfixed[[theta_group[i]]]
        theta_diag <- diag(param[[theta_group[i]]])
        theta_diag[free_diag] <- "1"
        diag(param[[theta_group[i]]]) <- theta_diag

        param[[logvars_group[i]]] <- setNames(
          rep(0, p),
          names(trans[[logvars_group[i]]])
        )
        param[[logvars_group[i]]][free_diag] <-
          trans[[logvars_group[i]]][free_diag]

      }

      param[[psi_group[i]]] <- trans[[psi_group[i]]]
      param[[psi_group[i]]][fixed[[psi_group[i]]]] <-
        fixed_values_list[[psi_group[i]]]

    }

    #### Sample covariance and thresholds ####

    if(control$free_S) {

      param[S_group[[i]]] <- cov_params[[i]][S_group[[i]]]

      for(j in seq_along(S_group[[i]])) {

        S_name <- S_group[[i]][j]
        statistic_labels <- rownames(NVCOV_cov[[i]][[j]])

        if(is.null(statistic_labels)) {
          stop("The sample covariance matrix must have parameter labels before it can be freed.")
        }

        free_statistics <- trans[[S_name]] %in% statistic_labels
        param[[S_name]][free_statistics] <-
          trans[[S_name]][free_statistics]

      }

    } else {
      param[S_group[[i]]] <- cov_params[[i]][S_group[[i]]]
    }

    if(control$free_taus) {
      param[taus_group[[i]]] <- trans[taus_group[[i]]]
    } else {
      param[taus_group[[i]]] <- cov_params[[i]][taus_group[[i]]]
    }

    if(dataList$cor == "poly") {
      param[S_group[[i]]] <- lapply(param[S_group[[i]]], FUN = \(x) {
        diag(x) <- 1
        return(x)
      })
    }

    #### Sample means ####

    if(control$free_M &&
       control$meanstructure &&
       dataList$cor == "pearson") {
      param[M_group[[i]]] <- trans[M_group[[i]]]
    } else {
      param[M_group[[i]]] <- means_params[[i]][M_group[[i]]]
    }

    #### Observed intercepts and latent means ####

    if(!control$meanstructure && !fiml_missing) {

      if(dataList$cor == "pearson") {
        sample_mean_name <- M_group[[i]][1L]
        param[[nu_group[i]]] <- means_params[[i]][[sample_mean_name]]
      } else {
        param[[nu_group[i]]] <- matrix(
          0, nrow = p, ncol = 1L,
          dimnames = dimnames(trans[[nu_group[i]]])
        )
      }

      param[[alpha_group[i]]] <- matrix(
        0, nrow = q, ncol = 1L,
        dimnames = dimnames(trans[[alpha_group[i]]])
      )

    } else {

      param[[nu_group[i]]] <- trans[[nu_group[i]]]
      param[[nu_group[i]]][fixed[[nu_group[i]]]] <-
        fixed_values_list[[nu_group[i]]]

      param[[alpha_group[i]]] <- trans[[alpha_group[i]]]
      param[[alpha_group[i]]][fixed[[alpha_group[i]]]] <-
        fixed_values_list[[alpha_group[i]]]

    }

    #### Unstandardized ordinal thresholds ####

    if(dataList$cor == "poly") {
      param[[kappa_group[i]]] <- trans[[kappa_group[i]]]
      param[[kappa_group[i]]][fixed[[kappa_group[i]]]] <-
        fixed_values_list[[kappa_group[i]]]
    }

    #### Delta parameterization ####

    if(control$deltaparam) {

      if(control$std.ov && !fiml_missing) {
        # param[[delta_group[i]]] <- matrix(
        #   1, nrow = p, ncol = 1L,
        #   dimnames = dimnames(trans[[delta_group[i]]])
        # )
        param[[delta_group[i]]] <- matrix(diag(param[[S_group[[i]]]]),
                                          nrow = p, ncol = 1L,
                                          dimnames = dimnames(trans[[delta_group[i]]]))
      } else {
        param[[delta_group[i]]] <- trans[[delta_group[i]]]
      }

    }

    #### Positive-definite constraint targets ####

    if(dataList$positive) {

      target_theta[[i]] <- matrix(
        0, nrow = p, ncol = p,
        dimnames = dimnames(trans[[theta_group[i]]])
      )
      target_theta[[i]][nonfixed[[theta_group[i]]]] <- 1

      target_psi[[i]] <- matrix(
        0, nrow = q, ncol = q,
        dimnames = dimnames(trans[[psi_group[i]]])
      )
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

  init_param <- vector("list", length = control$rstarts)

  for(rs in seq_len(control$rstarts)) {

    init_param[[rs]] <- list()

    for(i in seq_len(dataList$ngroups)) {

      p <- dataList$nitems[[i]]
      q <- dataList$nfactors[[i]]

      #### Factor loadings ####

      # init_param[[rs]][[lambda_group[i]]] <- rorth(p, q)
      # dimnames(init_param[[rs]][[lambda_group[i]]]) <-
      #   dimnames(trans[[lambda_group[i]]])
      # init_param[[rs]][[lambda_group[i]]][fixed[[lambda_group[i]]]] <-
      #   fixed_values_list[[lambda_group[i]]]

      S <- cov_params[[i]][[S_group[[i]][1L]]]
      smc <- diag(S)-1/diag(approx_Hinv(S))
      S_reduced <- S
      diag(S_reduced) <- smc
      eig <- eigen(S_reduced, symmetric = TRUE)
      values <- pmax(eig$values[seq_len(q)], .Machine$double.eps)
      lambda <- eig$vectors[, seq_len(q), drop = FALSE]%*%
        diag(sqrt(values), nrow = q)
      lambda <- lambda + matrix(rnorm(p*q, sd = 0.01), nrow = p, ncol = q)

      init_param[[rs]][[lambda_group[i]]] <- lambda
      dimnames(init_param[[rs]][[lambda_group[i]]]) <-
        dimnames(trans[[lambda_group[i]]])
      init_param[[rs]][[lambda_group[i]]][fixed[[lambda_group[i]]]] <-
        fixed_values_list[[lambda_group[i]]]

      #### Covariance matrices ####

      if(dataList$positive) {

        init_param[[rs]][[xtheta_group[i]]] <- rorth(p, p)
        dimnames(init_param[[rs]][[xtheta_group[i]]]) <-
          dimnames(trans[[xtheta_group[i]]])

        init_param[[rs]][[xpsi_group[i]]] <- rorth(q, q)
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

        lambda <- init_param[[rs]][[lambda_group[i]]]
        init_param[[rs]][[theta_group[i]]] <- diag(diag(S) - rowSums(lambda*lambda))
        dimnames(init_param[[rs]][[theta_group[i]]]) <-
          dimnames(trans[[theta_group[i]]])
        init_param[[rs]][[theta_group[i]]][fixed[[theta_group[i]]]] <-
          fixed_values_list[[theta_group[i]]]

        if(!control$deltaparam) {

          diag_indices <- seq.int(1L, p*p, by = p+1L)
          free_diag <- !(diag_indices %in% fixed[[theta_group[i]]])
          theta_diag <- diag(init_param[[rs]][[theta_group[i]]])
          minimum_variance <- pmax(0.05*diag(S),
                                   sqrt(.Machine$double.eps))
          theta_diag[free_diag] <- pmax(theta_diag[free_diag],
                                        minimum_variance[free_diag])
          diag(init_param[[rs]][[theta_group[i]]]) <- theta_diag

          init_param[[rs]][[logvars_group[i]]] <- setNames(
            rep(0, p),
            names(trans[[logvars_group[i]]])
          )
          init_param[[rs]][[logvars_group[i]]][free_diag] <-
            log(theta_diag[free_diag])

        }

        init_param[[rs]][[psi_group[i]]] <- diag(q)
        dimnames(init_param[[rs]][[psi_group[i]]]) <-
          dimnames(trans[[psi_group[i]]])
        init_param[[rs]][[psi_group[i]]][fixed[[psi_group[i]]]] <-
          fixed_values_list[[psi_group[i]]]

      }

      #### Observed intercepts and latent means ####

      if(dataList$cor == "pearson") {

        if(isTRUE(dataList$sample_stats_only)) {
          init_means <- dataList$sample.mean[[i]][dataList$item_label[[i]]]
        } else {
          init_means <- colMeans(dataList$data_per_group[[i]], na.rm = TRUE)
        }

      } else {

        init_means <- rep(0, p)

      }

      init_param[[rs]][[nu_group[i]]] <-
        matrix(init_means, ncol = 1L,
               dimnames = dimnames(trans[[nu_group[i]]]))
      init_param[[rs]][[nu_group[i]]][fixed[[nu_group[i]]]] <-
        fixed_values_list[[nu_group[i]]]

      init_param[[rs]][[alpha_group[i]]] <-
        matrix(0, nrow = q, ncol = 1L,
               dimnames = dimnames(trans[[alpha_group[i]]]))
      init_param[[rs]][[alpha_group[i]]][fixed[[alpha_group[i]]]] <-
        fixed_values_list[[alpha_group[i]]]

      #### Unstandardized ordinal thresholds ####

      if(dataList$cor == "poly") {

        initial_kappas <- lcfa_order_thresholds(
          threshold_blocks = cov_params[[i]][taus_group[[i]]],
          model_tau = dataList$model[[i]]$tau
        )

        init_param[[rs]][[kappa_group[i]]] <-
          matrix(initial_kappas,
                 nrow = nrow(trans[[kappa_group[i]]]),
                 ncol = ncol(trans[[kappa_group[i]]]),
                 dimnames = dimnames(trans[[kappa_group[i]]]))
        init_param[[rs]][[kappa_group[i]]][fixed[[kappa_group[i]]]] <-
          fixed_values_list[[kappa_group[i]]]

      }

      #### Optional free sample statistics ####

      for(nm in intersect(S_group[[i]], names(param))) {
        init_param[[rs]][[nm]] <- cov_params[[i]][[nm]]
      }

      for(nm in intersect(M_group[[i]], names(param))) {
        init_param[[rs]][[nm]] <- means_params[[i]][[nm]]
      }

      for(nm in intersect(taus_group[[i]], names(param))) {
        init_param[[rs]][[nm]] <- cov_params[[i]][[nm]]
      }

      #### Delta parameterization ####

      if(control$deltaparam) {
        init_param[[rs]][[delta_group[i]]] <-
          matrix(1, nrow = p, ncol = 1L,
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
    add_euclidean(alpha_group[i])
    add_euclidean(kappa_group[i])
    add_euclidean(delta_group[i])
    add_euclidean(logvars_group[i])
    add_euclidean(M_group[[i]])
    add_euclidean(S_group[[i]])
    add_euclidean(taus_group[[i]])

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

    if(control$deltaparam) {

      transforms[[k]] <- list(transform = "deltaparam",
                              parameters_in = c(delta_group[i],
                                                lambda_group[i],
                                                psi_group[i]),
                              parameters_out = list(diag(trans[[theta_group[i]]])),
                              extra = list(p = nrow(trans[[theta_group[i]]]),
                                           q = nrow(trans[[psi_group[i]]])))
      k <- k+1L

    } else if(!dataList$positive) {

      p <- dataList$nitems[[i]]
      diag_indices <- seq.int(1L, p*p, by = p+1L)
      theta_numeric <- suppressWarnings(
        as.numeric(dataList$model[[i]]$theta)
      )
      free_diag <- is.na(theta_numeric[diag_indices])

      if(any(free_diag)) {

        transforms[[k]] <- list(
          transform = "exponential",
          parameters_in = list(trans[[logvars_group[i]]][free_diag]),
          parameters_out = list(diag(trans[[theta_group[i]]])[free_diag])
        )
        k <- k+1L

      }

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

    #### Model-implied observed means ####

    transforms[[k]] <- list(transform = "meanstructure",
                            parameters_in = c(nu_group[i],
                                              lambda_group[i],
                                              alpha_group[i]),
                            parameters_out = meanshat_group[i],
                            extra = list(p = dataList$nitems[[i]],
                                         q = dataList$nfactors[[i]]))
    k <- k+1L

    #### Model-implied standardized ordinal thresholds ####

    if(dataList$cor == "poly") {

      threshold_items <- lcfa_threshold_items(
        model_tau = dataList$model[[i]]$tau,
        item_names = dataList$item_label[[i]]
      )-1L

      transforms[[k]] <- list(transform = "tau_param",
                              parameters_in = c(kappa_group[i],
                                                meanshat_group[i],
                                                model_group[i]),
                              parameters_out = tauhat_group[i],
                              extra = list(p = dataList$nitems[[i]],
                                           threshold_items = threshold_items))
      k <- k+1L

    }

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
  total_nobs <- sum(unlist(dataList$nobs))

  #### CFA discrepancy functions ####

  for(i in seq_len(dataList$ngroups)) {

    if(dataList$cor == "poly") {

      cfa_estimator <- switch(tolower(dataList$estimator),
                              uls = "cfa_dwls_poly",
                              dwls = "cfa_dwls_poly",
                              stop("Ordinal CFA currently supports uls or dwls."))

    } else {

      cfa_estimator <- switch(tolower(dataList$estimator),
                              uls = "cfa_dwls",
                              dwls = "cfa_dwls",
                              ml = "cfa_fml",
                              fml = "cfa_fml",
                              stop("Unknown estimator: ", dataList$estimator))

    }

    cov_params_i <- cov_params[[i]]
    cov_params_i <- cov_params_i[startsWith(names(cov_params_i), "S")]

    for(j in seq_along(cov_params_i)) {

      pick <- rownames(cov_params_i[[j]])
      S_group_ij <- S_group[[i]][[j]]
      M_group_ij <- M_group[[i]][[j]]
      p <- nrow(cov_params_i[[j]])

      #### Covariance weights ####

      if(control$estimator %in% c("uls", "ml", "fml")) {

        W_cov <- matrix(1, nrow = p, ncol = p)

      } else {

        sample_S_labels <- c(
          trans[[S_group_ij]][lower.tri(trans[[S_group_ij]],
                                        diag = !control$std.ov)]
        )
        covariance_weights <- lcfa_diagonal_weights(
          NVCOV = NVCOV_cov[[i]][[j]],
          labels = sample_S_labels,
          object_name = "sample covariance/correlation statistic"
        )

        W_cov <- matrix(0, nrow = p, ncol = p)
        W_cov[lower.tri(W_cov, diag = !control$std.ov)] <-
          covariance_weights
        W_cov[upper.tri(W_cov)] <- t(W_cov)[upper.tri(W_cov)]

      }

      model_covariance <- c(trans[[model_group[i]]][pick, pick])
      sample_covariance <- c(trans[[S_group_ij]])
      contribution_weight <- nobs_ij[[i]][[j]]/total_nobs

      if(dataList$cor == "poly") {

        model_thresholds <- c(trans[[tauhat_group[i]]])
        sample_thresholds <- lcfa_order_thresholds(
          threshold_blocks = trans[taus_group[[i]]],
          model_tau = dataList$model[[i]]$tau
        )

        if(length(model_thresholds) != length(sample_thresholds)) {
          stop("Model-implied and sample thresholds have different lengths.")
        }

        if(control$estimator == "uls") {

          w_thresholds <- rep(1, length(sample_thresholds))

        } else {

          w_thresholds <- lcfa_diagonal_weights(
            NVCOV = NVCOV_cov[[i]][[j]],
            labels = sample_thresholds,
            object_name = "sample threshold statistic"
          )

        }

        estimators[[k]] <- list(
          estimator = cfa_estimator,
          parameters = list(model_covariance,
                            sample_covariance,
                            model_thresholds,
                            sample_thresholds),
          extra = list(W = W_cov,
                       w = contribution_weight,
                       w_thresholds = w_thresholds,
                       q = nrow(trans[[psi_group[i]]]),
                       p = p,
                       n = nobs_ij[[i]][[j]],
                       group_index = i,
                       group_label = dataList$group_label[i],
                       pattern_index = j,
                       role = "cfa",
                       covariance_name = S_group_ij,
                       double_names = S_group_ij,
                       matrix_names = S_group_ij)
        )

      } else {

        model_means <- c(trans[[meanshat_group[i]]][pick, , drop = FALSE])
        sample_means <- c(trans[[M_group_ij]])

        if(control$estimator == "uls") {

          w_means <- rep(1, p)

        } else if(control$estimator == "dwls") {

          w_means <- lcfa_diagonal_weights(
            NVCOV = NVCOV_means[[i]][[j]],
            labels = sample_means,
            zero_allowed = TRUE,
            object_name = "sample mean statistic"
          )

        } else {

          w_means <- rep(1, p)

        }

        estimators[[k]] <- list(
          estimator = cfa_estimator,
          parameters = list(model_covariance,
                            sample_covariance,
                            model_means,
                            sample_means),
          extra = list(W = W_cov,
                       w = contribution_weight,
                       w_means = w_means,
                       q = nrow(trans[[psi_group[i]]]),
                       p = p,
                       n = nobs_ij[[i]][[j]],
                       group_index = i,
                       group_label = dataList$group_label[i],
                       pattern_index = j,
                       role = "cfa",
                       covariance_name = S_group_ij,
                       means_name = M_group_ij,
                       double_names = S_group_ij,
                       matrix_names = c(S_group_ij, M_group_ij))
        )

      }

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
                                           group_index = i,
                                           group_label = dataList$group_label[i],
                                           pattern_index = NA_integer_,
                                           role = "penalty",
                                           double_names = "logdetR psi"))
      k <- k+1L

      lower_indices <- which(lower.tri(trans[[theta_group[i]]], diag = TRUE))
      estimators[[k]] <- list(estimator = "logdetR",
                              parameters = theta_group[i],
                              extra = list(lower_indices = lower_indices-1L,
                                           p = nrow(trans[[theta_group[i]]]),
                                           logdetw = control$penalties$logdet$w,
                                           group_index = i,
                                           group_label = dataList$group_label[i],
                                           pattern_index = NA_integer_,
                                           role = "penalty",
                                           double_names = "logdetR theta"))
      k <- k+1L

    }

  }

  control_estimator <- create_estimators(estimators = estimators,
                                         structures = trans)

  #### Result ####

  return(control_estimator)

}
