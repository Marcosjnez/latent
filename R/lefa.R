# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 23/08/2026
#'
#' Exploratory Factor Analysis
#'
#' Fit an exploratory factor analysis model by first estimating an orthogonal
#' factor model with \code{lcfa()} and subsequently rotating the fitted model
#' with \code{lrotate()}.
#'
#' @usage
#' lefa(data = NULL, nfactors = 1L, estimator = "ml",
#'      projection = "oblq", rotation = "oblimin",
#'      model = NULL, ordered = FALSE, group = NULL,
#'      sample.cov = NULL, sample.mean = NULL, sample.nobs = NULL,
#'      positive = FALSE, penalties = TRUE,
#'      missing = "pairwise.complete.obs",
#'      std.lv = TRUE, std.ov = FALSE,
#'      meanstructure = TRUE,
#'      parameterization = NULL,
#'      likelihood = NULL, se = TRUE,
#'      message = FALSE, do.fit = TRUE,
#'      mimic = "latent", control = NULL,
#'      rotation.control = NULL, ...)
#'
#' @param data Optional data frame or matrix containing the observed variables.
#'   Alternatively, sample.cov can be supplied.
#' @param nfactors Integer. Number of factors used when \code{model = NULL}.
#' @param estimator Estimation method passed to \code{lcfa()}.
#' @param projection Rotation projection passed to \code{lrotate()}. Available
#'   options are \code{"orth"}, \code{"oblq"}, and \code{"poblq"}.
#' @param rotation Rotation criterion passed to \code{lrotate()}.
#' @param model Optional lavaan model syntax. If \code{NULL}, an exploratory
#'   lower-diagonal loading model is generated automatically.
#' @param ordered Logical value indicating whether indicators are ordinal. The
#'   character value \code{"yule"} requests Yule correlations.
#' @param group Optional character string identifying the grouping variable.
#' @param sample.cov Optional sample covariance matrix or list of covariance
#'   matrices passed to \code{lcfa()}.
#' @param sample.mean Optional sample mean vector or list of vectors passed to
#'   \code{lcfa()}.
#' @param sample.nobs Optional number of observations passed to \code{lcfa()}.
#' @param positive Logical. Request the positive-definite parameterization used
#'   by \code{lcfa()}.
#' @param penalties Logical value or list controlling regularization in
#'   \code{lcfa()}.
#' @param missing Missing-data method passed to \code{lcfa()}.
#' @param std.lv Logical. Standardize latent variables in the unrotated model.
#' @param std.ov Logical. Standardize observed variables in the unrotated model.
#' @param meanstructure Logical. Estimate the observed-variable mean structure.
#' @param parameterization Optional parameterization passed to \code{lcfa()}.
#' @param likelihood Character string controlling the likelihood convention in
#'   \code{lcfa()}.
#' @param se Logical or character controlling standard errors in \code{lcfa()}.
#' @param message Logical. Print progress messages during CFA estimation.
#' @param do.fit Logical. If \code{FALSE}, return the unrotated \code{lcfa}
#'   model specification without fitting or rotation.
#' @param mimic Retained for backward compatibility. Only \code{"latent"} is
#'   currently supported.
#' @param control Optional list of optimization controls passed to
#'   \code{lcfa()}.
#' @param rotation.control Optional list of optimization controls passed to
#'   \code{lrotate()}.
#' @param ... Additional arguments. CFA/lavaan arguments are passed to
#'   \code{lcfa()}; arguments required by the selected rotation criterion or
#'   projection are passed only to \code{lrotate()}.
#'
#' @return A fitted object of class \code{"lefa"}. The unrotated
#'   \code{lcfa} object is stored in its \code{extra} slot. If
#'   \code{do.fit = FALSE}, the unfitted \code{lcfa} specification is returned.
#'
#' @examples
#' \dontrun{
#' fit <- lefa(data = HolzingerSwineford1939,
#'             nfactors = 3L,
#'             rotation = "oblimin")
#' }
#'
#' @export
lefa <- function(data = NULL, nfactors = 1L, estimator = "ml",
                 projection = "oblq", rotation = "oblimin",
                 model = NULL, ordered = FALSE, group = NULL,
                 sample.cov = NULL, sample.mean = NULL, sample.nobs = NULL,
                 positive = FALSE, penalties = TRUE,
                 missing = "pairwise.complete.obs",
                 std.lv = TRUE, std.ov = FALSE,
                 meanstructure = TRUE,
                 parameterization = NULL,
                 likelihood = NULL, se = TRUE,
                 message = FALSE, do.fit = TRUE,
                 mimic = "latent", control = NULL,
                 rotation.control = NULL,
                 ...) {

  #### Check input arguments ####

  if(is.null(data) && is.null(sample.cov)) {
    stop("Either data or sample.cov must be provided")
  }

  if(!is.null(data) && !is.data.frame(data) && !is.matrix(data)) {
    stop("data must be NULL, a data.frame, or a matrix")
  }

  if(!is.null(control) && !is.list(control)) {
    stop("control must be NULL or a list")
  }

  if(!is.null(rotation.control) && !is.list(rotation.control)) {
    stop("rotation.control must be NULL or a list")
  }

  if(!is.null(model) &&
     (!is.character(model) || length(model) < 1L)) {
    stop("model must be NULL or lavaan model syntax")
  }

  if(!is.null(group) &&
     (!is.character(group) || length(group) != 1L || is.na(group))) {
    stop("group must be NULL or the name of one grouping variable")
  }

  if(length(positive) != 1L ||
     !is.logical(positive) ||
     is.na(positive)) {
    stop("positive must be TRUE or FALSE")
  }

  if(length(std.lv) != 1L ||
     !is.logical(std.lv) ||
     is.na(std.lv)) {
    stop("std.lv must be TRUE or FALSE")
  }

  if(length(std.ov) != 1L ||
     !is.logical(std.ov) ||
     is.na(std.ov)) {
    stop("std.ov must be TRUE or FALSE")
  }

  if(length(meanstructure) != 1L ||
     !is.logical(meanstructure) ||
     is.na(meanstructure)) {
    stop("meanstructure must be TRUE or FALSE")
  }

  if(length(message) != 1L ||
     !is.logical(message) ||
     is.na(message)) {
    stop("message must be TRUE or FALSE")
  }

  if(length(do.fit) != 1L ||
     !is.logical(do.fit) ||
     is.na(do.fit)) {
    stop("do.fit must be TRUE or FALSE")
  }

  estimator <- tolower(estimator)
  projection <- tolower(projection)
  rotation <- tolower(rotation)
  missing <- tolower(missing)

  supported_projection <- c("orth", "oblq", "poblq")
  supported_rotation <- c("cf", "geomin", "lclf", "oblimin",
                          "target", "varimax", "varimin", "xtarget")

  if(!(projection %in% supported_projection)) {
    stop("Unknown projection: ", projection)
  }

  if(!(rotation %in% supported_rotation)) {
    stop("Unknown rotation criterion: ", rotation)
  }

  if(length(mimic) != 1L ||
     !is.character(mimic) ||
     is.na(mimic) ||
     tolower(mimic) != "latent") {
    stop("Only mimic = 'latent' is currently supported")
  }

  if(is.null(model)) {

    if(length(nfactors) != 1L ||
       !is.numeric(nfactors) ||
       is.na(nfactors) ||
       nfactors < 1 ||
       nfactors != as.integer(nfactors)) {
      stop("nfactors must be a positive integer")
    }

    nfactors <- as.integer(nfactors)

  }

  #### Store original call ####

  mc <- match.call()

  #### Additional arguments ####

  dots <- list(...)

  if(length(dots) > 0L) {

    dot_names <- names(dots)

    if(is.null(dot_names) ||
       any(is.na(dot_names)) ||
       any(dot_names == "")) {
      stop("All arguments supplied through ... must be named")
    }

    if(anyDuplicated(dot_names)) {
      stop("Arguments supplied through ... must have unique names")
    }

  }

  dots_split <- split_dots_lefa(dots = dots,
                                projection = projection,
                                rotation = rotation)

  #### Create the EFA model ####

  efa_model <- create_lefa_model(data = data,
                                 sample.cov = sample.cov,
                                 nfactors = nfactors,
                                 model = model,
                                 group = group)

  #### Fit the orthogonal factor model ####

  efa_fit <- fit_lefa_cfa(data = data,
                          model = efa_model,
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
                          meanstructure = meanstructure,
                          parameterization = parameterization,
                          likelihood = likelihood,
                          se = se,
                          message = message,
                          do.fit = do.fit,
                          control = control,
                          dots = dots_split$cfa)

  #### Return the unfitted model ####

  if(!do.fit) {

    #### Result ####

    return(efa_fit)

  }

  #### Rotate the factor solution ####

  rotation_fit <- fit_lefa_rotation(fit = efa_fit,
                                    projection = projection,
                                    rotation = rotation,
                                    control = rotation.control,
                                    dots = dots_split$rotation)

  result <- as_lefa(rotation_fit, call = mc)

  #### Result ####

  return(result)

}

#### Function to create an lefa object ####

as_lefa <- function(fit, call) {

  fit@dataList$lefa_call <- call

  result <- new("lefa",
                version          = fit@version,
                call             = fit@call,
                timing           = fit@timing,
                dataList         = fit@dataList,
                modelInfo        = fit@modelInfo,
                Optim            = fit@Optim,
                parameters       = fit@parameters,
                transformed_pars = fit@transformed_pars,
                extra            = fit@extra)

  #### Result ####

  return(result)

}

#### Function to split CFA and rotation arguments ####

split_dots_lefa <- function(dots, projection, rotation) {

  rotation_names <- switch(
    rotation,
    cf      = "k",
    geomin  = "epsilon",
    lclf    = "epsilon",
    oblimin = "gamma",
    target  = c("target", "weight"),
    varimax = character(0L),
    varimin = character(0L),
    xtarget = c("target", "weight", "w", "psitarget", "psiweight"),
    character(0L)
  )

  if(projection == "poblq") {
    rotation_names <- c(rotation_names, "constraints")
  }

  rotation_names <- unique(rotation_names)

  rotation_dots <- dots[intersect(names(dots), rotation_names)]

  cfa_names <- setdiff(names(dots),
                       c(rotation_names, "orthogonal"))

  cfa_dots <- dots[cfa_names]

  result <- list(
    cfa = cfa_dots,
    rotation = rotation_dots
  )

  #### Result ####

  return(result)

}

#### Function to create the exploratory factor model ####

create_lefa_model <- function(data, sample.cov, nfactors, model, group) {

  if(!is.null(model)) {

    #### Result ####

    return(model)

  }

  if(is.null(data)) {

    if(is.list(sample.cov)) {
      model_data <- sample.cov[[1L]]
    } else {
      model_data <- sample.cov
    }

  } else {

    model_data <- data

  }

  if(!is.null(group) && !is.null(data)) {

    if(is.null(colnames(model_data)) ||
       !(group %in% colnames(model_data))) {
      stop("The grouping variable is not present in data")
    }

    model_columns <- setdiff(colnames(model_data), group)

    if(length(model_columns) == 0L) {
      stop("No observed variables remain after removing the grouping variable")
    }

    model_data <- model_data[, model_columns, drop = FALSE]

  }

  model <- make_lowerdiag_lavaan(data = model_data,
                                 nfactors = nfactors)

  #### Result ####

  return(model)

}

#### Function to fit the orthogonal CFA representation ####

fit_lefa_cfa <- function(data, model, estimator,
                         ordered, group,
                         sample.cov, sample.mean, sample.nobs,
                         positive, penalties,
                         missing, std.lv, std.ov,
                         meanstructure, parameterization,
                         likelihood, se, message,
                         do.fit, control,
                         dots) {

  args <- list(
    data = data,
    model = model,
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
    meanstructure = meanstructure,
    parameterization = parameterization,
    likelihood = likelihood,
    se = se,
    message = message,
    do.fit = do.fit,
    control = control,
    orthogonal = TRUE
  )

  args <- c(args, dots)

  result <- do.call(lcfa, args)

  #### Result ####

  return(result)

}

#### Function to determine whether rotation is required ####

rotate_lefa <- function(fit) {

  nfactors <- unlist(fit@dataList$nfactors,
                     use.names = FALSE)

  result <- any(nfactors > 1L)

  #### Result ####

  return(result)

}

#### Function to fit the rotation ####

fit_lefa_rotation <- function(fit, projection, rotation,
                              control, dots) {

  args <- list(
    fit = fit,
    projection = projection,
    rotation = rotation,
    do.fit = TRUE,
    control = control
  )

  args <- c(args, dots)

  result <- do.call(lrotate, args)

  #### Result ####

  return(result)

}
