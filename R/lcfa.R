# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 26/04/2026
#'
#' @title
#' Fit a Confirmatory Factor Analysis (CFA) model with lavaan syntax.
#'
#' @usage
#'
#' lcfa(data, model = NULL, estimator = "ml",
#' ordered = FALSE, group = NULL,
#' sample.cov = NULL, nobs = NULL,
#' positive = FALSE, penalties = TRUE,
#' missing = "pairwise.complete.obs",
#' std.lv = FALSE, do.fit = TRUE,
#' message = FALSE, mimic = 'latent',
#' control = NULL, ...)
#'
#' @param data data frame or matrix.
#' @param model lavaan's model syntax.
#' @param estimator Available estimators: "ml", "uls", and "dwls". Defaults to "ml".
#' @param ordered Logical. Defaults to TRUE.
#' @param group String. Name of the variable that splits the data in different groups.
#' @param sample.cov Covariance matrix between the items. Defaults to NULL.
#' @param nobs Number of observations. Defaults to NULL.
#' @param positive Force a positive-definite solution. Defaults to FALSE.
#' @param penalties list of penalty terms for the parameters.
#' @param missing Method to handle missing data.
#' @param std.lv Logical. Provide the parameters of the standardized model. Default is TRUE.
#' @param std.ov Logical. Standardize the observed variables before fitting. Default is FALSE.
#' @param acov String. "standard" or "robust". Default is "standard".
#' @param meanstructure Logical. Estimate the means of the variables. Default is FALSE.
#' @param do.fit TRUE to fit the model and FALSE to return only the model setup. Defaults to TRUE.
#' @param message Logical. Defaults to FALSE.
#' @param se Logical. Compute standard errors. Defaults to TRUE.
#' @param likelihood String. Use N (normal) or N-1 (wishart) in the denominator. Defaults to "normal" for ML and "wishart" otherwise.
#' @param control List of control parameters for the optimization algorithm. See 'details' for more information.
#' @param ... Additional lavaan arguments. See ?lavaan for more information.
#'
#' @details \code{lcfa} estimates confirmatory factor models.
#'
#' @return List with the following objects:
#' \item{version}{Version number of 'latent' when the model was estimated.}
#' \item{call}{Code used to estimate the model.}
#' \item{ModelInfo}{Model information.}
#' \item{Optim}{Output of the optimizer.}
#' \item{parameters}{Structure with all model parameters.}
#' \item{transparameters}{Structure with all transformed model parameters.}
#' \item{loglik}{Logarithm likelihood of the model.}
#' \item{penalized_loglik}{Logarithm likelihood + logarithm priors of the model.}
#'
#' @examples
#'
#' \dontrun{
#' # The famous Holzinger and Swineford (1939) example
#' HS.model <- ' visual  =~ x1 + x2 + x3
#'               textual =~ x4 + x5 + x6
#'               speed   =~ x7 + x8 + x9 '
#'
#' fit <- lcfa(model = HS.model, data = HolzingerSwineford1939)
#' summary(fit, digits = 3L)
#'}
#'
#' @export
lcfa <- function(data, model = NULL, estimator = "ml",
                 ordered = FALSE, group = NULL,
                 sample.cov = NULL, nobs = NULL,
                 positive = FALSE, penalties = FALSE,
                 missing = "pairwise.complete.obs",
                 std.lv = FALSE, std.ov = FALSE,
                 acov = "standard", meanstructure = TRUE,
                 parameterization = NULL,
                 likelihood = NULL, se = TRUE,
                 control = NULL, message = FALSE,
                 do.fit = TRUE,
                 ...) {

  ## store original call
  mc  <- match.call()

  estimator <- tolower(estimator)
  missing <- tolower(missing)

  if(isTRUE(ordered)) {

    cor <- "poly"
    std.ov <- TRUE

    # if(parameterization == "theta") {
    #   control$deltaparam <- FALSE
    # } else if(parameterization == "delta") {
    #   control$deltaparam <- TRUE
    # }

    if(positive) {
      control$deltaparam <- FALSE
      std.lv <- TRUE
    } else {
      control$deltaparam <- TRUE
    }

  } else if(ordered == "yule") {
    cor <- "yule"
    std.ov <- TRUE
    control$deltaparam <- TRUE
  } else {
    cor <- "pearson"
  }

  if(estimator == "ml" || is.null(likelihood)) {
    likelihood <- "normal"
  } else {
    likelihood <- "wishart"
  }

  if(missing == "fiml") {
    meanstructure <- TRUE
    std.ov <- FALSE
  }

  if(meanstructure) {
    if(estimator == "ml" || estimator == "fml") estimator <- "means_fml"
    if(estimator == "uls") estimator <- "means_uls"
    if(estimator == "dwls") estimator <- "means_dwls"
  }

  control$ordered <- ordered
  control$std.lv <- std.lv
  control$std.ov <- std.ov
  control$positive <- positive
  control$penalties <- penalties
  control$estimator <- tolower(estimator)
  control$meanstructure <- meanstructure
  control$missing <- missing
  control <- lcfa_control(control)

  args <- as.list(match.call(expand.dots = TRUE))[-1]

  #### Create the datalist ####

  dataList <- create_cfa_datalist(
    data = data,
    model = model,
    cor = cor,
    estimator = estimator,
    ordered = ordered,
    group = group,
    sample.cov = sample.cov,
    nobs = nobs,
    positive = positive,
    penalties = penalties,
    missing = missing,
    std.lv = std.lv,
    std.ov = std.ov,
    acov = acov,
    message = message,
    likelihood = likelihood,
    meanstructure = meanstructure,
    args = args,
    control = control,
    ...
  )

  #### Create the model ####

  full_model <- create_cfa_model(dataList = dataList,
                                 model = model,
                                 control = control)
  list2env(full_model, envir = environment())

  #### Create the modelInfo ####

  modelInfo <- create_cfa_modelInfo(dataList = dataList,
                                    full_model = full_model,
                                    control = control)

  #### Fit the model ####

  if(!do.fit) {

    lcfa_list <- new("latent",
                     version            = as.character( packageVersion('latent') ),
                     call               = mc,
                     timing             = numeric(),
                     dataList           = dataList,
                     modelInfo          = modelInfo,
                     Optim              = list(),
                     parameters         = list(),
                     transformed_pars   = list(),
                     loglik             = numeric(),
                     penalized_loglik   = numeric(),
                     loss               = numeric(),
                     penalized_loss     = numeric()
    )

    return(lcfa_list)

  }

  if(message) {
    msg <- "Fitting the model"
    w <- nchar(msg) + 4
    cat("\n", "+", strrep("-", w), "+\n",
        "|  ", msg, "  |\n",
        "+", strrep("-", w), "+\n\n", sep = "")
  }

  modelInfo$control_optimizer$cores <- min(modelInfo$control_optimizer$rstarts,
                                           modelInfo$control_optimizer$cores)
  # Fit the model:
  Optim <- optimizer(control_manifold = modelInfo$control_manifold,
                     control_transform = modelInfo$control_transform,
                     control_estimator = modelInfo$control_estimator,
                     control_optimizer = modelInfo$control_optimizer)
  names(Optim$parameters) <- modelInfo$parameters_labels
  names(Optim$transparameters) <- modelInfo$transparameters_labels

  # Collect all the information about the optimization:

  elapsed <- Optim$elapsed

  #### Estimated model structures ####

  # Create the structures of transformed parameters:
  transformed_pars <- fill_in(modelInfo$trans, Optim$transparameters)

  # Create the structures of untransformed parameters:
  parameters <- transformed_pars[names(modelInfo$param)]

  #### Process the fit information ####

  loss <- sum(unlist(lapply(Optim$outputs$estimators$doubles,
                            FUN = \(x) x[[1]])))
  loglik <- sum(unlist(lapply(Optim$outputs$estimators$doubles,
                              FUN = \(x) x[[2]])))
  penalty <- sum(unlist(lapply(Optim$outputs$estimators$doubles,
                               FUN = \(x) x[[5]])))
  penalized_loss <- loss + penalty
  penalized_loglik <- loglik + penalty

  #### latent object ####

  result <- new("lcfa",
                version            = as.character( packageVersion('latent') ),
                call               = mc,
                timing             = elapsed,
                dataList           = dataList,
                modelInfo          = modelInfo,
                Optim              = Optim,
                parameters         = parameters,
                transformed_pars   = transformed_pars,
                loglik             = loglik,
                penalized_loglik   = penalized_loglik,
                loss               = loss,
                penalized_loss     = penalized_loss
  )

  if(message) {
    msg <- "Computing standard errors"
    w <- nchar(msg) + 4
    cat("\n", "+", strrep("-", w), "+\n",
        "|  ", msg, "  |\n",
        "+", strrep("-", w), "+\n\n", sep = "")
  }

  #### Standard errors ####

  if(isTRUE(se)) {
    Optim$SE <- se(result, type = "standard", digits = 9)
  }

  result@Optim <- Optim

  # # Fit by group:
  # fit_by_group <- latInspect(result, what = "fit")
  # loss.group <- unlist(lapply(fit_by_group, FUN = \(x) x["loss"]))
  # penalized_loss.group <- unlist(lapply(fit_by_group, FUN = \(x) x["penalized_loss"]))
  # loglik.group <- unlist(lapply(fit_by_group, FUN = \(x) x["loglik"]))
  # penalized_loglik.group <- unlist(lapply(fit_by_group, FUN = \(x) x["penalized_loglik"]))
  # penalty.group <- unlist(lapply(fit_by_group, FUN = \(x) x["penalty"]))

  #### Return ####

  return(result)

}

create_cfa_datalist <- function(data, model = NULL, cor = "pearson",
                                estimator = "ml", ordered = FALSE,
                                group = NULL, sample.cov = NULL, nobs = NULL,
                                positive = FALSE, penalties = TRUE,
                                missing = "pairwise.complete.obs",
                                std.lv = TRUE, std.ov = FALSE,
                                acov = "standard", message = FALSE,
                                likelihood = NULL, meanstructure = TRUE,
                                args = NULL, control = NULL,
                                ...) {

  cor <- tolower(cor)
  estimator <- tolower(estimator)
  acov <- tolower(acov)
  missing <- tolower(missing)

  if(is.null(group)) {
    ngroups <- 1L
    group <- "group"
    group_label <- ""
    data$group <- group_label
  } else {
    group_label <- unique(data[[group]])
    ngroups <- length(group_label)
  }

  item_names <- extract_item_names_lavaan(model, ngroups = ngroups)

  # Remove group cases with all missing data for given item_names:
  keep <- rep(FALSE, nrow(data))
  for(i in seq_len(ngroups)) {
    group_i <- data[[group]] == group_label[i]
    items_i <- item_names[[i]]
    not_all_na_i <- !apply(is.na(data[group_i, items_i, drop = FALSE]), 1, all)
    keep[group_i] <- not_all_na_i
  }
  data <- data[keep, , drop = FALSE]

  get_group_data <- function(i) {
    if (ngroups > 1L) {
      data[data[[group]] == group_label[i], item_names[[i]], drop = FALSE]
    } else {
      data[, item_names[[i]], drop = FALSE]
    }
  }

  unwrap_single <- function(x) {
    if (length(x) == 0L || all(vapply(x, is.null, logical(1)))) {
      return(NULL)
    }
    if (ngroups == 1L) x[[1L]] else x
  }

  normalize_to_list <- function(x, ngroups) {
    if (is.null(x)) {
      return(vector("list", ngroups))
    }
    if (ngroups == 1L && !is.list(x)) {
      return(list(x))
    }
    x
  }

  sample.cov <- normalize_to_list(sample.cov, ngroups)

  if (!is.null(nobs) && ngroups == 1L && !is.list(nobs)) {
    nobs_in <- list(nobs)
  } else {
    nobs_in <- nobs
  }

  nobs_list <- vector("list", length = ngroups)
  X <- vector("list", length = ngroups)
  NACOV <- vector("list", length = ngroups)
  ACOV <- vector("list", length = ngroups)
  WLS.V <- vector("list", length = ngroups)
  thresholds <- vector("list", length = ngroups)
  fit_cov <- vector("list", length = ngroups)
  fit_means <- vector("list", length = ngroups)

  for (i in seq_len(ngroups)) {

    X[[i]] <- get_group_data(i)
    if(ngroups < 2) {
      control$subfix <- ""
    } else {
      control$subfix <- group_label[i]
    }

    # Item means estimation:
    fit_means[[i]] <- lmean(data = X[[i]],
                            std.ov = std.ov,
                            control = control,
                            do.fit = TRUE)

    # Covariance matrix estimation:
    if(cor == "pearson") {
      fit_cov[[i]] <- lpearson(data = X[[i]],
                               std.ov = std.ov,
                               acov = acov,
                               likelihood = likelihood,
                               missing = missing,
                               control = control,
                               do.fit = TRUE)
      if(missing == "fiml") {
        patterns <- split_by_missing_pattern(X[[i]])
        npatterns <- length(patterns)
        fit_means[[i]]@extra <- vector("list", length = npatterns)
        fit_cov[[i]]@extra <- vector("list", length = npatterns)
        subfix <- control$subfix
        for(j in seq_len(npatterns)) { # Run by missing data pattern
          control$subfix <- paste(subfix, ".pattern", j, sep = "")
          fit_means[[i]]@extra[[j]] <- lmean(data = patterns[[j]]$data,
                                             std.ov = std.ov,
                                             do.fit = TRUE,
                                             control = control,
                                             ...)
          fit_cov[[i]]@extra[[j]] <- lpearson(data = patterns[[j]]$data,
                                              model = model,
                                              std.ov = std.ov,
                                              acov = acov,
                                              likelihood = likelihood,
                                              missing = "pairwise.complete.obs",
                                              do.fit = TRUE,
                                              control = control,
                                              ...)
        }
      }
    } else if(cor == "poly") {
      fit_cov[[i]] <- lpoly(data = X[[i]],
                            method = "two-step",
                            control = control,
                            do.fit = TRUE)
    } else if(cor == "yule") {
      fit_cov[[i]] <- lyule(data = X[[i]],
                            control = control,
                            do.fit = TRUE)
    } else {
      stop("Unknown correlation type")
    }

    nobs_list[[i]] <- fit_cov[[i]]@dataList$nobs
    sample.cov[[i]] <- fit_cov[[i]]@transformed_pars$S
    NACOV[[i]] <- fit_cov[[i]]@Optim$SE$ACOV * fit_cov[[i]]@dataList$nobs
    ACOV[[i]] <- fit_cov[[i]]@Optim$SE$ACOV
    WLS.V[[i]] <- diag(ACOV[[i]])
    idx_taus <- startsWith(names(fit_cov[[i]]@transformed_pars), "taus")
    thresholds[[i]] <- fit_cov[[i]]@transformed_pars[idx_taus]

  }

  LAV <- lavaan::cfa(
    model = model,
    sample.cov = unwrap_single(sample.cov),
    sample.nobs = unwrap_single(nobs_list),
    group = group,
    NACOV = unwrap_single(NACOV),
    WLS.V = unwrap_single(WLS.V),
    ordered = ordered,
    std.lv = std.lv,
    std.ov = std.ov,
    meanstructure = meanstructure,
    do.fit = FALSE,
    warn = FALSE,
    ...
  )

  LAV@Options$positive <- positive

  lavmodel <- LAV@Model
  item_label <- LAV@Data@ov.names
  nobs_list <- LAV@Data@nobs
  factor_label <- replicate(ngroups, list(LAV@Model@dimNames[[1]][[2]]))

  model_out <- getmodel_fromlavaan(LAV)
  if (ngroups == 1L) {
    model_out <- list(model_out)
  }

  nitems <- as.list(lavmodel@nvar)
  npatterns <- lapply(nitems, function(p) 0.5 * p * (p + 1))
  nfactors <- lapply(model_out, function(x) ncol(x$lambda))

  if (is.null(args)) {
    args <- as.list(match.call(expand.dots = TRUE))[-1]
  }

  dataList <- list()
  dataList$ngroups <- ngroups
  dataList$data <- data
  dataList$data_per_group <- X
  dataList$nobs <- nobs_list
  dataList$nitems <- nitems
  dataList$npatterns <- npatterns
  dataList$nfactors <- nfactors
  dataList$positive <- positive
  dataList$estimator <- estimator
  dataList$cor <- cor
  dataList$group_label <- group_label
  dataList$item_label <- item_label
  dataList$factor_label <- factor_label
  dataList$LAV <- LAV
  dataList$args <- args
  dataList$model <- model_out
  dataList$sample.cov <- sample.cov
  dataList$NACOV <- NACOV
  dataList$WLS.V <- WLS.V
  dataList$thresholds <- thresholds
  dataList$fit_means <- fit_means
  dataList$fit_cov <- fit_cov

  return(dataList)

}

create_cfa_model <- function(dataList, model, control) {

  # Generate the model syntax and initial parameter values

  list2env(dataList, envir = environment())

  # Initialize the objects to store the initial parameters:
  trans <- vector("list")
  fixed <- nonfixed <- fixed_values_list <- vector("list")

  #### Model for the transformed parameters ####

  # Initialize the target matrices for positive-definite constraints:
  target_psi <- target_theta <- targets <- vector("list", length = ngroups)
  rest <- 0L

  lambda_group <- paste("lambda.", dataList$group_label, sep = "")
  psi_group <- paste("psi.", dataList$group_label, sep = "")
  theta_group <- paste("theta.", dataList$group_label, sep = "")
  xpsi_group <- paste("xpsi.", dataList$group_label, sep = "")
  xtheta_group <- paste("xtheta.", dataList$group_label, sep = "")
  model_group <- paste("model.", dataList$group_label, sep = "")
  nu_group <- paste("nu.", dataList$group_label, sep = "")
  delta_group <- paste("delta.", dataList$group_label, sep = "")
  tau_group <- paste("tau.", dataList$group_label, sep = "")
  S_group <- taus_group <- M_group <- vector("list", length = ngroups)
  means_params <- means_params_labels <- vector("list", length = ngroups)
  cov_params <- cov_params_labels <- vector("list", length = ngroups)
  acov_means <- vector("list", length = ngroups)
  acov_cov <- vector("list", length = ngroups)
  nobs_ij <- vector("list", length = ngroups)

  for(i in 1:ngroups) {

    if(control$missing == "fiml") {
      means_params[[i]] <- unlist(lapply(fit_means[[i]]@extra,
                                         FUN = \(x) x@parameters),
                                  recursive = FALSE)
      means_params_labels[[i]] <- unlist(lapply(fit_means[[i]]@extra,
                                                FUN = \(x) x@modelInfo$trans),
                                         recursive = FALSE)
      acov_means[[i]] <- lapply(fit_means[[i]]@extra,
                                  FUN = \(x) x@Optim$SE$ACOV)
      cov_params[[i]] <- unlist(lapply(fit_cov[[i]]@extra,
                                       FUN = \(x) x@parameters),
                                recursive = FALSE)
      cov_params_labels[[i]] <- unlist(lapply(fit_cov[[i]]@extra,
                                              FUN = \(x) x@modelInfo$trans),
                                       recursive = FALSE)
      acov_cov[[i]] <- lapply(fit_cov[[i]]@extra,
                                FUN = \(x) x@Optim$SE$ACOV)
      nobs_ij[[i]] <- lapply(fit_cov[[i]]@extra,
                                    FUN = \(x) x@dataList$nobs)
    } else {
      means_params[[i]] <- fit_means[[i]]@parameters
      means_params_labels[[i]] <- fit_means[[i]]@modelInfo$trans
      acov_means[[i]] <- list(fit_means[[i]]@Optim$SE$ACOV)
      cov_params[[i]] <- fit_cov[[i]]@parameters
      cov_params_labels[[i]] <- fit_cov[[i]]@modelInfo$trans
      acov_cov[[i]] <- list(fit_cov[[i]]@Optim$SE$ACOV)
      nobs_ij[[i]] <- fit_cov[[i]]@dataList$nobs
    }

    means_params_names <- names(means_params[[i]])
    M_group[[i]] <- means_params_names[startsWith(means_params_names, "means.")]
    cov_params_names <- names(cov_params[[i]])
    S_group[[i]] <- cov_params_names[startsWith(cov_params_names, "S.")]
    taus_group[[i]] <- cov_params_names[startsWith(cov_params_names, "taus.")]

  }

  data_param <- list(lambda_group = lambda_group,
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
                     cov_params = cov_params,
                     acov_means = acov_means,
                     acov_cov = acov_cov,
                     nobs_ij = nobs_ij)

  # Transformed parameters:
  list_struct <- vector("list")
  k <- 1L
  for(i in 1:ngroups) {

    # Lambda:
    list_struct[[k]] <- list(name = lambda_group[i],
                             type = "matrix",
                             dim = c(nitems[[i]], nfactors[[i]]),
                             rownames = item_label[[i]],
                             colnames = factor_label[[i]])
    k <- k+1L

    # Create additional parameters if there are positive-definite constraints:
    if(positive) {

      # Theta:
      list_struct[[k]] <- list(name = xtheta_group[i],
                               type = "matrix",
                               dim = c(nitems[[i]], nitems[[i]]),
                               rownames = item_label[[i]],
                               colnames = item_label[[i]])
      k <- k+1L

      # Psi:
      list_struct[[k]] <- list(name = xpsi_group[i],
                               type = "matrix",
                               dim = c(nfactors[[i]], nfactors[[i]]),
                               rownames = factor_label[[i]],
                               colnames = factor_label[[i]])
      k <- k+1L

    }

    # Theta:
    list_struct[[k]] <- list(name = theta_group[i],
                             type = "matrix",
                             dim = c(nitems[[i]], nitems[[i]]),
                             rownames = item_label[[i]],
                             colnames = item_label[[i]],
                             symmetric = TRUE)
    k <- k+1L

    # Psi:
    list_struct[[k]] <- list(name = psi_group[i],
                             type = "matrix",
                             dim = c(nfactors[[i]], nfactors[[i]]),
                             rownames = factor_label[[i]],
                             colnames = factor_label[[i]],
                             symmetric = TRUE)
    k <- k+1L

    # Model matrix:
    list_struct[[k]] <- list(name = model_group[i],
                             type = "matrix",
                             dim = c(nitems[[i]], nitems[[i]]),
                             rownames = item_label[[i]],
                             colnames = item_label[[i]],
                             symmetric = TRUE)
    k <- k+1L

    # Model means vector:
    if(control$meanstructure) {
      list_struct[[k]] <- list(name = nu_group[i],
                               type = "matrix",
                               dim = c(nitems[[i]], 1),
                               rownames = item_label[[i]],
                               colnames = "intrcp")
      k <- k+1L
    }

    # latent variances:
    if(control$deltaparam) {
      list_struct[[k]] <- list(name = delta_group[i],
                               type = "matrix",
                               dim = c(nitems[[i]], 1L),
                               rownames = item_label[[i]],
                               colnames = "latent.var")
      k <- k+1L
    }

  }

  trans <- create_parameters(list_struct)
  if(control$meanstructure) {
    trans <- c(trans, unlist(means_params_labels, recursive = FALSE))
  }
  trans <- c(trans, unlist(cov_params_labels, recursive = FALSE))

  #### Replace latent labels by lavaan labels ####

  for(i in 1:ngroups) {

    # Get the positions of parameters and fixed values:

    # Ensure the same label ordering than in model:
    group_i <- c(lambda_group[i], theta_group[i], psi_group[i], nu_group[i])

    nonfixed[group_i] <- lapply(model[[i]][1:4], FUN = \(x) {
      which(is.na(suppressWarnings(as.numeric(x))))
    })

    fixed[group_i] <- lapply(model[[i]][1:4], FUN = \(x) {
      which(!is.na(suppressWarnings(as.numeric(x))))
    })

    fixed_values_list[group_i] <- lapply(model[[i]][1:4], FUN = \(x) {
      numerals <- suppressWarnings(as.numeric(x))
      inds <- which(!is.na(numerals))
      return(numerals[inds])
    })

    trans[[lambda_group[i]]][nonfixed[[lambda_group[i]]]] <- model[[i]]$lambda[nonfixed[[lambda_group[i]]]]
    trans[[theta_group[i]]][nonfixed[[theta_group[i]]]] <- model[[i]]$theta[nonfixed[[theta_group[i]]]]
    trans[[psi_group[i]]][nonfixed[[psi_group[i]]]] <- model[[i]]$psi[nonfixed[[psi_group[i]]]]
    if(control$meanstructure) {
      trans[[nu_group[i]]][nonfixed[[nu_group[i]]]] <- model[[i]]$nu[nonfixed[[nu_group[i]]]]
    }

  }

  #### Model for the parameters ####

  param <- list()
  item_label <- dataList$item_label

  for(i in 1:ngroups) {

    param[[lambda_group[i]]] <- trans[[lambda_group[i]]]
    # Insert fixed values in the model:
    param[[lambda_group[i]]][fixed[[lambda_group[i]]]] <-
      model[[i]]$lambda[fixed[[lambda_group[i]]]]

    if(positive) {

      # Theta:
      param[[xtheta_group[i]]] <- trans[[xtheta_group[i]]]
      # Psi:
      param[[xpsi_group[i]]] <- trans[[xpsi_group[i]]]

    } else {

      # Theta:
      param[[theta_group[i]]] <- trans[[theta_group[i]]]
      # Insert fixed values in the model:
      param[[theta_group[i]]][fixed[[theta_group[i]]]] <-
        model[[i]]$theta[fixed[[theta_group[i]]]]

      if(control$deltaparam) {
        diag(param[[theta_group[i]]]) <- "1"
      }

      # Psi:
      param[[psi_group[i]]] <- trans[[psi_group[i]]]
      # Insert fixed values in the model:
      param[[psi_group[i]]][fixed[[psi_group[i]]]] <-
        model[[i]]$psi[fixed[[psi_group[i]]]]

    }

    # Fix the sample covariance matrix:
    if(control$free_S) {
      param[unlist(S_group[[i]])] <- trans[unlist(S_group[[i]])]
    } else {
      param[S_group[[i]]] <- cov_params[[i]][S_group[[i]]]
    }

    if(cor %in% c("poly", "polys", "polychoric", "polychorics")) {
      fix_diag <- TRUE
    } else {
      fix_diag <- FALSE
    }

    # Fix the diagonal of the sample covariance matrix:
    if(fix_diag) {
      param[S_group[[i]]] <- lapply(param[S_group[[i]]],
                                    FUN = \(x) {
                                      diag(x) <- 1; return(x)
                                    })
    }

    param[taus_group[[i]]] <- cov_params[[i]][taus_group[[i]]]

    if(control$meanstructure) {

      if(control$free_M) {
        param[M_group[[i]]] <- trans[M_group[[i]]]
      } else {
        param[M_group[[i]]] <- means_params[[i]][M_group[[i]]]
      }

      if(control$std.ov) {
        param[[nu_group[[i]]]] <- matrix(rep(0, nitems[[i]]), ncol = 1L,
                                         dimnames = list(item_label[[i]],
                                                         "intrcp"))
      } else {
        param[[nu_group[i]]] <- trans[[nu_group[i]]]
      }

    }

    # # Insert fixed values in the model:
    # param[[nu_group[i]]][fixed[[nu_group[i]]]] <- model[[i]]$nu[fixed[[nu_group[i]]]]

    if(control$deltaparam) {

      if(control$std.lv) {
        param[[delta_group[i]]] <- trans[[delta_group[i]]]
      } else {
        param[[delta_group[i]]] <- matrix(rep(1, nitems[[i]]), ncol = 1L,
                                          dimnames = list(item_label[[i]],
                                                          "latent.vars"))
      }

    }

    # Create the target matrices for positive-definite constraints:
    if(positive) {

      target_theta[[i]] <- matrix(0, nrow = nitems[[i]], ncol = nitems[[i]])
      target_theta[[i]][nonfixed[[theta_group[i]]]] <- 1

      target_psi[[i]] <- matrix(0, nrow = nfactors[[i]], ncol = nfactors[[i]])
      target_psi[[i]][nonfixed[[psi_group[i]]]] <- 1

      # UGLY FIX THIS TO COUNT THE DEGREES OF FREEDOM AUTOMATICALLY
      q <- nfactors[[i]]
      p <- nitems[[i]]
      lower_theta <- lower.tri(diag(p), diag = TRUE)
      lower_psi <- lower.tri(diag(q), diag = TRUE)
      nconstraints <- sum(unlist(c(target_theta[[i]][lower_theta],
                                   target_psi[[i]][lower_psi])) == 0)
      rest <- rest + 0.5*q*(q-1) + 0.5*p*(p-1) + nconstraints

    }

  }

  # param <- param[names(trans) %in% names(param)]

  #### Fixed parameters ####

  # # Replace the parameters by custom values, if available:
  # if(!is.null(model)) {
  #
  #   # Replace the parameters by custom values:
  #
  #   nm <- intersect(names(model), names(param))
  #   nm <- nm[!vapply(model[nm], is.null, logical(1))]
  #   param[nm] <- model[nm]
  #
  # }

  #### Create the initial values for the parameters ####

  # Collect the unique nontransformed parameters and the unique transformed parameters:

  init_param <- vector("list", length = control$rstarts)

  for(rs in 1:control$rstarts) {

    init_param[[rs]] <- vector("list")

    for(i in 1:ngroups) {

      init_param[[rs]][[lambda_group[i]]] <- rorth(nitems[[i]], nfactors[[i]])
      init_param[[rs]][[lambda_group[i]]][fixed[[lambda_group[i]]]] <- fixed_values_list[[lambda_group[i]]]

      if(positive) {

        # init_param[[rs]][[xtheta_group[i]]] <- rpoblq(nitems[[i]], nitems[[i]], constraints = target_theta[[i]])
        # init_param[[rs]][[xpsi_group[i]]] <- rpoblq(nfactors[[i]], nfactors[[i]], constraints = target_psi[[i]])
        init_param[[rs]][[xtheta_group[i]]] <- rorth(nitems[[i]], nitems[[i]])
        init_param[[rs]][[xpsi_group[i]]] <- rorth(nfactors[[i]], nfactors[[i]])

        init_param[[rs]][[theta_group[i]]] <- crossprod(init_param[[rs]][[xtheta_group[i]]])
        init_param[[rs]][[psi_group[i]]] <- crossprod(init_param[[rs]][[xpsi_group[i]]])

      } else {

        init_param[[rs]][[theta_group[i]]] <- diag(runif(nitems[[i]]))
        init_param[[rs]][[theta_group[i]]][fixed[[theta_group[i]]]] <- fixed_values_list[[theta_group[i]]]

        init_param[[rs]][[psi_group[i]]] <- diag(nfactors[[i]])
        init_param[[rs]][[psi_group[i]]][fixed[[psi_group[i]]]] <- fixed_values_list[[psi_group[i]]]

      }

      Lambda <- init_param[[rs]][[lambda_group[i]]]
      Theta <- init_param[[rs]][[theta_group[i]]]
      Psi <- init_param[[rs]][[psi_group[i]]]
      # init_param[[rs]][[model_group[i]]] <- Lambda %*% Psi %*% t(Lambda) + Theta

      if(control$meanstructure) {
        init_param[[rs]][[nu_group[i]]] <- matrix(colMeans(dataList$data_per_group[[i]],
                                                           na.rm = TRUE), ncol = 1L)
      }

      # init_param[[rs]][S_group[[i]]] <- cov_params[[i]][S_group[[i]]]
      # if(control$meanstructure) {
      #   init_param[[rs]][M_group[[i]]] <- means_params[[i]][M_group[[i]]]
      # }

      if(control$deltaparam) {

        init_param[[rs]][[delta_group[i]]] <- matrix(1, nitems[[i]], ncol = 1L)
        rownames(init_param[[rs]][[delta_group[i]]]) <- item_label[[i]]

      }

    }

  }

  #### Custom initial values ####

  # Replace initial starting values by custom starting values:

  if(!is.null(control$start)) {

    nm <- names(control$start)
    nm <- nm[!vapply(control$start, is.null, logical(1))]

    for (i in seq_len(control$rstarts)) {
      common_nm <- intersect(nm, names(init_param[[i]]))
      for (j in common_nm) {
        init_param[[i]][[j]] <- insert_object(init_param[[i]][[j]],
                                              control$start[[j]])
      }
    }

  }

  #### Return ####

  result <- list(param = param,
                 trans = trans,
                 init_param = init_param,
                 target_psi = target_psi,
                 target_theta = target_theta,
                 rest = rest,
                 data_param = data_param)

  return(result)

}

create_cfa_modelInfo <- function(dataList, full_model, control) {

  # Generate control_manifold, control_transform, and control_estimator

  list2env(dataList, envir = environment())
  list2env(full_model, envir = environment())
  list2env(data_param, envir = environment())

  #### Manifolds ####

  manifolds <- list()
  k <- 1L

  for(i in 1:ngroups) {

    manifolds[[k]] <- list(manifold = "euclidean",
                           parameters = lambda_group[i])
    k <- k+1L

    manifolds[[k]] <- list(manifold = "euclidean",
                           parameters = nu_group[i])
    k <- k+1L

    manifolds[[k]] <- list(manifold = "euclidean",
                           parameters = delta_group[i])
    k <- k+1L

    manifolds[[k]] <- list(manifold = "euclidean",
                           parameters = M_group[[i]])
    k <- k+1L

    manifolds[[k]] <- list(manifold = "euclidean",
                           parameters = S_group[[i]])
    k <- k+1L

    if(positive) {

      manifolds[[k]] <- list(manifold = "poblq",
                             parameters = xpsi_group[i],
                             extra = list(
                               p = nfactors[[i]],
                               q = nfactors[[i]],
                               constraints = target_psi[[i]]))
      k <- k+1L

      manifolds[[k]] <- list(manifold = "poblq",
                             parameters = xtheta_group[i],
                             extra = list(
                               p = nitems[[i]],
                               q = nitems[[i]],
                               constraints = target_theta[[i]]))
      k <- k+1L

    } else {

      manifolds[[k]] <- list(manifold = "euclidean",
                             parameters = psi_group[i])
      k <- k+1L

      manifolds[[k]] <- list(manifold = "euclidean",
                             parameters = theta_group[i])
      k <- k+1L

    }

  }

  control_manifold <- create_manifolds(manifolds = manifolds,
                                       structures = param)

  #### Transformations ####

  transforms <- list()
  dots <- list()
  k <- 1L

  for(i in 1:ngroups) {

    if(positive) {

      lower_psi <- lower.tri(trans[[psi_group[i]]], diag = TRUE)
      lower_theta <- lower.tri(trans[[theta_group[i]]], diag = TRUE)

      dots$p <- nrow(trans[[psi_group[i]]])
      transforms[[k]] <- list(transform = "crossprod",
                              parameters_in = xpsi_group[i],
                              parameters_out = psi_group[i],
                              extra = dots)
      k <- k+1L

      dots$p <- nrow(trans[[theta_group[i]]])
      transforms[[k]] <- list(transform = "crossprod",
                              parameters_in = xtheta_group[i],
                              parameters_out = theta_group[i],
                              extra = dots)
      k <- k+1L

    }

    if(control$deltaparam) {

      dots$p <- nrow(trans[[theta_group[i]]])
      dots$q <- nrow(trans[[psi_group[i]]])
      transforms[[k]] <- list(transform = "deltaparam",
                              parameters_in = c(delta_group[i],
                                                lambda_group[i],
                                                psi_group[i]),
                              parameters_out = list(diag(trans[[theta_group[i]]])),
                              extra = dots)
      k <- k+1L

    }

    # Model matrix correlation:
    lower_psi <- lower.tri(trans[[psi_group[i]]], diag = TRUE)
    lower_theta <- lower.tri(trans[[theta_group[i]]], diag = TRUE)
    lower_diag <- lower.tri(trans[[model_group[i]]], diag = TRUE)

    dots$p <- nitems[[i]]
    dots$q <- nfactors[[i]]
    transforms[[k]] <- list(transform = "factor_cor",
                            parameters_in = c(lambda_group[i],
                                              psi_group[i],
                                              theta_group[i]),
                            parameters_out = model_group[i],
                            extra = dots)
    k <- k+1L

  }

  control_transform <- create_transforms(transforms = transforms,
                                         structures = trans)

  #### Estimators ####

  estimators <- list()
  k <- 1L

  for(i in 1:ngroups) {

    estimator <- tolower(estimator)
    cfa_estimator <- switch(estimator,
                            uls = "cfa_dwls",
                            dwls  = "cfa_dwls",
                            ml = "cfa_fml",
                            fml  = "cfa_fml",
                            means_fml = "cfa_means_fml",
                            means_dwls = "cfa_means_dwls",
                            means_uls = "cfa_means_dwls",
                            stop("Unknown estimator: ", estimator)
    )

    idx <- startsWith(names(cov_params[[i]]), "S.")
    cov_params[[i]] <- cov_params[[i]][idx]
    for(j in seq_len(length(cov_params[[i]]))) {

      pick <- rownames(cov_params[[i]][[j]])
      S_group_ij <- S_group[[i]][[j]]
      M_group_ij <- M_group[[i]][[j]]

      p <- nrow(cov_params[[i]][[j]])
      if(control$estimator %in% c("uls", "means_uls", "ml", "fml", "means_fml")) {
        W_cov <- matrix(1, nrow = p, ncol = p)
      } else {
        idx <- startsWith(rownames(acov_cov[[i]][[j]]), "S.")
        W_cov <- matrix(NA_real_, nrow = p, ncol = p)
        W_cov[lower.tri(W_cov, diag = !control$std.ov)] <-
          diag(acov_cov[[i]][[j]][idx, idx]) / nobs_ij[[i]][[j]]
        W_cov[upper.tri(W_cov)] <- t(W_cov)[upper.tri(W_cov)]
        W_cov <- 1 / W_cov
        # if(control$std.ov) {
        #   diag(W_cov) <- 1
        # } else {
        #   diag(W_cov) <- 0
        # }
        diag(W_cov) <- 100000
        # if(control$std.ov) diag(W_cov) <- 0
      }
      w_means <- diag(acov_means[[i]][[j]])

      estimators[[k]] <- list(estimator = cfa_estimator,
                              parameters = list(c(trans[[model_group[i]]][pick, pick]),
                                                c(trans[[S_group_ij]]),
                                                c(trans[[nu_group[i]]][pick, ]),
                                                c(trans[[M_group_ij]])),
                              extra = list(W = W_cov,
                                           w = nobs_ij[[i]][[j]] /
                                             sum(unlist(dataList$nobs)),
                                           w_means = w_means,
                                           q = nrow(trans[[psi_group[i]]]),
                                           p = p,
                                           n = nobs_ij[[i]][[j]]))

      k <- k + 1L

    }

  }

  if(control$reg) {

    for(i in 1:ngroups) {

      # For the psi matrix:

      lower_indices <- which(lower.tri(trans[[psi_group[i]]], diag = TRUE))
      logdetw <- control$penalties$logdet$w #* nfactors[[i]]
      estimators[[k]] <- list(estimator = "logdetR",
                              parameters = psi_group[i],
                              extra = list(lower_indices = lower_indices-1L,
                                           p = nrow(trans[[psi_group[i]]]),
                                           logdetw = logdetw))
      k <- k+1L

      # For the theta matrix:

      lower_indices <- which(lower.tri(trans[[theta_group[i]]], diag = TRUE))
      logdetw <- control$penalties$logdet$w * nitems[[i]]
      estimators[[k]] <- list(estimator = "logdetR",
                              parameters = theta_group[i],
                              extra = list(lower_indices = lower_indices-1L,
                                           p = nrow(trans[[theta_group[i]]]),
                                           logdetw = control$penalties$logdet$w))
      k <- k+1L

    }

  }

  control_estimator <- create_estimators(estimators = estimators,
                                         structures = trans)

  #### Pass the initial values to vectors ####

  idx_transformed <- unlist(lapply(control_transform,
                                   FUN = \(x) unlist(x$indices_out)+1L))
  inits <- create_init(trans, param, init_param,
                       idx_transformed = idx_transformed, control)

  parameters <- inits$parameters
  parameters_labels <- names(parameters[[1]])
  nparam <- length(parameters_labels)

  transparameters <- inits$transparameters
  transparameters_labels <- names(transparameters[[1]])
  ntrans <- length(transparameters_labels)

  trans2param <- match(parameters_labels, transparameters_labels)

  #### Set up the optimizer ####

  # Create defaults for the control of the optimizer:
  control_optimizer <- control
  control_optimizer$parameters <- parameters
  control_optimizer$transparameters <- transparameters
  control_optimizer$init_param <- init_param
  control_optimizer$transparam2param <- trans2param-1L

  #### Collect all the model information ####

  modelInfo <- list(param = param,
                    trans = trans,
                    nparam = nparam - rest,
                    ntrans = ntrans,
                    parameters_labels = parameters_labels,
                    transparameters_labels = transparameters_labels,
                    dof = sum(unlist(npatterns)) - nparam + rest,
                    control_manifold = control_manifold,
                    control_transform = control_transform,
                    control_estimator = control_estimator,
                    control_optimizer = control_optimizer)

  #### Return ####

  return(modelInfo)

}

lcfa_control <- function(control) {

  # Auxiliary function for lcfa.R

  # Control input

  if(is.null(control$opt)) {
    if(control$positive) {
      control$opt <- "grad"
      if(is.null(control$rstarts)) {
        control$rstarts <- 10L
      }
    } else {
      control$opt <- "lbfgs"
    }
  }

  if(!control$positive) {

    control$penalties <- FALSE

  }

  if(isFALSE(control$penalties)) {

    control$reg <- FALSE

  } else if(isTRUE(control$penalties)) {

    control$reg <- TRUE

    control$penalties <- list(
      logdet = list(w = 1e-03)
    )

  } else if(is.list(control$penalties)) {

    if(control$penalties$logdet$w <= 0) {
      stop("The penalty w must be positive")
    }
    control$reg <- TRUE

  } else {

    stop("penalties should be TRUE, FALSE, or a list")

  }

  if(is.null(control$free_S)) {
    control$free_S <- FALSE
  }

  if(is.null(control$free_M)) {
    control$free_M <- FALSE
  }

  if(is.null(control$deltaparam)) {
    control$deltaparam <- FALSE
  }

  if(is.null(control$step_maxit)) {
    control$step_maxit <- 30L
  }

  if(is.null(control$c1)) {
    control$c1 <- 0.5
  }

  if(is.null(control$c2)) {
    control$c2 <- 0.5
  }

  if(is.null(control$step_eps)) {
    control$step_eps <- 1e-09
  }

  if(is.null(control$df_eps)) {
    control$df_eps <- 1e-09
  }

  if(is.null(control$M)) {
    control$M <- 100L
  }

  if(is.null(control$eps)) {
    control$eps <- 1e-06
  }

  if(is.null(control$ss_fac)) {
    control$ss_fac <- 2
  }

  if(is.null(control$maxit)) {
    control$maxit <- 1000L
  }

  if(is.null(control$rstarts)) {
    control$rstarts <- 1L
  }

  if(is.null(control$cores)) {
    control$cores <- 1L
  }

  if(is.null(control$tcg_maxit)) {
    control$tcg_maxit <- 10L
  }

  return(control)

}
