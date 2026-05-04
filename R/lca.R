# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 04/05/2026
#'
#' @title
#' Latent Class Analysis.
#' @description
#'
#' Estimate latent class models with gaussian and multinomial item models,
#' with or without covariates to predict class membership.
#'
#' @usage
#'
#' lca(data, nclasses = 2L, item = rep("gaussian", ncol(data)),
#'     X = NULL, penalties = TRUE, model = NULL,
#'     start = NULL, do.fit = TRUE, verbose = TRUE, control = NULL)
#'
#' @param data data frame or matrix.
#' @param nclasses Number of latent classes.
#' @param item Character vector with the model for each item (i.e., "gaussian" or "multinomial"). Defaults to "gaussian" for all the items.
#' @param X Matrix of covariates.
#' @param penalties Boolean or list of penalty terms for the parameters.
#' @param model List of parameter labels. See 'details' for more information.
#' @param start List of starting values for the parameters. See 'details' for more information.
#' @param do.fit TRUE to fit the model and FALSE to return only the model setup. Defaults to TRUE.
#' @param control List of control parameters for the optimization algorithm. See 'details' for more information.
#' @param verbose Print information of model estimation. Defaults to FALSE.
#'
#' @details \code{lca} estimates models with categorical and continuous data.
#'
#' @return List with the following objects:
#' \item{version}{Version number of 'latent' when the model was estimated.}
#' \item{call}{Code used to estimate the model.}
#' \item{ModelInfo}{Model information.}
#' \item{Optim}{Output of the optimizer.}
#' \item{user_model}{Structure with most relevant parameters.}
#' \item{parameters}{Structure with all model parameters.}
#' \item{transparameters}{Structure with all transformed model parameters.}
#' \item{posterior}{Posterior probability of class membership.}
#' \item{state}{Class with highest posterior probability.}
#' \item{loglik}{Logarithm likelihood of the model.}
#' \item{penalized_loglik}{Logarithm likelihood + logarithm priors of the model.}
#' \item{loglik_case}{Logarithm likelihood of each pattern.}
#' \item{summary_table}{Table of summaries of the fitted model. Useful when all the items are 'multinomial'.}
#' \item{ClassConditional}{Parameters that are conditional on the class memberships.}
#' \item{ResponseConditional}{Probability of class memberships conditional on item response. Only available when all the items have a 'multinomial' likelihood.}
#' \item{probCat}{Marginal probability of an item response. Only available when all the items have a 'multinomial' likelihood.}
#'
#' @examples
#'
#' \dontrun{
#'
#' fit <- lca(data = gss82, nclasses = 3L,
#'            item = rep("multinomial", ncol(gss82)),
#'            penalties = TRUE, do.fit = TRUE)
#' fit@timing
#' fit@loglik # -2754.643
#' fit@penalized_loglik # -2759.507
#' fit@Optim$opt$iterations
#'
#' # Plot model fit info:
#' fit
#'
#' # Get fit indices:
#' getfit(fit)
#'
#' # Get a summary:
#' summary(fit)
#'
#' # Inspect model objects:
#' latInspect(fit, what = "coefs", digits = 3)
#' latInspect(fit, what = "classes", digits = 3)
#' latInspect(fit, what = "profile", digits = 3)
#' latInspect(fit, what = "posterior", digits = 3)
#' latInspect(fit, what = "table", digits = 3)
#' latInspect(fit, what = "pattern", digits = 3)
#'
#' # Get standard errors:
#' SE <- se(fit, type = "standard", digits = 4)
#' SE$table
#'
#' # Get confidence intervals:
#' CI <- ci(fit, type = "standard", confidence = 0.95, digits = 2)
#' CI$table
#' }
#'
#' @references
#'
#' None yet.
#'
#' @export
lca <- function(data, nclasses = 1L,
                gaussian = NULL,
                multinomial = NULL,
                X = NULL,
                penalties = TRUE,
                model = NULL,
                start = NULL,
                control = NULL,
                do.fit = TRUE,
                verbose = TRUE) {

  ## store original call
  mc  <- match.call()
  args <- as.list(match.call(expand.dots = TRUE))[-1]

  if(is.character(X)) {
    X <- data[, X]
  }

  ngaussian <- length(gaussian)
  nmultinomial <- length(multinomial)
  model_labels <- rep(c("gaussian", "multinomial"),
                      times = c(ngaussian, nmultinomial))
  model_vector <- c(gaussian, multinomial)
  idx <- match(colnames(data), model_vector); idx <- idx[!is.na(idx)]
  data <- data[, model_vector[idx]]
  item <- model_labels[idx]

  original_model <- model
  original_X <- X

  # Fix the measurement part of the model if a previous llca fit is available:
  # This is usually done for two-step estimators of covariate coefficients.
  if(class(model) == "llca") {

    model <- model@parameters
    model$beta <- NULL

  }

  # Check control parameters:
  control$penalties <- penalties
  control$start <- start
  control <- lca_control(control)

  #### Initial input checks ####

  # Number of subjects:
  nobs <- nrow(data)
  # Number of items:
  nitems <- ncol(data)

  # Check that data is either a data.frame or a matrix:
  if(!is.data.frame(data) & !is.matrix(data)) {
    stop("data must be a matrix or data.frame")
  }

  # Check that item is a character vector with a string for each column of data:
  if(!is.character(item) || length(item) != nitems) {

    stop("item must be a character vector with as many elements as columns in data")

  }

  # Process the covariates:
  if(is.null(X)) {

    X <- matrix(1, nrow = nobs, ncol = 1L)
    colnames(X) <- "(Intercept)"

  } else {

    # Check that X is either a data.frame or a matrix:
    if(!is.data.frame(X) & !is.matrix(X)) {
      stop("data must be a matrix or data.frame")
    }

    if(nrow(X) != nobs) {
      stop("Number of cases in the data and covariates do not match")
    }

    # Transform characters into factors:
    X_df <- as.data.frame(X)
    X_df[] <- lapply(X_df, function(x) if (is.character(x)) factor(x) else x)

    # Create the design matrix:
    X <- model.matrix(~ . + 1, X_df)
    # Center the variables:
    # X[, -1] <- apply(X[, -1], MARGIN = 2, FUN = \(x) x-mean(x))

    # Put an underscore between the variable names and their level names:
    for (v in names(X_df)[sapply(X_df, is.factor)]) {
      i <- which(startsWith(colnames(X), v))
      colnames(X)[i] <- paste0(v, "_", make.names(levels(X_df[[v]]))[seq_along(i)])
    }

  }

  #### Process the data ####

  # Check if any column in "data" is modeled with the multinomial model and
  # transform it to a factor:

  condition <- item == "multinomial"
  factor_indices <- which(condition)
  # Transform into factors:
  data[, factor_indices] <- lapply(data[, factor_indices, drop = FALSE],
                                   FUN = factor)
  # Save category names for data modeled with the multinomial model:
  factor_names <- lapply(data[, factor_indices, drop = FALSE], levels)
  data[, factor_indices] <- lapply(data[, factor_indices, drop = FALSE],
                                   FUN = function(col) {
                                     if (is.factor(col)) as.integer(col) - 1L else col
                                   })

  # Get the number of possible response patterns:
  if(all(condition)) { # If all the items are multinomial...

    # Count the number of categories:
    Ks <- apply(data[, factor_indices, drop = FALSE], MARGIN = 2,
                FUN = \(x) length(unique(x[!is.na(x)])))
    # If all the items are multinomial, the number of possible patterns is:
    npossible_patterns <- min(prod(Ks)-1, nobs)
    # npossible_patterns <- min(prod(Ks), npatterns)

  } else { # If any item is not multinomial...

    npossible_patterns <- nobs

  }

  # Combine the data and the covariates:
  dt <- cbind(X, data)
  pcov <- ncol(X) # Number of covariates

  # Convert data to a data.table object:
  dt <- data.table::as.data.table(dt)

  ## Collect some information from the data ##
  counts_dt <- dt[, .(index = .I[1], count = .N), by = names(dt)]
  # Data matrix with the unique response patterns:
  patterns <- as.matrix(counts_dt[, names(dt), with = FALSE])[, -(1:pcov) , drop = FALSE]
  # Covariates with unique response patterns:
  cov_patterns2 <- as.matrix(counts_dt[, names(dt), with = FALSE])[, 1:pcov , drop = FALSE]
  # Number of unique response patterns:
  npatterns <- nrow(counts_dt)
  # Counts of each response pattern:
  weights <- counts_dt$count
  # Indices to map the original data to the matrix of unique patterns:
  full2short <- counts_dt$index
  # Indices to map the matrix of unique patterns to the original data:
  short2full <- match(do.call(paste, dt),
                      do.call(paste, counts_dt[, -c("index", "count"), with = FALSE]))

  ## Collect some information from the covariates ##
  cov_dt <- data.table::as.data.table(X)
  counts_cov_dt <- cov_dt[, .(index = .I[1], count = .N), by = names(cov_dt)]
  # Covariates with unique patterns:
  cov_patterns <- as.matrix(counts_cov_dt[, names(cov_dt), with = FALSE])
  # Number of unique covariate patterns:
  ncov_patterns <- nrow(cov_patterns)
  # Indices to map the original data to the matrix of unique patterns:
  cov_full2short <- counts_cov_dt$index
  # Indices to map the matrix of unique patterns to the original data:
  cov_short2full <- match(do.call(paste, cov_dt),
                          do.call(paste, counts_cov_dt[, -c("index", "count"),
                                                       with = FALSE]))

  # Put in a list the objects generated form the data:
  dataList <- vector("list")
  dataList$data <- data
  dataList$X <- X
  dataList$item <- item
  dataList$nobs <- nobs
  dataList$patterns <- patterns
  dataList$cov_patterns2 <- cov_patterns2
  dataList$npatterns <- npatterns
  dataList$npossible_patterns <- npossible_patterns
  dataList$nitems <- nitems
  dataList$weights <- weights
  dataList$full2short <- full2short
  dataList$short2full <- short2full
  dataList$item_names <- colnames(data)
  dataList$factor_indices <- factor_indices
  dataList$factor_names <- factor_names
  dataList$cov_patterns <- cov_patterns
  dataList$ncov_patterns <- ncov_patterns
  dataList$cov_full2short <- cov_full2short
  dataList$cov_short2full <- cov_short2full
  dataList$original_X <- original_X
  dataList$args <- args
  dataList$original_model <- original_model
  dataList$nclasses <- nclasses

  #### Initialize objects to store all the models ####

  model0 <- model
  NCLASSES <- nclasses
  nmodels <- length(NCLASSES)

  llca_list <- list()

  # Loop to create and fit the models for different nclasses:
  for(NK in 1:nmodels) {

    # NK <- 1L
    if(verbose) {

      print(paste0("Model nclasses = ", NCLASSES[NK]) )

    }

    # Ensure that nclasses is an integer:
    nclasses <- as.integer(NCLASSES[NK])
    # Check that nclasses is positive:
    if(nclasses < 1L) {
      stop("nclasses must be a positive integer")
    }

    #### Create the model ####

    # Get the model specification:
    full_model <- create_lca_model(dataList = dataList, nclasses = nclasses,
                                   item = item, model = model,
                                   control = control)
    list2env(full_model, envir = environment())

    #### Create the structures ####

    # Generate the structures for optimization:
    modelInfo <- create_lca_modelInfo(dataList = dataList,
                                      full_model = full_model,
                                      control = control)
    modelInfo$original_model <- original_model

    #### Fit the model ####

    result <- vector("list") # Initialize the object to be returned

    # Fit the model or just get the model specification:
    if(!do.fit) {

      llca_list[[NK]] <- new("llca",
                             version            = as.character(packageVersion('latent')),
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
                             penalized_loss     = numeric(),
                             extra              = list()
      )

      next # Go to the next model with a different value of "nclasses"

    }

    modelInfo$control_optimizer$cores <- min(modelInfo$control_optimizer$rstarts,
                                             modelInfo$control_optimizer$cores)
    # Perform the optimization (fit the model):
    Optim <- optimizer(control_manifold = modelInfo$control_manifold,
                       control_transform = modelInfo$control_transform,
                       control_estimator = modelInfo$control_estimator,
                       control_optimizer = modelInfo$control_optimizer)

    # Collect all the information about the optimization:

    names(Optim$parameters) <- modelInfo$parameters_labels
    names(Optim$transparameters) <- modelInfo$transparameters_labels

    #### Process the outputs ####

    # Create the transformed parameters:
    transformed_pars <- fill_in(modelInfo$trans,
                                Optim$transparameters)

    parameters <- transformed_pars[names(modelInfo$param)]

    #### Process the fit information ####

    loss <- sum(unlist(lapply(Optim$outputs$estimators$doubles,
                              FUN = \(x) x[[1]])))
    loglik <- sum(unlist(lapply(Optim$outputs$estimators$doubles,
                                FUN = \(x) x[[2]])))
    penalty <- sum(unlist(lapply(Optim$outputs$estimators$doubles,
                                 FUN = \(x) x[[5]])))
    penalized_loss <- loss + penalty
    penalized_loglik <- loglik - penalty

    #### latent object ####

    elapsed <- Optim$elapsed

    llca_list[[NK]] <- new("llca",
                           version            = as.character(packageVersion('latent')),
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
                           penalized_loss     = penalized_loss,
                           extra              = list()
    )


  }

  #### Return ####

  class(llca_list) <- "llcalist"

  if(nmodels == 1) {
    return(llca_list[[1]])
  } else{
    return(llca_list)
  }

}

create_lca_model <- function(dataList, nclasses, item,
                             model = NULL, control) {

  # Generate the model syntax and initial parameter values

  list2env(dataList, envir = environment())

  class_names <- paste("Class", 1:nclasses, sep = "")
  pattern_names <- paste("pattern", 1:npatterns, sep = "")

  #### Model for the transformed parameters ####

  p <- ncol(X) # Number of predictors
  pred_names <- colnames(X) # Names of predictors

  list_struct <- vector("list")
  k <- 1L

  # Model for the betas:

  list_struct[[k]] <- list(name = "beta",
                           type = "matrix",
                           dim = c(p, nclasses),
                           rownames = pred_names,
                           colnames = class_names)
  k <- k+1L

  list_struct[[k]] <- list(name = "theta",
                           type = "matrix",
                           dim = c(npatterns, nclasses),
                           rownames = pattern_names,
                           colnames = class_names)
  k <- k+1L

  list_struct[[k]] <- list(name = "class",
                           type = "matrix",
                           dim = c(npatterns, nclasses),
                           rownames = pattern_names,
                           colnames = class_names)
  k <- k+1L

  list_struct[[k]] <- list(name = "loglik",
                           type = "array",
                           dim = c(npatterns, nitems, nclasses),
                           dimnames = list(pattern_names,
                                           item_names,
                                           class_names))
  k <- k+1L

  # Model for gaussian items:
  if(any(item == "gaussian")) {

    gauss <- which(item == "gaussian")
    Jgauss <- length(gauss) # Number of gaussian items

    for(j in 1:Jgauss) {

      list_struct[[k]] <- list(name = item_names[gauss[j]],
                               type = "matrix",
                               dim = c(3, nclasses),
                               rownames = c("mean",
                                            "stdv",
                                            "log(stdv)"),
                               colnames = class_names)
      k <- k+1L

    }

  }

  # Model for multinomial items:
  if(any(item == "multinomial")) {

    multinom <- which(item == "multinomial")
    Jmulti <- length(multinom) # Number of multinomial items

    K <- unlist(lapply(factor_names, FUN = length))
    if(any(K > 30)) {
      stop("You cannot use a multinomial likelihood to model a variable with more than 30 unique categories")
    }

    for(j in 1:Jmulti) {

      log_probs <- paste("log", item_names[multinom[j]], sep = "_")
      list_struct[[k]] <- list(name = log_probs,
                               type = "matrix",
                               dim = c(K[j], nclasses),
                               rownames = factor_names[[j]],
                               colnames = class_names)
      k <- k+1L

      list_struct[[k]] <- list(name = item_names[multinom[j]],
                               type = "matrix",
                               dim = c(K[j], nclasses),
                               rownames = factor_names[[j]],
                               colnames = class_names)
      k <- k+1L

    }

  }

  # Create the full transparameter structure:
  trans <- create_parameters(list_struct)

  #### Model for the parameters ####

  param <- list()

  # Model for the betas:
  param$beta <- trans$beta
  param$beta[, 1] <- "0"

  # Model for gaussian items:
  if(any(item == "gaussian")) {

    cl <- length(param)
    idx_gauss <- (cl + 1L):(cl + Jgauss)
    param[idx_gauss] <- trans[item_names[gauss]]
    param[idx_gauss] <- lapply(trans[item_names[gauss]],
                               FUN = \(x) {
                                 x[2, ] <- "1"
                                 return(x)
                               })
    names(param)[idx_gauss] <- item_names[gauss]

  }

  # Model for multinomial items:
  if(any(item == "multinomial")) {

    cl <- length(param)
    idx_multinom <- (cl + 1L):(cl + Jmulti)
    item_names_multinom <- paste("log", item_names[multinom], sep = "_")
    param[idx_multinom] <- trans[item_names_multinom]
    names(param)[idx_multinom] <- item_names_multinom

    param[idx_multinom] <- lapply(param[idx_multinom],
                                  FUN = \(x) {
                                    x[1, ] <- "0"; return(x)
                                  })

  }

  #### Fixed parameters ####

  lca_all <- trans

  # Replace the transformed parameter by custom values, if available:
  if(!is.null(model)) {

    # Replace the parameters by custom values:

    nm <- intersect(names(model), names(param))
    nm <- nm[!vapply(model[nm], is.null, logical(1))]
    param[nm] <- model[nm]

    nm <- intersect(names(model), names(lca_all))
    nm <- nm[!vapply(model[nm], is.null, logical(1))]
    lca_all[nm] <- model[nm]

  }

  lca_all[names(param)] <- param

  #### Create the initial values for the parameters ####

  # Initial values for the log coefficients (betas):
  # betas are sampled from a normal distribution:
  init_param <- vector("list", length = control$rstarts)

  if(any(item == "gaussian")) {

    init_mean <- init_sd <- init_logsd <- list()
    for(j in 1:Jgauss) {

      init_mean[[j]] <- rep(mean(data[, gauss[j]], na.rm = TRUE),
                            times = nclasses)
      init_sd[[j]] <- rep(sd(data[, gauss[j]], na.rm = TRUE), times = nclasses)
      init_logsd[[j]] <- log(init_sd[[j]])

    }

  }

  pi_hat_list <- list()
  if(any(item == "multinomial")) {

    eta_hat_list <- list()
    for(j in 1:Jmulti) {

      int_vector <- data[, multinom[j]]
      int_vector <- int_vector[!is.na(int_vector)] # Remove missing values
      nsize <- length(int_vector)
      props <- count(int_vector, nsize, K[j]) / nsize
      pi_hat_list[[j]] <- props %*% t(rep(1, nclasses))
      log_props <- log(props)
      eta_hat_list[[j]] <- (log_props-log_props[1]) %*% t(rep(1, nclasses))

    }

    pi_hat_vector <- unlist(pi_hat_list)
    vars <- (1-pi_hat_vector)/(nobs*pi_hat_vector)
    sds <- sqrt(vars)
    Ks <- length(unlist(eta_hat_list))

  }

  for(i in 1:control$rstarts) {

    init_param[[i]] <- vector("list", length = length(param))
    names(init_param[[i]]) <- names(param)

    # Initial values for betas:
    init_beta <- matrix(rnorm(p*nclasses), nrow = p, ncol = nclasses)
    init_param[[i]][["beta"]] <- init_beta
    init_param[[i]]$beta[, 1] <- 0

    # Initial values for gaussian items:
    if(any(item == "gaussian")) {

      for(j in 1:Jgauss) {

        rmean <- init_mean[[j]] + rnorm(nclasses,
                                        mean = 0,
                                        sd = init_sd[[j]]/sqrt(nobs))

        init_param[[i]][[item_names[gauss[j]]]] <- rbind(rmean,
                                                         init_logsd[[j]],
                                                         init_sd[[j]])
        dimnames(init_param[[i]][[item_names[gauss[j]]]]) <-
          dimnames(param[[item_names[gauss[j]]]])

      }

    }

    # Initial values for multinomial items:
    if(any(item == "multinomial")) {

      item_names_multinom <- paste("log", item_names[multinom], sep = "_")

      for(j in 1:Jmulti) {

        sds <- (1-pi_hat_list[[j]])/(nobs*pi_hat_list[[j]])

        init_param[[i]][[item_names_multinom[j]]] <- eta_hat_list[[j]] +
          rnorm(length(eta_hat_list[[j]]), 0, sds)
        init_param[[i]][[item_names_multinom[j]]][1, ] <- 0

        dimnames(init_param[[i]][[item_names[item_names_multinom[j]]]]) <-
          dimnames(param[[item_names[item_names_multinom[j]]]])

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
                 lca_all = lca_all,
                 init_param = init_param,
                 pi_hat_list = pi_hat_list)

  return(result)

}

create_lca_modelInfo <- function(dataList, full_model, control) {

  # Generate control_manifold, control_transform, and control_estimator

  list2env(dataList, envir = environment())
  list2env(full_model, envir = environment())

  nclasses <- ncol(trans$class)

  #### Manifolds ####

  manifolds <- list(
    list(manifold = "euclidean", parameters = names(param))
  )

  control_manifold <- create_manifolds(manifolds = manifolds,
                                       structures = param)

  #### Transformations ####

  transforms <- list()
  k <- 1L

  # betas to thetas:
  transforms[[k]] <- list(transform = "column_space",
                          parameters_in = "beta",
                          parameters_out = "theta",
                          extra = list(X = cov_patterns2))
  k <- k+1L

  # thetas to classes:
  for(s in 1:npatterns) {

    transforms[[k]] <- list(transform = "softmax",
                            parameters_in = list(trans$theta[s, ]),
                            parameters_out = list(trans$class[s, ]))
    k <- k + 1L

  }

  # Conditional item likelihoods (gaussian):
  if(any(item == "gaussian")) {

    gauss <- which(item == "gaussian")
    Jgauss <- length(gauss) # Number of gaussian items

    means <- c(do.call(rbind, lapply(trans[item_names[gauss]],
                                     FUN = \(x) x[1, ])))
    stdv <- c(do.call(rbind, lapply(trans[item_names[gauss]],
                                    FUN = \(x) x[2, ])))
    log_stdv <- c(do.call(rbind, lapply(trans[item_names[gauss]],
                                        FUN = \(x) x[3, ])))

    transforms[[k]] <- list(transform = "exponential",
                            parameters_in = list(log_stdv),
                            parameters_out = list(stdv))
    k <- k+1L

    y <- as.matrix(patterns[, gauss])
    transforms[[k]] <- list(transform = "normal",
                            parameters_in = list(means, stdv),
                            parameters_out = list(trans$loglik[, gauss, ]),
                            extra = list(y = y, S = npatterns, J = Jgauss,
                                         I = nclasses))
    k <- k+1L

  }

  # Conditional item likelihoods (multinomial):
  if(any(item == "multinomial")) {

    # All of this is redundant in get_full_lca_model:
    multinom <- which(item == "multinomial")
    Jmulti <- length(multinom) # Number of multinomial items

    # Softmax transformations:
    for(j in 1:Jmulti) {

      log_probs <- paste("log", item_names[multinom][j], sep = "_")
      probs <- item_names[multinom[j]]

      for(i in 1:nclasses) {

        transforms[[k]] <- list(transform = "softmax",
                                parameters_in = list(trans[[log_probs]][, i]),
                                parameters_out = list(trans[[probs]][, i]))
        k <- k+1L

      }
    }

    # Multinomial transformation:

    y <- as.matrix(patterns[, multinom])
    K <- unlist(lapply(dataList$factor_names, FUN = length))
    transforms[[k]] <- list(transform = "multinomial",
                            parameters_in = list(trans[item_names[multinom]]),
                            parameters_out = list(trans$loglik[, multinom, ]),
                            extra = list(y = y, S = npatterns, J = Jmulti,
                                         I = nclasses, K = K))

  }

  control_transform <- create_transforms(transforms = transforms,
                                         structures = trans)

  #### Estimators ####

  estimators <- list()
  G <- 1L

  estimators[[G]] <- list(estimator = "lca",
                          parameters = c("class", "loglik"),
                          extra = list(S = npatterns,
                                       J = nitems,
                                       I = nclasses,
                                       weights = weights))
  G <- G + 1L

  # Choose whether using Bayes constants:
  if(control$reg) {

    # Bayes Constant for class probabilities:
    alpha <- control$penalties$class$alpha
    if(alpha != 0) {

      # Get the indices corresponding to the unique covariate patterns:
      dt_uniq_X <- data.table::as.data.table(cov_patterns2)
      counts_X <- dt_uniq_X[, .(index = .I[1], count = .N), by = names(dt_uniq_X)]
      uniques <- counts_X$index
      U <- length(uniques)

      for(i in uniques) {

        estimators[[G]] <- list(estimator = "bayesconst1",
                                parameters = list(trans$class[i, ]),
                                extra = list(K = nclasses,
                                             alpha = alpha,
                                             U = U,
                                             N = nobs))
        G <- G+1L

      }
    }

    # Bayes Constant for standard deviations:
    alpha <- control$penalties$sd$alpha
    if(any(item == "gaussian") & alpha != 0) {

      Y <- data[, gauss, drop = FALSE]
      # sigma_class <- split(trans$sigma, rep(1:nclasses, each = Jgauss))
      varshat <- apply(Y, MARGIN = 2, FUN = var, na.rm = TRUE)*(nobs-1)/nobs

      for(i in 1:nclasses) {

        stdv_by_class <- unlist(lapply(trans[item_names[gauss]],
                                       FUN = \(x) x[2, i]))
        estimators[[G]] <- list(estimator = "bayesconst3",
                                parameters = list(stdv_by_class),
                                extra = list(K = nclasses,
                                             varshat = varshat,
                                             alpha = alpha,
                                             N = nobs))
        G <- G+1L

      }

    }

    # Bayes Constant for multinomial probabilities:
    alpha <- control$penalties$prob$alpha
    if(any(item == "multinomial") & alpha != 0) {

      for(j in 1:Jmulti) {
        for(i in 1:nclasses) {

          pihat <- pi_hat_list[[j]][, i]
          probs_by_class <- unlist(lapply(trans[item_names[multinom[j]]],
                                          FUN = \(x) x[, i]))

          estimators[[G]] <- list(estimator = "bayesconst2",
                                  # parameters = list(trans$peta[[j]][, i]),
                                  parameters = list(probs_by_class),
                                  extra = list(K = nclasses,
                                               pihat = pihat,
                                               alpha = alpha,
                                               N = nobs))
          G <- G+1L

        }
      }

    }

    # Gaussian regularization for coefficients:
    alpha <- control$penalties$beta$alpha
    if(alpha != 0) {

      p <- nrow(trans$beta)-1L
      q <- ncol(trans$beta)
      means <- matrix(0, nrow = p, ncol = q)
      sds <- apply(cov_patterns2[, -1], MARGIN = 2, sd, na.rm = TRUE)
      sds <- matrix(sds, nrow = p, ncol = q) / alpha

      estimators[[G]] <- list(estimator = "gaussian_loglik",
                              parameters = list(trans$beta[-1, ]),
                              extra = list(means = means,
                                           sds = sds,
                                           alpha = alpha,
                                           N = nobs))
      G <- G+1L


    }

    # Ridge regularization for coefficients:
    lambda <- control$penalties$beta$lambda
    power <- control$penalties$beta$power
    if(!is.null(lambda) & !is.null(power)) {

      if(lambda != 0 && power != 0) {

        estimators[[G]] <- list(estimator = "ridge",
                                parameters = list(trans$beta[-1, ]),
                                extra = list(lambda = lambda,
                                             power = power,
                                             N = nobs))
        G <- G+1L

      }
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

  #### Return ####

  modelInfo <- list(param = param,
                    trans = trans,
                    nparam = nparam,
                    ntrans = ntrans,
                    parameters_labels = parameters_labels,
                    transparameters_labels = transparameters_labels,
                    dof = npatterns - nparam,
                    control_manifold = control_manifold,
                    control_transform = control_transform,
                    control_estimator = control_estimator,
                    control_optimizer = control_optimizer)

  return(modelInfo)

}

lca_control <- function(control) {

  # Auxiliary function for lca.R

  # Control input

  if(isFALSE(control$penalties)) {

    control$reg <- FALSE

  } else if(isTRUE(control$penalties)) {

    control$reg <- TRUE

    control$penalties <- list(
      # beta  = list(alpha = 0, lambda = 0, power = 0),
      beta  = list(alpha = 0),
      class = list(alpha = 1),
      prob  = list(alpha = 1),
      sd    = list(alpha = 1)
    )

  } else if(is.list(control$penalties)) {

    control$reg <- TRUE

  } else {

    stop("penalties should be TRUE, FALSE, or a list")

  }

  if(is.null(control$step_maxit)) {
    control$step_maxit <- 30L
  } else if(control$step_maxit < 1L) {
    stop("step_maxit must be an integer greater than 0")
  }

  if(is.null(control$c1)) {
    control$c1 <- 0.5
  } else if(control$c1 < 0L) {
    stop("c1 must be a positive number, preferable lower than c2")
  }

  if(is.null(control$c2)) {
    control$c2 <- 0.5
  } else if(control$c2 < 0L) {
    stop("c2 must be a positive number, preferable larger than c1")
  }

  if(is.null(control$step_eps)) {
    control$step_eps <- 1e-09
  } else if(control$step_eps < 0) {
    stop("step_eps must be a positive number, preferable close to 0")
  }

  if(is.null(control$df_eps)) {
    control$df_eps <- 1e-09
  } else if(control$df_eps < 0) {
    stop("df_eps must be a positive number, preferable close to 0")
  }

  if(is.null(control$eps)) {
    control$eps <- 1e-05
  } else if(control$eps < 0) {
    stop("eps must be a positive number, preferable close to 0")
  }

  if(is.null(control$M)) {
    control$M <- 100L
  } else if(control$M < 0L) {
    stop("M must be a positive integer")
  }

  if(is.null(control$ss_fac)) {
    control$ss_fac <- 2
  } else if(control$ss_fac <= 1) {
    stop("ss_fac must be a positive integer larger than 1")
  }

  if(is.null(control$maxit)) {
    control$maxit <- 1000L
  } else if(control$maxit < 0L) {
    stop("maxit must be a positive integer")
  }

  if(is.null(control$cores)) {
    control$cores <- parallel::detectCores()-1L
  } else if(control$cores < 0L) {
    stop("cores must be a positive integer")
  }

  if(is.null(control$tcg_maxit)) {
    control$tcg_maxit <- 10L
  } else if(control$tcg_maxit < 0L) {
    stop("tcg_maxit must be a positive integer")
  }

  if(is.null(control$opt)) {
    control$opt <- "lbfgs"
  }

  if(is.null(control$rstarts)) {
    control$rstarts <- 16L
  } else if(control$rstarts < 0L) {
    stop("rstarts must be a positive integer")
  }

  return(control)

}
