# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 28/03/2026

get_full_lca_model <- function(data_list, nclasses, item,
                               model = NULL, control = NULL) {

  # Generate the model syntax and initial parameter values

  list2env(data_list, envir = environment())

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
  lca_all <- create_parameters(list_struct)

  transparameters_labels <- unname(unique(unlist(lca_all)))
  ntrans <- length(transparameters_labels)

  #### Model for the parameters ####

  param <- list()

  # Model for the betas:
  param$beta <- lca_all$beta
  param$beta[, 1] <- "0"

  # Model for gaussian items:
  if(any(item == "gaussian")) {

    cl <- length(param)
    idx_gauss <- (cl + 1L):(cl + Jgauss)
    param[idx_gauss] <- lca_all[item_names[gauss]]
    param[idx_gauss] <- lapply(lca_all[item_names[gauss]],
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
    param[idx_multinom] <- lca_all[item_names_multinom]
    names(param)[idx_multinom] <- item_names_multinom

    param[idx_multinom] <- lapply(param[idx_multinom],
                                      FUN = \(x) {
                                        x[1, ] <- "0"; return(x)
                                        })

  }

  #### Fixed parameters ####

  trans <- lca_all

  # Replace the transformed parameter by custom values, if available:
  if(!is.null(model)) {

    # Replace the parameters by custom values:

    nm <- intersect(names(model), names(param))
    nm <- nm[!vapply(model[nm], is.null, logical(1))]
    param[nm] <- model[nm]

    nm <- intersect(names(model), names(trans))
    nm <- nm[!vapply(model[nm], is.null, logical(1))]
    trans[nm] <- model[nm]

  }

  trans[names(param)] <- param

  #### Arrange labels ####

  # Collect the parameter labels:
  vector_param <- unname(unlist(param))
  indicator <- is.na(suppressWarnings(as.numeric(vector_param)))
  fixed_indices <- which(!indicator)
  fixed_values <- as.numeric(vector_param[fixed_indices])
  parameters_indices <- which(indicator)
  parameters_labels <- unique(vector_param[parameters_indices])
  nparam <- length(parameters_labels)

  # Indices of fixed parameters:
  alltrans <- suppressWarnings(as.numeric(unname(unlist(trans))))
  fixed_indices <- which(!is.na(alltrans))
  fixed_values <- alltrans[fixed_indices]

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
    init_beta <- replicate(nclasses, rnorm(p))
    init_param[[i]][["beta"]] <- fill_list_with_vector(param$beta,
                                                       init_beta)
    init_param[[i]]$beta[, 1] <- "0"

    # Initial values for gaussian items:
    if(any(item == "gaussian")) {

      for(j in 1:Jgauss) {

        rmean <- init_mean[[j]] + rnorm(nclasses,
                                        mean = 0,
                                        sd = init_sd[[j]]/sqrt(nobs))
        rlogsd <- init_logsd[[j]]

        values <- vector(length = 3*nclasses)
        values[seq(1, 3*nclasses, by = 3)] <- rmean
        values[seq(2, 3*nclasses, by = 3)] <- init_sd[[j]]
        values[seq(3, 3*nclasses, by = 3)] <- rlogsd

        init_param[[i]][item_names[gauss[j]]] <- fill_list_with_vector(param[item_names[gauss[j]]],
                                                                       values)

      }

    }

    # Initial values for multinomial items:
    if(any(item == "multinomial")) {

      item_names_multinom <- paste("log", item_names[multinom], sep = "_")

      values <- unlist(eta_hat_list) + rnorm(sum(Ks), 0, sds)
      # values <- unlist(eta_hat_list) + rnorm(sum(Ks), 0, 0.01)
      init_param[[i]][item_names_multinom] <- fill_list_with_vector(param[item_names_multinom],
                                                                    values)

      for(j in 1:Jmulti) {
        init_param[[i]][[item_names_multinom[j]]][1, ] <- "0"
      }

    }

  }

  # Transform all the initial parameters in numerical values:
  init_param <- allnumeric(init_param)

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

  #### Pass the initial values to vectors ####

  parameters <- vector("list", length = control$rstarts)
  select_params <- match(parameters_labels, vector_param)
  for(i in 1:control$rstarts) {
    parameters[[i]] <- unlist(init_param[[i]])[select_params]
    names(parameters[[i]]) <- parameters_labels
  }

  # Initialize the vector of transformed parameters:
  # Relate the transformed parameters to the parameters:
  trans2param <- match(parameters_labels, transparameters_labels)
  transparameters <- vector("list", length = control$rstarts)
  for(i in 1:control$rstarts) { # For each random start i...
    transparameters[[i]] <- vector(length = ntrans)
    # Copy the parameter values to the transparameters:
    transparameters[[i]][fixed_indices] <- fixed_values # Fixed parameters
    transparameters[[i]][trans2param] <- parameters[[i]]
    names(transparameters[[i]]) <- transparameters_labels
  }

  #### Set up the optimizer ####

  # Create defaults for the control of the optimizer:
  control$init_param <- init_param
  control$parameters <- parameters
  control$transparameters <- transparameters
  control$transparam2param <- trans2param-1L

  #### Return ####

  result <- list(parameters_labels = parameters_labels,
                 nparam = nparam,
                 transparameters_labels = transparameters_labels,
                 ntrans = ntrans,
                 param = param,
                 trans = trans,
                 lca_all = lca_all,
                 pi_hat_list = pi_hat_list,
                 init_param = init_param,
                 control = control)

  return(result)

}

get_short_lca_model <- function(data_list, nclasses, item,
                                lca_all, param, model = NULL) {

  # This function displays the reduced LCA model in logarithm and probability
  # scales. This is a short version of the full model syntax created with
  # get_full_lca_covariate_model.

  # Create a short summary model:

  data <- data_list$data
  gauss <- which(item == "gaussian")
  multinom <- which(item == "multinomial")
  K <- lapply(data_list$factor_names, FUN = length)
  nitems <- ncol(data)
  items <- vector("list", length = nitems)
  names(items) <- colnames(data)

  item_names <- data_list$item_names

  prob_model <- list()
  prob_model$beta <- lca_all$beta

  #### Transformed parameters model ####

  if(any(item == "gaussian")) {

    prob_model[item_names[gauss]] <- lca_all[data_list$item_names[gauss]]

  }

  if(any(item == "multinomial")) {

    prob_model[item_names[multinom]] <- lca_all[data_list$item_names[multinom]]

  }

  #### Parameters model ####

  log_model <- list()
  log_model$beta <- lca_all$beta

  if(any(item == "gaussian")) {

    log_model[item_names[gauss]] <- lapply(lca_all[data_list$item_names[gauss]],
                                            FUN = \(x) x[-2, ])

  }

  if(any(item == "multinomial")) {

    log_probs <- paste("log", data_list$item_names[multinom], sep = "_")
    log_model[item_names[multinom]] <- lca_all[log_probs]

  }

  # rownames(log_model$beta) <- colnames(data_list$cov_patterns2)

  #### Return ####

  # Return the model in the logarithm and probability scales:
  result <- list(prob_model = prob_model, log_model = log_model)

  return(result)

}

create_lca_modelInfo <- function(data_list, full_model, control) {

  # Generate control_manifold, control_transform, and control_estimator

  list2env(data_list, envir = environment())
  list2env(full_model, envir = environment())

  nclasses <- ncol(lca_all$class)

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
                            parameters_in = list(lca_all$theta[s, ]),
                            parameters_out = list(lca_all$class[s, ]))
    k <- k + 1L

  }

  # Conditional item likelihoods (gaussian):
  if(any(item == "gaussian")) {

    gauss <- which(item == "gaussian")
    Jgauss <- length(gauss) # Number of gaussian items

    means <- c(do.call(rbind, lapply(lca_all[item_names[gauss]],
                                     FUN = \(x) x[1, ])))
    stdv <- c(do.call(rbind, lapply(lca_all[item_names[gauss]],
                                    FUN = \(x) x[2, ])))
    log_stdv <- c(do.call(rbind, lapply(lca_all[item_names[gauss]],
                                        FUN = \(x) x[3, ])))

    transforms[[k]] <- list(transform = "exponential",
                            parameters_in = list(log_stdv),
                            parameters_out = list(stdv))
    k <- k+1L

    y <- as.matrix(patterns[, gauss])
    transforms[[k]] <- list(transform = "normal",
                            parameters_in = list(means, stdv),
                            parameters_out = list(lca_all$loglik[, gauss, ]),
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
                                parameters_in = list(lca_all[[log_probs]][, i]),
                                parameters_out = list(lca_all[[probs]][, i]))
        k <- k+1L

      }
    }

    # Multinomial transformation:

    y <- as.matrix(patterns[, multinom])
    K <- unlist(lapply(data_list$factor_names, FUN = length))
    transforms[[k]] <- list(transform = "multinomial",
                            parameters_in = list(lca_all[item_names[multinom]]),
                            parameters_out = list(lca_all$loglik[, multinom, ]),
                            extra = list(y = y, S = npatterns, J = Jmulti,
                                         I = nclasses, K = K))

  }

  control_transform <- create_transforms(transforms = transforms,
                                         structures = lca_all)

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
                                parameters = list(lca_all$class[i, ]),
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
      # sigma_class <- split(lca_all$sigma, rep(1:nclasses, each = Jgauss))
      varshat <- apply(Y, MARGIN = 2, FUN = var, na.rm = TRUE)*(nobs-1)/nobs

      for(i in 1:nclasses) {

        stdv_by_class <- unlist(lapply(lca_all[item_names[gauss]],
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
          probs_by_class <- unlist(lapply(lca_all[item_names[multinom[j]]],
                                         FUN = \(x) x[, i]))

          estimators[[G]] <- list(estimator = "bayesconst2",
                                  # parameters = list(lca_all$peta[[j]][, i]),
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

      p <- nrow(lca_all$beta)-1L
      q <- ncol(lca_all$beta)
      means <- matrix(0, nrow = p, ncol = q)
      sds <- apply(cov_patterns2[, -1], MARGIN = 2, sd, na.rm = TRUE)
      sds <- matrix(sds, nrow = p, ncol = q) / alpha

      estimators[[G]] <- list(estimator = "gaussian_loglik",
                              parameters = list(lca_all$beta[-1, ]),
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
                                parameters = list(lca_all$beta[-1, ]),
                                extra = list(lambda = lambda,
                                             power = power,
                                             N = nobs))
        G <- G+1L

      }
    }

  }

  control_estimator <- create_estimators(estimators = estimators,
                                         structures = lca_all)

  #### Return ####

  modelInfo <- list(nparam = nparam,
                    dof = npossible_patterns - nparam,
                    ntrans = ntrans,
                    parameters_labels = parameters_labels,
                    transparameters_labels = transparameters_labels,
                    param = param,
                    trans = trans,
                    lca_all = lca_all,
                    control_manifold = control_manifold,
                    control_transform = control_transform,
                    control_estimator = control_estimator,
                    control_optimizer = control)

  return(modelInfo)

}
