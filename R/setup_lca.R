# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 01/04/2026

create_lca_model <- function(data_list, nclasses, item,
                             model = NULL, control) {

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

get_short_lca_model <- function(data_list, nclasses, item,
                                trans, param, model = NULL) {

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
  prob_model$beta <- trans$beta

  #### Transformed parameters model ####

  if(any(item == "gaussian")) {

    prob_model[item_names[gauss]] <- trans[data_list$item_names[gauss]]

  }

  if(any(item == "multinomial")) {

    prob_model[item_names[multinom]] <- trans[data_list$item_names[multinom]]

  }

  #### Parameters model ####

  log_model <- list()
  log_model$beta <- trans$beta

  if(any(item == "gaussian")) {

    log_model[item_names[gauss]] <- lapply(trans[data_list$item_names[gauss]],
                                           FUN = \(x) x[-2, ])

  }

  if(any(item == "multinomial")) {

    log_probs <- paste("log", data_list$item_names[multinom], sep = "_")
    log_model[item_names[multinom]] <- trans[log_probs]

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
    K <- unlist(lapply(data_list$factor_names, FUN = length))
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
