# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 12/12/2025

get_full_lca_model <- function(data_list, nclasses, item,
                               model = NULL, control = NULL) {

  # Generate the model syntax and initial parameter values

  list2env(data_list, envir = environment())

  class_names <- paste("Class", 1:nclasses, sep = "")
  pattern_names <- paste("pattern", 1:npatterns, sep = "")

  # Initialize the objects to store the initial parameters:
  lca_param <- list()
  lca_all <- list()

  #### Model for the nontransformed parameters ####

  # Model for the betas:
  p <- ncol(X)
  labels <- paste("b", rep(1:p, times = nclasses), "|",
                  rep(1:nclasses, each = p), sep = "")
  lca_param$beta <- matrix(labels, nrow = p, ncol = nclasses)
  rownames(lca_param$beta) <- colnames(X)
  colnames(lca_param$beta) <- class_names
  lca_all$beta <- lca_param$beta
  lca_param$beta[, 1] <- "0"

  # Initial values for the log coefficients (betas):
  # betas are sampled from a normal distribution:
  init_param <- vector("list", length = control$rstarts)
  for(i in 1:control$rstarts) {
    for(k in 1:nclasses) {
      init_beta <- replicate(nclasses, rnorm(p))
      init_param[[i]]$beta <- fill_list_with_vector(lca_param$beta, init_beta)
      init_param[[i]]$beta[, 1] <- "0"
    }
  }

  # Model for gaussian items:
  if(any(item == "gaussian")) {

    gauss <- which(item == "gaussian")
    Jgauss <- length(gauss) # Number of gaussian items
    repitems <- rep(gauss, times = nclasses)
    repclasses <- rep(1:nclasses, each = Jgauss)
    mu <- paste("mu[", repitems, "|", repclasses, "]", sep = "")
    s <- paste("s[", repitems, "|", repclasses, "]", sep = "")
    lca_param$mu <- matrix(mu, nrow = Jgauss, ncol = nclasses)
    lca_param$s <- matrix(s, nrow = Jgauss, ncol = nclasses)
    colnames(lca_param$mu) <- colnames(lca_param$s) <-
      paste("Class", 1:nclasses, sep = "")
    rownames(lca_param$mu) <- rownames(lca_param$s) <-
      colnames(data)[gauss]
    lca_all$mu <- lca_param$mu
    lca_all$s <- lca_param$s

    # Initial values for mu and s:
    # For mu, they will be the mean of the items
    # For s, they will be the sd of the items
    init_mu <- rep(colMeans(data[, gauss, drop = FALSE], na.rm = TRUE),
                   times = nclasses)
    init_sd <- rep(apply(data[, gauss, drop = FALSE], MARGIN = 2,
                         FUN = sd, na.rm = TRUE), times = nclasses)
    init_s <- log(init_sd)
    for(i in 1:control$rstarts) {
      init_mui <- rnorm(Jgauss*nclasses, init_mu, init_sd/sqrt(nobs))
      init_param[[i]]$mu <- fill_list_with_vector(lca_param$mu, init_mui)
      init_param[[i]]$s <- fill_list_with_vector(lca_param$s, init_s)
    }

  }

  # Model for multinomial items:
  pi_hat_list <- list()
  if(any(item == "multinomial")) {

    multinom <- which(item == "multinomial")
    Jmulti <- length(multinom) # Number of multinomial items
    # Collect the number of response categories:
    K <- unlist(lapply(factor_names, FUN = length))
    if(any(K > 30)) {
      stop("You cannot model with a multinomial likelihood a variable with more than 30 unique categories")
    }
    eta <- vector("list", length = Jmulti)
    for(j in 1:Jmulti) {
      eta[[j]] <- matrix(NA, nrow = K[j], ncol = nclasses)
      eta[[j]][] <- paste("eta[", 1:K[j], "|", j, "|" ,
                          rep(1:nclasses, each = K[j]), "]", sep = "")
      colnames(eta[[j]]) <- class_names
      rownames(eta[[j]]) <- factor_names[[j]]
    }
    names(eta) <- colnames(data)[multinom]
    lca_all$eta <- eta
    lca_param$eta <- lapply(eta, FUN = \(x) {x[1,] <- 0; return(x)})

    # Initial values for eta:

    pi_hat_list <- eta_hat_list <- vector("list", length = Jmulti)

    # pi_hat_list <- apply(Y, 2, \(x) table(x)/nobs)
    for(j in 1:Jmulti) { # Much faster
      int_vector <- data[, multinom[j]]
      int_vector <- int_vector[!is.na(int_vector)] # Remove missing values
      nsize <- length(int_vector)
      props <- count(int_vector, nsize, K[j]) / nsize
      pi_hat_list[[j]] <- props %*% t(rep(1, nclasses))
      log_props <- log(props)
      eta_hat_list[[j]] <- (log_props-log_props[1]) %*% t(rep(1, nclasses))
    }

    # eta is sampled from a normal distribution:
    Ks <- length(unlist(eta_hat_list))
    for(i in 1:control$rstarts) {
      # init_eta <- rnorm(sum(Ks))
      init_eta <- unlist(eta_hat_list) + rnorm(sum(Ks), 0, 0.1)
      init_param[[i]]$eta <- fill_list_with_vector(lca_param$eta, init_eta)
      for(j in 1:Jmulti) {
        init_param[[i]]$eta[[j]][1, ] <- "0"
      }
    }

  }

  #### Replace the parameter by custom values ####

  if(!is.null(model$beta)) {
    lca_param$beta <- model$beta
  }

  if(!is.null(model$mu)) {
    lca_param$mu <- model$mu
  }

  if(!is.null(model$s)) {
    lca_param$s <- model$s
  }

  if(!is.null(model$eta)) {
    lca_param$eta <- model$eta
  }

  # Transform all the initial parameters in numerical values:
  init_param <- allnumeric(init_param)

  # Replace initial starting values by custom starting values:

  if(!is.null(control$start$beta)) {
    for(i in 1:control$rstarts) {
      init_param[[i]]$beta <- insert_object(init_param[[i]]$beta,
                                            control$start$beta)
    }
  }

  if(!is.null(control$start$mu)) {
    for(i in 1:control$rstarts) {
      init_param[[i]]$mu <- insert_object(init_param[[i]]$mu,
                                          control$start$mu)
    }
  }

  if(!is.null(control$start$s)) {
    for(i in 1:control$rstarts) {
      init_param[[i]]$s <- insert_object(init_param[[i]]$s,
                                         control$start$s)
    }
  }

  if(!is.null(control$start$eta)) {
    for(i in 1:control$rstarts) {
      init_param[[i]]$eta <- insert_object(init_param[[i]]$eta,
                                           control$start$eta)
    }
  }

  #### Arrange parameter labels ####

  vector_param <- unname(unlist(lca_param))
  indicator <- is.na(suppressWarnings(as.numeric(vector_param)))
  fixed_indices <- which(!indicator)
  fixed_values <- as.numeric(vector_param[fixed_indices])
  parameters_indices <- which(indicator)
  parameters_labels <- unique(vector_param[parameters_indices])
  nparam <- length(parameters_labels)

  #### Create the initial values for the parameters ####

  parameters <- vector("list", length = control$rstarts)
  select_params <- match(parameters_labels, vector_param)
  for(i in 1:control$rstarts) {
    parameters[[i]] <- unlist(init_param[[i]])[select_params]
    names(parameters[[i]]) <- parameters_labels
  }

  #### Model for the transformed parameters ####

  lca_all$theta <- matrix(NA, nrow = npatterns, ncol = nclasses)
  lca_all$class <- matrix(NA, nrow = npatterns, ncol = nclasses)
  colnames(lca_all$class) <- colnames(lca_all$theta) <-
    class_names

  for(s in 1:npatterns) {

    lca_all$theta[s, ] <- paste("theta", 1:nclasses, "[", s, "]", sep = "")
    lca_all$class[s, ] <- paste("class", 1:nclasses, "[", s, "]", sep = "")

  }

  rownames(lca_all$theta) <- rownames(lca_all$class) <- pattern_names

  lca_all$loglik <- array(NA, dim = c(npatterns, nitems, nclasses),
                            dimnames = list(pattern_names,
                                            item_names,
                                            class_names))

  # Model for gaussian items:
  if(any(item == "gaussian")) {

    # Model for the transformed parameters:
    repitems <- rep(gauss, times = nclasses)
    repclasses <- rep(1:nclasses, each = Jgauss)
    Srepitems <- rep(repitems, each = npatterns)
    Srepclasses <- rep(repclasses, each = npatterns)
    repS <- rep(1:npatterns, times  = Jgauss*nclasses)
    sigma <- paste("sigma[", repitems, "|", repclasses, "]", sep = "")
    loglik_gauss <- paste("loglik[", repS, ",", Srepitems, ",", Srepclasses, "]",
                          sep = "")
    loglik_gauss <- array(loglik_gauss, dim = c(npatterns, Jgauss, nclasses))
    lca_all$sigma <- matrix(sigma, nrow = Jgauss, ncol = nclasses)
    # lca_all$loglik_gauss <- loglik_gauss
    lca_all$loglik[, gauss, ] <- loglik_gauss

  }

  # Model for multinomial items:
  if(any(item == "multinomial")) {

    # Model for the transformed parameters:
    repitems <- rep(multinom, times = nclasses)
    repclasses <- rep(1:nclasses, each = Jmulti)
    Srepitems <- rep(repitems, each = npatterns)
    Srepclasses <- rep(repclasses, each = npatterns)
    repS <- rep(1:npatterns, times  = Jmulti*nclasses)
    loglik_multinom <- paste("loglik[", repS, ",", Srepitems, ",", Srepclasses, "]", sep = "")
    loglik_multinom <- array(loglik_multinom, dim = c(npatterns, Jmulti, nclasses))
    peta <- vector("list", length = Jmulti)
    for(j in 1:Jmulti) {
      peta[[j]] <- matrix(NA, nrow = K[j], ncol = nclasses)
      peta[[j]][] <- paste("peta[", 1:K[j], "|", j, "|" ,
                           rep(1:nclasses, each = K[j]), "]", sep = "")
      colnames(peta[[j]]) <- class_names
      rownames(peta[[j]]) <- data_list$factor_names[[j]]
    }
    names(peta) <- colnames(data)[multinom]
    lca_all$peta <- peta
    lca_all$loglik[, multinom, ] <- loglik_multinom

  }

  #### Create the initial values for the transformed parameters ####

  # # Arrange transparameter labels:
  # vector_trans <- unname(unlist(lca_all))
  # fixed_labels <- vector_trans[1:length(vector_param)][!indicator]
  # transparameters_labels <- unique(c(vector_trans, parameters_labels))
  # ntrans <- length(transparameters_labels)
  #
  # # Indices of free parameters:
  # free_indices <- match(parameters_labels, transparameters_labels)
  # # Indices of fixed parameters:
  # fixed_indices <- match(fixed_labels, transparameters_labels)
  #
  # # Initialize the vector of transformed parameters:
  # transparameters <- vector("list", length = control$rstarts)
  # for(i in 1:control$rstarts) { # For each random start i...
  #   transparameters[[i]] <- vector(length = ntrans)
  #   # Copy the parameter values to the transparameters:
  #   transparameters[[i]][free_indices] <- parameters[[i]]
  #   transparameters[[i]][fixed_indices] <- fixed_values # Fixed parameters
  # }

  # Replace the transformed parameter by custom values:

  lca_trans <- lca_all

  if(!is.null(model$beta)) {
    lca_trans$beta <- model$beta
  }

  if(!is.null(model$mu)) {
    lca_trans$mu <- model$mu
  }

  if(!is.null(model$s)) {
    lca_trans$s <- model$s
  }

  if(!is.null(model$eta)) {
    lca_trans$eta <- model$eta
  }

  if(!is.null(model$theta)) {
    lca_trans$theta <- model$theta
  }

  if(!is.null(model$class)) {
    lca_trans$class <- model$class
  }

  if(!is.null(model$peta)) {
    lca_trans$peta <- model$peta
  }

  if(!is.null(model$sigma)) {
    lca_trans$sigma <- model$sigma
  }

  if(!is.null(model$loglik)) {
    lca_trans$loglik <- model$loglik
  }

  # Arrange transparameter labels:
  transparameters_labels <- unname(unlist(lca_all))
  ntrans <- length(transparameters_labels)

  # Indices of fixed parameters:
  alltrans <- suppressWarnings(as.numeric(unname(unlist(lca_trans))))
  fixed_indices <- which(!is.na(alltrans))
  fixed_values <- alltrans[fixed_indices]

  # Initialize the vector of transformed parameters:
  transparameters <- vector("list", length = control$rstarts)
  for(i in 1:control$rstarts) { # For each random start i...
    transparameters[[i]] <- vector(length = ntrans)
    # Copy the parameter values to the transparameters:
    transparameters[[i]][fixed_indices] <- fixed_values # Fixed parameters
  }

  #### Relate the transformed parameters to the parameters ####

  # Relate the parameters to the transformed parameters:
  param2trans <- match(transparameters_labels, parameters_labels)
  param2trans <- param2trans[!is.na(param2trans)] # Remove NAs
  # Relate the transformed parameters to the parameters:
  trans2param <- match(parameters_labels, transparameters_labels)

  #### Set up the optimizer ####

  # Create defaults for the control of the optimizer:
  control$init_param <- init_param
  control$parameters <- parameters
  control$transparameters <- transparameters
  control$param2transparam <- param2trans-1L
  control$transparam2param <- trans2param-1L

  #### Return ####

  result <- list(parameters_labels = parameters_labels,
                 nparam = nparam,
                 transparameters_labels = transparameters_labels,
                 ntrans = ntrans,
                 lca_param = lca_param,
                 lca_trans = lca_trans,
                 lca_all = lca_all,
                 pi_hat_list = pi_hat_list,
                 init_param = init_param,
                 control = control)

  return(result)

}

get_short_lca_model <- function(data_list, nclasses, item,
                                lca_all, model = NULL) {

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

  #### Transformed parameters model ####

  if(any(item == "gaussian")) {

    for(j in gauss) {
      items[[j]] <- matrix(NA, nrow = 2, ncol = nclasses)
      rownames(items[[j]]) <- c("Means", "Stds")
      colnames(items[[j]]) <- paste("Class", 1:nclasses, sep = "")
    }
    # Fill the gaussian model:
    mu_sorted <- sort(lca_all$mu)
    sigma_sorted <- sort(lca_all$sigma)
    out <- c(rbind(mu_sorted, sigma_sorted))
    items[gauss] <- fill_list_with_vector(items[gauss], out)

  }

  if(any(item == "multinomial")) {

    i <- 1L
    for(j in multinom) {
      items[[j]] <- matrix(NA, nrow = K[[i]], ncol = nclasses)
      # rownames(items[[j]]) <- paste("Category", 1:K[i], sep = "")
      rownames(items[[j]]) <- data_list$factor_names[[i]]
      colnames(items[[j]]) <- paste("Class", 1:nclasses, sep = "")
      i <- i+1L
    }

    # Fill the multinomial model:
    items[multinom] <- fill_list_with_vector(items[multinom],
                                             unlist(lca_all$peta))

  }

  prob_model <- list()
  prob_model$beta <- lca_all$beta
  prob_model$items <- items

  #### Parameters model ####

  if(any(item == "gaussian")) {

    for(j in gauss) {
      rownames(items[[j]]) <- c("Means", "logStds")
    }

    # Fill the gaussian model:
    s_sorted <- sort(lca_all$s)
    out <- c(rbind(mu_sorted, s_sorted))
    items[gauss] <- fill_list_with_vector(items[gauss], out)

  }

  if(any(item == "multinomial")) {

    # Fill the multinomial model:
    items[multinom] <- fill_list_with_vector(items[multinom],
                                             unlist(lca_all$eta))

  }

  log_model <- list()
  log_model$beta <- lca_all$beta
  log_model$items <- items
  rownames(log_model$beta) <- colnames(data_list$cov_patterns2)

  #### Return ####

  # Return the model in the logarithm and probability scales:
  result <- list(prob_model = prob_model, log_model = log_model)

  return(result)

}

get_lca_structures <- function(data_list, full_model, control) {

  # Generate control_manifold, control_transform, and control_estimator

  list2env(data_list, envir = environment())
  list2env(full_model, envir = environment())

  nclasses <- ncol(lca_all$class)

  #### Manifolds ####

  mani_and_labs <- list(
    list("euclidean", lca_param)
    )
  control_manifold <- create_manifolds(manifolds_and_labels = mani_and_labs,
                                       param_structures = lca_param)

  #### Transformations ####

  control_transform <- list()
  k <- 1L

  # betas to thetas:
  labels_in <- lca_all$beta # px1
  labels_out <- lca_all$theta # X*beta = Sxp by px1
  indices_in <- match(labels_in, transparameters_labels)-1L
  indices_out <- match(labels_out, transparameters_labels)-1L
  control_transform[[k]] <- list(transform = "column_space",
                                 # labels_in = labels_in,
                                 # labels_out = labels_out,
                                 indices_in = indices_in,
                                 indices_out = indices_out,
                                 X = cov_patterns2)
  k <- k+1L

  # thetas to classes:
  for(s in 1:npatterns) {

    labels_in <- lca_all$theta[s, ]
    labels_out <- lca_all$class[s, ]
    indices_in <- match(labels_in, transparameters_labels)-1L
    indices_out <- match(labels_out, transparameters_labels)-1L
    control_transform[[k]] <- list(transform = "softmax",
                                   # labels_in = labels_in,
                                   # labels_out = labels_out,
                                   indices_in = indices_in,
                                   indices_out = indices_out)
    k <- k + 1L

  }

  # Conditional item likelihoods (gaussian):
  if(any(item == "gaussian")) {

    gauss <- which(item == "gaussian")
    Jgauss <- length(gauss) # Number of gaussian items
    labels_in <- lca_all$s
    labels_out <- lca_all$sigma
    indices_in <- match(labels_in, transparameters_labels)-1L
    indices_out <- match(labels_out, transparameters_labels)-1L
    control_transform[[k]] <- list(transform = "exponential",
                                   # labels_in = labels_in,
                                   # labels_out = labels_out,
                                   indices_in = indices_in,
                                   indices_out = indices_out)
    k <- k+1L

    labels_in <- c(lca_all$mu, lca_all$sigma)
    indices_in_mu <- match(lca_all$mu, transparameters_labels)-1L
    indices_in_sigma <- match(lca_all$sigma, transparameters_labels)-1L
    labels_out <- lca_all$loglik[, gauss, ]
    indices_out <- match(labels_out, transparameters_labels)-1L

    y <- as.matrix(patterns[, gauss])

    control_transform[[k]] <- list(transform = "normal",
                                   # labels_out = labels_out,
                                   indices_mu = indices_in_mu,
                                   indices_sigma = indices_in_sigma,
                                   indices_out = indices_out,
                                   y = y,
                                   S = npatterns,
                                   J = Jgauss,
                                   I = nclasses)
    k <- k+1L

  }

  # Conditional item likelihoods (multinomial):
  if(any(item == "multinomial")) {

    # All of this is redundant in get_full_lca_model:
    multinom <- which(item == "multinomial")
    Jmulti <- length(multinom) # Number of multinomial items

    # Softmax transformations:
    for(j in 1:Jmulti) {
      for(i in 1:nclasses) {

        labels_in <- lca_all$eta[[j]][, i]
        labels_out <- lca_all$peta[[j]][, i]
        indices_in <- match(labels_in, transparameters_labels)-1L
        indices_out <- match(labels_out, transparameters_labels)-1L
        control_transform[[k]] <- list(transform = "softmax",
                                       # labels_in = labels_in,
                                       # labels_out = labels_out,
                                       indices_in = indices_in,
                                       indices_out = indices_out)
        k <- k+1L

      }
    }

    # Multinomial transformation:
    labels_in <- unname(unlist(lca_all$peta))
    labels_out <- c(lca_all$loglik[, multinom, ])
    indices_in <- match(labels_in, transparameters_labels)-1L
    indices_out <- match(labels_out, transparameters_labels)-1L

    y <- as.matrix(patterns[, multinom])
    K <- unlist(lapply(data_list$factor_names, FUN = length))

    control_transform[[k]] <- list(transform = "multinomial",
                                   # labels_out = labels_out,
                                   # labels_in = labels_in,
                                   indices_in = indices_in,
                                   indices_out = indices_out,
                                   y = y,
                                   S = npatterns,
                                   J = Jmulti,
                                   I = nclasses,
                                   K = K)

  }

  #### Estimators ####

  control_estimator <- list()
  k <- 0L

  item_loglik <- lca_all$loglik
  indices_classes <- match(lca_all$class, transparameters_labels)-1L
  indices_cubeloglik <- match(item_loglik, transparameters_labels)-1L

  control_estimator[[1]] <- list(estimator = "lca",
                                 # labels = labels,
                                 indices_classes = indices_classes,
                                 indices_cubeloglik = indices_cubeloglik,
                                 S = npatterns,
                                 J = nitems,
                                 I = nclasses,
                                 weights = weights)

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
      G <- 2

      for(i in uniques) {

        labels <- lca_all$class[i, ]
        indices <- match(labels, transparameters_labels)-1L
        control_estimator[[G]] <- list(estimator = "bayesconst1",
                                       # labels = labels,
                                       indices = indices,
                                       K = nclasses,
                                       alpha = alpha,
                                       U = U,
                                       N = data_list$nobs)
        G <- G+1L

      }
    }

    G <- length(control_estimator) + 1L

    # Bayes Constant for standard deviations:
    alpha <- control$penalties$sd$alpha
    if(any(item == "gaussian") & alpha != 0) {

      Y <- data[, gauss, drop = FALSE]
      sigma_class <- split(lca_all$sigma, rep(1:nclasses, each = Jgauss))
      varshat <- apply(Y, MARGIN = 2, FUN = var, na.rm = TRUE)*(nobs-1)/nobs
      for(i in 1:nclasses) {

        labels <- sigma_class[[i]]
        indices <- match(labels, transparameters_labels)-1L
        control_estimator[[G]] <- list(estimator = "bayesconst3",
                                       # labels = labels,
                                       indices = indices,
                                       K = nclasses,
                                       varshat = varshat,
                                       alpha = alpha,
                                       N = data_list$nobs)
        G <- G+1L

      }

    }

    # Bayes Constant for multinomial probabilities:
    alpha <- control$penalties$prob$alpha
    if(any(item == "multinomial") & alpha != 0) {

      for(j in 1:Jmulti) {
        for(i in 1:nclasses) {

          pihat <- pi_hat_list[[j]][, i]
          labels <- lca_all$peta[[j]][, i]
          indices <- match(labels, transparameters_labels)-1L
          control_estimator[[G]] <- list(estimator = "bayesconst2",
                                         # labels = labels,
                                         indices = indices,
                                         K = nclasses,
                                         pihat = pihat,
                                         alpha = control$penalties$prob$alpha,
                                         N = data_list$nobs)
          G <- G+1L

        }
      }

    }

    # Gaussian regularization for coefficients:
    alpha <- control$penalties$beta$alpha
    if(alpha != 0) {

      # means <- control$means
      # sds <- control$sds
      #
      # if(is.null(means) || is.null(sds)) {
      #   stop("Please, provide means and sds for the beta regularization")
      # }

      p <- nrow(lca_all$beta)-1L
      q <- ncol(lca_all$beta)
      means <- matrix(0, nrow = p, ncol = q)
      sds <- apply(cov_patterns2[, -1], MARGIN = 2, sd, na.rm = TRUE)
      sds <- matrix(sds, nrow = p, ncol = q) / alpha

      labels <- lca_all$beta[-1, ] # Remove the intercept
      indices <- match(labels, transparameters_labels)-1L
      control_estimator[[G]] <- list(estimator = "gaussian_loglik",
                                     # labels = labels,
                                     indices = indices,
                                     means = means,
                                     sds = sds,
                                     alpha = alpha,
                                     N = data_list$nobs)
      G <- G+1L


    }

    # Ridge regularization for coefficients:
    lambda <- control$penalties$beta$lambda
    power <- control$penalties$beta$power
    if(!is.null(lambda) & !is.null(power)) {

      if(lambda != 0 && power != 0) {

        labels <- lca_all$beta[-1, ] # Remove the intercept
        indices <- match(labels, transparameters_labels)-1L
        control_estimator[[G]] <- list(estimator = "ridge",
                                       # labels = labels,
                                       indices = indices,
                                       lambda = lambda,
                                       power = power,
                                       N = data_list$nobs)
        G <- G+1L

      }
    }

  }

  #### Return ####

  result <- list(control_manifold = control_manifold,
                 control_transform = control_transform,
                 control_estimator = control_estimator)

  return(result)

}
