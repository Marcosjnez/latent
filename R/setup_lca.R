# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 08/09/2025

fill_list_with_vector <- function(lst, values) {

# Insert a vector of values in an arbitrary list lst

i <- 1

assign_recursive <- function(x) {
  if (is.list(x)) {
    lapply(x, assign_recursive)
  } else if (is.matrix(x)) {
    dims <- dim(x)
    n <- prod(dims)
    x[] <- values[i:(i + n - 1)]
    i <<- i + n
    x
  } else if (is.atomic(x)) {
    n <- length(x)
    x[] <- values[i:(i + n - 1)]
    i <<- i + n
    x
  } else {
    stop("Unsupported type")
  }
}

assign_recursive(lst)

}

allnumeric <- function(lst) {

  # Transform all elements of a list into numeric values
  lst <- rapply(
    lst,
    function(x) {
      if (is.matrix(x)) { storage.mode(x) <- "double"; x }           # keep as matrix
      else if (is.factor(x)) as.numeric(as.character(x))             # factors -> numeric values
      else if (is.atomic(x)) as.numeric(x)                           # vectors -> numeric
      else x
    },
    how = "replace",
    classes = c("matrix","array","factor","numeric","integer","logical","character")
  )

  return(lst)

}

get_full_lca_model <- function(data_list, nclasses, item, model = NULL,
                               control = NULL) {

  # Generate the model syntax and initial parameter values

  list2env(data_list, envir = environment())

  # Initialize the objects to store the initial parameters:
  lca_param <- list()
  lca_trans <- list()

  #### Model for the nontransformed parameters ####

  theta <- paste("theta", 1:nclasses, sep = "")
  lca_trans$theta <- theta
  theta[1] <- "0"
  lca_param$theta <- theta

  # Initial values for the log class probabilities (theta):
  # theta are sampled from a normal distribution:
  init_param <- vector("list", length = control$rstarts)
  for(i in 1:control$rstarts) {
    init_theta <- rnorm(nclasses)
    init_param[[i]] <- init_theta
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
    lca_trans$mu <- lca_param$mu
    lca_trans$s <- lca_param$s

    # Initial values for mu and s:
    # For mu, they will be the mean of the items
    # For s, they will be the sd of the items
    init_mu <- rep(colMeans(data[, gauss, drop = FALSE], na.rm = TRUE), times = nclasses)
    init_sd <- rep(apply(data[, gauss, drop = FALSE], MARGIN = 2, FUN = sd, na.rm = TRUE),
                  times = nclasses)
    init_s <- log(init_sd)
    for(i in 1:control$rstarts) {
      init_mui <- rnorm(Jgauss*nclasses, init_mu, init_sd/sqrt(nobs))
      init_param[[i]] <- c(init_param[[i]], init_mui, init_s)
    }

  }

  pi_hat_list <- list()
  # Model for multinomial items:
  if(any(item == "multinomial")) {

    multinom <- which(item == "multinomial")
    Jmulti <- length(multinom) # Number of multinomial items
    # Collect the number of response categories:
    K <- unlist(lapply(data_list$factor_names, FUN = length))
    eta <- vector("list", length = Jmulti)
    for(j in 1:Jmulti) {
      eta[[j]] <- matrix(NA, nrow = K[j], ncol = nclasses)
      eta[[j]][] <- paste("eta[", j, "|" , rep(1:nclasses, each = K[j]),
                          "|", 1:K[j], "]", sep = "")
    }
    lca_trans$eta <- eta
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
      init_param[[i]] <- c(init_param[[i]], init_eta)
    }

  }

  #### Model for the transformed parameters ####

  class <- paste("class", 1:nclasses, sep = "")
  lca_trans$class <- class
  lca_trans$loglik <- array(NA, dim = c(npatterns, nitems, nclasses))

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
    lca_trans$sigma <- matrix(sigma, nrow = Jgauss, ncol = nclasses)
    # lca_trans$loglik_gauss <- loglik_gauss
    lca_trans$loglik[, gauss, ] <- loglik_gauss

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
      peta[[j]][] <- paste("peta[", j, "|" , rep(1:nclasses, each = K[j]),
                           "|", 1:K[j], "]", sep = "")
    }
    lca_trans$peta <- peta
    lca_trans$loglik[, multinom, ] <- loglik_multinom

  }

  #### Arrange labels ####

  # Arrange parameter labels:
  vector_param <- unname(unlist(lca_param))
  indicator <- is.na(suppressWarnings(as.numeric(vector_param)))
  fixed_indices <- which(!indicator)
  fixed_values <- as.numeric(vector_param[fixed_indices])
  parameters_indices <- which(indicator)
  parameters_labels <- unique(vector_param[parameters_indices])
  nparam <- length(parameters_labels)

  # Arrange transparameter labels:
  vector_trans <- unname(unlist(lca_trans))
  fixed_labels <- vector_trans[1:length(vector_param)][!indicator]
  transparameters_labels <- unique(c(parameters_labels, vector_trans))
  ntrans <- length(transparameters_labels)

  #### Create the initial values for the parameters ####

  parameters <- vector("list", length = control$rstarts)
  select_params <- match(parameters_labels, vector_param)
  for(i in 1:control$rstarts) {
    parameters[[i]] <- init_param[[i]][select_params]
    names(parameters[[i]]) <- parameters_labels
  }

  #### Create the initial values for the transformed parameters ####

  free_indices <- match(parameters_labels, transparameters_labels)
  fixed_indices <- match(fixed_labels, transparameters_labels)
  transparameters <- vector("list", length = control$rstarts)
  for(i in 1:control$rstarts) {
    transparameters[[i]] <- vector(length = ntrans)
    transparameters[[i]][free_indices] <- parameters[[i]]
    transparameters[[i]][fixed_indices] <- fixed_values
  }

  #### Relate the transformed parameters to the parameters ####

  param2trans <- match(transparameters_labels, parameters_labels)
  param2trans <- param2trans[!is.na(param2trans)]
  # Relate the parameters to the transformed parameters:
  trans2param <- match(parameters_labels, transparameters_labels)

  #### Set up the optimizer ####

  # Create defaults for the control of the optimizer:
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
                 pi_hat_list = pi_hat_list,
                 control = control)

  return(result)

}

get_short_lca_model <- function(data_list, nclasses, item, lca_trans,
                                model = NULL) {

  # This function displays the reduced LCA model in logarithm and probability
  # scales. This is a short version of the full model syntax created with
  # get_full_lca_model.

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
      rownames(items[[j]]) <- c("Means", "Sigma")
      colnames(items[[j]]) <- paste("Class", 1:nclasses, sep = "")
    }
    # Fill the gaussian model:
    mu_sorted <- sort(lca_trans$mu)
    sigma_sorted <- sort(lca_trans$sigma)
    out <- c(rbind(mu_sorted, sigma_sorted))
    items[gauss] <- fill_list_with_vector(items[gauss], out)

  }

  if(any(item == "multinomial")) {

    i <- 1L
    for(j in multinom) {
      items[[j]] <- matrix(NA, nrow = K[[i]], ncol = nclasses)
      # rownames(items[[j]]) <- paste("Category", 1:K[i], sep = "")
      rownames(items[[j]]) <- data_list$factor_names[[i]]
      colnames(items[[j]]) <- lca_trans$class
      i <- i+1L
    }

    # Fill the multinomial model:
    items[multinom] <- fill_list_with_vector(items[multinom],
                                             unlist(lca_trans$peta))

  }

  prob_model <- list()
  prob_model$classes <- lca_trans$class
  prob_model$items <- items

  #### Parameters model ####

  if(any(item == "gaussian")) {

    for(j in gauss) {
      rownames(items[[j]]) <- c("Means", "logSigma")
    }

    # Fill the gaussian model:
    s_sorted <- sort(lca_trans$s)
    out <- c(rbind(mu_sorted, s_sorted))
    items[gauss] <- fill_list_with_vector(items[gauss], out)

  }

  if(any(item == "multinomial")) {

    # Fill the multinomial model:
    items[multinom] <- fill_list_with_vector(items[multinom],
                                             unlist(lca_trans$eta))

  }

  log_model <- list()
  log_model$classes <- lca_trans$theta
  log_model$items <- items

  #### Return ####

  # Return the model in the logarithm and probability scales:
  result <- list(prob_model = prob_model, log_model = log_model)

  return(result)

}

get_lca_structures <- function(data_list, full_model, control) {

  # Generate control_manifold, control_transform, and control_estimator

  list2env(data_list, envir = environment())
  list2env(full_model, envir = environment())

  nclasses <- length(lca_trans$class)

  #### Manifolds ####

  control_manifold <- list()
  indices <- 1:nparam
  control_manifold[[1]] <- list(manifold = "euclidean",
                                labels = parameters_labels,
                                indices = list(indices-1L))
  #### Transformations ####

  control_transform <- list()
  labels_in <- lca_trans$theta
  labels_out <- lca_trans$class
  indices_in <- match(labels_in, transparameters_labels)
  indices_out <- match(labels_out, transparameters_labels)
  control_transform[[1]] <- list(transform = "softmax",
                                 labels_in = labels_in,
                                 indices_in = list(indices_in-1L),
                                 labels_out = labels_out,
                                 indices_out = list(indices_out-1L))
  k <- 2L

  if(any(item == "gaussian")) {

    gauss <- which(item == "gaussian")
    Jgauss <- length(gauss) # Number of gaussian items
    labels_in <- lca_trans$s
    labels_out <- lca_trans$sigma
    indices_in <- match(labels_in, transparameters_labels)
    indices_out <- match(labels_out, transparameters_labels)
    control_transform[[k]] <- list(transform = "exponential",
                                   labels_in = labels_in,
                                   indices_in = list(indices_in-1L),
                                   labels_out = labels_out,
                                   indices_out = list(indices_out-1L))
    k <- k+1L

    labels_in <- c(lca_trans$mu, lca_trans$sigma)
    indices_in <- match(labels_in, transparameters_labels)
    indices_in_mu <- match(lca_trans$mu, transparameters_labels)
    indices_in_sigma <- match(lca_trans$sigma, transparameters_labels)

    labels_out <- lca_trans$loglik[, gauss, ]
    indices_out <- match(labels_out, transparameters_labels)

    indices_in <- list(indices_in-1L, indices_in_mu-1L, indices_in_sigma-1L)
    indices_out <- list(indices_out-1L)

    y <- as.matrix(patterns[, gauss])

    control_transform[[k]] <- list(transform = "normal",
                                   labels_in = labels_in,
                                   indices_in = indices_in,
                                   labels_out = labels_out,
                                   indices_out = indices_out,
                                   y = y,
                                   S = npatterns,
                                   J = Jgauss,
                                   I = nclasses)
    k <- k+1L

  }

  if(any(item == "multinomial")) {

    # All of this is redundant in get_full_lca_model:
    multinom <- which(item == "multinomial")
    Jmulti <- length(multinom) # Number of multinomial items

    # Softmax transformations:
    for(j in 1:Jmulti) {
      for(i in 1:nclasses) {

        labels_in <- lca_trans$eta[[j]][, i]
        labels_out <- lca_trans$peta[[j]][, i]
        indices_in <- match(labels_in, transparameters_labels)
        indices_out <- match(labels_out, transparameters_labels)
        control_transform[[k]] <- list(transform = "softmax",
                                       labels_in = labels_in,
                                       indices_in = list(indices_in-1L),
                                       labels_out = labels_out,
                                       indices_out = list(indices_out-1L))
        k <- k+1L

      }
    }

    # Multinomial transformation:
    labels_in <- unname(unlist(lca_trans$peta))
    labels_out <- c(lca_trans$loglik[, multinom, ])
    indices_in <- match(labels_in, transparameters_labels)
    indices_out <- match(labels_out, transparameters_labels)
    indices_eta <- match(unname(unlist(lca_trans$eta)),
                         transparameters_labels)

    indices_in <- list(indices_in-1L, indices_eta-1L)
    indices_out <- list(indices_out-1L)
    y <- as.matrix(patterns[, multinom])
    K <- unlist(lapply(data_list$factor_names, FUN = length))

    control_transform[[k]] <- list(transform = "multinomial",
                                   labels_in = labels_in,
                                   indices_in = indices_in,
                                   labels_out = labels_out,
                                   indices_out = indices_out,
                                   y = y,
                                   S = npatterns,
                                   J = Jmulti,
                                   I = nclasses,
                                   K = K)

  }

  #### Estimators ####

  item_loglik <- c(lca_trans$loglik)

  labels <- c(lca_trans$class, item_loglik)
  all_indices <- match(labels, transparameters_labels)
  indices_classes <- match(lca_trans$class, transparameters_labels)
  indices_items <- match(item_loglik, transparameters_labels)
  indices_theta <- match(lca_trans$theta, transparameters_labels)
  indices <- list(all_indices-1L, indices_classes-1L, indices_items-1L,
                  indices_theta-1L)
  SJ <- npatterns*nitems
  hess_indices <- lapply(0:(nclasses - 1), function(i) {
    nclasses + seq(1 + i * SJ, (i + 1) * SJ)-1L })

  control_estimator <- list()
  control_estimator[[1]] <- list(estimator = "lca",
                                 labels = labels,
                                 indices = indices,
                                 S = npatterns,
                                 J = nitems,
                                 I = nclasses,
                                 weights = weights,
                                 hess_indices = hess_indices)

  # Choose whether using Bayes constants:
  if(control$reg) {

    # Bayes Constant for class probabilities:
    alpha <- control$penalties$class$alpha
    if(alpha != 0) {

      labels <- lca_trans$class
      indices <- match(labels, transparameters_labels)
      control_estimator[[2]] <- list(estimator = "bayesconst1",
                                     labels = labels,
                                     indices = list(indices-1L),
                                     K = nclasses,
                                     alpha = alpha)
    }

    G <- length(control_estimator) + 1L

    # Bayes Constant for standard deviations:
    alpha <- control$penalties$sd$alpha
    if(any(item == "gaussian") & alpha != 0) {

      Y <- data[, gauss, drop = FALSE]
      sigma_class <- split(lca_trans$sigma, rep(1:nclasses, each = Jgauss))
      varshat <- apply(Y, MARGIN = 2, FUN = var, na.rm = TRUE)*(nobs-1)/nobs
      for(i in 1:nclasses) {

        labels <- sigma_class[[i]]
        indices <- match(labels, transparameters_labels)
        control_estimator[[G]] <- list(estimator = "bayesconst3",
                                       labels = labels,
                                       indices = list(indices-1L),
                                       K = nclasses,
                                       varshat = varshat,
                                       alpha = alpha)
        G <- G+1L

      }

    }

    # Bayes Constant for multinomial probabilities:
    alpha <- control$penalties$prob$alpha
    if(any(item == "multinomial") & alpha != 0) {

      for(j in 1:Jmulti) {
        for(i in 1:nclasses) {

          pihat <- pi_hat_list[[j]][, i]
          labels <- lca_trans$peta[[j]][, i]
          indices <- match(labels, transparameters_labels)
          control_estimator[[G]] <- list(estimator = "bayesconst2",
                                         labels = labels,
                                         indices = list(indices-1L),
                                         K = nclasses,
                                         pihat = pihat,
                                         alpha = control$penalties$prob$alpha)
          G <- G+1L

        }
      }

    }
  }

  #### Return ####

  result <- list(control_manifold = control_manifold,
                 control_transform = control_transform,
                 control_estimator = control_estimator)

  return(result)

}

get_lca_covariate_model <- function(fit, X) {

  p <- ncol(X)

  thetas <- fit@modelInfo$model$classes


}

# ng <- obj@Data@ngroups
# tech <- lavaan::lavTech(obj)
# par_mat <- lavaan::lavMatrixRepresentation(lavaan::partable(obj),
#                                            representation = "LISREL")
