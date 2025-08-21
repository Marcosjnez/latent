# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 17/08/2025
#'
#' @title
#' Get the default model for Latent Class Analysis.
#' @description
#'
#' Get the default model for Latent Class Analysis.
#'
#' @usage
#'
#' setup_lca(data, item, nclasses, model, control)
#'
#' @param data data.frame or matrix of response.
#' @param item Character vector with the model for each item.
#' @param nclasses Number of latent classes.
#' @param model List of parameter labels. See 'details' for more information.
#' @param constraints Should the model be checked for identification? Defaults to TRUE.
#'
#' @details \code{get_lca} generates the model for the probability of belonging
#' to the classes and the conditional response probabilities. These models may
#' be modified by the user to set equality constraints or to fix parameters.
#'
#' @return List with the following objects:
#' \item{none}{.}
#' \item{none}{.}
#'
#' @references
#'
#' None yet.
#'
#' @export
setup_lca <- function(data, nclasses = 2L, item = rep("gaussian", ncol(data)),
                      model = NULL, control = NULL) {

  #### Initial input checks ####

  # Check that data is either a data.frame or a matrix:
  if(!is.data.frame(data) & !is.matrix(data)) {
    stop("data must be a matrix or data.frame")
  }

  # Number of items:
  nitems <- ncol(data)

  # Check that item is a character vector with a string for each column of data:
  if(!is.character(item) || length(item) != nitems) {

    stop("item must be a character vector with as many elements as columns in data")

  }

  #### Process the data ####

  # Number of subjects:
  nobs <- nrow(data)
  # Convert data to a data.table object:
  dt <- data.table::as.data.table(data)
  ## Collect some information from the data ##
  counts_dt <- dt[, .(index = .I[1], count = .N), by = names(dt)]
  # Data matrix with the unique response patterns:
  patterns <- as.matrix(counts_dt[, names(dt), with = FALSE])
  # Number of unique response patterns:
  npatterns <- nrow(counts_dt)
  # Counts of each response pattern:
  weights <- counts_dt$count
  # Indices to map the original data to the matrix of unique patterns:
  full2short <- counts_dt$index
  # Indices to map the matrix of unique patterns to the original data:
  short2full <- match(do.call(paste, dt),
                      do.call(paste, counts_dt[, -c("index", "count"), with = FALSE]))

  # Put in a list the objects generated form the data:
  data_list <- vector("list")
  data_list$dt <- dt
  data_list$patterns <- patterns
  data_list$npatterns <- npatterns
  data_list$nitems <- nitems
  data_list$weights <- weights
  data_list$full2short <- full2short
  data_list$short2full <- short2full

  #### Create the model ####

  # Get the model specification:
  full_model <- get_full_lca_model(data_list = data_list, item = item,
                                   nclasses = nclasses, lca_trans = lca_trans,
                                   model = model, control = control)
  list2env(full_model, envir = environment())

  # Get the short model specification (in logarithm and probability scale) with
  # labels for each parameter:
  short_model <- get_short_lca_model(data = data, item = item, nclasses = nclasses,
                                     lca_trans = lca_trans, model = model)

  #### Create the structures ####

  # Generate the structures for optimization:
  structures <- get_lca_structures(data_list = data_list,
                                   full_model = full_model)

  #### Return ####

  result <- list(full_model = full_model,
                 short_model = short_model,
                 structures = structures)

  return(result)

}

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

  # Transform all elements in a list into numeric values
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
    lca_param$mu <- mu
    lca_param$s <- s
    lca_trans$mu <- mu
    lca_trans$s <- s

    # Initial values for mu and s:
    # For mu, they will be the mean of the items
    # For s, they will be the sd of the items
    init_mu <- rep(colMeans(data[, gauss], na.rm = TRUE), times = nclasses)
    init_sd <- rep(apply(data[, gauss], MARGIN = 2, FUN = sd, na.rm = TRUE),
                  times = nclasses)
    init_s <- log(init_sd)
    for(i in 1:control$rstarts) {
      init_mui <- rnorm(Jgauss*nclasses, init_mu, init_sd/sqrt(nobs))
      init_param[[i]] <- c(init_param[[i]], init_mui, init_s)
    }

  }

  # Model for multinomial items:
  if(any(item == "multinomial")) {

    multinom <- which(item == "multinomial")
    Jmulti <- length(multinom) # Number of multinomial items
    # Collect the number of response categories:
    K <- apply(data[, multinom], MARGIN = 2, FUN = \(x) length(unique(x)))
    Ks <- rep(K, times = nclasses)
    repitems <- rep(multinom, times = nclasses)
    repclasses <- rep(1:nclasses, each = Jmulti)
    repitemsK <- rep(repitems, times = Ks)
    repclassesK <- rep(repclasses, times = Ks)
    KsJ <- lapply(Ks, FUN = \(x) 1:x)
    eta <- paste("eta[", repitemsK, "|", repclassesK, "|", unlist(KsJ), "]",
                 sep = "")
    eta <- split(eta, rep(seq_along(Ks), Ks))
    lca_trans$eta <- eta
    eta <- lapply(eta, function(x) { x[1] <- 0; x })
    lca_param$eta <- eta

    # Initial values for eta:
    # eta is sampled from a normal distribution:
    for(i in 1:control$rstarts) {
      init_eta <- rnorm(sum(Ks))
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
    lca_trans$sigma <- sigma
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
    peta <- paste("peta[", repitemsK, "|", repclassesK, "|", unlist(KsJ), "]", sep = "")
    peta <- split(peta, rep(seq_along(Ks), Ks))
    lca_trans$peta <- peta
    # lca_trans$loglik_multinom <- loglik_multinom
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

  # Create the initial values for the parameters:
  parameters <- vector("list", length = control$rstarts)
  select_params <- match(parameters_labels, vector_param)
  for(i in 1:control$rstarts) {
    parameters[[i]] <- init_param[[i]][select_params]
    names(parameters[[i]]) <- parameters_labels
  }

  # Create the initial values for the transformed parameters:
  free_indices <- match(parameters_labels, transparameters_labels)
  fixed_indices <- match(fixed_labels, transparameters_labels)
  transparameters <- vector("list", length = control$rstarts)
  for(i in 1:control$rstarts) {
    transparameters[[i]] <- vector(length = ntrans)
    transparameters[[i]][free_indices] <- parameters[[i]]
    transparameters[[i]][fixed_indices] <- fixed_values
  }

  # Relate the transformed parameters to the parameters:
  param2trans <- match(transparameters_labels, parameters_labels)
  param2trans <- param2trans[!is.na(param2trans)]
  # Relate the parameters to the transformed parameters:
  trans2param <- match(parameters_labels, transparameters_labels)

  #### Set up the optimizer ####

  # Create defaults for the control of the optimizer:
  control <- lca_control(control)
  control$parameters <- parameters
  control$transparameters <- transparameters
  control$param2transparam <- param2trans-1L
  control$transparam2param <- trans2param-1L

  #### Return ####

  result <- list(parameters_labels = parameters_labels,
                 # parameters = parameters,
                 nparam = nparam,
                 transparameters_labels = transparameters_labels,
                 # transparameters = transparameters,
                 ntrans = ntrans,
                 lca_param = lca_param,
                 lca_trans = lca_trans,
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
  K <- apply(data[, multinom], MARGIN = 2, FUN = \(x) length(unique(x)))
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
      items[[j]] <- matrix(NA, nrow = K[i], ncol = nclasses)
      # rownames(items[[j]]) <- paste("Category", 1:K[i], sep = "")
      rownames(items[[j]]) <- data_list$factor_names[[j]]
      colnames(items[[j]]) <- lca_trans$class
      i <- i+1L
    }

    # Fill the multinomial model:
    peta_sorted <- sort(unlist(lca_trans$peta))
    items[multinom] <- fill_list_with_vector(items[multinom], peta_sorted)

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
    eta_sorted <- sort(unlist(lca_trans$eta))
    items[multinom] <- fill_list_with_vector(items[multinom], eta_sorted)

  }

  log_model <- list()
  log_model$classes <- lca_trans$theta
  log_model$items <- items

  #### Return ####

  # Return the model in the logarithm and probability scales:
  result <- list(prob_model = prob_model, log_model = log_model)

  return(result)

}

get_lca_structures <- function(data_list, full_model) {

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
    labels_out <- lca_trans$loglik[, gauss, ]
    indices_in <- match(labels_in, transparameters_labels)
    indices_out <- match(labels_out, transparameters_labels)
    mu_indices <- match(lca_trans$mu, transparameters_labels[indices_in])
    mu_indices <- rep(mu_indices, each = npatterns)
    sigma_indices <- match(lca_trans$sigma, transparameters_labels[indices_in])
    sigma_indices <- rep(sigma_indices, each = npatterns)
    indices_in <- list(indices_in-1L, mu_indices-1L, sigma_indices-1L)
    y <- rep(c(patterns[, gauss]), times = nclasses)
    control_transform[[k]] <- list(transform = "normal",
                                   labels_in = labels_in,
                                   indices_in = indices_in,
                                   labels_out = labels_out,
                                   indices_out = list(indices_out-1L),
                                   y = y)
    k <- k+1L

  }

  if(any(item == "multinomial")) {

    # All of this is redundant in get_full_lca_model:
    multinom <- which(item == "multinomial")
    Jmulti <- length(multinom) # Number of multinomial items
    K <- apply(data[, multinom], MARGIN = 2, FUN = \(x) length(unique(x)))
    Ks <- rep(K, times = nclasses)

    for(i in 1:(Jmulti*nclasses)) {

      labels_in <- lca_trans$eta[[i]]
      labels_out <- lca_trans$peta[[i]]
      indices_in <- match(labels_in, transparameters_labels)
      indices_out <- match(labels_out, transparameters_labels)
      control_transform[[k]] <- list(transform = "softmax",
                                     labels_in = labels_in,
                                     indices_in = list(indices_in-1L),
                                     labels_out = labels_out,
                                     indices_out = list(indices_out-1L))
      k <- k+1L

    }

    labels_in <- unname(unlist(lca_trans$peta))
    labels_out <- c(lca_trans$loglik[, multinom, ])
    indices_in <- match(labels_in, transparameters_labels) # JxIxK
    indices_out <- match(labels_out, transparameters_labels) # SxJxI
    p_indices <- match(labels_in, transparameters_labels[indices_in])
    cumsumKs <- c(0, cumsum(Ks))
    # For each datapoint, create a vector of possible peta contributions to the likelihood:
    peta_indices <- c()
    for(ks in 1:length(Ks)) {
      selection <- (cumsumKs[ks]+1L):cumsumKs[ks+1L]
      peta_indices <- c(peta_indices, rep(p_indices[selection], times = npatterns))
    }
    df <- as.data.frame(patterns[, multinom])
    # Create a dummy vector indicating which peta contributes to the likelihood:
    dummy <- c()
    for(j in 1:Jmulti) {
      dummy <- c(dummy, c(t(model.matrix(~ factor(df[,j])-1))))
    }
    dummy <- rep(dummy, times = nclasses)
    dummy_indices <- which(as.logical(dummy))
    # dummy <- unname(cumsum(rep(K, each = npatterns)) - K[1] + c(patterns[, multinom]))
    # dummy_indices <- rep(dummy, Jmulti) + rep(seq(0L, Jmulti-1L) * (max(dummy) + 1L),
    #                                          each = npatterns*Jmulti)
    # petas contributing to the likelihood for each (by class):
    peta_indices <- peta_indices[dummy_indices]
    indices_in <- list(indices_in-1L, peta_indices-1L)

    control_transform[[k]] <- list(transform = "multinomial",
                                   labels_in = labels_in,
                                   indices_in = indices_in,
                                   labels_out = labels_out,
                                   indices_out = list(indices_out-1L))

  }

  #### Estimators ####

  indices <- list()
  full_loglik <- c(lca_trans$loglik)
  all_indices <- match(c(lca_trans$class, full_loglik), transparameters_labels)
  indices_classes <- match(lca_trans$class, transparameters_labels[all_indices])
  indices_items <- match(full_loglik, transparameters_labels[all_indices])
  indices[[1]] <- all_indices-1L
  indices[[2]] <- indices_classes-1L
  indices[[3]] <- indices_items-1L
  labels <- transparameters_labels[all_indices]
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

  #### Return ####

  result <- list(control_manifold = control_manifold,
                 control_transform = control_transform,
                 control_estimator = control_estimator)

}

# indices <- control_transform[[11]]$indices_in[[2]]+1L
# cbind(control_transform[[11]]$labels_in[indices],
# control_transform[[11]]$labels_out)

