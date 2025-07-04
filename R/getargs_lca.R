# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 22/05/2025

getargs_lca <- function (data_list, item, model, control)
{
  Y <- data_list$Y
  n <- data_list$n
  uniq_indices <- data_list$uniq_indices
  map2full <- data_list$map2full
  dt <- data_list$dt
  estimators <- unique(item)
  estimators[estimators == "gaussian"] <- "lca_gaussian"
  estimators[estimators == "multinomial"] <- "lca_multinomial"
  nestimators <- length(estimators)
  arguments <- vector("list", length = nestimators)
  transform_setup <- estimator_setup <- list()
  k <- 1
  ep1 <- extract_param(model$log_model)
  parameter_vector <- ep1$parameter_vector
  indices_param_vector <- ep1$indices_param_vector
  indices_param_vector2 <- ep1$indices_param_vector2
  full_parameter_vector <- ep1$full_parameter_vector
  full_fixed_vector <- ep1$full_fixed_vector
  fixed_vector <- ep1$fixed_vector
  indices_full_param_vector <- ep1$indices_full_param_vector
  indices_full_fixed_vector <- ep1$indices_full_fixed_vector
  nparam <- length(parameter_vector)
  ep2 <- extract_param(model$prob_model)
  transparameter_vector <- ep2$parameter_vector
  transfixed_vector <- ep2$fixed_vector
  full_transparameter_vector <- ep2$full_parameter_vector
  full_transfixed_vector <- ep2$full_fixed_vector
  indices_full_transparam_vector <- ep2$indices_full_param_vector
  indices_full_transfixed_vector <- ep2$indices_full_fixed_vector
  indices_transparam_vector2 <- ep2$indices_param_vector2
  ntransparam <- length(transparameter_vector)
  gauss <- "gaussian" %in% item
  categ <- "multinomial" %in% item
  classes <- model$prob_model$classes
  conditionals <- model$prob_model$conditionals
  id_classes <- model$log_model$classes
  id_conditionals <- model$log_model$conditionals
  nclasses <- length(classes)
  init <- vector("list", length = control$rstarts)
  for (i in 1:control$rstarts) {
    init[[i]] <- model$log_model
    init[[i]]$classes <- rnorm(nclasses)
  }
  indices <- 1:nparam
  manifold_setup <- list()
  manifold_setup[[1]] <- list(manifold = "euclidean", indices = indices -
                                1L)
  classes_indices <- which(parameter_vector %in% id_classes)
  classes_transindices <- which(transparameter_vector %in%
                                  classes)
  X <- suppressWarnings(as.numeric(id_classes))
  vector_indices <- which(is.na(X))
  transform_setup[[1]] <- list(transform = "softmax", indices = classes_indices -
                                 1L, target_indices = classes_transindices - 1L, vector_indices = vector_indices -
                                 1L, X = X)
  if (gauss) {
    indices <- which(item == "gaussian")
    normal_data <- Y[, item == "gaussian", drop = FALSE]
    estimator_setup[[k]] <- lca_gaussian(normal_data, n,
                                         uniq_indices, map2full, classes, conditionals[indices],
                                         transparameter_vector, transfixed_vector)
    mean_labels <- unlist(lapply(conditionals[indices],
                                 FUN = function(x) x[1, ]))
    id_mean_labels <- unlist(lapply(id_conditionals[indices],
                                    FUN = function(x) x[1, ]))
    id_mean_indices <- which(parameter_vector %in% id_mean_labels)
    mean_transindices <- which(transparameter_vector %in%
                                 mean_labels)
    X <- suppressWarnings(as.numeric(id_mean_labels))
    vector_indices <- which(is.na(X))
    transform_setup[[2]] <- list(transform = "identity",
                                 indices = id_mean_indices - 1L, target_indices = mean_transindices -
                                   1L, vector_indices = vector_indices - 1L, X = X)
    sd_labels <- unlist(lapply(conditionals[indices], FUN = function(x) x[2,
    ]))
    id_sd_labels <- unlist(lapply(id_conditionals[indices],
                                  FUN = function(x) x[2, ]))
    id_sd_indices <- which(parameter_vector %in% id_sd_labels)
    sd_transindices <- which(transparameter_vector %in%
                               sd_labels)
    X <- suppressWarnings(as.numeric(id_sd_labels))
    vector_indices <- which(is.na(X))
    transform_setup[[3]] <- list(transform = "exponential",
                                 indices = id_sd_indices - 1L, target_indices = sd_transindices -
                                   1L, vector_indices = vector_indices - 1L, X = X)
    means <- rep(colMeans(normal_data, na.rm = TRUE), each = nclasses)
    if (control$opt == "em") {
      sds <- rep(apply(normal_data, MARGIN = 2, var, na.rm = TRUE),
                 each = nclasses)
    }
    else {
      sds <- rep(log(apply(normal_data, MARGIN = 2, sd, na.rm = TRUE)),
                 each = nclasses)
    }
    init_gaussian_conditionals <- c(rbind(means, sds))
    gaussian_conditionals <- fill_list_with_vector(conditionals[indices],
                                                   init_gaussian_conditionals)
    for (i in 1:control$rstarts) {
      init[[i]]$conditionals[indices] <- convert_all_to_numeric(gaussian_conditionals)
    }
    k <- k + 1
  }
  if (categ) {
    indices <- which(item == "multinomial")
    categorical_data <- Y[, item == "multinomial", drop = FALSE]
    estimator_setup[[k]] <- lca_multinomial(categorical_data,
                                            n, uniq_indices, map2full, classes, conditionals[indices],
                                            transparameter_vector, transfixed_vector)
    j <- length(transform_setup)
    for (i in 1:length(indices)) {
      id_res <- apply(id_conditionals[indices][[i]], MARGIN = 2,
                      FUN = function(x) which(parameter_vector %in%
                                                x), simplify = "list")
      transres <- apply(conditionals[indices][[i]], MARGIN = 2,
                        FUN = function(x) which(transparameter_vector %in%
                                                  x), simplify = "list")
      for (g in 1:length(id_res)) {
        j <- j + 1
        X <- suppressWarnings(as.numeric(id_conditionals[indices][[i]][,
                                                                       g]))
        vector_indices <- which(is.na(X))
        transform_setup[[j]] <- list(transform = "softmax",
                                     indices = id_res[[g]] - 1L, target_indices = transres[[g]] -
                                       1L, vector_indices = vector_indices - 1L,
                                     X = X)
      }
    }
    nelem_multinomial <- sum(unlist(lapply(id_conditionals[indices],
                                           FUN = function(x) prod(dim(x)))))
    for (i in 1:control$rstarts) {
      init_multinomial_conditionals <- rnorm(nelem_multinomial)
      multinomial_conditionals <- fill_list_with_vector(conditionals[indices],
                                                        init_multinomial_conditionals)
      init[[i]]$conditionals[indices] <- convert_all_to_numeric(multinomial_conditionals)
    }
    k <- k + 1
  }
  duplications <- duplicated(transform_setup)
  if (any(duplications)) {
    remove <- -which(duplications)
    transform_setup <- transform_setup[remove]
  }
  estimator_setup[[k]] <- list()
  estimator_setup[[k]]$estimator <- "latentloglik_combination"
  estimator_setup[[k]]$S <- estimator_setup[[k - 1]]$S
  estimator_setup[[k]]$n <- n
  estimator_setup[[k]]$indices <- which(transparameter_vector %in%
                                          classes) - 1L
  estimator_setup[[k]]$nclasses <- nclasses
  estimator_setup[[k]]$classes <- estimator_setup[[k - 1]]$classes
  estimator_setup[[k]]$indices_classes <- estimator_setup[[k -
                                                             1]]$indices_classes
  estimator_setup[[k]]$indices_target_classes <- estimator_setup[[k -
                                                                    1]]$indices_target_classes
  parameters <- transparameters <- posterior <- vector("list",
                                                       length = control$rstarts)
  S <- nrow(Y)
  init_vector <- unlist(init)
  init_vector[indices_full_fixed_vector] <- full_fixed_vector
  init <- fill_list_with_vector(init, init_vector)
  init <- convert_all_to_numeric(init)
  nid_classes <- length(classes_indices)
  init_vector <- unlist(init)[indices_param_vector]
  if ("parameters" %in% names(control)) {
    parameters <- control$parameters
  }
  else {
    for (i in 1:control$rstarts) {
      init_vector[1:nid_classes] <- rnorm(nid_classes)
      parameters[[i]] <- init_vector
    }
  }
  if ("transparameters" %in% names(control)) {
    transparameters <- control$transparameters
  }
  else {
    for (i in 1:control$rstarts) {
      transparameters[[i]] <- vector(length = ntransparam)
    }
  }
  if ("posterior" %in% names(control)) {
    posterior <- control$posterior
  }
  else {
    for (i in 1:control$rstarts) {
      posterior[[i]] <- matrix(runif(S * nclasses), nrow = S,
                               ncol = nclasses)
      posterior[[i]] <- posterior[[i]]/rowSums(posterior[[i]])
    }
  }
  control_setup <- control
  control_setup$parameters <- parameters
  control_setup$transparameters <- transparameters
  control_setup$posterior <- posterior
  control_setup$S <- estimator_setup[[1]]$S
  control_setup$nlatent <- length(estimator_setup[[1]]$classes)
  result <- list(manifold_setup = manifold_setup, transform_setup = transform_setup,
                 estimator_setup = estimator_setup, control_setup = control_setup,
                 parameter_vector = parameter_vector, transparameter_vector = transparameter_vector,
                 parameters = parameters, transparameters = transparameters,
                 posterior = posterior, nparam = nparam, ntransparam = ntransparam,
                 indices_param_vector = indices_param_vector, indices_param_vector2 = indices_param_vector2,
                 indices_full_param_vector = indices_full_param_vector,
                 indices_full_fixed_vector = indices_full_fixed_vector,
                 full_fixed_vector = full_fixed_vector, full_transfixed_vector = full_transfixed_vector,
                 indices_full_transparam_vector = indices_full_transparam_vector,
                 indices_full_transfixed_vector = indices_full_transfixed_vector,
                 indices_transparam_vector2 = indices_transparam_vector2)
  return(result)
}

convert_all_to_numeric <- function(x) {

  if (is.list(x)) {
    lapply(x, convert_all_to_numeric)
  } else if (is.matrix(x)) {
    matrix(as.numeric(x), nrow = nrow(x), dimnames = dimnames(x))
  } else if (is.atomic(x)) {
    as.numeric(x)
  } else {
    stop("Unsupported type")
  }

}
