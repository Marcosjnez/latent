# Auxiliary function for lca.R

getargs_lca <- function(data_list, model, classes, conditionals, control) {

  # This function generates all the arguments that are necessary to execute the
  # estimators for either normal and categorical data

  Y <- data_list$Y
  n <- data_list$n
  uniq_indices <- data_list$uniq_indices
  map2full <- data_list$map2full

  # Create the vector of unique estimators:
  estimators <- unique(model)
  estimators[estimators == "gaussian"] <- "lca_gaussian"
  estimators[estimators == "multinomial"] <- "lca_multinomial"

  nestimators <- length(estimators)
  arguments <- vector("list", length = nestimators)

  estimator_setup <- list()
  k <- 1

  ep <- extract_param(classes, conditionals)
  parameter_vector <- ep$parameter_vector
  fixed_vector <- ep$fixed_vector
  opt <- control$opt
  do.fit <- control$do.fit

  # Count the number of class and conditional parameters:
  nclasses <- length(classes)
  ncond <- sum(grepl("[A-Za-z]", unique(unlist(c(classes, conditionals)))))
  cond <- vector(length = ncond)

  # Generate the arguments that are necessary to execute the estimators for
  # either gaussian and multinomial data

  gauss <- "gaussian" %in% model
  categ <- "multinomial" %in% model

  # Initial values for the gaussian conditional parameters:
  if(gauss) {
    # For gaussian data:
    indices <- which(model == "gaussian")
    normal_data <- Y[, model == "gaussian", drop = FALSE]
    estimator_setup[[k]] <- lca_gaussian(normal_data, n, uniq_indices, map2full, classes,
                                         conditionals[indices],
                                         parameter_vector, fixed_vector, opt)
    means <- rep(colMeans(normal_data), each = nclasses)
    if(opt == "em") {
      sds <- rep(apply(normal_data, MARGIN = 2, var), each = nclasses)
    } else {
      sds <- rep(log(apply(normal_data, MARGIN = 2, sd)), each = nclasses)
    }
    indices <- estimator_setup[[1]]$indices + 1
    indices_conditionals2 <- indices[estimator_setup[[1]]$indices_conditionals2 + 1]
    indices_target_conditionals2 <- estimator_setup[[1]]$indices_target_conditionals2 + 1
    cond[indices_conditionals2] <- c(rbind(means, sds))[indices_target_conditionals2]
    k <- k+1 # In case that there is a multinomial model
  }

  # Initial values for the multinomial conditional parameters:
  if(categ) {
    # For multinomial data:
    indices <- which(model == "multinomial")
    categorical_data <- Y[, model == "multinomial", drop = FALSE]
    estimator_setup[[k]] <- lca_multinomial(categorical_data, n, uniq_indices, map2full,
                                            classes,
                                            conditionals[indices],
                                            parameter_vector, fixed_vector, opt)
    indices <- estimator_setup[[k]]$indices + 1
    cond_indices <- indices[estimator_setup[[k]]$indices_conditionals2 + 1]
    cond[cond_indices] <- rnorm(length(estimator_setup[[k]]$indices_conditionals2))
    # GET BETTER STARTING VALUES
    k <- k+1
  }

  if(control$softmax) {
    # TRUE for derivative-based algorithms:
    estimator_setup[[k]] <- list()
    estimator_setup[[k]]$estimator <- "latentloglik_combination_softmax"
    estimator_setup[[k]]$S <- estimator_setup[[k-1]]$S
    estimator_setup[[k]]$n <- estimator_setup[[k-1]]$n
    estimator_setup[[k]]$indices <- 1:(nclasses-1) - 1
    estimator_setup[[k]]$nclasses <- nclasses
    estimator_setup[[k]]$classes <- estimator_setup[[k-1]]$classes
    estimator_setup[[k]]$indices_classes <- estimator_setup[[k-1]]$indices_classes
    estimator_setup[[k]]$indices_target_classes <- estimator_setup[[k-1]]$indices_target_classes
  } else {
    # TRUE for EM:
    estimator_setup[[k]] <- list()
    estimator_setup[[k]]$estimator <- "latentloglik_combination"
    estimator_setup[[k]]$S <- estimator_setup[[k-1]]$S
    estimator_setup[[k]]$n <- estimator_setup[[k-1]]$n
    estimator_setup[[k]]$indices <- 1:(nclasses-1) - 1
    estimator_setup[[k]]$nclasses <- nclasses
    estimator_setup[[k]]$classes <- estimator_setup[[k-1]]$classes
    estimator_setup[[k]]$indices_classes <- estimator_setup[[k-1]]$indices_classes
    estimator_setup[[k]]$indices_target_classes <- estimator_setup[[k-1]]$indices_target_classes
  }

  # Initial values for the parameters and posterior:
  parameters <- posterior <- vector("list", length = control$rstarts)
  class_indices <- estimator_setup[[1]]$indices_classes + 1
  S <- nrow(Y)

  if("parameter" %in% names(control)) {
    parameters <- control$parameters
  } else {
    for(i in 1:control$rstarts) {
      parameters[[i]] <- cond
      parameters[[i]][class_indices] <- rnorm(length(class_indices))
    }
  }

  if("posterior" %in% names(control)) {
    posterior <- control$posterior
  } else {
    for(i in 1:control$rstarts) {
      posterior[[i]] <- matrix(runif(S * nclasses), nrow = S, ncol = nclasses)
      posterior[[i]] <- posterior[[i]] / rowSums(posterior[[i]])
    }
  }

  result <- list(estimator_setup = estimator_setup,
                 parameters = parameters,
                 posterior = posterior)

  # classes: model for the classes
  # conditionals: list with the model for the conditional parameters of each item
  # estimator_setup: list of objects for evaluating the estimators in C++

  return(result)

}
