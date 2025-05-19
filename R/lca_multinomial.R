# Auxiliary function for lca.R

lca_multinomial <- function(Y, n, uniq_indices, map2full, classes, conditionals,
                            parameter_vector, fixed_vector) {

  S <- nrow(Y) # Number of different response patterns
  J <- ncol(Y) # Number of items
  K <- vector(length = J) # Number of categories per item

  nclasses <- length(classes) # Number of classes in each group

  indices_target_classes <- indices_target_conditionals <- indices_classes <-
    indices_conditionals <- c()

  classes_hat <- vector(length = nclasses)
  conditionals_hat <- vector("list", length = J)
  # Find which elements in the targets correspond to a parameter:
  indices_target_classes <- which(classes %in% parameter_vector) # Which classes are estimated
  indices_target_conditionals2 <- which(unlist(conditionals) %in% parameter_vector) # Which conditionals

  # Relate the parameters in the targets to the parameters in the parameter vector:
  if(length(indices_target_classes) == 0) {
    indices_classes <- logical(0)
  } else {
    indices_classes <- match(classes[indices_target_classes], parameter_vector) # Which classes are estimated in group i
  }
  if(length(indices_target_conditionals2) == 0) {
    indices_conditionals2 <- logical(0)
  } else {
    indices_conditionals2 <- match(unlist(conditionals)[indices_target_conditionals2], parameter_vector) # Which conditionals are estimated in group i
  }

  # indices <- match(classes[indices_target_classes], parameter_vector)
  duplicated_indices <- c(indices_classes, indices_conditionals2)
  indices <- unique(duplicated_indices)

  # Find which elements in the targets correspond to a fixed value:
  indices_fixtarget_classes <- which(classes %in% fixed_vector) # Which classes are fixed in group i
  indices_fixtarget_conditionals <- which(unlist(conditionals) %in% fixed_vector) # Which conditionals are fixed in group i

  # Relate the elements in the targets to the fixed values in the fixed vector:
  if(length(indices_fixtarget_classes) == 0) {
    indices_fix_classes <- logical(0)
  } else {
    indices_fix_classes <- match(classes[indices_fixtarget_classes], fixed_vector) # Which classes are fixed in group i
  }
  if(length(indices_fixtarget_conditionals) == 0) {
    indices_fix_conditionals <- logical(0)
  } else {
    indices_fix_conditionals <- match(unlist(conditionals)[indices_fixtarget_conditionals], fixed_vector) # Which classes are fixed in group i
  }

  # Non-specified elements in classes are set to 0:
  classes_hat <- rep(0, nclasses)
  classes_hat[indices_fixtarget_classes] <- fixed_vector[indices_fix_classes]
  class(classes_hat) <- "numeric"
  if(any(classes_hat < 0)) {
    constant_class <- 1
    fix_logit <- TRUE
  } else {
    constant_class <- 1-sum(classes_hat)
    fix_logit <- FALSE
  }

  # Non-specified elements in conditionals are set to 0:
  for(j in 1:J) {
    K[j] <- nrow(conditionals[[j]])
    conditionals_hat[[j]] <- matrix(0, nrow = K[j], ncol = nclasses)
  }
  flat_vector <- unlist(conditionals_hat)
  flat_vector[indices_fixtarget_conditionals] <- fixed_vector[indices_fix_conditionals]

  flat_vector_logical <- rep(TRUE, times = length(flat_vector))
  flat_vector_logical[indices_fixtarget_conditionals] <- FALSE
  conditionals_hat_logical <- conditionals_hat

  indices_conditionals <- indices_target_conditionals <-
    constant_cond <- vector("list", length = J)

  dims <- lapply(conditionals_hat, dim)
  start <- 1
  for (j in 1:J) {
    # Calculate the number of elements in the current matrix
    num_elements <- prod(dims[[j]])
    # Extract the relevant portion of the flat vector for this matrix
    conditionals_hat[[j]] <- matrix(flat_vector[start:(start + num_elements - 1)],
                                    nrow = dims[[j]][1])
    class(conditionals_hat[[j]]) <- "numeric"

    fill <- flat_vector_logical[start:(start + num_elements - 1)]
    conditionals_hat_logical[[j]] <- matrix(fill, nrow = dims[[j]][1])

    # Update the start position for the next matrix
    start <- start + num_elements

    for(i in 1:nclasses) {
      logical_vector <- conditionals_hat_logical[[j]][, i]
      targets <- which(logical_vector)
      indices_target_conditionals[[j]][[i]] <- targets-1
      indices_conditionals[[j]][[i]] <- match(conditionals[[j]][targets, i],
                                              parameter_vector[indices])-1 # Which conditionals are estimated in group i
      constant_cond[[j]][i] <- 1-sum(conditionals_hat[[j]][, i])
    }
  }

  nparam <- length(indices)
  # indices <- 1:nparam-1

  # logclasses <- log(classes_hat)
  # logclasses[is.infinite(logclasses)] <- 0

  # For each item, set the minimum score to 0
  for(j in 1:ncol(Y)) {
    minimum <- min(Y[, j])
    Y[, j] <- Y[, j] - minimum
  }

  export <- list(estimator = "lca_multinomial",
                 Y = Y,
                 S = S,
                 J = J,
                 K = K,
                 n = n,
                 nclasses = nclasses,
                 uniq_indices = uniq_indices,
                 map2full = map2full,
                 classes = classes_hat,
                 conditionals = conditionals_hat,
                 indices_classes = indices_classes-1,
                 indices_target_classes = indices_target_classes-1,
                 indices_conditionals = indices_conditionals,
                 indices_target_conditionals = indices_target_conditionals,
                 constant_class = constant_class,
                 fix_logit = fix_logit,
                 constant_cond = constant_cond,
                 nparam = nparam,
                 indices = indices-1,
                 indices_conditionals2 = match(indices_conditionals2, indices)-1,
                 indices_target_conditionals2 = indices_target_conditionals2-1)

  return(export)

}
