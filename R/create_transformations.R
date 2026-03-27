# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 25/03/2026

create_transforms <- function(transforms, structures) {

  # Collect the unique parameter labels that are not fixed values:
  param_structures_vector <- unname(unique(c(unlist(structures))))

  # Handle the transformations:
  transforms_and_labels <- transforms

  # Check that a list was provided:
  if(!is.list(transforms_and_labels)) {
    stop("transforms_and_labels should be a list")
  }
  ntransforms <- length(transforms_and_labels) # Number of transforms
  # Pick the transforms:
  transforms <- lapply(transforms_and_labels, FUN = \(x) x$transform)
  # Pick the in parameter labels going to each transform:
  labels_in <- lapply(transforms_and_labels, FUN = \(x) x$parameters_in)
  # Pick the out parameter labels going to each transform:
  labels_out <- lapply(transforms_and_labels, FUN = \(x) x$parameters_out)
  # Pick the extra objects going to each manifold:
  dots <- lapply(transforms_and_labels, FUN = \(x) x$extra)

  # loop for each transform:
  result <- vector("list", length = ntransforms)
  k <- 1L
  for(i in seq_len(ntransforms)) {

    # Pick the transform label:
    transform <- transforms[[i]]

    # Choose which extra objects from 'dots' should be kept for this transform:
    transform_objects <- switch(transform,
                                XYt = c("p", "q"),
                                XY      = c("p", "q"),
                                softmax   = character(0),
                                normal      = c("y", "S", "J", "I"),
                                multinomial = c("y", "K", "S", "J", "I"),
                                matrix_inverse = c("p"),
                                logarithm   = character(0),
                                identity   = character(0),
                                factor_cor   = c("p", "q"),
                                exponential   = character(0),
                                crossprod   = c("p"),
                                column_space   = c("X"),
                                deltaparam   = c("p", "q"),
                                stop("Unknown transform: ", transform)
    )

    # Pick only those objects from 'dots':
    extra <- dots[[i]][transform_objects]

    # Ensure the required extras are present and named:
    if (length(transform_objects) > 0L) {
      missing <- setdiff(transform_objects, names(extra))
      if (length(missing) > 0L) {
        stop("Missing required object(s) for transform '", transform, "': ",
             paste(missing, collapse = ", "))
      }
    }

    #### labels_in ####

    indices_in <- vector("list", length = length(labels_in[[i]]))
    for(j in 1:length(labels_in[[i]])) {

      # Collect the labels_in:
      if(is.list(labels_in[[i]][j])) {
        labels_in_vector <- unname(c(unlist(labels_in[[i]][[j]])))
      } else {
        labels_in_vector <- unname(c(unlist(structures[labels_in[[i]][[j]]])))
      }

      # Get the indices of the param_structures_vector that are in the labels_in_vector:
      m_in <- match(labels_in_vector, param_structures_vector)
      if (anyNA(m_in)) { # Check for wrong parameter labels_in
        stop("Some labels_in were not found in param_structures_vector: ",
             paste(labels_in_vector[is.na(m_in)], collapse = ", "))
      }

      indices_in[[j]] <- m_in-1L # C++ indexing starts at 0

    }

    #### labels_out ####

    indices_out <- vector("list", length = length(labels_out[[i]]))
    for(j in 1:length(labels_out[[i]])) {

      # Collect the labels_out:
      if(is.list(labels_out[[i]][j])) {
        labels_out_vector <- unname(c(unlist(labels_out[[i]][[j]])))
      } else {
        labels_out_vector <- unname(c(unlist(structures[labels_out[[i]][[j]]])))
      }

      # Get the indices of the param_structures_vector that are in the labels_out_vector:
      m_out <- match(labels_out_vector, param_structures_vector)
      if (anyNA(m_out)) { # Check for wrong parameter labels_out
        stop("Some labels_out were not found in param_structures_vector: ",
             paste(labels_out_vector[is.na(m_out)], collapse = ", "))
      }

      indices_out[[j]] <- m_out-1L # C++ indexing starts at 0

    }

    #### result ####

    result[[k]] <- c(
      list(transform = transform,
           indices_in  = indices_in,
           indices_out = indices_out),
      extra
    )
    k <- k+1L

  }

  return(result)

}
