# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 14/02/2026

extra_transforms <- function(transform, labels_in, labels_out, dots) {

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
                              stop("Unknown transform: ", transform)
  )

  # Pick only those objects from 'dots':
  extra <- dots[transform_objects]

  # Ensure the required extras are present and named:
  if (length(transform_objects) > 0L) {
    missing <- setdiff(transform_objects, names(dots))
    if (length(missing) > 0L) {
      stop("Missing required object(s) for transform '", transform, "': ",
           paste(missing, collapse = ", "))
    }
  }

  # Define the manifold:
  result <- c(
    list(transform, labels_in, labels_out),
    extra
  )

  return(result)

}

create_transforms <- function(transforms_and_labels, param_structures) {

  if(!is.list(transforms_and_labels)) {
    stop("transforms_and_labels should be a list")
  }
  ntransforms <- length(transforms_and_labels) # Number of transforms
  # Pick the transforms:
  transforms <- lapply(transforms_and_labels, FUN = \(x) x[[1]])
  # Pick the in parameter labels going to each transform:
  labels_in <- lapply(transforms_and_labels, FUN = \(x) x[[2]])
  # Pick the out parameter labels going to each transform:
  labels_out <- lapply(transforms_and_labels, FUN = \(x) x[[3]])

  # Collect the unique parameter labels that are not fixed values:
  param_structures_unique <- unname(unique(unlist(param_structures)))
  nonfixed_param_structures <- which(is.na(suppressWarnings(as.numeric(param_structures_unique))))
  param_structures_vector <- param_structures_unique[nonfixed_param_structures]

  # loop for each transform:
  result <- vector("list", length = ntransforms)
  for(i in seq_len(ntransforms)) {

    # Pick the transform label:
    transform <- transforms[[i]]

    #### labels_in ####

    # Collect the unique labels_in that are not fixed values:
    labels_in_unique <- unname(unique(unlist(labels_in[[i]])))
    nonfixed_labels_in <- which(is.na(suppressWarnings(as.numeric(labels_in_unique))))
    labels_in_vector <- labels_in_unique[nonfixed_labels_in]

    # Get the indices of the param_structures_vector that are in the labels_in_vector:
    m_in <- match(labels_in_vector, param_structures_vector)
    if (anyNA(m_in)) { # Check for wrong parameter labels_in
      stop("Some labels_in were not found in param_structures_vector: ",
           paste(labels_in_vector[is.na(m_in)], collapse = ", "))
    }

    indices_in <- list(m_in - 1L) # C++ indexing starts at 0

    #### labels_out ####

    # Collect the unique labels_out that are not fixed values:
    labels_out_unique <- unname(unique(unlist(labels_out[[i]])))
    nonfixed_labels_out <- which(is.na(suppressWarnings(as.numeric(labels_out_unique))))
    labels_out_vector <- labels_out_unique[nonfixed_labels_out]

    # Get the indices of the param_structures_vector that are in the labels_out_vector:
    m_out <- match(labels_out_vector, param_structures_vector)
    if (anyNA(m_out)) { # Check for wrong parameter labels_out
      stop("Some labels_out were not found in param_structures_vector: ",
           paste(labels_out_vector[is.na(m_out)], collapse = ", "))
    }

    indices_out <- list(m_out - 1L) # C++ indexing starts at 0

    #### Extra arguments ####

    # Additional arguments for the specific transforms:
    extra <- transforms_and_labels[[i]]
    if (length(extra) > 3L) {
      extra <- extra[4:length(extra)]
    } else {
      extra <- list()
    }

    result[[i]] <- c(
      list(transform = transform,
           indices_in  = indices_in,
           indices_out = indices_out),
      extra
    )

  }

  return(result)

}
