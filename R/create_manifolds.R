# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 14/02/2026

extra_manifolds <- function(manifold, labels, dots = NULL) {

  # Choose which extra objects from 'dots' should be kept for this manifold:
  mani_objects <- switch(manifold,
                         euclidean = character(0),
                         unit      = character(0),
                         simplex   = character(0),
                         orth      = c("p", "q"),
                         oblq      = c("p", "q"),
                         poblq     = c("p", "q", "constraints"),
                         stop("Unknown manifold: ", manifold)
  )

  # Pick only those objects from 'dots':
  extra <- dots[mani_objects]

  # Ensure the required extras are present and named:
  if (length(mani_objects) > 0L) {
    missing <- setdiff(mani_objects, names(dots))
    if (length(missing) > 0L) {
      stop("Missing required object(s) for manifold '", manifold, "': ",
           paste(missing, collapse = ", "))
    }
  }

  # Define the manifold:
  result <- c(
    list(manifold, labels),
    extra
  )

  return(result)

}

create_manifolds <- function(manifolds_and_labels, param_structures) {

  if(!is.list(manifolds_and_labels)) {
    stop("manifolds_and_labels should be a list")
  }
  nmanifolds <- length(manifolds_and_labels) # Number of manifolds
  # Pick the manifolds:
  manifolds <- lapply(manifolds_and_labels, FUN = \(x) x[[1]])
  # Pick the parameter labels going to each manifold:
  labels <- lapply(manifolds_and_labels, FUN = \(x) x[[2]])

  # Collect the unique parameter labels that are not fixed values:
  param_structures_unique <- unname(unique(unlist(param_structures)))
  nonfixed_param_structures <- which(is.na(suppressWarnings(as.numeric(param_structures_unique))))
  param_structures_vector <- param_structures_unique[nonfixed_param_structures]

  # loop for each manifold:
  result <- vector("list", length = nmanifolds)
  for(i in seq_len(nmanifolds)) {

    # Pick the manifold label:
    manifold <- manifolds[[i]]

    # Collect the unique labels that are not fixed values:
    labels_unique <- unname(unique(unlist(labels[[i]])))
    nonfixed_labels <- which(is.na(suppressWarnings(as.numeric(labels_unique))))
    labels_vector <- labels_unique[nonfixed_labels]

    # Get the indices of the param_structures_vector that are in the labels_vector:
    m <- match(labels_vector, param_structures_vector)
    if (anyNA(m)) { # Check for wrong parameter labels
      stop("Some labels were not found in param_structures_vector: ",
           paste(labels_vector[is.na(m)], collapse = ", "))
    }

    indices <- m - 1L # C++ indexing starts at 0

    # Additional objects for the specific manifolds:
    extra <- manifolds_and_labels[[i]]
    if (length(extra) > 2L) {
      extra <- extra[3:length(extra)]
    } else {
      extra <- list()
    }

    result[[i]] <- c(
      list(manifold = manifold,
           indices  = indices),
      extra
    )

  }

  return(result)

}

get_manifold <- function(manifolds,
                         structures) {

  # Collect the unique parameter labels that are not fixed values:
  unique_str <- unname(unique(unlist(structures)))
  nonfixed_str <- which(is.na(suppressWarnings(as.numeric(unique_str))))
  vector_structures <- unique_str[nonfixed_str]

  # Handle the manifolds:
  manifolds_and_labels <- manifolds

  # Check that a list was provided:
  if(!is.list(manifolds_and_labels)) {
    stop("manifolds should be a list")
  }
  nmanifolds <- length(manifolds_and_labels) # Number of manifolds
  # Pick the manifolds:
  manifolds <- lapply(manifolds_and_labels, FUN = \(x) x$manifold)
  # Pick the parameter labels going to each manifold:
  inputs <- lapply(manifolds_and_labels, FUN = \(x) x$parameters)
  # Pick the extra objects going to each manifold:
  dots <- lapply(manifolds_and_labels, FUN = \(x) x$extra)

  # For each manifold:
  result <- vector("list", length = nmanifolds)
  k <- 1L
  for(i in 1:nmanifolds) {

    manifold <- manifolds[[i]]

    # Choose which extra objects from 'dots' should be kept for this manifold:
    mani_objects <- switch(manifold,
                           euclidean = character(0),
                           unit      = character(0),
                           simplex   = character(0),
                           orth      = c("p", "q"),
                           oblq      = c("p", "q"),
                           poblq     = c("p", "q", "constraints"),
                           stop("Unknown manifold: ", manifold)
    )

    # Pick only those objects from 'dots':
    extra <- dots[[i]][mani_objects]

    # Ensure the required extras are present and named:
    if (length(mani_objects) > 0L) {
      missing <- setdiff(mani_objects, names(dots[[i]]))
      if (length(missing) > 0L) {
        stop("Missing required object(s) for manifold '", manifold, "': ",
             paste(missing, collapse = ", "))
      }
    }

    ninputs <- length(inputs[[i]]) # Number of manifolds

    for(j in 1:ninputs) {

      # Collect the unique subset of parameter labels that are not fixed values:
      labels_unique <- unname(unique(unlist(structures[inputs[[i]][j]])))
      nonfixed_labels <- which(is.na(suppressWarnings(as.numeric(labels_unique))))
      labels_vector <- labels_unique[nonfixed_labels]

      # Get the indices of the vector_structures that are in the labels_vector:
      m <- match(labels_vector, vector_structures)
      if (anyNA(m)) { # Check for wrong parameter labels
        stop("Some parameters were not found in structures: ",
             paste(labels_vector[is.na(m)], collapse = ", "))
      }

      indices <- m - 1L # C++ indexing starts at 0

      result[[k]] <- c(
        list(manifold = manifold,
             indices  = indices),
        extra
      )
      k <- k+1L

    }

  }

  return(result)

}
