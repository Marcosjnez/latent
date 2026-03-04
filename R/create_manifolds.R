# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 02/03/2026

get_manifold <- function(manifolds, structures) {

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
      missing <- setdiff(mani_objects, names(extra))
      if (length(missing) > 0L) {
        stop("Missing required object(s) for manifold '", manifold, "': ",
             paste(missing, collapse = ", "))
      }
    }

    ninputs <- length(inputs[[i]]) # Number of manifolds
    for(j in 1:ninputs) {

      # Collect the unique subset of parameter labels that are not fixed values:
      if(is.list(inputs[[i]][j])) {
        labels_unique <- unname(unique(c(unlist(inputs[[i]][[j]]))))
      } else {
        labels_unique <- unname(unique(c(unlist(structures[inputs[[i]][[j]]]))))
      }
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
