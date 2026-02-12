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

    indices <- list(m - 1L) # C++ indexing starts at 0

    # Additional arguments for the specific manifolds:
    extra <- manifolds_and_labels[[i]]
    if (length(extra) > 2L) {
      extra <- extra[3:length(extra)]
    } else {
      extra <- list()
    }

    result[[i]] <- list(manifold = manifold,
                        indices = indices,
                        extra = extra)

  }

  return(result)

}
