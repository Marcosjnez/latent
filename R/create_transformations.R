# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 22/08/2026

create_transforms <- function(transforms, structures) {

  #### Check inputs ####

  if(!is.list(transforms)) {
    stop("transforms should be a list")
  }

  param_structures_vector <- unname(unique(c(unlist(structures))))
  transform_names <- lapply(transforms, FUN = \(x) x$transform)
  labels_in <- lapply(transforms, FUN = \(x) x$parameters_in)
  labels_out <- lapply(transforms, FUN = \(x) x$parameters_out)
  dots <- lapply(transforms, FUN = \(x) x$extra)

  #### Create transformation control objects ####

  result <- vector("list", length = length(transforms))

  for(i in seq_along(transforms)) {

    transform <- transform_names[[i]]

    transform_objects <- switch(
      transform,
      XYt = c("p", "q"),
      XY = c("p", "q"),
      softmax = character(0L),
      matrix_inverse = "p",
      logarithm = character(0L),
      identity = character(0L),
      factor_cor = c("p", "q"),
      meanstructure = c("p", "q"),
      tau_param = c("p", "threshold_items"),
      exponential = character(0L),
      crossprod = "p",
      column_space = "X",
      deltaparam = c("p", "q"),
      multinomial = c("y", "K", "S", "J", "I"),
      normal = c("y", "S", "J", "I"),
      mvnormal = c("y", "S", "J", "I"),
      sum_vectors = character(0L),
      sqrt_vector = character(0L),
      pos_incrsng = character(0L),
      stop("Unknown transform: ", transform)
    )

    extra_i <- dots[[i]]

    if(is.null(extra_i)) {
      extra_i <- list()
    }

    extra <- extra_i[transform_objects]

    if(length(transform_objects) > 0L) {

      missing_objects <- setdiff(transform_objects, names(extra))

      if(length(missing_objects) > 0L) {
        stop("Missing required object(s) for transform '", transform, "': ",
             paste(missing_objects, collapse = ", "))
      }

    }

    #### Input indices ####

    indices_in <- vector("list", length = length(labels_in[[i]]))

    for(j in seq_along(labels_in[[i]])) {

      if(is.list(labels_in[[i]][j])) {
        labels_in_vector <- unname(c(unlist(labels_in[[i]][[j]])))
      } else {
        labels_in_vector <- unname(c(unlist(structures[labels_in[[i]][[j]]])))
      }

      match_in <- match(labels_in_vector, param_structures_vector)

      if(anyNA(match_in)) {
        stop("Some labels_in were not found in param_structures_vector: ",
             paste(labels_in_vector[is.na(match_in)], collapse = ", "))
      }

      indices_in[[j]] <- match_in-1L

    }

    #### Output indices ####

    indices_out <- vector("list", length = length(labels_out[[i]]))

    for(j in seq_along(labels_out[[i]])) {

      if(is.list(labels_out[[i]][j])) {
        labels_out_vector <- unname(c(unlist(labels_out[[i]][[j]])))
      } else {
        labels_out_vector <- unname(c(unlist(structures[labels_out[[i]][[j]]])))
      }

      match_out <- match(labels_out_vector, param_structures_vector)

      if(anyNA(match_out)) {
        stop("Some labels_out were not found in param_structures_vector: ",
             paste(labels_out_vector[is.na(match_out)], collapse = ", "))
      }

      indices_out[[j]] <- match_out-1L

    }

    result[[i]] <- c(list(transform = transform,
                          indices_in = indices_in,
                          indices_out = indices_out),
                     extra)

  }

  #### Result ####

  return(result)

}
