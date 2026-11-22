# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 13/09/2026

# Match all labels at once, so long models only build the lookup table once.
.parameter_indices <- function(specs, structures, labels_all) {
  blocks <- lapply(specs, function(spec) {
    if(is.list(spec)) return(lapply(unname(spec), unlist, use.names = FALSE))
    unknown <- setdiff(spec, names(structures))
    if(length(unknown)) stop("Unknown parameter block: ", paste(unknown, collapse = ", "))
    lapply(unname(structures[spec]), unlist, use.names = FALSE)
  })
  labels <- unlist(blocks, use.names = FALSE)
  indices <- match(labels, labels_all) - 1L
  if(anyNA(indices)) stop("Unknown parameter label: ", labels[which(is.na(indices))[1L]])
  offset <- 0L
  lapply(blocks, function(block) lapply(block, function(x) {
    index <- offset + seq_along(x)
    offset <<- offset + length(x)
    indices[index]
  }))
}

create_transforms <- function(transforms, structures) {

  #### Check inputs ####

  if(!is.list(transforms)) {
    stop("transforms should be a list")
  }

  param_structures_vector <- unname(unique(c(unlist(structures))))
  transform_names <- lapply(transforms, FUN = \(x) x$transform)
  labels_in <- lapply(transforms, FUN = \(x) x$parameters_in)
  labels_out <- lapply(transforms, FUN = \(x) x$parameters_out)
  inputs <- .parameter_indices(labels_in, structures, param_structures_vector)
  outputs <- .parameter_indices(labels_out, structures, param_structures_vector)
  dots <- lapply(transforms, FUN = \(x) x$extra)

  #### Create transformation control objects ####

  result <- vector("list", length = length(transforms))

  for(i in seq_along(transforms)) {

    transform <- transform_names[[i]]

    transform_objects <- switch(
      transform,
      XYt = c("p", "q"),
      XY = c("p", "q"),
      matrix_vector = c("p", "q"),
      XtY = c("p", "q", "r"),
      XtYX = c("p", "q"),
      XYXt = c("p", "q"),
      XYz = c("p", "q", "r"),
      invAx = c("p"),
      invA_X_invAt = c("p"),
      softmax = character(0L),
      matrix_inverse = "p",
      inverse_of_I_minusX = "p",
      logarithm = character(0L),
      identity = character(0L),
      sem_cov_model = c("p", "q"),
      sem_means_model = c("p", "q"),
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
      square_vector = character(0L),
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

    result[[i]] <- c(list(transform = transform,
                          indices_in = inputs[[i]],
                          indices_out = outputs[[i]]),
                     extra)

  }

  #### Result ####

  return(result)

}
