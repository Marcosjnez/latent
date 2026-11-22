# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 06/07/2026

create_parameters <- function(parameters) {

  type <- lapply(parameters, FUN = \(x) x$type)
  name <- lapply(parameters, FUN = \(x) x$name)
  dims <- lapply(parameters, FUN = \(x) x$dim)
  rnames <- lapply(parameters, FUN = \(x) x$rownames)
  cnames <- lapply(parameters, FUN = \(x) x$colnames)
  symmetric <- lapply(parameters, FUN = \(x) isTRUE(x$symmetric))
  array_names <- lapply(parameters, FUN = \(x) x$dimnames)
  labels <- lapply(parameters, FUN = \(x) x$labels)
  nlist <- length(type)

  result <- vector("list", length = nlist)

  for(i in seq_len(nlist)) {

    if(type[[i]] == "scalar") {

      if(!is.null(labels[[i]])) {
        result[[i]] <- labels[[i]]
      } else {
        result[[i]] <- name[[i]]
      }

    } else if(type[[i]] == "vector") {

      dim1 <- dims[[i]][1]

      if(!is.null(labels[[i]])) {

        result[[i]] <- labels[[i]]

        if(length(result[[i]]) != dim1) {
          stop("Custom labels for parameter block '", name[[i]],
               "' do not match the requested dimensions.")
        }

      } else {
        result[[i]] <- paste(name[[i]], seq_len(dim1), sep = "")
      }

      if(!is.null(rnames[[i]])) names(result[[i]]) <- rnames[[i]]

    } else if(type[[i]] == "matrix") {

      dim1 <- dims[[i]][1]
      dim2 <- dims[[i]][2]

      if(!is.null(labels[[i]])) {

        result[[i]] <- labels[[i]]

        if(!identical(dim(result[[i]]), as.integer(c(dim1, dim2)))) {
          stop("Custom labels for parameter block '", name[[i]],
               "' do not match the requested dimensions.")
        }

      } else {

        vector_names <- paste(name[[i]], "[", rep(seq_len(dim1), times = dim2),
                              ",", rep(seq_len(dim2), each = dim1), "]",
                              sep = "")
        result[[i]] <- matrix(vector_names, nrow = dim1, ncol = dim2)

      }

      if(!is.null(rnames[[i]])) rownames(result[[i]]) <- rnames[[i]]
      if(!is.null(cnames[[i]])) colnames(result[[i]]) <- cnames[[i]]

      if(symmetric[[i]]) {
        upper_indices <- upper.tri(result[[i]])
        result[[i]][upper_indices] <- t(result[[i]])[upper_indices]
      }

    } else if(type[[i]] == "array") {

      if(!is.null(labels[[i]])) {

        result[[i]] <- labels[[i]]

        if(!identical(dim(result[[i]]), as.integer(dims[[i]]))) {
          stop("Custom labels for parameter block '", name[[i]],
               "' do not match the requested dimensions.")
        }

        if(!is.null(array_names[[i]])) dimnames(result[[i]]) <- array_names[[i]]

      } else {

        indices <- arrayInd(seq_len(prod(dims[[i]])), dims[[i]])
        vector_names <- if(nrow(indices)) {
          paste0(name[[i]], "[", do.call(paste, c(as.data.frame(indices), sep = ",")), "]")
        } else character(0)

        result[[i]] <- array(vector_names,
                             dim = dims[[i]],
                             dimnames = array_names[[i]])

      }

    } else {
      stop("Unknown type for parameter block '", name[[i]], "'.")
    }

  }

  names(result) <- unlist(name)

  return(result)

}
