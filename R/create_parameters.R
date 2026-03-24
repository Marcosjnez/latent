# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 16/03/2026

create_parameters <- function(parameters) {

  type <- lapply(parameters, FUN = \(x) x$type)
  name <- lapply(parameters, FUN = \(x) x$name)
  dims <- lapply(parameters, FUN = \(x) x$dim)
  rnames <- lapply(parameters, FUN = \(x) x$rownames)
  cnames <- lapply(parameters, FUN = \(x) x$colnames)
  symmetric <- lapply(parameters, FUN = \(x) isTRUE(x$symmetric))
  array_names <- lapply(parameters, FUN = \(x) x$dimnames)
  nlist <- length(type)

  result <- vector("list", length = nlist)

  for(i in 1:nlist) {

    if(type[[i]] == "scalar") {
      result[[i]] <- name[[i]]
    } else if(type[[i]] == "vector") {
      result[[i]] <- paste(name[[i]], 1:dim, sep = "")
    } else if(type[[i]] == "matrix") {
      dim1 <- dims[[i]][1]
      dim2 <- dims[[i]][2]
      vector_names <- paste(name[[i]], "[", rep(1:dim1, times = dim2),
                             ",", rep(1:dim2, each = dim1), "]", sep = "")
      result[[i]] <- matrix(vector_names, nrow = dim1, ncol = dim2)
      if(!is.null(rnames[[i]])) rownames(result[[i]]) <- rnames[[i]]
      if(!is.null(cnames[[i]])) colnames(result[[i]]) <- cnames[[i]]
      if(symmetric[[i]]) {
        upper_indices <- upper.tri(result[[i]])
        result[[i]][upper_indices] <- t(result[[i]])[upper_indices]
      }
    } else if(type[[i]] == "array") {

      dim1 <- dims[[i]][1]
      dim2 <- dims[[i]][2]
      dim3 <- dims[[i]][3]

      dim1_expanded <- rep(1:dim1, times  = dim2*dim3)
      dim2_expanded <- rep(rep(1:dim2, times = dim3), each = dim1)
      dim3_expanded <- rep(rep(1:dim3, each = dim2), each = dim1)

      vector_names <- paste("loglik[", dim1_expanded, ",",
                            dim2_expanded, ",",
                            dim3_expanded, "]", sep = "")

      result[[i]] <- array(vector_names,
                           dim = dims[[i]],
                           dimnames = array_names[[i]])

    }

  }

  names(result) <- unlist(name)

  return(result)

}
