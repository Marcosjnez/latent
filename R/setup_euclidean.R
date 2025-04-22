setup_euclidean <- function(...) {

  list2env(..., envir = environment())
  # attach(list(...))

  # indices
  result <- list(manifold = "euclidean", indices = indices)

  return(result)

}
