setup_unit <- function(...) {

  list2env(..., envir = environment())
  # attach(list(...))

  # indices
  result <- list(manifold = "unit", indices = indices)

  return(result)

}
