setup_oblq <- function(...) {

  list2env(..., envir = environment())
  # attach(list(...))

  # q, indices
  result <- list(manifold = "oblq", q = q,
                 indices = indices)

  return(result)

}
