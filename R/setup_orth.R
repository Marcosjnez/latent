setup_orth <- function(...) {

  list2env(..., envir = environment())
  # attach(list(...))

  # q, indices
  result <- list(manifold = "orth", q = q,
                 indices = indices)

  return(result)

}
