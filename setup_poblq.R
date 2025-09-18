setup_poblq <- function(...) {

  list2env(..., envir = environment())
  # attach(list(...))

  # q, indices, PhiTarget
  result <- list(manifold = "poblq", q = q, PhiTarget = PhiTarget,
                 indices = indices)

  return(result)

}
