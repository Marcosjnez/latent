setup_all_manifolds <- function(manifolds, arguments) {

  result <- list()

  for(i in 1:length(manifolds)) {

    manifold <- manifolds[i]
    args <- arguments[[i]]
    result[[i]] <- setup_manifold(manifold, args)

  }

  return(result)

}
