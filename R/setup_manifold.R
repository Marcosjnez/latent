setup_manifold <- function(manifold, ...) {

  if(manifold == "euclidean") {

    result <- setup_euclidean(...)

  } else if(manifold == "unit") {

    result <- setup_unit(...)

  } else if(manifold == "orth"){

    result <- setup_orth(...)

  } else if(manifold == "oblq") {

    result <- setup_oblq(...)

  } else if(manifold == "poblq") {

    result <- setup_poblq(...)

  }

  return(result)

}
