extract_param <- function(model) {

  modelvector <- unlist(model)

  # Get the parameters (those elements that are not digits):
  z <- suppressWarnings(as.numeric(modelvector))
  full_parameter_vector <- modelvector[which(is.na(z))]
  # Get the fixed values (those elements that are digits):
  full_fixed_vector <- modelvector[which(!is.na(z))]

  indices_full_param_vector <- which(modelvector %in% full_parameter_vector)
  indices_full_fixed_vector <- which(modelvector %in% full_fixed_vector)

  # Find the unique elements in all the targets:
  uniques <- unique(modelvector)
  # Get the parameters (those elements that are not digits):
  z <- suppressWarnings(as.numeric(uniques))
  parameter_vector <- uniques[which(is.na(z))]
  # Get the fixed values (those elements that are digits):
  fixed_vector <- uniques[which(!is.na(z))]

  indices_param_vector <- match(parameter_vector, modelvector)

  result <- list(parameter_vector = parameter_vector, # Unique parameter vector
                 fixed_vector = fixed_vector,         # Unique fixed vector
                 indices_param_vector = indices_param_vector,
                 full_parameter_vector = full_parameter_vector,
                 full_fixed_vector = full_fixed_vector,
                 indices_full_param_vector = indices_full_param_vector,
                 indices_full_fixed_vector = indices_full_fixed_vector)

  return(result)

}
