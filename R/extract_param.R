# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 22/05/2025

extract_param <- function(model) {

  ## Extract labels and associated indices ##

  ## Separate parameter labels from fixed values ##
  modelvector <- unlist(model)
  z <- suppressWarnings(as.numeric(modelvector))
  # Get the parameter labels:
  full_parameter_vector <- modelvector[which(is.na(z))]
  # Get the fixed values (digits):
  full_fixed_vector <- modelvector[which(!is.na(z))]

  # Indices that map the parameter labels to the full model vector:
  indices_full_param_vector <- which(modelvector %in% full_parameter_vector)
  # Indices that map the fixed values to the full model vector:
  indices_full_fixed_vector <- which(modelvector %in% full_fixed_vector)

  ## Separate the unique parameter labels from the unique fixed values ##
  uniques <- unique(modelvector) # Unique elements in the full model vector
  z <- suppressWarnings(as.numeric(uniques))
  # Get the unique parameter labels:
  parameter_vector <- uniques[which(is.na(z))]
  # Get the unique fixed values (digits):
  fixed_vector <- uniques[which(!is.na(z))]

  # Indices that map the unique parameter labels to the full model vector:
  indices_param_vector <- match(parameter_vector, modelvector)
  # Indices that map the full model vector to the unique parameter labels :
  indices_param_vector2 <- match(modelvector, parameter_vector)

  result <- list(full_parameter_vector = full_parameter_vector, # Full parameter labels vector
                 full_fixed_vector = full_fixed_vector, # Full fixed vector
                 indices_full_param_vector = indices_full_param_vector,
                 indices_full_fixed_vector = indices_full_fixed_vector,
                 parameter_vector = parameter_vector, # Unique parameter labels vector
                 fixed_vector = fixed_vector, # Unique fixed vector
                 indices_param_vector = indices_param_vector,
                 indices_param_vector2 = indices_param_vector2)

  return(result)

}
