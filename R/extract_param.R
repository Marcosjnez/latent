extract_param <- function(classes, conditionals) {

  # Find the unique elements in all the targets:
  uniques <- unique(unlist(c(classes, conditionals)))
  # Get the parameters (those elements that are not digits):
  z <- suppressWarnings(as.numeric(uniques))
  parameter_vector <- uniques[which(is.na(z))]
  # Get the fixed values (those elements that are digits):
  fixed_vector <- uniques[which(!is.na(z))]

  result <- list(parameter_vector = parameter_vector, fixed_vector = fixed_vector)

}
