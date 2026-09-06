# Author: Marcos Jimenez
# Modification date: 06/09/2026

# Moment-estimation controls are independent of the CFA/EFA/rotation controls.
# Only shared input validation belongs here; the selected moment estimator
# supplies defaults and validates its remaining estimator-specific options.
normalize_control_moments <- function(control.moments) {

  if(is.null(control.moments)) {
    control.moments <- list()
  }

  if(!is.list(control.moments) || is.data.frame(control.moments)) {
    stop("control.moments must be NULL or a named list")
  }

  if(length(control.moments) > 0L) {

    control_names <- names(control.moments)
    if(is.null(control_names) || anyNA(control_names) ||
       any(control_names == "") || anyDuplicated(control_names)) {
      stop("control.moments must have unique, non-empty names")
    }

  }

  cores <- control.moments$cores
  if(!is.null(cores)) {

    if(!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) ||
       cores < 1 || cores > .Machine$integer.max || cores != floor(cores)) {
      stop("control.moments$cores must be a positive integer")
    }
    control.moments$cores <- as.integer(cores)

  }

  #### Result ####

  return(control.moments)

}
