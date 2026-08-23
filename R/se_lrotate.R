# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 23/08/2026

#### Standard errors for rotations of fitted lcfa models ####

se_lrotate_multistep <- function(fit, parameters = NULL, digits = 4L) {

  #### Check inputs ####

  if(!inherits(fit, "multistep")) {
    stop("fit must inherit from class 'multistep'.")
  }

  if(length(fit@Optim$transparameters) == 0L) {
    stop("The rotation model has not been fitted.")
  }

  source_models <- fit@extra[
    vapply(fit@extra,
           FUN = \(x) inherits(x, "lcfa"),
           FUN.VALUE = logical(1L))
  ]

  if(length(source_models) != 1L) {
    stop("The rotation must contain exactly one fitted lcfa source model.")
  }

  source <- source_models[[1L]]

  if(is.null(parameters)) {
    parameters <- fit@modelInfo$trans[names(fit@modelInfo$param)]
  } else {
    selected_parameters <- unique(unlist(parameters))

    if(!all(selected_parameters %in%
            fit@modelInfo$transparameters_labels)) {
      stop("Unknown parameters.")
    }
  }

  if(!is.null(digits) &&
     (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits < 0L || digits != as.integer(digits))) {
    stop("digits must be NULL or a non-negative integer.")
  }

  #### Unrestricted rotation model ####

  current_labels <- fit@modelInfo$parameters_labels
  fit_full <- unrestricted_model_multistep(fit)
  full_labels <- fit_full@modelInfo$parameters_labels

  missing_current <- setdiff(current_labels, full_labels)

  if(length(missing_current) > 0L) {
    stop("The unrestricted rotation model does not contain current-step parameter(s): ",
         paste(missing_current, collapse = ", "))
  }

  previous_labels <- full_labels[
    !(full_labels %in% current_labels)
  ]

  #### Covariance of the fitted lcfa parameters ####

  source_labels <- source@modelInfo$step_labels

  if(is.null(source_labels)) {
    source_labels <- source@modelInfo$parameters_labels
  }

  missing_previous <- setdiff(previous_labels, source_labels)

  if(length(missing_previous) > 0L) {
    stop("The fitted lcfa covariance matrix does not contain parameter(s): ",
         paste(missing_previous, collapse = ", "))
  }

  A <- stored_vcov_multistep(source)
  A <- A[previous_labels, previous_labels, drop = FALSE]
  A <- (A+t(A))/2

  #### Hessian of the rotation parameters ####

  H2 <- hessian(fit)
  H2 <- H2[current_labels, current_labels, drop = FALSE]
  H2 <- (H2+t(H2))/2

  H2_inv <- invert_hessian_multistep(
    H2,
    object_name = "rotation step"
  )

  #### Cross derivatives between lcfa and rotation parameters ####

  H_full <- hessian(fit_full)
  H_full <- H_full[full_labels, full_labels, drop = FALSE]
  H_full <- (H_full+t(H_full))/2

  if(length(previous_labels) > 0L) {

    C <- H_full[previous_labels, current_labels, drop = FALSE]

    B <- t(C)%*%A%*%C
    B <- (B+t(B))/2

    v12 <- -A%*%C%*%H2_inv
    v22 <- H2_inv%*%B%*%H2_inv
    v22 <- (v22+t(v22))/2

  } else {

    C <- matrix(0,
                nrow = 0L,
                ncol = length(current_labels),
                dimnames = list(character(0L), current_labels))

    B <- matrix(0,
                nrow = length(current_labels),
                ncol = length(current_labels),
                dimnames = list(current_labels, current_labels))

    v12 <- matrix(0,
                  nrow = 0L,
                  ncol = length(current_labels),
                  dimnames = list(character(0L), current_labels))

    v22 <- B

  }

  #### Joint covariance of lcfa and rotation parameters ####

  joint_labels <- c(previous_labels, current_labels)

  v <- matrix(0,
              nrow = length(joint_labels),
              ncol = length(joint_labels),
              dimnames = list(joint_labels, joint_labels))

  if(length(previous_labels) > 0L) {
    v[previous_labels, previous_labels] <- A
    v[previous_labels, current_labels] <- v12
    v[current_labels, previous_labels] <- t(v12)
  }

  v[current_labels, current_labels] <- v22
  v <- (v+t(v))/2

  #### Covariance of rotated transformed parameters ####

  VCOV <- vcov(fit_full,
               v = v,
               parameters = parameters)

  est <- fill_in(parameters,
                 fit_full@Optim$transparameters,
                 miss = NA)

  table_se <- fill_in(parameters,
                      VCOV$se,
                      miss = NA)

  table <- combine_est_se(est, table_se, digits = digits)

  free_VCOV <- v[current_labels, current_labels, drop = FALSE]

  sample_se <- source@dataList$se_type

  if(is.null(sample_se)) {
    sample_se <- NA_character_
  }

  #### Result ####

  result <- list(table = table,
                 table_se = table_se,
                 se = c(VCOV$se),
                 vcov = VCOV$vcov,
                 VCOV = VCOV$vcov,
                 free_VCOV = free_VCOV,
                 jacob = VCOV$jacob,
                 H = H2,
                 B = B,
                 A = A,
                 C = C,
                 joint_vcov = v,
                 type = "multistep",
                 sample_se = sample_se)

  return(result)

}
