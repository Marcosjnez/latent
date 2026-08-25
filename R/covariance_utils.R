# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 25/08/2026

#### Validate and symmetrize a covariance matrix ####

validate_covariance_matrix <- function(V, labels = NULL,
                                       object_name = "variance-covariance matrix") {

  if(!is.matrix(V)) {
    V <- as.matrix(V)
  }

  if(length(dim(V)) != 2L ||
     nrow(V) != ncol(V) ||
     any(!is.finite(V))) {
    stop(object_name, " must be a finite square numeric matrix.")
  }

  if(!isSymmetric(V, tol = sqrt(.Machine$double.eps))) {
    stop(object_name, " must be symmetric.")
  }

  if(!is.null(labels)) {

    if(length(labels) != nrow(V)) {
      stop("The dimensions of ", object_name,
           " do not match the parameter labels.")
    }

    if(!is.null(rownames(V)) || !is.null(colnames(V))) {

      if(is.null(rownames(V)) || is.null(colnames(V)) ||
         !all(labels %in% rownames(V)) ||
         !all(labels %in% colnames(V))) {
        stop("The names of ", object_name,
             " do not match the parameter labels.")
      }

      V <- V[labels, labels, drop = FALSE]

    }

    rownames(V) <- colnames(V) <- labels

  }

  V <- (V+t(V))/2

  #### Result ####

  return(V)

}

#### Stable inverse of a symmetric information matrix ####

invert_information_matrix <- function(H, labels = NULL,
                                      object_name = "information matrix",
                                      tolerance = sqrt(.Machine$double.eps),
                                      warn = TRUE) {

  H <- validate_covariance_matrix(H, labels = labels,
                                  object_name = object_name)

  if(nrow(H) == 0L) {

    #### Result ####

    return(H)

  }

  eigenvalues <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  scale <- max(1, max(abs(eigenvalues)))
  threshold <- tolerance*scale
  rank <- sum(eigenvalues > threshold)
  condition <- if(all(eigenvalues > threshold)) {
    max(eigenvalues)/min(eigenvalues)
  } else {
    Inf
  }

  if(warn && rank < nrow(H)) {
    warning("The ", object_name, " is singular or not positive definite; ",
            "an approximate symmetric inverse was used.")
  }

  result <- approx_Hinv(H)

  if(!is.matrix(result)) {
    result <- as.matrix(result)
  }

  if(any(!is.finite(result))) {
    stop("The ", object_name, " could not be inverted.")
  }

  result <- (result+t(result))/2
  dimnames(result) <- dimnames(H)
  attr(result, "rank") <- rank
  attr(result, "condition") <- condition
  attr(result, "minimum_eigenvalue") <- min(eigenvalues)

  #### Result ####

  return(result)

}

#### Hessian or KKT inverse operator ####

invert_hessian_latent <- function(fit, H,
                                  labels = fit@modelInfo$parameters_labels,
                                  object_name = "Hessian") {

  H <- validate_covariance_matrix(H,
                                  labels = labels,
                                  object_name = object_name)

  if(nrow(H) == 0L) {

    #### Result ####

    return(H)

  }

  se_method <- fit@modelInfo$control_optimizer$se_method

  if(is.null(se_method)) {
    se_method <- "Hessian"
  }

  if(!is.character(se_method) ||
     length(se_method) != 1L ||
     is.na(se_method)) {
    stop("se_method must be 'Hessian' or 'KKT'.")
  }

  se_method <- toupper(se_method)

  if(!(se_method %in% c("HESSIAN", "KKT"))) {
    stop("se_method must be 'Hessian' or 'KKT'.")
  }

  if(se_method == "KKT") {

    derivatives <- constraints_derivs(
      fit,
      parameters = labels
    )

    dconstr <- derivatives$dconstr
    d2constr <- derivatives$d2constr

    empty_constraints <- is.null(dconstr) ||
      is.null(d2constr) ||
      is.null(dim(dconstr)) ||
      is.null(dim(d2constr)) ||
      ncol(dconstr) == 0L ||
      nrow(d2constr) == 0L ||
      ncol(d2constr) == 0L

  } else {

    empty_constraints <- TRUE

  }

  if(empty_constraints) {

    result <- invert_information_matrix(
      H,
      labels = labels,
      object_name = object_name
    )

    #### Result ####

    return(result)

  }

  #### Active constraint derivatives ####

  nparam <- length(labels)
  nconstraints <- ncol(dconstr)
  expected_dimension <- nparam*nconstraints

  if(nrow(dconstr) != nparam ||
     nrow(d2constr) != expected_dimension ||
     ncol(d2constr) != expected_dimension) {
    stop("The first and second constraint derivatives are not aligned.")
  }

  dconstr <- as.matrix(dconstr)

  if(any(!is.finite(dconstr))) {
    stop("The first constraint derivatives contain non-finite values.")
  }

  constraint_idx <- which(colSums(abs(dconstr)) > 0)

  if(length(constraint_idx) == 0L) {

    result <- invert_information_matrix(
      H,
      labels = labels,
      object_name = object_name
    )

    #### Result ####

    return(result)

  }

  dconstr <- dconstr[, constraint_idx, drop = FALSE]
  nconstraints <- ncol(dconstr)

  #### Euclidean gradient and Lagrange multipliers ####

  control_optimizer <- fit@modelInfo$control_optimizer
  control_optimizer$parameters[[1L]] <- fit@Optim$parameters
  control_optimizer$transparameters[[1L]] <- fit@Optim$transparameters

  gradient <- get_grad(
    control_manifold = fit@modelInfo$control_manifold,
    control_transform = fit@modelInfo$control_transform,
    control_estimator = fit@modelInfo$control_estimator,
    control_optimizer = control_optimizer
  )$g
  gradient <- c(gradient)

  free_labels <- fit@modelInfo$parameters_labels

  if(length(gradient) != length(free_labels) ||
     any(!is.finite(gradient)) ||
     !all(labels %in% free_labels)) {
    stop("The gradient does not match the KKT parameter coordinates.")
  }

  names(gradient) <- free_labels
  gradient <- gradient[labels]

  lambda <- -approx_Hinv(crossprod(dconstr)) %*%
    crossprod(dconstr, gradient)
  lambda <- c(lambda)

  if(length(lambda) != nconstraints ||
     any(!is.finite(lambda))) {
    stop("The KKT Lagrange multipliers could not be computed.")
  }

  #### Hessian of the Lagrangian ####

  H_L <- H

  for(j in seq_len(nconstraints)) {

    constraint <- constraint_idx[j]
    idx <- (constraint-1L)*nparam+seq_len(nparam)

    H_constraint <- as.matrix(
      d2constr[idx, idx, drop = FALSE]
    )

    if(any(!is.finite(H_constraint))) {
      stop("A constraint Hessian contains non-finite values.")
    }

    H_constraint <- (H_constraint+t(H_constraint))/2
    H_L <- H_L+lambda[j]*H_constraint

  }

  H_L <- (H_L+t(H_L))/2

  #### KKT inverse ####

  KKT <- rbind(
    cbind(H_L, dconstr),
    cbind(t(dconstr), matrix(0,
                             nrow = nconstraints,
                             ncol = nconstraints))
  )
  KKT <- (KKT+t(KKT))/2

  KKT_inv <- approx_inv(KKT) # KKT_inv may not be symmetric

  if(!is.matrix(KKT_inv)) {
    KKT_inv <- as.matrix(KKT_inv)
  }

  if(any(!is.finite(KKT_inv))) {
    stop("The KKT matrix could not be inverted.")
  }

  parameter_idx <- seq_len(nparam)
  result <- KKT_inv[parameter_idx, parameter_idx, drop = FALSE]
  result <- (result+t(result))/2
  rownames(result) <- colnames(result) <- labels

  #### Result ####

  return(result)

}

#### Standard errors from a covariance matrix ####

standard_errors_from_vcov <- function(V,
                                      tolerance = sqrt(.Machine$double.eps),
                                      object_name = "variance-covariance matrix") {

  V <- validate_covariance_matrix(V, object_name = object_name)
  variances <- diag(V)

  if(length(variances) == 0L) {

    result <- numeric(0L)

    #### Result ####

    return(result)

  }

  scale <- max(1, max(abs(variances)))
  negligible <- variances < 0 & variances >= -tolerance*scale
  variances[negligible] <- 0

  invalid <- which(variances < 0)

  if(length(invalid) > 0L) {

    labels <- rownames(V)

    if(is.null(labels)) {
      labels <- as.character(invalid)
    } else {
      labels <- labels[invalid]
    }

    # stop("Negative parameter variance(s) were found in ", object_name, ": ",
    #      paste(labels, collapse = ", "))
    warning("Negative parameter variance(s) were found in ", object_name, ": ",
         paste(labels, collapse = ", "))

  }

  result <- sqrt(variances)

  if(!is.null(rownames(V))) {
    names(result) <- rownames(V)
  }

  #### Result ####

  return(result)

}
