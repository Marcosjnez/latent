# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 21/08/2026

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
