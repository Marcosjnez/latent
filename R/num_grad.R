num_grad <- function(f, X, ..., cor = FALSE, cov = FALSE, eps = 1e-09) {

  # Compute the derivative of a scalar or vector-valued function
  # X should be a scalar, vector or matrix of parameters
  # Set cor = TRUE if your parameter X is a correlation matrix
  # Set cov = TRUE if your parameter X is a covariance matrix

  fv <- f(X, ...)
  k <- length(c(fv))

  if(!is.matrix(X)) { # If X is a vector

    p <- length(X) # Dimensionality
    dfdX <- matrix(NA, k, p)
    for(i in 1:p) {
      epsilon <- rep(0, p)
      epsilon[i] <- eps
      dfdX[, i] <- c(f(X + epsilon, ...) - f(X - epsilon, ...)) / (2*eps)
    }

    return(dfdX)

  }

  p <- nrow(X)
  q <- ncol(X)

  if(cov) {

    if(p != q) stop("Not a covariance matrix")
    qcov <- q*(q+1)/2 # Dimensionality
    dfdX <- matrix(NA, k, qcov)
    ij <- 0
    for(j in 1:q) {
      for(i in j:q) {
        ij <- ij+1
        epsilon <- matrix(0, q, q)
        epsilon[i, j] <- eps
        epsilon[j, i] <- eps
        dfdX[, ij] <- c(f(X + epsilon, ...) - f(X - epsilon, ...)) / (2*eps)
      }
    }

  } else if(cor) {

    if(p != q) stop("Not a correlation matrix")
    qcor <- q*(q-1)/2 # Dimensionality
    dfdX <- matrix(NA, k, qcor)
    ij <- 0
    for(j in 1:(q-1)) {
      for(i in (j+1):q) {
        ij <- ij+1
        epsilon <- matrix(0, q, q)
        epsilon[i, j] <- eps
        epsilon[j, i] <- eps
        dfdX[, ij] <- c(f(X + epsilon, ...) - f(X - epsilon, ...)) / (2*eps)
      }
    }

  } else { # X is a matrix but not a correlation nor covariance matrix

    pq <- p*q # Dimensionality
    dfdX <- matrix(NA, k, pq)
    ij <- 0
    for(i in 1:pq) {
      ij <- ij+1
      epsilon <- matrix(0, p, q)
      epsilon[ij] <- eps
      dfdX[, ij] <- c(f(X + epsilon, ...) - f(X - epsilon, ...)) / (2*eps)
    }

  }

  if(k == 1) {
    if(cov) { # Is this correct?
      gtemp <- matrix(0, q, q)
      gtemp[lower.tri(gtemp, diag = TRUE)] <- c(dfdX)
      gtemp[upper.tri(gtemp)] <- t(gtemp)[upper.tri(gtemp)]
      dfdX <- gtemp
    } else if(cor) {
      gtemp <- matrix(0, q, q)
      gtemp[lower.tri(gtemp)] <- c(dfdX)
      dfdX <- t(gtemp) + gtemp
      diag(dfdX) <- 0
    } else{
      dfdX <- matrix(c(dfdX), p, q)
    }
  }

  return(dfdX)

}
num_dgrad <- function(g, X, dX, ..., eps = 1e-09) {

  # Compute the differential along the direction dX

  dg <- (g(X + eps*dX, ...) - g(X - eps*dX, ...)) / (2*eps)
  return(dg)

}
