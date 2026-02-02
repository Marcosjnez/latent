# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 29/01/2026

effects_coding <- function(coeffs, vcov) {

  # Transform from zeros constraints to effects coding constraints

  # intercept_index is the position of the intercept in the coeffs matrix
  intercept_index <- which(apply(coeffs, MARGIN = 2, FUN = \(x) all(x == 0)))

  p <- nrow(coeffs)
  q <- ncol(coeffs)

  # Centering matrix:
  C <- diag(q) - matrix(1 / q, q, q)

  # effects coding parameterization:
  beta_new <- coeffs %*% C

  # Jacobian of coeffs %*% C:
  J <- kronecker(t(C[-intercept_index, ]), diag(p))

  # standard errors under effects coding parameterization:
  vcov_new <- J %*% vcov %*% t(J)
  se_new <- sqrt(diag(vcov_new))

  # Return:
  result <- list(se = se, vcov = vcov,
                 beta_new = beta_new,
                 se_new = se_new,
                 vcov_new = vcov_new)

  return(result)

}

move_intercept <- function(coeffs, vcov, column = 1L) {

  p <- nrow(coeffs)
  q <- ncol(coeffs)

  Theta <- beta[, -1, drop = FALSE]

  beta_new <- matrix(0, p, q, dimnames = dimnames(beta))
  beta_new[, q] <- 0
  beta_new[, 1] <- -Theta[, q - 1]
  if (q > 2) {
    beta_new[, 2:(q - 1)] <- Theta[, 1:(q - 2), drop=FALSE] - Theta[, q - 1]
  }

  R <- matrix(0, nrow = q - 1, ncol = q - 1)
  R[q - 1, 1] <- -1
  if (q > 2) {
    for (n in 2:(q - 1)) {
      R[n - 1, n] <- 1
      R[ - 1, n] <- -1
    }
  }

  A <- kronecker(t(R), diag(p))
  vcov_new <- A %*% vcov %*% t(A)
  se_new <- sqrt(diag(vcov_new))

  # Return:
  result <- list(se = se, vcov = vcov,
                 beta_new = beta_new,
                 se_new = se_new,
                 vcov_new = vcov_new)

  return(result)

}
