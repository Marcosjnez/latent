# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 07/06/2025

getargs_cfa_ml <- function(parameter_vector, fixed_vector,
                           data, lambda, phi, psi,
                           missing, positive = FALSE) {

  # Initialize the vector of initial parameter estimates:
  init <- vector(length = length(parameter_vector))

  p <- nrow(lambda) # Number of variables
  q <- ncol(lambda) # Number of factors
  free_indices_phi <- free_indices_psi <- integer(0)

  indexes_lambda <- indexes_phi <- indexes_psi <-
    indexes_target <- indexes_targetphi <- indexes_targetpsi <-
    lambda_hat <- phi_hat <- psi_hat <- uniquenesses_hat <-
    free_indices_phi <- free_indices_psi <- list()
  indexes_factorvars <- indexes_uniquenesses <- c()
  original_targetphi <- phi
  original_targetpsi <- psi

  # Find which elements in the targets correspond to a parameter:
  indexes_target <- which(lambda %in% parameter_vector) # Which lambdas are estimated
  if(positive) {
    indexes_targetphi <- which(phi %in% parameter_vector) # Which phis are estimated
    free_indices_phi <- which(is.na(suppressWarnings(as.numeric(diag(original_targetphi)))))
    indexes_targetpsi <- which(psi %in% parameter_vector) # Which psis are estimated
    free_indices_psi <- which(is.na(suppressWarnings(as.numeric(diag(original_targetpsi)))))
  } else {
    indexes_targetphi <- which(phi %in% parameter_vector & lower.tri(phi, diag = TRUE)) # Which phis are estimated
    indexes_targetpsi <- which(psi %in% parameter_vector & lower.tri(psi, diag = TRUE)) # Which psis are estimated
  }

  # Get the indexes for the factor variances:
  indexes_factorvars <- c(indexes_factorvars,
                          which(parameter_vector %in% diag(phi)))

  # Get the indexes for the uniquenesses:
  indexes_uniquenesses <- c(indexes_uniquenesses,
                            which(parameter_vector %in% diag(psi)))

  # Relate the parameters in the targets to the parameters in the parameter vector:
  if(length(indexes_target) == 0) {
    indexes_lambda <- logical(0)
  } else {
    indexes_lambda <- match(lambda[indexes_target], parameter_vector) # Which lambdas are estimated
  }
  if(length(indexes_targetphi) == 0) {
    indexes_phi <- logical(0)
  } else {
    indexes_phi <- match(phi[indexes_targetphi], parameter_vector) # Which phis are estimated
  }
  if(length(indexes_targetpsi) == 0) {
    indexes_psi <- logical(0)
  } else {
    indexes_psi <- match(psi[indexes_targetpsi], parameter_vector) # Which psis are estimated
  }

  # Find which elements in the targets for correspond to a fixed value:
  indexes_fixtarget <- which(lambda %in% fixed_vector) # Which lambdas are fixed
  indexes_fixtargetphi <- which(phi %in% fixed_vector) # Which phis are fixed
  indexes_fixtargetpsi <- which(psi %in% fixed_vector) # Which psis are fixed

  # Relate the elements in the targets to the fixed values in the fixed vector:
  if(length(indexes_fixtarget) == 0) {
    indexes_fixlambda <- logical(0)
  } else {
    indexes_fixlambda <- match(lambda[indexes_fixtarget], fixed_vector) # Which lambdas are fixed
  }
  if(length(indexes_fixtargetphi) == 0) {
    indexes_fixphi <- logical(0)
  } else {
    indexes_fixphi <- match(phi[indexes_fixtargetphi], fixed_vector) # Which phis are fixed
  }
  if(length(indexes_fixtargetpsi) == 0) {
    indexes_fixpsi <- logical(0)
  } else {
    indexes_fixpsi <- match(psi[indexes_fixtargetpsi], fixed_vector) # Which psis are fixed
  }

  # non-specified elements in lambda are 0:
  lambda_hat <- matrix(0, p, q)
  lambda_hat[indexes_fixtarget] <- fixed_vector[indexes_fixlambda]
  # non-specified elements in phi are zero if off-diagonal and 1 if diagonal:
  phi_hat <- matrix(0, q, q); #diag(phi_hat) <- 1
  psi_hat <- matrix(0, p, p); #diag(psi_hat) <- 1
  if(positive) {
    phi_hat[fill_Phi_Target] <- 1
    psi_hat[fill_Psi_Target] <- 1
  }
  # non-specified elements in phi and psi are zero if off-diagonal and estimated if diagonal:
  phi_hat[indexes_fixtargetphi] <- fixed_vector[indexes_fixphi]
  psi_hat[indexes_fixtargetpsi] <- fixed_vector[indexes_fixpsi]

  class(lambda_hat) <- "numeric"
  class(phi_hat) <- "numeric"
  class(psi_hat) <- "numeric"

  # Initial lambda and uniqueness values based on the eigendecomposition of
  # the reduced correlation matrix:
  if(nrow(data) == ncol(data)) {
    S <- data
  } else {
    S <- cor(data, use = missing)
  }
  R <- S
  u <- 1/diag(solve(S))
  diag(S) <- u
  e <- eigen(S)
  D <- matrix(0, q, q)
  diag(D) <- sqrt(e$values[1:q])
  V <- e$vectors[, 1:q]
  VD <- V %*% D
  VV <- VD %*% t(VD)
  # init[indexes_lambda] <- VD[indexes_target]
  # init[indexes_phi] <- diag(q)[indexes_targetphi]
  # init[indexes_psi] <- diag(u)[indexes_targetpsi]
  matrices <- list(lambda = VD, phi = diag(q), psi = diag(u))

  # indices <- sort(unique(c(indexes_lambda, indexes_phi, indexes_psi)))
  duplicated_indices <- c(indexes_lambda, indexes_phi, indexes_psi)
  indices <- unique(duplicated_indices)
  nparam <- length(indices)
  logdetR <- log(det(R))

  result <- list(estimator = "cfa_ml",
                 # data = data,
                 # nfactors = q,
                 # nobs = nobs,
                 # free_indices_phi = free_indices_phi,
                 # free_indices_psi = free_indices_psi,
                 R = R,
                 logdetR = logdetR,
                 lambda = lambda_hat,
                 phi = phi_hat,
                 psi = psi_hat,
                 p = nrow(R),
                 nparam = nparam,
                 lambda_indexes = match(indexes_lambda, indices)-1,
                 phi_indexes = match(indexes_phi, indices)-1,
                 psi_indexes = match(indexes_psi, indices)-1,
                 target_indexes = indexes_target-1,
                 targetphi_indexes = indexes_targetphi-1,
                 targetpsi_indexes = indexes_targetpsi-1,
                 indices = indices-1,
                 duplicated_indices = duplicated_indices-1,
                 matrices = matrices)

  return(result)

}
