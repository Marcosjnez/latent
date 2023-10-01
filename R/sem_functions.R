MGf_minres <- function(x, S, ldetS, Inv_W, p, q, Lambda, Phi, Psi,
                       indexes_lambda, indexes_phi, indexes_psi,
                       indexes_target, indexes_targetphi, indexes_targetpsi) {

  # Objective function for minres

  n <- length(S) # Number of groups
  f <- 0 # Objective value (to be incremented)

  # q is not required to be constant across the groups
  # p is not required to be constant across the groups

  for(i in 1:n) {

    Lambda[[i]][indexes_target[[i]]] <- x[indexes_lambda[[i]]]
    Phi[[i]][indexes_targetphi[[i]]] <- x[indexes_phi[[i]]]
    Psi[[i]][indexes_targetpsi[[i]]] <- x[indexes_psi[[i]]]
    Rhat <- Lambda[[i]] %*% Phi[[i]] %*% t(Lambda[[i]]) + Psi[[i]]
    res <- S[[i]] - Rhat
    f <- f + 0.5*sum(res*res)

  }

  return(f)

}
MGg_minres <- function(x, S, ldetS, Inv_W, p, q, Lambda, Phi, Psi,
                       indexes_lambda, indexes_phi, indexes_psi,
                       indexes_target, indexes_targetphi, indexes_targetpsi) {

  # Gradient for minres

  n <- length(S) # Number of groups
  g <- rep(0, times = length(x)) # gradient (to be incremented)

  # q is not required to be constant across the groups
  # p is not required to be constant across the groups

  for(i in 1:n) {

    Lambda[[i]][indexes_target[[i]]] <- x[indexes_lambda[[i]]]
    Phi[[i]][indexes_targetphi[[i]]] <- x[indexes_phi[[i]]]
    Psi[[i]][indexes_targetpsi[[i]]] <- x[indexes_psi[[i]]]
    Rhat <- Lambda[[i]] %*% Phi[[i]] %*% t(Lambda[[i]]) + Psi[[i]]
    res <- S[[i]] - Rhat

    g1 <- (res %*% Lambda[[i]] %*% Phi[[i]])[indexes_target[[i]]]
    g2 <- (t(Lambda[[i]]) %*% res %*% Lambda[[i]])[indexes_targetphi[[i]]]
    res2 <- res
    res2[lower.tri(res2)] <- 2*res[lower.tri(res)]
    res2[upper.tri(res2)] <- 2*res[upper.tri(res)]
    g3 <- 0.5*res2[indexes_targetpsi[[i]]]

    g[indexes_lambda[[i]]] <- g[indexes_lambda[[i]]] -2*c(g1)
    g[indexes_phi[[i]]] <- g[indexes_phi[[i]]] -2*c(g2)
    g[indexes_psi[[i]]] <- g[indexes_psi[[i]]] -2*c(g3)

  }

  return(g)

}
MGf_ml <- function(x, S, ldetS, Inv_W, p, q, Lambda, Phi, Psi,
                   indexes_lambda, indexes_phi, indexes_psi,
                   indexes_target, indexes_targetphi, indexes_targetpsi) {

  # Objective function for ml

  n <- length(S) # Number of groups
  f <- 0 # Objective value (to be incremented)

  # q is not required to be constant across the groups
  # p is not required to be constant across the groups

  for(i in 1:n) {

    Lambda[[i]][indexes_target[[i]]] <- x[indexes_lambda[[i]]]
    Phi[[i]][indexes_targetphi[[i]]] <- x[indexes_phi[[i]]]
    Psi[[i]][indexes_targetpsi[[i]]] <- x[indexes_psi[[i]]]
    Rhat <- Lambda[[i]] %*% Phi[[i]] %*% t(Lambda[[i]]) + Psi[[i]]
    detRhat <- det(Rhat)
    ldetRhat <- log(detRhat)

    f <- f + ldetRhat - ldetS[[i]] + sum(S[[i]]*solve(Rhat)) - p[i]

  }

  return(f)

}
MGg_ml <- function(x, S, ldetS, Inv_W, p, q, Lambda, Phi, Psi,
                   indexes_lambda, indexes_phi, indexes_psi,
                   indexes_target, indexes_targetphi, indexes_targetpsi) {

  # Gradient for minres

  n <- length(S) # Number of groups
  g <- rep(0, times = length(x)) # gradient (to be incremented)

  # q is not required to be constant across the groups
  # p is not required to be constant across the groups

  for(i in 1:n) {

    Lambda[[i]][indexes_target[[i]]] <- x[indexes_lambda[[i]]]
    Phi[[i]][indexes_targetphi[[i]]] <- x[indexes_phi[[i]]]
    Psi[[i]][indexes_targetpsi[[i]]] <- x[indexes_psi[[i]]]
    Rhat <- Lambda[[i]] %*% Phi[[i]] %*% t(Lambda[[i]]) + Psi[[i]]

    Rhat_inv <- solve(Rhat)
    Ri_res_Ri <- 2*Rhat_inv %*% (Rhat - S[[i]]) %*% Rhat_inv
    dlogdetRhatdU <- 2*Rhat_inv
    diag(dlogdetRhatdU) <- 0.5*diag(dlogdetRhatdU)
    dRhat_invdU <- 2*(-Rhat_inv %*% S[[i]] %*% Rhat_inv)
    diag(dRhat_invdU) <- 0.5*diag(dRhat_invdU)

    g1 <- c(Ri_res_Ri %*% Lambda[[i]] %*% Phi[[i]])[indexes_target[[i]]]
    g2 <- c(t(Lambda[[i]]) %*% Ri_res_Ri %*% Lambda[[i]])[indexes_targetphi[[i]]]
    g3 <- c(dRhat_invdU + dlogdetRhatdU)[indexes_targetpsi[[i]]]

    g[indexes_lambda[[i]]] <- g[indexes_lambda[[i]]] + c(g1)
    g[indexes_phi[[i]]] <- g[indexes_phi[[i]]] + c(g2)
    g[indexes_psi[[i]]] <- g[indexes_psi[[i]]] + c(g3)

  }


  return(g)

}
MGf_dwls <- function(x, S, ldetS, Inv_W, p, q, Lambda, Phi, Psi,
                       indexes_lambda, indexes_phi, indexes_psi,
                       indexes_target, indexes_targetphi, indexes_targetpsi) {

  # Objective function for minres

  n <- length(S) # Number of groups
  f <- 0 # Objective value (to be incremented)

  # q is not required to be constant across the groups
  # p is not required to be constant across the groups

  for(i in 1:n) {

    Lambda[[i]][indexes_target[[i]]] <- x[indexes_lambda[[i]]]
    Phi[[i]][indexes_targetphi[[i]]] <- x[indexes_phi[[i]]]
    Psi[[i]][indexes_targetpsi[[i]]] <- x[indexes_psi[[i]]]
    Rhat <- Lambda[[i]] %*% Phi[[i]] %*% t(Lambda[[i]]) + Psi[[i]]
    res <- S[[i]] - Rhat
    f <- f + 0.5*sum(res*res * Inv_W[[i]])

  }

  return(f)

}
MGg_dwls <- function(x, S, ldetS, Inv_W, p, q, Lambda, Phi, Psi,
                       indexes_lambda, indexes_phi, indexes_psi,
                       indexes_target, indexes_targetphi, indexes_targetpsi) {

  # Gradient for minres

  n <- length(S) # Number of groups
  g <- rep(0, times = length(x)) # gradient (to be incremented)

  # q is not required to be constant across the groups
  # p is not required to be constant across the groups

  for(i in 1:n) {

    Lambda[[i]][indexes_target[[i]]] <- x[indexes_lambda[[i]]]
    Phi[[i]][indexes_targetphi[[i]]] <- x[indexes_phi[[i]]]
    Psi[[i]][indexes_targetpsi[[i]]] <- x[indexes_psi[[i]]]
    Rhat <- Lambda[[i]] %*% Phi[[i]] %*% t(Lambda[[i]]) + Psi[[i]]
    res <- S[[i]] - Rhat

    g1 <- ((Inv_W[[i]] * res) %*% Lambda[[i]] %*% Phi[[i]])[indexes_target[[i]]]
    g2 <- (t(Lambda[[i]]) %*% (Inv_W[[i]] * res) %*% Lambda[[i]])[indexes_targetphi[[i]]]
    res2 <- res
    res2[lower.tri(res2)] <- 2*res[lower.tri(res)]
    res2[upper.tri(res2)] <- 2*res[upper.tri(res)]
    g3 <- 0.5*res2[indexes_targetpsi[[i]]] * Inv_W[[i]][indexes_targetpsi[[i]]]

    g[indexes_lambda[[i]]] <- g[indexes_lambda[[i]]] -2*c(g1)
    g[indexes_phi[[i]]] <- g[indexes_phi[[i]]] -2*c(g2)
    g[indexes_psi[[i]]] <- g[indexes_psi[[i]]] -2*c(g3)

  }

  return(g)

}
MGCFA <- function(S, target, targetphi, targetpsi,
                  method = "minres", W = NULL) {

  # S is a list of correlation matrices (one for each group)
  # target is a list of matrices for the loadings indicating the parameters
  # (as characters) and fixed values (numeric)
  # targetphi is a list of matrices for the factor correlations indicating the
  # parameters (as characters) and fixed values (numeric)
  # targetpsi is a list of matrices for the uniquenesses and correlated errors
  # indicating the parameters (as characters) and fixed values (numeric)

  # The targets must be fully specified, that is, every entry in a target
  # must contain either a character or a numeric value

  n <- length(S) # Number of groups
  p <- q <- vector(length = n) # Number of variables (p) and factors (q) in each group
  Lambda <- Phi <- Psi <-
    indexes_lambda <- indexes_phi <- indexes_psi <-
    indexes_target <- indexes_targetphi <- indexes_targetpsi <-
    lambda_hat <- phi_hat <- psi_hat <- uniquenesses_hat <-
    S_hat <- residuals <- list()
  indexes_factorvars <- indexes_uniquenesses <- c()

  # Find the unique elements in all the targets:
  uniques <- unique(unlist(c(target, targetphi, targetpsi)))
  # Get the parameters (those elements that are not digits):
  z <- suppressWarnings(as.numeric(uniques))
  parameter_vector <- uniques[which(is.na(z))]
  # Get the fixed values (those elements that are digits):
  fixed_vector <- uniques[which(!is.na(z))]
  # Initialize the vector of initial parameter estimates:
  init <- vector(length = length(parameter_vector))

  for(i in 1:n) { # For each group...

    p[i] <- nrow(target[[i]]) # Number of variables in group i
    q[i] <- ncol(target[[i]]) # Number of factors in group i

    # Find which elements in the targets correspond to a parameter:
    indexes_target[[i]] <- which(target[[i]] %in% parameter_vector) # Which lambdas are estimated in group i
    indexes_targetphi[[i]] <- which(targetphi[[i]] %in% parameter_vector) # Which phis are estimated in group i
    indexes_targetpsi[[i]] <- which(targetpsi[[i]] %in% parameter_vector) # Which psis are estimated in group i

    # Get the indexes for the factor variances:
    indexes_factorvars <- c(indexes_factorvars,
                              which(parameter_vector %in% diag(targetphi[[i]])))

    # Get the indexes for the uniquenesses:
    indexes_uniquenesses <- c(indexes_uniquenesses,
                              which(parameter_vector %in% diag(targetpsi[[i]])))

    # Relate the parameters in the targets to the parameters in the parameter vector:
    if(length(indexes_target[[i]]) == 0) {
      indexes_lambda[[i]] <- logical(0)
    } else {
      indexes_lambda[[i]] <- match(target[[i]][indexes_target[[i]]], parameter_vector) # Which lambdas are estimated in group i
    }
    if(length(indexes_targetphi[[i]]) == 0) {
      indexes_phi[[i]] <- logical(0)
    } else {
      indexes_phi[[i]] <- match(targetphi[[i]][indexes_targetphi[[i]]], parameter_vector) # Which phis are estimated in group i
    }
    if(length(indexes_targetpsi[[i]]) == 0) {
      indexes_psi[[i]] <- logical(0)
    } else {
      indexes_psi[[i]] <- match(targetpsi[[i]][indexes_targetpsi[[i]]], parameter_vector) # Which psis are estimated in group i
    }

    # Find which elements in the targets for correspond to a fixed value:
    indexes_fixtarget <- which(target[[i]] %in% fixed_vector) # Which lambdas are fixed in group i
    indexes_fixtargetphi <- which(targetphi[[i]] %in% fixed_vector) # Which phis are fixed in group i
    indexes_fixtargetpsi <- which(targetpsi[[i]] %in% fixed_vector) # Which psis are fixed in group i

    # Relate the elements in the targets to the fixed values in the fixed vector:
    if(length(indexes_fixtarget) == 0) {
      indexes_fixlambda <- logical(0)
    } else {
      indexes_fixlambda <- match(target[[i]][indexes_fixtarget], fixed_vector) # Which lambdas are fixed in group i
    }
    if(length(indexes_fixtargetphi) == 0) {
      indexes_fixphi <- logical(0)
    } else {
      indexes_fixphi <- match(targetphi[[i]][indexes_fixtargetphi], fixed_vector) # Which phis are fixed in group i
    }
    if(length(indexes_fixtargetpsi) == 0) {
      indexes_fixpsi <- logical(0)
    } else {
      indexes_fixpsi <- match(targetpsi[[i]][indexes_fixtargetpsi], fixed_vector) # Which psis are fixed in group i
    }

    # non-specified elements in lambda are 0:
    lambda_hat[[i]] <- matrix(0, p[i], q[i])
    lambda_hat[[i]][indexes_fixtarget] <- fixed_vector[indexes_fixlambda]
    # non-specified elements in phi are zero if off-diagonal and 1 if diagonal:
    phi_hat[[i]] <- matrix(0, q[i], q[i]); diag(phi_hat[[i]]) <- 1
    phi_hat[[i]][indexes_fixtargetphi] <- fixed_vector[indexes_fixphi]
    # non-specified elements in psi are zero if off-diagonal and estimated if diagonal:
    psi_hat[[i]] <- matrix(0, p[i], p[i])
    psi_hat[[i]][indexes_fixtargetpsi] <- fixed_vector[indexes_fixpsi]

    class(lambda_hat[[i]]) <- "numeric"
    class(phi_hat[[i]]) <- "numeric"
    class(psi_hat[[i]]) <- "numeric"

    # Initial lambda and uniqueness values based on the eigendecomposition of
    # the reduced correlation matrix:
    Si <- S[[i]]
    u <- 1/diag(solve(Si))
    diag(Si) <- u
    e <- eigen(Si)
    D <- matrix(0, q[i], q[i])
    diag(D) <- sqrt(e$values[1:q[i]])
    V <- e$vectors[, 1:q[i]]
    VD <- V %*% D
    VV <- VD %*% t(VD)
    init[indexes_lambda[[i]]] <- VV[indexes_target[[i]]]
    init[indexes_phi[[i]]] <- diag(q[i])[indexes_targetphi[[i]]]
    init[indexes_psi[[i]]] <- diag(u)[indexes_targetpsi[[i]]]

  }

  indexes_factorvars <- unique(indexes_factorvars)
  indexes_uniquenesses <- unique(indexes_uniquenesses)
  lambda_p <- length(unique(unlist(indexes_lambda)))
  phi_p <- length(unique(unlist(indexes_phi)))
  psi_p <- length(unique(unlist(indexes_psi)))

  lower_psi <- rep(-0.995, psi_p) # Lower bounds for correlated residuals
  upper_psi <- rep(0.995, psi_p) # Upper bounds for correlated residuals
  lower <- c(rep(-Inf, lambda_p), rep(-0.995, phi_p), lower_psi)
  upper <- c(rep(Inf, lambda_p), rep(0.995, phi_p), upper_psi)
  lower[indexes_uniquenesses] <- 0.001
  upper[indexes_uniquenesses] <- 0.999
  lower[indexes_factorvars] <- 0.001
  upper[indexes_factorvars] <- Inf

  # lower_psi <- rep(-Inf, psi_p) # Lower bounds for correlated residuals
  # upper_psi <- rep(Inf, psi_p) # Upper bounds for correlated residuals
  # lower <- c(rep(-Inf, lambda_p), rep(-Inf, phi_p), lower_psi)
  # upper <- c(rep(Inf, lambda_p), rep(Inf, phi_p), upper_psi)
  # lower[indexes_uniquenesses] <- 0.005

  # x <- init
  x <- c(stats::runif(lambda_p), rep(0, phi_p), rep(0.05, psi_p))

  if(method == "minres") {

    ldetS <- NULL
    Inv_W <- NULL
    f <- MGf_minres
    g <- MGg_minres

  } else if(method == "ml") {

    ldetS <- list()
    Inv_W <- NULL
    for(i in 1:n) ldetS[[i]] <- log(det(S[[i]]))
    f <- MGf_ml
    g <- MGg_ml

  } else if(method == "dwls") {

    ldetS <- NULL
    Inv_W <- list()
    for(i in 1:n) Inv_W[[i]] <- 1/W[[i]]
    f <- MGf_dwls
    g <- MGg_dwls


  }

  cfa <- stats::nlminb(x, objective = f, gradient = g,
                       lower = lower, upper = upper,
                       S = S, ldetS = ldetS, p = p, q = q,
                       Lambda = lambda_hat, Phi = phi_hat, Psi = psi_hat,
                       indexes_target = indexes_target,
                       indexes_targetphi = indexes_targetphi,
                       indexes_targetpsi = indexes_targetpsi,
                       indexes_lambda = indexes_lambda,
                       indexes_phi = indexes_phi, indexes_psi = indexes_psi,
                       Inv_W = Inv_W,
                       control = list(iter.max = 1e4, eval.max = 1e4))

  # Arrange the parameter estimates in the lambda, phi, and psi matrices:
  for(i in 1:n) {

    # Arrange lambda parameter estimates:
    lambda_hat[[i]][indexes_target[[i]]] <- cfa$par[indexes_lambda[[i]]]

    # Arrange phi parameter estimates:
    phi_hat[[i]][indexes_targetphi[[i]]] <- cfa$par[indexes_phi[[i]]]

    # Arrange psi parameter estimates:
    psi_hat[[i]][indexes_targetpsi[[i]]] <- cfa$par[indexes_psi[[i]]]

    # Model matrix:
    S_hat[[i]] <- lambda_hat[[i]] %*% phi_hat[[i]] %*% t(lambda_hat[[i]]) + psi_hat[[i]]
    uniquenesses_hat[[i]] <- diag(psi_hat[[i]])
    diag(S_hat[[i]]) <- diag(S[[i]]) # Fix rounding errors from the optimization
    residuals[[i]] <- S[[i]] - S_hat[[i]]

  }

  # Degrees of freedom:
  df <- sum(p*(p+1)/2) - length(parameter_vector)

  results <- list(f = cfa$objective, convergence = cfa$convergence,
                  iterations = cfa$iterations, df = df,
                  lambda = lambda_hat, phi = phi_hat,
                  psi = psi_hat, uniquenesses = uniquenesses_hat,
                  model = S_hat, residuals = residuals,
                  parameter_vector = cfa$par)

  return(results)

}
