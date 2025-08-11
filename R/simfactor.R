# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 10/08/2025

grid_search <- function(delta, G) {

  n <- 1000
  x <- seq(-1e3, 1e3, length.out = n)
  y <- vector(length = n)
  for(i in 1:n) y[i] <- root_ml(x[i], delta, G)

  index <- which.min(y)
  x <- seq(x[index-1], x[index+1], length.out = n)
  for(i in 1:n) y[i] <- root_ml(x[i], delta, G)

  index <- which.min(y)
  x <- seq(x[index-1], x[index+1], length.out = n)
  for(i in 1:n) y[i] <- root_ml(x[i], delta, G)

  return(x[which.min(y)])

}
opt <- function(x, delta, G) {

  # x <- stats::runif(1, 0, 1)
  # det(diag(nrow(G)) + x*G)
  root <- stats::optim(x, fn = root_ml, gr = groot_ml, estimator = "L-BFGS-B",
                       lower = -Inf, upper = Inf, G = G, delta = delta)
  k <- root$par

  return(k)

}
opt_error <- function(x, delta, G) {

  x <- tryCatch({opt(x, delta, G)}, error = return(x))
  return(x)

}
f_minres <- function(x, S, ldetS, Inv_W, q, indexes_lambda, indexes_phi, indexes_psi) {

  p <- nrow(S)
  lambda_p <- length(indexes_lambda)
  Lambda <- matrix(0, p, q)
  Lambda[indexes_lambda] <- x[1:lambda_p]
  phi_p <- length(indexes_phi)
  Phi <- matrix(0, q, q)
  Phi[indexes_phi] <- x[(lambda_p+1):(lambda_p + phi_p)]
  Phi <- t(Phi) + Phi
  diag(Phi) <- 1
  Psi <- matrix(0, p, p)
  Psi[indexes_psi] <- x[-(1:(lambda_p + phi_p))]
  Psi[upper.tri(Psi)] <- t(Psi)[upper.tri(Psi)]
  Rhat <- Lambda %*% Phi %*% t(Lambda) + Psi
  res <- S - Rhat
  f <- 0.5*sum(res*res)

  return(f)

}
g_minres <- function(x, S, ldetS, Inv_W, q, indexes_lambda, indexes_phi, indexes_psi) {

  p <- nrow(S)
  lambda_p <- length(indexes_lambda)
  Lambda <- matrix(0, p, q)
  Lambda[indexes_lambda] <- x[1:lambda_p]
  phi_p <- length(indexes_phi)
  Phi <- matrix(0, q, q)
  Phi[indexes_phi] <- x[(lambda_p+1):(lambda_p + phi_p)]
  Phi <- t(Phi) + Phi
  diag(Phi) <- 1
  Psi <- matrix(0, p, p)
  Psi[indexes_psi] <- x[-(1:(lambda_p + phi_p))]
  Psi[upper.tri(Psi)] <- t(Psi)[upper.tri(Psi)]
  Rhat <- Lambda %*% Phi %*% t(Lambda) + Psi
  res <- S - Rhat

  g1 <- (res %*% Lambda %*% Phi)[indexes_lambda]
  g2 <- (t(Lambda) %*% res %*% Lambda)[indexes_phi]
  # g <- -2*c(g1, g2, 0.5*diag(res))
  res2 <- res
  res2[lower.tri(res2)] <- 2*res[lower.tri(res)]
  g <- -2*c(g1, g2, 0.5*res2[indexes_psi])

  return(g)

}
f_ml <- function(x, S, ldetS, Inv_W, q, indexes_lambda, indexes_phi, indexes_psi) {

  p <- nrow(S)
  lambda_p <- length(indexes_lambda)
  Lambda <- matrix(0, p, q)
  Lambda[indexes_lambda] <- x[1:lambda_p]
  phi_p <- length(indexes_phi)
  Phi <- matrix(0, q, q)
  Phi[indexes_phi] <- x[(lambda_p+1):(lambda_p + phi_p)]
  Phi <- t(Phi) + Phi
  diag(Phi) <- 1
  Psi <- matrix(0, p, p)
  Psi[indexes_psi] <- x[-(1:(lambda_p + phi_p))]
  Psi[upper.tri(Psi)] <- t(Psi)[upper.tri(Psi)]
  Rhat <- Lambda %*% Phi %*% t(Lambda) + Psi
  f <- log(det(Rhat)) - ldetS + sum(S*solve(Rhat)) - p

  return(f)

}
g_ml <- function(x, S, ldetS, Inv_W, q, indexes_lambda, indexes_phi, indexes_psi) {

  p <- nrow(S)
  lambda_p <- length(indexes_lambda)
  Lambda <- matrix(0, p, q)
  Lambda[indexes_lambda] <- x[1:lambda_p]
  phi_p <- length(indexes_phi)
  Phi <- matrix(0, q, q)
  Phi[indexes_phi] <- x[(lambda_p+1):(lambda_p + phi_p)]
  Phi <- t(Phi) + Phi
  diag(Phi) <- 1
  Psi <- matrix(0, p, p)
  Psi[indexes_psi] <- x[-(1:(lambda_p + phi_p))]
  Psi[upper.tri(Psi)] <- t(Psi)[upper.tri(Psi)]

  Rhat <- Lambda %*% Phi %*% t(Lambda) + Psi
  Rhat_inv <- solve(Rhat)
  Ri_res_Ri <- 2*Rhat_inv %*% (Rhat - S) %*% Rhat_inv
  Ri_res_Ri2 <- Ri_res_Ri
  Ri_res_Ri2[lower.tri(Ri_res_Ri2)] <- 2*Ri_res_Ri[lower.tri(Ri_res_Ri)]

  # Joreskog (page 10; 1965) Testing a simple structure in factor analysis
  # g <- c(c(Ri_res_Ri %*% Lambda %*% Phi)[indexes_lambda],
  #        c(t(Lambda) %*% Ri_res_Ri %*% Lambda)[indexes_phi],
  #        diag(Ri_res_Ri)*0.5)
  g <- c(c(Ri_res_Ri %*% Lambda %*% Phi)[indexes_lambda],
         c(t(Lambda) %*% Ri_res_Ri %*% Lambda)[indexes_phi],
         Ri_res_Ri2[indexes_psi]*0.5)

  return(g)

}
f_dwls <- function(x, S, ldetS, Inv_W, q, indexes_lambda, indexes_phi, indexes_psi) {

  p <- nrow(S)
  lambda_p <- length(indexes_lambda)
  Lambda <- matrix(0, p, q)
  Lambda[indexes_lambda] <- x[1:lambda_p]
  phi_p <- length(indexes_phi)
  Phi <- matrix(0, q, q)
  Phi[indexes_phi] <- x[(lambda_p+1):(lambda_p + phi_p)]
  Phi <- t(Phi) + Phi
  diag(Phi) <- 1
  Psi <- matrix(0, p, p)
  Psi[indexes_psi] <- x[-(1:(lambda_p + phi_p))]
  Psi[upper.tri(Psi)] <- t(Psi)[upper.tri(Psi)]
  Rhat <- Lambda %*% Phi %*% t(Lambda) + Psi
  res <- S - Rhat
  f <- 0.5*sum(res*res*Inv_W)

  return(f)

}
g_dwls <- function(x, S, ldetS, Inv_W, q, indexes_lambda, indexes_phi, indexes_psi) {

  p <- nrow(S)
  lambda_p <- length(indexes_lambda)
  Lambda <- matrix(0, p, q)
  Lambda[indexes_lambda] <- x[1:lambda_p]
  phi_p <- length(indexes_phi)
  Phi <- matrix(0, q, q)
  Phi[indexes_phi] <- x[(lambda_p+1):(lambda_p + phi_p)]
  Phi <- t(Phi) + Phi
  diag(Phi) <- 1
  Psi <- matrix(0, p, p)
  Psi[indexes_psi] <- x[-(1:(lambda_p + phi_p))]
  Psi[upper.tri(Psi)] <- t(Psi)[upper.tri(Psi)]
  Rhat <- Lambda %*% Phi %*% t(Lambda) + Psi
  res <- S - Rhat

  g1 <- ((Inv_W * res) %*% Lambda %*% Phi)[indexes_lambda]
  g2 <- (t(Lambda) %*% (Inv_W*res) %*% Lambda)[indexes_phi]
  # g <- -2*c(g1, g2, 0.5*diag(res))
  res2 <- res
  res2[lower.tri(res2)] <- 2*res[lower.tri(res)]
  g <- -2*c(g1, g2, 0.5*res2[indexes_psi] * Inv_W[indexes_psi])

  return(g)

}

CFA <- function(S, target, targetphi, targetpsi = diag(nrow(target)),
                estimator = "minres", W = NULL) {

  p <- nrow(target)
  q <- ncol(target)
  indexes_lambda <- which(target != 0) # Which lambdas are estimated
  indexes_phi <- which(targetphi != 0 & lower.tri(targetphi)) # Which phis are estimated
  indexes_psi <- which(targetpsi != 0 & lower.tri(targetpsi, diag = TRUE)) # Which psies are estimated
  lambda_p <- length(indexes_lambda) # Number of lambda parameters
  phi_p <- length(indexes_phi) # Number of phi parameters
  psi_p <- length(indexes_psi) # Number of psi parameters

  init_diag_psi <- 1/diag(solve(S)) # Initial diagonal psi parameter values
  init_psi <- rep(0, times = psi_p)
  diag_indexes <- (p+1)*0:(p-1)+1 # Indexes for the diagonal of Psi
  offdiag_indexes <- which(targetpsi != 0 & lower.tri(targetpsi)) # Indexes for the off-diagonal of Psi
  cor_res_indexes <- which(indexes_psi %in% offdiag_indexes) # Indexes for correlated residuals
  # Allocate init_diag_psi in the positions of the vector corresponding to the diagonal of Psi:
  init_psi[-cor_res_indexes] <- init_diag_psi

  lower_psi <- rep(0.005, psi_p) # Lower bounds for the uniquenessess
  lower_psi[cor_res_indexes] <- -0.995 # Lower bounds for correlated residuals
  upper_psi <- rep(0.995, psi_p) # Upper bounds for correlated residuals
  lower <- c(rep(-Inf, lambda_p), rep(-1, phi_p), lower_psi)
  upper <- c(rep(Inf, lambda_p), rep(1, phi_p), upper_psi)

  x <- c(stats::runif(lambda_p), rep(0, phi_p), init_psi)

  if(estimator == "minres") {

    ldetS <- NULL
    Inv_W <- NULL
    f <- f_minres
    g <- g_minres

  } else if(estimator == "ml") {

    ldetS <- log(det(S))
    Inv_W <- NULL
    f <- f_ml
    g <- g_ml

  } else if(estimator == "dwls") {

    ldetS <- NULL
    Inv_W <- 1/W
    Inv_W[is.infinite(Inv_W)] <- 0
    f <- f_dwls
    g <- g_dwls

  }

  cfa <- stats::nlminb(x, objective = f, gradient = g,
                       lower = lower, upper = upper,
                       S = S, ldetS = ldetS, q = q, indexes_lambda = indexes_lambda,
                       indexes_phi = indexes_phi, indexes_psi = indexes_psi,
                       Inv_W = Inv_W,
                       control = list(iter.max = 1e4, eval.max = 1e4))

  # Arrange lambda parameter estimates:
  lambda_hat <- matrix(0, p, q)
  lambda_hat[indexes_lambda] <- cfa$par[1:lambda_p]

  # Arrange phi parameter estimates:
  phi_hat <- matrix(0, q, q)
  phi_hat[indexes_phi] <- cfa$par[(lambda_p+1):(lambda_p + phi_p)]
  phi_hat <- t(phi_hat) + phi_hat
  diag(phi_hat) <- 1

  # Arrange psi parameter estimates:
  psi_hat <- matrix(0, p, p)
  psi_hat[indexes_psi] <- cfa$par[-(1:(lambda_p + phi_p))]
  psi_hat[upper.tri(psi_hat)] <- t(psi_hat)[upper.tri(psi_hat)]

  # Model matrix:
  S_hat <- lambda_hat %*% phi_hat %*% t(lambda_hat) + psi_hat
  uniquenesses_hat <- 1 - diag(lambda_hat %*% phi_hat %*% t(lambda_hat))
  diag(S_hat) <- 1 # Fix rounding errors from the optimization
  residuals <- S - S_hat

  # Degrees of freedom:
  df <- p*(p+1)/2 - (lambda_p + phi_p + psi_p)

  results <- list(f = cfa$objective, convergence = cfa$convergence,
                  iterations = cfa$iterations, df = df,
                  lambda = lambda_hat, phi = phi_hat,
                  psi = psi_hat, uniquenesses = uniquenesses_hat,
                  model = S_hat, residuals = residuals, pars = cfa$par)

  return(results)

}
root_ml <- function(x, delta, G) {

  I <- diag(nrow(G))
  f <- x*sum(diag(G)) - log(det(I + x*G)) - delta

  return(f^2)

}
groot_ml <- function(x, delta, G) {

  I <- diag(nrow(G))
  f <- x*sum(diag(G)) - log(det(I + x*G)) - delta
  g <- sum(diag(G)) - sum(diag(solve(I + x*G) %*% G))
  g <- 2*g*f

  return(g)

}
dxt <- function(X) {

  # derivative wrt transpose (just a permutation matrix)

  p <- nrow(X)
  q <- ncol(X)
  pq <- p*q

  res <- array(0, dim = c(pq, pq))
  null <- matrix(0, p, q)

  for(i in 1:pq) {
    temp <- null
    temp[i] <- 1
    res[, i] <- c(t(temp))
  }

  return(res)

}
gLRhat <- function(Lambda, Phi) {

  # derivative of Lambda wrt Rhat

  p <- nrow(Lambda)
  g1 <- (Lambda %*% Phi) %x% diag(p)
  g21 <- diag(p) %x% (Lambda %*% Phi)
  g2 <- g21 %*% dxt(Lambda)
  g <- g1 + g2

  return(g)

}
gPRhat <- function(Lambda, Phi) {

  g1 <- Lambda %x% Lambda
  g2 <- g1 %*% dxt(Phi)
  g <- g1 + g2
  g <- g[, which(lower.tri(Phi)), drop = FALSE]

  return(g)

}
guRhat <- function(p) {

  gu <- matrix(0, p*p, p)

  for(i in 1:p) {

    index <- (i-1)*p + i
    gu[index, i] <- 1

  }

  return(gu)

}
gURhat <- function(p) {

  pcov <- p*(p+1)*0.5

  Psi <- diag(p)
  gPsi <- diag(p) %x% diag(p)
  gPsi <- gPsi + dxt(Psi) %*% gPsi
  gPsi <- gPsi[, lower.tri(Psi, diag = TRUE), drop = FALSE]
  gPsi[gPsi != 0] <- 1

  return(gPsi)

}
cudeck <- function(R, lambda, Phi, Psi,
                   fit = "rmsr", misfit = "close",
                   estimator = "minres", efa = FALSE) {

  # Cudeck, R., & Browne, M. W. (1992). Constructing a covariance matrix that
  # yields a specified minimizer and a specified minimum discrepancy function value.
  # Psychometrika, 57(3), 357–369. https://doi.org/10.1007/BF02295424

  p <- nrow(lambda)
  q <- ncol(lambda)
  uniquenesses <- diag(Psi)

  # Count the number of parameters
  nlambda <- sum(lambda != 0)
  nphi <- sum(Phi[lower.tri(Phi)] != 0)
  npsi <- sum(Psi[lower.tri(Psi, diag = TRUE)] != 0)
  npars <- nlambda + nphi + npsi
  df <- p*(p+1)/2 - npars # Degrees of freedom

  if(nlambda + nphi > p*q - 0.5*q*(q-1)) {
    warning("The model is not identified. There exists infinite solutions for the model parameters.")
  }

  if(nlambda + nphi + npsi > p*(p+1)/2) {
    warning("The true model has negative degrees of freedom.")
  }

  # Create the matrix of derivatives wrt the correlation model:

  if(efa) {

    dS_dL <- gLRhat(lambda, Phi)
    dS_dU <- gURhat(p)[, which(Psi[lower.tri(Psi, diag = TRUE)] != 0)]
    gS <- cbind(dS_dL, dS_dU)

  } else {

    dS_dL <- gLRhat(lambda, Phi)[, which(lambda != 0)]
    dS_dP <- gPRhat(lambda, Phi)[, which(Phi[lower.tri(Phi)] != 0)]
    dS_dU <- gURhat(p)[, which(Psi[lower.tri(Psi, diag = TRUE)] != 0)]
    gS <- cbind(dS_dL, dS_dP, dS_dU)

  }

  if(estimator == "minres" || estimator == "uls") {

    # Select the nonduplicated elements of the correlation matrix wrt each parameter
    B <- -gS[lower.tri(R, diag = TRUE), ]

  } else if(estimator == "ml") {

    # K <- transition(p)
    # MP_inv <- solve(t(K) %*% K) %*% t(K)
    # D <- MP_inv %*% t(MP_inv)
    indexes <- vector(length = p)
    indexes[1] <- 1
    for(i in 2:p) {
      increment <- i
      indexes[i] <- indexes[i-1]+increment
    }
    D <- matrix(0, p*(p+1)/2, p*(p+1)/2)
    diag(D) <- 2
    diag(D)[indexes] <- 1
    R_inv <- solve(R)
    # vecs <- apply(gS, 2, FUN = function(x) -t(R_inv %*% matrix(x, p, p) %*% R_inv))
    # B <- t(vecs[which(upper.tri(R, diag = TRUE)), ]) %*% D
    # B <- t(B)
    B <- -apply(gS, 2, FUN = function(x) t((R_inv %*% matrix(x, p, p) %*% R_inv)[which(upper.tri(R, diag = TRUE))]) %*% D)
    # The error must be orthogonal to the derivative of each parameter wrt correlation model

  }

  # Generate random error:

  m <- p+1
  U <- replicate(p, stats::runif(m, -1, 1))
  A1 <- t(U) %*% U
  sq <- diag(1/sqrt(diag(A1)))
  A2 <- sq %*% A1 %*% sq
  diag_u <- diag(sqrt(uniquenesses))
  y <- diag_u %*% A2 %*% diag_u
  y <- y[lower.tri(y, diag = TRUE)]
  # y <- A2[lower.tri(A2, diag = TRUE)]
  # v <- MASS::ginv(t(B) %*% B) %*% t(B) %*% y
  # e <- y - B %*% v # equation 7 from Cudeck and Browne (1992)
  Q <- qr.Q(qr(B))
  e <- y - Q %*% t(Q) %*% y

  # Adjust the error to satisfy the desired amount of misfit:

  if(estimator == "minres" || estimator == "uls") {

    E <- matrix(0, p, p)
    E[lower.tri(E, diag = TRUE)] <- e
    E <- t(E) + E
    # diag(E) <- diag(E)/2
    diag(E) <- 0

    if(fit == "rmsr") {
      if(misfit == "close") {
        r2 <- mean(1-uniquenesses)
        misfit <- 0.05*r2
      } else if(misfit == "acceptable") {
        r2 <- mean(1-uniquenesses)
        misfit <- 0.10*r2
      }
      delta <- misfit^2*0.5*p*(p-1)
      # delta <- (1-misfit2)*(0.5*(sum(R_error^2) - p))
    } else if(fit == "cfi") {
      null_f <- 0.5*(sum(R^2) - p)
      delta <- (1-misfit)*null_f
    } else if(fit == "rmsea") {
      delta <- misfit^2 * df
    } else if(fit == "raw") {
      delta <- misfit
    }

    k <- sqrt(2*delta/sum(E*E))
    E <- k*E

  } else if(estimator == "ml") {

    E <- matrix(0, p, p)
    E[upper.tri(R, diag = TRUE)] <- e
    E <- t(E) + E
    diag(E) <- 0

    if(fit == "rmsr") {
      if(misfit == "close") {
        r2 <- mean(1-uniquenesses)
        misfit <- 0.05*r2
      } else if(misfit == "acceptable") {
        r2 <- mean(1-uniquenesses)
        misfit <- 0.10*r2
      }
      delta <- "A given RMSR is compatible with multiple maximum likelihood discrepancy values and is not provided"
    } else if(fit == "cfi") {
      null_f <- -log(det(R))
      delta <- (1-misfit)*null_f
    } else if(fit == "rmsea") {
      delta <- misfit^2 * df
    } else if(fit == "raw") {
      delta <- misfit
    }

    if(fit == "rmsr") {

      k <- sqrt((0.5*p*(p-1))*2*misfit^2/sum(E*E))
      E <- k*E

    } else {

      constant <- 1e-04 / sqrt(mean(E*E))
      E <- constant*E # Fix this to avoid NAs
      R_inv <- solve(R)
      G <- R_inv %*% E
      x <- suppressWarnings(grid_search(delta, G))
      # x <- sqrt(2*delta/sum(G*G)) # Initial value suggested by Cudeck
      k <- opt(x, delta, G)
      # limits <- c(-1e05, 1e05)
      # k <- GSS(delta, G, limits)
      # k <- grad_descend(delta, G)
      E <- k*E

    }
  }

  R_error <- R + E

  if(estimator == "ml" & fit == "rmsr") {
    delta <- log(det(R)) - log(det(R_error)) + sum(R_error*solve(R)) - nrow(R)
  }

  # Check for positiveness:
  minimum_eigval <- min(eigen(R_error, symmetric = TRUE, only.values = TRUE)$values)
  positive <- TRUE
  if(minimum_eigval <= 0) {
    # warning("The matrix was not positive-definite. The amount of misfit may be too big.")
    # positive <- FALSE
    return(cudeck(R, lambda, Phi, Psi, fit, misfit, estimator))
  }

  return(list(R_error = R_error, fit = fit, delta = delta, misfit = misfit, positive = positive))

}
yuan <- function(R, lambda, Phi, Psi,
                 fit = "rmsr", misfit = "close",
                 estimator = "minres") {

  # Yuan, K.-H., & Hayashi, K. (2003). Bootstrap approach to inference and power
  # analysis based on three test statistics for covariance structure models.
  # British Journal of Mathematical and Statistical Psychology, 56(1), 93–110.
  # https://doi.org/10.1348/000711003321645368

  p <- nrow(R)
  q <- ncol(lambda)
  uniquenesses <- diag(Psi)

  # Count the number of parameters
  nlambda <- sum(lambda != 0)
  nphi <- sum(Phi[lower.tri(Phi)] != 0)
  npsi <- sum(Psi[lower.tri(Psi, diag = TRUE)] != 0)
  npars <- nlambda + nphi + npsi
  df <- p*(p+1)/2 - npars # Degrees of freedom

  if(nlambda + nphi > p*q - 0.5*q*(q-1)) {
    warning("The population model is not identified. There exists infinite solutions for the model parameters.")
  }

  if(nlambda + nphi + npsi > p*(p+1)/2) {
    warning("The population model has negative degrees of freedom.")
  }

  # Add an small error to the population parameters
  # lambda_error <- lambda - 1e-04
  # Rerror <- lambda_error %*% Phi %*% t(lambda_error) + Psi; diag(Rerror) <- 1
  Rerror <- R
  Rerror[lower.tri(R)] <- Rerror[lower.tri(R)] + stats::runif(0.5*p*(p-1), -1e-06, 1e-06)
  Rerror[upper.tri(R)] <- t(Rerror)[upper.tri(R)]

  # Create the FA model
  target <- ifelse(lambda != 0, 1, 0)
  targetphi <- ifelse(Phi != 0, 1, 0)
  targetpsi <- ifelse(Psi != 0, 1, 0)
  cfa <- CFA(Rerror, target, targetphi, targetpsi, estimator = estimator, W = NULL)
  Phat <- cfa$model

  # Get the error matrix:
  E <- Rerror - Phat
  # Hopefully, the error is orthogonal to the derivative of each parameter wrt correlation model

  # Adjust the error to satisfy the desired amount of misfit:

  if(estimator == "minres" || estimator == "uls") {

    if(fit == "rmsr") {
      if(misfit == "close") {
        r2 <- mean(1-uniquenesses)
        misfit <- 0.05*r2
      } else if(misfit == "acceptable") {
        r2 <- mean(1-uniquenesses)
        misfit <- 0.10*r2
      }
      delta <- misfit^2*0.5*p*(p-1)
      # delta <- (1-misfit2)*(0.5*(sum(R_error^2) - p))
    } else if(fit == "cfi") {
      null_f <- 0.5*(sum(R^2) - p)
      delta <- (1-misfit)*null_f
    } else if(fit == "rmsea") {
      delta <- misfit^2 * df
    } else if(fit == "raw") {
      delta <- misfit
    }

    k <- sqrt(2*delta/sum(E*E))
    E <- k*E

  } else if(estimator == "ml") {

    if(fit == "rmsr") {
      delta <- "A given RMSR is compatible with multiple maximum likelihood discrepancy values and is not provided"
    } else if(fit == "cfi") {
      null_f <- -log(det(R))
      delta <- (1-misfit)*null_f
    } else if(fit == "rmsea") {
      delta <- misfit^2 * df
    } else if(fit == "raw") {
      delta <- misfit
    }

    if(fit == "rmsr") {

      k <- sqrt((0.5*p*(p-1))*2*misfit^2/sum(E*E))
      E <- k*E

    } else {

      constant <- 1e-04 / sqrt(mean(E*E))
      E <- constant*E # Fix this to avoid NAs
      R_inv <- solve(R)
      G <- R_inv %*% E
      x <- suppressWarnings(grid_search(delta, G))
      # x <- sqrt(2*delta/sum(G*G)) # Initial value suggested by Cudeck
      k <- opt(x, delta, G)
      # limits <- c(-1e05, 1e05)
      # k <- GSS(delta, G, limits)
      # k <- grad_descend(delta, G)
      E <- k*E

    }
  }

  R_error <- Phat + E

  if(estimator == "ml" & fit == "rmsr") {
    delta <- log(det(R)) - log(det(R_error)) + sum(R_error*solve(R)) - nrow(R)
  }

  # Check for positiveness:
  minimum_eigval <- min(eigen(R_error, symmetric = TRUE, only.values = TRUE)$values)
  positive <- TRUE
  if(minimum_eigval <= 0) {
    # warning("The matrix was not positive-definite. The amount of misfit may be too big.")
    # positive <- FALSE
    return(yuan(R, lambda, Phi, Psi, fit, misfit, estimator))
  }

  return(list(R_error = R_error, fit = fit, delta = delta, misfit = misfit, positive = positive))

}

#' @title
#' Simulate factor structures with misspecification errors.
#' @description
#'
#' Simulate factor and bifactor structures with crossloadings, correlated factors, and more.
#'
#' @usage
#'
#' simfactor(nfactors = 5, nitems = 6, loadings = "medium",
#' crossloadings = 0, correlations = 0,
#' estimator = "minres", fit = "rmsr", misfit = 0,
#' error_method = "cudeck", efa = FALSE,
#' ngenerals = 0, loadings_g = "medium", correlations_g = 0,
#' pure = FALSE,
#' lambda = NULL, Phi = NULL, Psi = NULL)
#'
#' @param nfactors Number of factors.
#' @param nitems Number of items per factor.
#' @param loadings Loadings' magnitude on the factors: "low", "medium" or "high". Defaults to "medium".
#' @param crossloadings Magnitude of the cross-loadings among the group factors. Defaults to 0.
#' @param correlations Correlation among the factors. Defaults to 0.
#' @param estimator estimator used to generate population error: "minres" or "ml".
#' @param fit Fit index to control the population error.
#' @param misfit Misfit value to generate population error.
#' @param error_method Method used to control population error: c("yuan", "cudeck"). Defaults to "cudeck".
#' @param efa Reproduce the error with EFA or CFA. Defaults to FALSE (CFA).
#' @param ngenerals Number of general factors.
#' @param loadings_g Loadings' magnitude on the general factors: "low", "medium" or "high". Defaults to "medium".
#' @param correlations_g Correlation among the general factors. Defaults to 0.
#' @param pure Fix a pure item on each general factor. Defaults to FALSE.
#' @param lambda Custom loading matrix. If Phi is NULL, then all the factors will be correlated at the value given in correlations.
#' @param Phi Custom Phi matrix. If lambda is NULL, then Phi should be conformable to the loading matrix specified with the above arguments.
#' @param Psi Custom Psi matrix.
#'
#' @details \code{simfactor} generates bi-factor and generalized bifactor patterns with cross-loadings, pure items and
#' correlations among the general and group factors. When \code{crossloading} is different than 0, one cross-loading
#' is introduced for an item pertaining to each group factor. When \code{pure} is TRUE, one item loading of each group
#' factor is removed so that the item loads entirely on the general factor. To maintain the item communalities
#' constant upon these modifications, the item loading on the other factors may shrunk (if adding cross-loadings)
#' or increase (if setting pure items).
#'
#' Loading magnitudes may range between 0.3-0.5 ("low"), 0.4-0.6 ("medium") and 0.5-0.7 ("high"). Custom ranges can be supplied as vectors (i.e., c(0.2, 0.5))
#'
#' @return List with the following objects:
#' \item{lambda}{Population loading matrix.}
#' \item{Phi}{Population factor correlation matrix.}
#' \item{Psi}{Population covariance matrix between the errors.}
#' \item{R}{Model correlation matrix.}
#' \item{R_error}{Model correlation matrix with misspecification errors.}
#' \item{uniquenesses}{Population uniquenesses.}
#' \item{delta}{Minimum of the loss function that correspond to the misfit value.}
#'
#' @references
#' Cudeck, R., & Browne, M. W. (1992). Constructing a covariance matrix that yields a specified minimizer and a specified minimum discrepancy function value. \emph{Psychometrika, 57}(3), 357–369. \doi{10.1007/BF02295424}
#'
#' Jiménez, M., Abad, F. J., Garcia-Garzon, E., & Garrido, L. E. (2023). Exploratory Bi-factor Analysis with Multiple General Factors. \emph{Multivariate Behavioral Research, 58}(6), 1072–1089. \doi{10.1080/00273171.2023.2189571}
#'
#' Jiménez, M., Abad, F. J., Garcia-Garzon, E., Golino, H., Christensen, A. P., & Garrido, L. E. (2023). Dimensionality assessment in bifactor structures with multiple general factors: A network psychometrics approach. \emph{Psychological Methods}. Advance online publication. \doi{10.1037/met0000590}
#'
#' Yuan, K.-H., & Hayashi, K. (2003). Bootstrap approach to inference and power analysis based on three test statistics for covariance structure models. \emph{British Journal of Mathematical and Statistical Psychology, 56}(1), 93–110. \doi{10.1348/000711003321645368}
#'
#' @examplesIf requireNamespace("MASS", quietly = TRUE)
#' # Simulate data:
#' sim <- simfactor(nfactors = 3, nitems = 4, correlations = 0.40,
#'                  crossloadings = 0.30)
#' sim$lambda
#' sim$Phi
#' scores <- MASS::mvrnorm(1e3, rep(0, nrow(sim$R_error)), Sigma = sim$R_error)
#' s <- cor(scores)
#'
#' @export
simfactor <- function(nfactors = 3, nitems = 4,
                      loadings = "medium", crossloadings = 0,
                      correlations_g = 0.3, correlations = 0,
                      estimator = "uls", fit = "rmsr", misfit = "close",
                      error_method = "cudeck", efa = FALSE,
                      ngenerals = 0, loadings_g = "medium",
                      pure = FALSE,
                      lambda = NULL, Phi = NULL, Psi = NULL) {

  ng <- ngenerals # Save the number of general factors for recursive iteration
  condition <- ngenerals == 0
  # Create a general factor by default, and then remove it upon the condition:
  if(condition) ngenerals <- 1

  lambda_null <- is.null(lambda)
  Phi_null <- is.null(Phi)
  Psi_null <- is.null(Psi)

  if(!lambda_null & Phi_null) {

    # If lambda is provided but not Phi, then create Phi correlating all
    # the factors at the value given by correlations

    nfactors <- ncol(lambda)
    Phi <- matrix(correlations, nfactors, nfactors)
    diag(Phi) <- 1

  }

  if(lambda_null) {

    # Configure loadings:

    if(crossloadings > 0.4) {

      stop("Crossloadings are too large")

    }

    if(is.numeric(loadings_g)) {
      loadings_g. = loadings_g # custom loadings on the general factors
    } else {
      if(loadings_g == "low") {
        loadings_g. = c(.3, .5)
      } else if(loadings_g == "medium") {
        loadings_g. = c(.4, .6)
      } else if(loadings_g == "high") {
        loadings_g. = c(.5, .7)
      }
    }

    if(is.numeric(loadings)) {
      loadings. = loadings # custom loadings on the group factors
    } else {
      if(loadings == "low") {
        loadings. = c(.3, .5)
      } else if(loadings == "medium") {
        loadings. = c(.4, .6)
      } else if(loadings == "high") {
        loadings. = c(.5, .7)
      }
    }

    # Total number of group factors:
    n_groups <- ngenerals * nfactors

    # Total number of items:
    n_items <- n_groups * nitems

    # Number of items per general factor:
    items_per_general <- n_items / ngenerals

    # Total number of factors:
    n_factors <- ngenerals + n_groups

    # Initialize the population loading matrix:
    lambda <- matrix(NA, nrow = n_items, ncol = n_factors)

    # Item loadings on the group factors:

    sequen <- seq(loadings.[2], loadings.[1], length.out = nitems)

    for(i in 0:(n_groups-1)) {

      start_row <- 1 + i*nitems
      end_row <- start_row + nitems - 1
      # Loadings are created by equal increments:
      # lambda[start_row:end_row , 1+i+ngenerals] <- sequen
      # Loadings are simulated from a uniform distribution:
      lambda[start_row:end_row , 1+i+ngenerals] <- stats::runif(nitems, loadings.[1], loadings.[2])
      # lambda[start_row:end_row , 1+i+ngenerals] <- mean(loadings.)

    }

    # Simulate item loadings on the general factors from a uniform distribution:

    for(i in 0:(ngenerals-1)) {

      start_row <- 1 + i*items_per_general
      end_row <- start_row + items_per_general - 1
      lambda[start_row:end_row , i+1] <- stats::runif(items_per_general, loadings_g.[1], loadings_g.[2])
      # lambda[start_row:end_row , i+1] <- mean(loadings_g.)

    }

    colnames(lambda) <- c(paste("G", 1:ngenerals, sep = ""), paste("S", 1:n_groups, sep = ""))

    # Pure items on the general factors:

    if(pure) {

      # row_indexes <- seq(from = 1, to = n_items, by = nitems)
      value <- sequen[floor(nitems/2 + 1)]
      row_indexes <- unlist(apply(lambda, 2, FUN = function(x) which(x == value)))
      column_indexes <- apply(lambda[row_indexes, ], 1, FUN = function(x) which(x > 0))
      n <- n_groups * ngenerals
      indexes <- which(!is.na(lambda[row_indexes, 1:ngenerals]))
      m <- sqrt(lambda[row_indexes, 1:ngenerals][indexes]^2 +
                  lambda[row_indexes, ][which(lambda[row_indexes, ] == value)]^2)
      lambda[row_indexes, 1:ngenerals][indexes] <- m
      lambda[row_indexes, ][which(lambda[row_indexes, ] == value)] <- 0.01

    }

    # Cross-loadings on the group factors:

    if(crossloadings != 0 & nfactors > 1) {

      ratio <- nfactors
      row_index <- seq(nitems+1, n_items, by = nitems)
      col_index <- seq(ngenerals+1, n_factors-1, by = 1)

      if(ratio < length((row_index))) {
        delete <- seq(ratio, length(row_index), by = ratio)
        row_index <- row_index[-delete]
        col_index <- col_index[-delete]
      }

      # Insert cross-loadings and then recalibrate the loadings on the general and group factors to maintain the previous communality:

      for(i in 1:length(row_index)) {

        row_indexes <- row_index[i]:(row_index[i])
        col_index_2 <- which(lambda[row_indexes[1], ] > 0)
        lambda[row_indexes, col_index[i]] <- crossloadings
        # lambda[row_indexes, col_index_2] <- sqrt(lambda[row_indexes, col_index_2]^2 - crossloadings^2/2)
        # Skip the recalibration of the items

      }

      for(i in 1:ngenerals) {

        row_index <- items_per_general*(i-1)+1
        row_indexes <- row_index:(row_index)
        col_index_2 <- which(lambda[row_indexes[1], ] > 0)
        lambda[row_indexes, ngenerals+i*ratio] <- crossloadings
        # lambda[row_indexes, col_index_2] <- sqrt(lambda[row_indexes, col_index_2]^2 - crossloadings^2/2)
        # Skip the recalibration of the items

      }

    }

    lambda[is.na(lambda)] <- 0
    rownames(lambda) <- paste("item", 1:nrow(lambda), sep = "")

  }

  # if ngenerals == 0, remove the general factor:
  if(condition & lambda_null) lambda <- lambda[, -1, drop = FALSE]

  if(Phi_null) {

    # Factor correlations:

    Phi <- matrix(0, n_factors, n_factors)
    Phi[1:ngenerals, 1:ngenerals] <- correlations_g
    Phi[-(1:ngenerals), -(1:ngenerals)] <- correlations
    diag(Phi) <- 1

    # if ngenerals == 0, remove the general factor:
    if(condition) Phi <- Phi[-1, , drop = FALSE][, -1, drop = FALSE]

  } else {

    if(nrow(Phi) != ncol(lambda) | ncol(Phi) != ncol(lambda)) {
      stop("The Phi matrix that was provided is not conformable with the lambda matrix generated by the argument specifications")
    }

  }

  if(!Psi_null) {
    Psi <- 0.5*(t(Psi) + Psi) # Make Psi symmetric
    diag(Psi) <- 0 # The actual uniquenesses depend on Lambda and Phi
  } else {
    p <- nrow(lambda)
    Psi <- matrix(0, p, p)
  }

  # Population model correlation matrix:

  R <- lambda %*% Phi %*% t(lambda) + Psi
  uniquenesses <- 1 - diag(R)
  diag(Psi) <- uniquenesses
  diag(R) <- 1
  R_error <- R
  delta <- 0

  # Execute sim_factor recursively until no communality is greater than 1:

  if( any(uniquenesses < 0) ) {

    if(!lambda_null | !Phi_null | !Psi_null) {
      stop("The provided lambda, Phi or Psi matrices produced negative variances")
    }

    warning("At least one communality greater than 1 was found, probably due to a high misfit value \n Resampling...")

    sim <- simfactor(ngenerals = ng, nfactors = nfactors,
                      nitems = nitems,
                      loadings_g = loadings_g, loadings = loadings,
                      crossloadings = crossloadings, pure = pure,
                      correlations_g = correlations_g, correlations = correlations,
                      estimator = estimator,
                      fit = fit, misfit = misfit, error_method = error_method,
                      lambda = lambda, Phi = Phi, Psi = Psi)

    return(sim)

  }

  # Add population error to the population model correlation matrix:

  if(misfit != 0 & misfit != "zero") {

    if(error_method == "cudeck") {

      cudeck_ <- cudeck(R = R, lambda = lambda, Phi = Phi, Psi = Psi,
                        fit = fit, misfit = misfit, efa = efa,
                        estimator = estimator)
      R_error <- cudeck_$R_error
      delta <- cudeck_$delta
      misfit <- cudeck_$misfit

    } else if(error_method == "yuan") {

      yuan_ <- yuan(R = R, lambda = lambda, Phi = Phi, Psi = Psi,
                    fit = fit, misfit = misfit,
                    estimator = estimator)
      R_error <- yuan_$R_error
      delta <- yuan_$delta
      misfit <- yuan_$misfit

    }

  } else {

    R_error <- R

  }

  return( list(lambda = lambda, Phi = Phi, Psi = Psi,
               uniquenesses = uniquenesses,
               R = R, R_error = R_error, fit = fit, delta = delta,
               misfit = misfit) )

}
