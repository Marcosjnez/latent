# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 10/08/2025

.positive_ml_root <- function(values, discrepancy) {
  # Find the root before the perturbed covariance reaches its singular boundary.
  magnitude <- max(abs(values))
  values <- values/magnitude
  limit <- if(min(values) < 0) (1-.Machine$double.eps)/-min(values) else Inf
  upper <- min(1, limit)
  f <- function(k) discrepancy(k*values)
  while(f(upper) < 0) {
    if(upper >= limit) stop("The target exceeds the numerically positive-definite range.")
    upper <- min(2*upper, limit)
  }
  stats::uniroot(f, c(0, upper), tol = 1e-13)$root/magnitude
}

grid_search <- function(delta, G) {

  if(delta == 0) return(0)
  values <- eigen(G, only.values = TRUE)$values
  if(any(abs(Im(values)) > 1e-8 * max(Mod(values)))) stop("Expected a covariance error direction with real eigenvalues.")
  values <- Re(values)
  if(!any(values != 0) || !is.finite(delta) || delta < 0) stop("Invalid ML discrepancy or error direction.")

  .positive_ml_root(values, function(z) sum(z-log1p(z))-delta)
}

opt <- function(x, delta, G) {
  grid_search(delta, G)
}

opt_error <- function(x, delta, G) {
  tryCatch(opt(x, delta, G), error = function(e) x)
}

f_minres <- function(x, S, ldetS, Inv_W, q, indexes_lambda, indexes_phi, indexes_psi) {

  p <- nrow(S)
  lambda_p <- length(indexes_lambda)
  Lambda <- matrix(0, p, q)
  Lambda[indexes_lambda] <- x[seq_len(lambda_p)]
  phi_p <- length(indexes_phi)
  Phi <- matrix(0, q, q)
  Phi[indexes_phi] <- x[lambda_p + seq_len(phi_p)]
  Phi <- t(Phi) + Phi
  diag(Phi) <- 1
  Psi <- matrix(0, p, p)
  Psi[indexes_psi] <- x[lambda_p + phi_p + seq_along(indexes_psi)]
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
  Lambda[indexes_lambda] <- x[seq_len(lambda_p)]
  phi_p <- length(indexes_phi)
  Phi <- matrix(0, q, q)
  Phi[indexes_phi] <- x[lambda_p + seq_len(phi_p)]
  Phi <- t(Phi) + Phi
  diag(Phi) <- 1
  Psi <- matrix(0, p, p)
  Psi[indexes_psi] <- x[lambda_p + phi_p + seq_along(indexes_psi)]
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
  Lambda[indexes_lambda] <- x[seq_len(lambda_p)]
  phi_p <- length(indexes_phi)
  Phi <- matrix(0, q, q)
  Phi[indexes_phi] <- x[lambda_p + seq_len(phi_p)]
  Phi <- t(Phi) + Phi
  diag(Phi) <- 1
  Psi <- matrix(0, p, p)
  Psi[indexes_psi] <- x[lambda_p + phi_p + seq_along(indexes_psi)]
  Psi[upper.tri(Psi)] <- t(Psi)[upper.tri(Psi)]
  Rhat <- Lambda %*% Phi %*% t(Lambda) + Psi
  U <- tryCatch(chol(Rhat), error = function(e) NULL)
  if(is.null(U)) return(Inf)
  f <- 2*sum(log(diag(U))) - ldetS + sum(S*chol2inv(U)) - p

  return(f)

}
g_ml <- function(x, S, ldetS, Inv_W, q, indexes_lambda, indexes_phi, indexes_psi) {

  p <- nrow(S)
  lambda_p <- length(indexes_lambda)
  Lambda <- matrix(0, p, q)
  Lambda[indexes_lambda] <- x[seq_len(lambda_p)]
  phi_p <- length(indexes_phi)
  Phi <- matrix(0, q, q)
  Phi[indexes_phi] <- x[lambda_p + seq_len(phi_p)]
  Phi <- t(Phi) + Phi
  diag(Phi) <- 1
  Psi <- matrix(0, p, p)
  Psi[indexes_psi] <- x[lambda_p + phi_p + seq_along(indexes_psi)]
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
  Lambda[indexes_lambda] <- x[seq_len(lambda_p)]
  phi_p <- length(indexes_phi)
  Phi <- matrix(0, q, q)
  Phi[indexes_phi] <- x[lambda_p + seq_len(phi_p)]
  Phi <- t(Phi) + Phi
  diag(Phi) <- 1
  Psi <- matrix(0, p, p)
  Psi[indexes_psi] <- x[lambda_p + phi_p + seq_along(indexes_psi)]
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
  Lambda[indexes_lambda] <- x[seq_len(lambda_p)]
  phi_p <- length(indexes_phi)
  Phi <- matrix(0, q, q)
  Phi[indexes_phi] <- x[lambda_p + seq_len(phi_p)]
  Phi <- t(Phi) + Phi
  diag(Phi) <- 1
  Psi <- matrix(0, p, p)
  Psi[indexes_psi] <- x[lambda_p + phi_p + seq_along(indexes_psi)]
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
                estimator = "minres", W = NULL, start = NULL) {

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
  diagonal <- match(indexes_psi, diag_indexes)
  init_psi[!is.na(diagonal)] <- init_diag_psi[diagonal[!is.na(diagonal)]]

  lower_psi <- rep(0.005, psi_p) # Lower bounds for the uniquenessess
  lower_psi[cor_res_indexes] <- -0.995 # Lower bounds for correlated residuals
  upper_psi <- rep(0.995, psi_p) # Upper bounds for correlated residuals
  lower <- c(rep(-Inf, lambda_p), rep(-1, phi_p), lower_psi)
  upper <- c(rep(Inf, lambda_p), rep(1, phi_p), upper_psi)

  x <- if(is.null(start)) c(stats::runif(lambda_p), rep(0, phi_p), init_psi) else start

  if(estimator %in% c("minres", "uls")) {

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
  lambda_hat[indexes_lambda] <- cfa$par[seq_len(lambda_p)]

  # Arrange phi parameter estimates:
  phi_hat <- matrix(0, q, q)
  phi_hat[indexes_phi] <- cfa$par[lambda_p + seq_len(phi_p)]
  phi_hat <- t(phi_hat) + phi_hat
  diag(phi_hat) <- 1

  # Arrange psi parameter estimates:
  psi_hat <- matrix(0, p, p)
  psi_hat[indexes_psi] <- cfa$par[lambda_p + phi_p + seq_along(indexes_psi)]
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
# dxt <- function(p, q) {
#
#   # derivative wrt transpose (just a permutation matrix)
#
#   pq <- p*q
#
#   res <- array(0, dim = c(pq, pq))
#   null <- matrix(0, p, q)
#
#   for(i in 1:pq) {
#     temp <- null
#     temp[i] <- 1
#     res[, i] <- c(t(temp))
#   }
#
#   return(res)
#
# }
gLRhat <- function(lambda, phi) {
  p <- nrow(lambda)
  g <- (lambda %*% phi) %x% diag(p)
  transpose <- as.vector(t(matrix(seq_len(p*p), p)))
  g + g[transpose, , drop = FALSE]
}

gPRhat <- function(lambda, q) {
  pairs <- which(lower.tri(diag(q)), arr.ind = TRUE)
  g <- lambda %x% lambda
  g[, pairs[, 1] + q*(pairs[, 2]-1L), drop = FALSE] +
    g[, pairs[, 2] + q*(pairs[, 1]-1L), drop = FALSE]
}

guRhat <- function(p) {
  gURhat(p, which(diag(p) != 0, arr.ind = TRUE))
}

gURhat <- function(p, pairs = which(lower.tri(diag(p), diag = TRUE), arr.ind = TRUE)) {
  g <- matrix(0, p*p, nrow(pairs))
  columns <- seq_len(nrow(pairs))
  g[cbind(pairs[, 1] + p*(pairs[, 2]-1L), columns)] <- 1
  g[cbind(pairs[, 2] + p*(pairs[, 1]-1L), columns)] <- 1
  g
}

.scale_population_error <- function(R, E, fit, misfit, estimator, df, uniquenesses) {

  p <- nrow(R)
  off <- lower.tri(R)
  if(fit == "rmsr" && is.character(misfit)) {
    fraction <- switch(misfit, close = 0.05, acceptable = 0.10,
                       stop("Unknown RMSR target."))
    misfit <- fraction*mean(1-uniquenesses)
  }
  if(!is.numeric(misfit) || length(misfit) != 1L || !is.finite(misfit) || misfit < 0)
    stop("Supply a non-negative numeric misfit for the requested fit index.")
  if(fit == "cfi" && misfit > 1) stop("CFI must be between 0 and 1.")

  if(fit == "rmsr") {
    delta <- misfit^2 * sum(off)
  } else if(fit == "raw") {
    delta <- misfit
  } else if(fit == "rmsea") {
    if(df <= 0) stop("RMSEA requires positive degrees of freedom.")
    delta <- misfit^2 * df
  } else {
    # CFI uses the independence discrepancy of the generated correlation matrix.
    # At the desired CFI, the model discrepancy is (1-CFI) times the null discrepancy.
    if(estimator == "ml") {
      G <- solve(R, E)
      values <- Re(eigen(G, only.values = TRUE)$values)
      null <- -2*sum(log(diag(chol(R))))
      k <- .positive_ml_root(values, function(z) sum(z) - misfit*sum(log1p(z)) - (1-misfit)*null)
    } else {
      r <- R[off]
      e <- E[off]
      aa <- misfit*sum(e^2)
      bb <- (1-misfit)*sum(r*e)
      cc <- (1-misfit)*sum(r^2)
      if(aa == 0 && bb >= 0) return(NULL)
      k <- if(aa == 0) -cc/(2*bb) else (bb + sqrt(bb^2 + aa*cc))/aa
    }
  }

  if(fit != "cfi") {
    if(delta == 0) {
      k <- 0
    } else if(fit == "rmsr" || estimator != "ml") {
      k <- sqrt(delta/sum(E[off]^2))
    } else {
      k <- grid_search(delta, solve(R, E))
    }
  }
  R_error <- R + k*E
  U <- tryCatch(chol(R_error), error = function(e) NULL)
  if(is.null(U)) return(NULL)
  if(fit == "cfi" || (fit == "rmsr" && estimator == "ml")) {
    delta <- if(estimator == "ml") {
      2*sum(log(diag(chol(R))) - log(diag(U))) + k*sum(diag(solve(R, E)))
    } else sum((k*E[off])^2)
  }
  list(R_error = R_error, fit = fit, delta = delta, misfit = misfit, positive = TRUE)
}

.population_error_space <- function(R, lambda, Phi, Psi, estimator, efa) {

  p <- nrow(lambda)
  q <- ncol(lambda)
  off <- lower.tri(R)
  dL <- gLRhat(lambda, Phi)
  dU <- gURhat(p, which(Psi != 0 & lower.tri(Psi, diag = TRUE), arr.ind = TRUE))
  if(efa) {
    derivatives <- cbind(dL, dU)
  } else {
    derivatives <- dL[, as.vector(lambda != 0), drop = FALSE]
    free_phi <- Phi[lower.tri(Phi)] != 0
    if(any(free_phi)) derivatives <- cbind(derivatives, gPRhat(lambda, q)[, free_phi, drop = FALSE])
    derivatives <- cbind(derivatives, dU)
  }
  if(estimator == "ml") {
    inverse <- solve(R)
    B <- vapply(seq_len(ncol(derivatives)), function(j) {
      (inverse %*% matrix(derivatives[, j], p, p) %*% inverse)[off]
    }, numeric(sum(off)))
    dim(B) <- c(sum(off), ncol(derivatives))
  } else {
    B <- derivatives[as.vector(off), , drop = FALSE]
  }
  # Project only the off-diagonal errors, since the diagonal must stay fixed.
  basis <- qr(B)
  if(basis$rank >= sum(off)) stop("The model leaves no direction for population error.")
  df <- p*(p+1)/2 - (if(efa) qr(derivatives)$rank else ncol(derivatives))
  list(basis = basis, df = df)
}

cudeck <- function(R, lambda, Phi, Psi, fit = "rmsr", misfit = "close",
                   estimator = "minres", efa = FALSE) {
  p <- nrow(R)
  off <- lower.tri(R)
  space <- .population_error_space(R, lambda, Phi, Psi, estimator, efa)
  for(attempt in seq_len(100L)) {
    U <- matrix(stats::runif((p+1)*p, -1, 1), p+1, p)
    residual <- stats::cov2cor(crossprod(U)) * sqrt(outer(diag(Psi), diag(Psi)))
    E <- matrix(0, p, p)
    E[off] <- qr.resid(space$basis, residual[off])
    E <- E + t(E)
    if(!any(E != 0)) next
    result <- .scale_population_error(R, E, fit, misfit, estimator, space$df, diag(Psi))
    if(!is.null(result)) return(result)
  }
  stop("Could not generate a positive-definite matrix with the requested misfit.")
}

yuan <- function(R, lambda, Phi, Psi, fit = "rmsr", misfit = "close",
                 estimator = "minres", efa = FALSE) {
  p <- nrow(R)
  off <- lower.tri(R)
  target <- lambda != 0
  targetphi <- Phi != 0
  targetpsi <- Psi != 0
  start <- c(lambda[target], Phi[targetphi & lower.tri(Phi)], Psi[targetpsi & lower.tri(Psi, diag = TRUE)])
  space <- .population_error_space(R, lambda, Phi, Psi, estimator, efa)
  for(attempt in seq_len(100L)) {
    Rerror <- R
    Rerror[off] <- R[off] + stats::runif(sum(off), -1e-6, 1e-6)
    Rerror[upper.tri(R)] <- t(Rerror)[upper.tri(R)]
    cfa <- CFA(Rerror, target, targetphi, targetpsi, estimator = estimator, start = start)
    if(cfa$convergence != 0) next
    # Remove the residual score before enlarging the small fitted residual.
    E <- matrix(0, p, p)
    E[off] <- qr.resid(space$basis, (Rerror-cfa$model)[off])
    E <- E + t(E)
    if(!any(E != 0)) next
    result <- .scale_population_error(R, E, fit, misfit, estimator, space$df, diag(Psi))
    if(!is.null(result)) return(result)
  }
  stop("Could not generate a positive-definite matrix with the requested misfit.")
}

#' @title
#' Simulate factor structures with misspecification errors
#' @description
#'
#' Simulate factor and bifactor structures with crossloadings, correlated factors, and more.
#'
#' @usage
#'
#' simfactor(nfactors = 3, nitems = 4, loadings = "medium",
#' crossloadings = 0, correlations_g = 0.3, correlations = 0,
#' estimator = "uls", fit = "rmsr", misfit = "close",
#' error_method = "cudeck", efa = FALSE,
#' ngenerals = 0, loadings_g = "medium",
#' pure = FALSE,
#' lambda = NULL, Phi = NULL, Psi = NULL)
#'
#' @param nfactors Number of factors.
#' @param nitems Number of items per factor.
#' @param loadings Loadings' magnitude on the factors: "low", "medium" or "high". Defaults to "medium".
#' @param crossloadings Magnitude of the cross-loadings among the group factors. Defaults to 0.
#' @param correlations Correlation among the factors. Defaults to 0.
#' @param estimator Estimator used to generate population error. Defaults to "uls".
#' @param fit Fit index used to control the population error. One of \code{"rmsr"},
#'   \code{"rmsea"}, \code{"cfi"}, or \code{"raw"} for the estimator discrepancy.
#' @param misfit Target value of the selected fit index. For RMSR, \code{"close"} and
#'   \code{"acceptable"} use 5 and 10 percent of the mean communality. A numeric
#'   zero or \code{"zero"} requests the model without population error.
#' @param error_method Method used to control population error: c("yuan", "cudeck"). Defaults to "cudeck".
#' @param efa Reproduce the error with EFA or CFA. Defaults to FALSE (CFA).
#' @param ngenerals Number of general factors.
#' @param loadings_g Loadings' magnitude on the general factors: "low", "medium" or "high". Defaults to "medium".
#' @param correlations_g Correlation among the general factors. Defaults to 0.3.
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

  estimator <- match.arg(estimator, c("uls", "minres", "ml"))
  fit <- match.arg(fit, c("rmsr", "raw", "cfi", "rmsea"))
  error_method <- match.arg(error_method, c("cudeck", "yuan"))
  generated <- is.null(lambda)
  if(generated) {
    counts <- c(nfactors, nitems, ngenerals)
    if(length(counts) != 3L || any(!is.finite(counts) | counts != floor(counts)) ||
       min(nfactors, nitems) < 1 || ngenerals < 0) stop("Invalid factor or item count.")
    loading_range <- function(x) {
      if(is.character(x)) x <- switch(x, low = c(.3, .5), medium = c(.4, .6), high = c(.5, .7))
      if(length(x) == 1L && is.numeric(x)) x <- rep(x, 2L)
      if(!is.numeric(x) || length(x) != 2L || any(!is.finite(x)) || x[1] > x[2])
        stop("Supply low, medium, high or an ordered numeric loading range.")
      x
    }
    group_range <- loading_range(loadings)
    general_range <- loading_range(loadings_g)
    if(abs(crossloadings) > .4) stop("Crossloadings are too large")
    if(pure && ngenerals == 0) stop("Pure general-factor items require a general factor.")
    if(pure && nitems == 1L && crossloadings != 0) stop("Pure and cross-loaded items must be distinct.")
    n_groups <- max(1L, ngenerals)*nfactors
    p <- n_groups*nitems
    q <- ngenerals+n_groups
    group <- rep(seq_len(n_groups), each = nitems)
    general <- rep(seq_len(max(1L, ngenerals)), each = nfactors*nitems)
  } else {
    if(!is.matrix(lambda) || !is.numeric(lambda) || !length(lambda) || any(!is.finite(lambda)))
      stop("lambda must be a finite numeric matrix.")
    p <- nrow(lambda)
    q <- ncol(lambda)
  }

  if(is.null(Phi)) {
    Phi <- matrix(correlations, q, q)
    if(generated && ngenerals > 0) {
      generals <- seq_len(ngenerals)
      Phi[generals, generals] <- correlations_g
      Phi[generals, -generals] <- Phi[-generals, generals] <- 0
    }
    diag(Phi) <- 1
  }
  if(!is.matrix(Phi) || !all(dim(Phi) == q) || !isSymmetric(unname(Phi)) ||
     any(!is.finite(Phi)) || any(diag(Phi) != 1) ||
     min(eigen(Phi, symmetric = TRUE, only.values = TRUE)$values) < -1e-10)
    stop("Phi must be a conformable factor correlation matrix.")

  diagonal_psi <- is.null(Psi)
  if(diagonal_psi) Psi <- matrix(0, p, p)
  if(!is.matrix(Psi) || !all(dim(Psi) == p) || any(!is.finite(Psi)))
    stop("Psi must be a conformable residual covariance matrix.")
  Psi <- (Psi + t(Psi))/2
  for(attempt in seq_len(100L)) {
    if(generated) {
      lambda <- matrix(0, p, q, dimnames = list(paste0("item", seq_len(p)),
                       c(if(ngenerals) paste0("G", seq_len(ngenerals)), paste0("S", seq_len(n_groups)))))
      lambda[cbind(seq_len(p), ngenerals+group)] <- stats::runif(p, group_range[1], group_range[2])
      if(ngenerals) {
        lambda[cbind(seq_len(p), general)] <- stats::runif(p, general_range[1], general_range[2])
      }
      if(pure) {
        rows <- (seq_len(n_groups)-1L)*nitems + floor(nitems/2) + 1L
        s <- cbind(rows, ngenerals+group[rows])
        g <- cbind(rows, general[rows])
        covariance <- Phi[cbind(g[, 2], s[, 2])]
        h2 <- lambda[g]^2 + lambda[s]^2 + 2*lambda[g]*lambda[s]*covariance
        lambda[g] <- ifelse(lambda[g] < 0, -1, 1)*sqrt(pmax(h2, 0))
        lambda[s] <- 0
      }
      if(crossloadings != 0 && nfactors > 1L) {
        rows <- seq.int(1L, p, by = nitems)
        # First item of each group also loads on the previous group of the same general factor.
        previous <- group[rows] - 1L
        first <- previous %% nfactors == 0L
        previous[first] <- previous[first] + nfactors
        columns <- ngenerals + previous
        main <- lambda[rows, , drop = FALSE]
        covariance <- main %*% Phi
        h2 <- rowSums(main * covariance)
        # Scale the original loadings to retain each communality after adding the cross-loading.
        linear <- crossloadings * covariance[cbind(seq_along(rows), columns)]
        discriminant <- linear^2 + h2*(h2-crossloadings^2)
        scale <- (sqrt(pmax(discriminant, 0))-linear)/h2
        if(any(discriminant < -1e-14 | !is.finite(scale) | scale < 0)) {
          if(attempt == 100L) stop("The cross-loadings exceed the available item communalities.")
          next
        }
        lambda[rows, ] <- main * scale
        lambda[cbind(rows, columns)] <- crossloadings
      }
    }
    common <- tcrossprod(lambda %*% Phi, lambda)
    uniquenesses <- 1 - diag(common)
    diag(Psi) <- uniquenesses
    if(all(uniquenesses >= 0) && (diagonal_psi ||
       min(eigen(Psi, symmetric = TRUE, only.values = TRUE)$values) >= -1e-10)) break
    if(!generated || attempt == 100L) stop("The loadings and covariances produce invalid residual variances or covariances.")
  }
  R <- common + Psi
  diag(R) <- 1
  delta <- 0

  # Add population error to the population model correlation matrix:

  exact <- misfit %in% c(0, "zero") || (fit == "cfi" && identical(misfit, 1))
  if(!exact) {

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
                    estimator = estimator, efa = efa)
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
