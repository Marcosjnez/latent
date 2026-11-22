library(latent)

# A scalar quadratic with its optimum at zero.
controls <- function(start, opt = "lbfgs", step = "armijo") {
  list(control_manifold = list(list(manifold = "euclidean", indices = 0L)),
       control_transform = list(),
       control_estimator = list(list(estimator = "target", indices = list(0L),
                                     target = matrix(0), weight = matrix(1))),
       control_optimizer = list(parameters = list(start), transparameters = list(start),
                                transparam2param = 0L, cores = 1L, rstarts = 1L,
                                opt = opt, step = step, ss = 0.5, maxit = 1L, eps = 1e-9))
}

for(opt in c("grad", "lbfgs", "newton")) {
  for(step in if(opt == "newton") c("armijo", "wolfe", "trust") else c("armijo", "wolfe")) {
    fit <- do.call(optimizer, controls(0, opt, step))
    stopifnot(fit$convergence, fit$iterations == 0L, fit$f == 0)

    # The final allowed step can reach the optimum.
    fit <- do.call(optimizer, controls(0.5, opt, step))
    stopifnot(fit$convergence, fit$iterations == 1L, abs(fit$f) < 1e-14)
  }
}

# Armijo must restore the accepted point when its only trial increases the loss.
conf <- controls(0.05, "grad")
conf$control_estimator[[1L]]$weight[] <- 10
conf$control_optimizer$step_maxit <- 1L
fit <- do.call(optimizer, conf)
stopifnot(!fit$convergence, drop(fit$parameters) == 0.05, abs(fit$f - 0.125) < 1e-14)

# EM uses the observed gradient to distinguish a stalled M step from a solution.
conf$control_optimizer$opt <- "em"
conf$control_optimizer$mopt <- "grad"
conf$control_optimizer$mstep_maxit <- 1L
fit <- do.call(optimizer, conf)
stopifnot(!fit$convergence, drop(fit$parameters) == 0.05)
conf$control_estimator[[1L]]$weight[] <- 1
fit <- do.call(optimizer, conf)
stopifnot(fit$convergence, fit$iterations == 1L, fit$f < 1e-16)

# Fixed coordinates can differ between starts and must survive postestimation.
conf <- controls(0)
conf$control_optimizer$cores <- conf$control_optimizer$rstarts <- 2L
conf$control_optimizer$parameters <- list(0, 0)
conf$control_optimizer$transparameters <- list(c(0, 9), c(0, 1))
conf$control_estimator[[1L]]$indices <- list(1L)
fit <- do.call(optimizer, conf)
object <- methods::new("latent", modelInfo = conf, Optim = fit)
stopifnot(fit$start == 2L, fit$f == 0.5, get_f(object)$f == fit$f)

# A singular starting matrix may fail even though the other start is valid.
conf$control_transform <- list(list(transform = "inverse_of_I_minusX", p = 1L,
                                    indices_in = list(0L), indices_out = list(1L)))
conf$control_optimizer$parameters <- list(1, 0)
conf$control_optimizer$transparameters <- list(c(1, 0), c(0, 0))
conf$control_estimator[[1L]]$target[] <- 1
fit <- do.call(optimizer, conf)
stopifnot(fit$start == 2L, fit$convergence, fit$f == 0)

# The loss and derivatives must agree for both kinds of coefficient penalty.
conf <- controls(2)
for(penalty in list(list(estimator = "ridge", lambda = 2, power = 2, N = 4),
                    list(estimator = "gaussian_loglik", means = matrix(1),
                         sds = matrix(.5), alpha = 1, N = 4))) {
  penalty$indices <- list(0L)
  conf$control_estimator <- list(penalty)
  object <- methods::new("latent", modelInfo = conf,
                         Optim = list(parameters = 2, transparameters = 2))
  set.seed(12)
  seed <- .Random.seed
  obj <- get_hess(object)
  stopifnot(identical(.Random.seed, seed), obj$f == get_f(object)$f)
  evaluate <- function(theta) {
    object@Optim$parameters <- theta
    get_f(object)$f
  }
  numeric_gradient <- (evaluate(2 + 1e-5) - evaluate(2 - 1e-5))/2e-5
  stopifnot(abs(drop(get_grad(object)$g) - numeric_gradient) < 1e-8)
}

# Retraction preserves an orthonormal point and applies a simplex step only once.
conf <- controls(c(0.6, 0.8))
conf$control_optimizer$transparam2param <- 0:1
conf$control_manifold <- list(list(manifold = "orth", indices = 0:1, p = 2L, q = 1L))
conf$control_estimator <- list(list(estimator = "target", indices = list(0:1),
                                    target = matrix(c(0.6, 0.8)), weight = matrix(1, 2L)))
fit <- do.call(optimizer, conf)
stopifnot(fit$convergence, fit$iterations == 0L, fit$f < 1e-25)

conf$control_optimizer$opt <- "grad"
conf$control_optimizer$parameters <- conf$control_optimizer$transparameters <- list(c(0.5, 0.5))
conf$control_manifold <- list(list(manifold = "simplex", indices = 0:1))
conf$control_estimator[[1L]]$target[] <- c(0.7, 0.3)
fit <- do.call(optimizer, conf)
stopifnot(max(abs(drop(fit$parameters) - c(0.6, 0.4))) < 1e-14)

# Check the matrix-product chain rule for both a full and a symmetric output.
X <- c(0.7, -0.1, 0.09, 0.6)
S <- c(0.6, 0.1, 0.1, 0.7)
for(name in c("crossprod", "XYXt", "XtYX", "sem_cov_model")) {
  inputs <- switch(name, crossprod = list(X), sem_cov_model = list(X, S, c(diag(2L))),
                   list(X, S))
  pars <- unlist(inputs)
  n <- length(pars)
  indices <- lapply(seq_along(inputs), function(i) seq_along(inputs[[i]]) - 1L + 4L*(i-1L))
  for(output in list(0:3, c(0L, 1L, 1L, 2L))) {
    output <- n + output
    selected <- unique(output)
    conf <- controls(pars)
    conf$control_optimizer$transparam2param <- seq_len(n) - 1L
    conf$control_optimizer$transparameters <- list(c(pars, numeric(length(selected))))
    conf$control_manifold[[1L]]$indices <- seq_len(n) - 1L
    conf$control_transform <- list(list(transform = name, p = 2L, q = 2L,
                                        indices_in = indices, indices_out = list(output)))
    conf$control_estimator <- list(list(estimator = "target", indices = list(selected),
                                        target = matrix(seq_along(selected)/5),
                                        weight = matrix(1, length(selected))))
    object <- methods::new("latent", modelInfo = conf,
                           Optim = list(parameters = pars,
                                        transparameters = conf$control_optimizer$transparameters[[1L]]))
    direction <- seq_along(pars)/n
    evaluate <- function(h) {
      object@Optim$parameters <- pars + h*direction
      get_grad(object)
    }
    plus <- evaluate(1e-5)
    minus <- evaluate(-1e-5)
    numerical <- (plus$f-minus$f)/2e-5
    result <- get_dgrad(object, direction)
    jacobian <- as.matrix(get_jacob(object)$jacob)[selected+1L, seq_len(n), drop = FALSE]
    stopifnot(abs(drop(crossprod(direction, result$g)) - numerical) < 1e-8,
              max(abs(result$dg - (plus$g-minus$g)/2e-5)) < 1e-8,
              max(abs(jacobian %*% direction -
                        (plus$transparameters-minus$transparameters)[selected+1L]/2e-5)) < 1e-8)
  }
}
