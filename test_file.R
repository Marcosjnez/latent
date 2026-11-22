# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 22/11/2026

#### Store a dataset ####

# usethis::use_data(dataset, overwrite = TRUE)

#### Build the package ####

# roxygen2::roxygenise()

#### Install latent ####

# pak::pak("marcosjnez/latent")
# devtools::install_github("marcosjnez/latent", force = TRUE)

# Install a previous version:
# remotes::install_github("marcosjnez/latent", ref = "6e3603e")

#### Check derivatives ####

library(latent)
model <- 'visual =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6'
fit <- lcfa(data = HolzingerSwineford1939, model = model, positive = FALSE)

# Check the free coordinates slightly away from the optimum.
point <- fit
parameters <- c(fit@Optim$parameters)
direction <- seq_along(parameters)/length(parameters)
point@Optim$parameters <- parameters + 0.01*direction
point@modelInfo$control_optimizer$dparameters <- direction
x <- grad_comp(point, compute = "all", eps = 1e-5)
stopifnot(max(abs(x$g - x$numg)) < 1e-6,
          max(abs(x$dg - x$numdg)) < 1e-6)

# Calculate the Hessian matrix using numerical approximations:
G <- function(parameters) {

  point@Optim$parameters <- parameters
  g <- c(get_grad(point)$g)

  return(g)

}

H <- numDeriv::jacobian(func = G, x = point@Optim$parameters)
stopifnot(max(abs(H - get_hess(point)$h)) < 1e-6)

# Propagate the free-parameter covariance at the fitted solution.
v <- information(fit)$VCOV
VCOV <- get_vcov(fit, vcov = v)
VCOV$vcov
Matrix::diag(VCOV$vcov)

J <- get_jacob(fit)$jacob[, match(fit@modelInfo$parameters_labels,
                                   fit@modelInfo$transparameters_labels), drop = FALSE]
stopifnot(max(abs(VCOV$vcov - J %*% v %*% Matrix::t(J))) < 1e-8)

#### To-do ####
# ML three-step SEs still need to propagate uncertainty in the estimated
# classification-error matrix for both modal and proportional assignment.
