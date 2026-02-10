# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 13/12/2025

#### Install latent ####

# devtools::install_github("marcosjnez/latent", force = TRUE)

#### Build the package ####

# roxygen2::roxygenise()

#### LCA (multinomial) ####

library(latent)
set.seed(2026)
fit <- lca(data = gss82, nclasses = 3,
           # multinomial = c("X1", "X2"),
           # poisson = ,
           # beta = ,
           # mimic = "LT",
           # model = formula(X1 ~ 1 + cluster1,
           #                 X2 ~ 1 + cluster2),
           item = rep("multinomial", ncol(gss82)),
           penalties = TRUE,
           control = list(opt = "newton",
                          step_maxit = 100,
                          tcg_maxit = 100),
           do.fit = TRUE)
fit@loglik # -2754.643
fit@penalized_loglik # -2759.507
fit@Optim$iterations # 67
fit@Optim$ng
max(fit@Optim$rg)
fit@Optim$convergence
fit@timing

# Plot model fit info:
fit

# Get fit indices:
getfit(fit)

# Get a summary:
summary(fit)

# Inspect model objects:
latInspect(fit, what = "pattern", digits = 3)
latInspect(fit, what = "coefs", digits = 3)
latInspect(fit, what = "classes", digits = 3)
latInspect(fit, what = "profile", digits = 3)
latInspect(fit, what = "posterior", digits = 3)
latInspect(fit, what = "table", digits = 3)

# Get standard errors:
SE <- se(fit, type = "standard", model = "model", digits = 4)
SE$table

# Get confidence intervals:
CI <- ci(fit, type = "standard", model = "model",
         confidence = 0.95, digits = 2)
CI$table

#### LCA (gaussian) ####

library(latent)
set.seed(2026)
fit <- lca(data = empathy[, 1:6], nclasses = 4L,
           item = rep("gaussian", ncol(empathy[, 1:6])),
           penalties = TRUE, do.fit = TRUE,
           control = list(opt = "lbfgs",
                          step_maxit = 100,
                          tcg_maxit = 100))

fit@loglik # -1841.336
fit@penalized_loglik # -1844.333
fit@Optim$iterations # 52
fit@Optim$ng
max(fit@Optim$rg)
fit@Optim$convergence # TRUE
fit@timing

# Plot model fit info:
fit

# Get fit indices:
getfit(fit)

# Inspect model objects:
latInspect(fit, what = "coefs", digits = 3)
latInspect(fit, what = "classes", digits = 3)
latInspect(fit, what = "profile", digits = 3)
latInspect(fit, what = "posterior", digits = 3)

# Get standard errors:
SE <- se(fit, type = "standard", model = "user", digits = 4)
SE$table

# Get confidence intervals:
CI <- ci(fit, type = "standard", model = "user",
         confidence = 0.95, digits = 2)
CI$table

#### Mixed LCA (multinomial and gaussian) ####

library(latent)
set.seed(2026)
fit <- lca(data = cancer[, 1:6], nclasses = 3L,
           item = c("gaussian", "gaussian",
                    "multinomial", "multinomial",
                    "gaussian", "gaussian"),
           control = list(opt = "newton",
                          eps = 1e-06,
                          step_maxit = 100,
                          tcg_maxit = 100),
           do.fit = TRUE)
fit@loglik # -5784.701
fit@penalized_loglik # -5795.573
fit@Optim$iterations
fit@Optim$ng
max(fit@Optim$rg)
fit@Optim$convergence
fit@timing

# Plot model fit info:
fit

# Get fit indices:
getfit(fit)

# Inspect model objects:
latInspect(fit, what = "coefs", digits = 3)
latInspect(fit, what = "classes", digits = 3)
latInspect(fit, what = "profile", digits = 3)
latInspect(fit, what = "posterior", digits = 3)

# Get standard errors:
SE <- se(fit, type = "standard", model = "user", digits = 4)
SE$table

# Get confidence intervals:
CI <- ci(fit, type = "standard", model = "user",
         confidence = 0.95, digits = 2)
CI$table

# hypothesis(fit, "b1|2 - b1|3 = 0")

#### LCA with covariates (gaussian) ####

library(latent)
data <- empathy[, 1:6]
X <- as.matrix(empathy[, 7:8]) # Covariates

set.seed(2026)
fit0 <- lca(data = data, X = NULL, model = NULL,
            item = rep("gaussian", ncol(data)),
            nclasses = 4L, penalties = TRUE,
            control = list(opt = "newton",
                           step_maxit = 100,
                           tcg_maxit = 100),
            do.fit = TRUE)
fit0@timing
fit0@loglik # -1841.336
fit0@penalized_loglik # -1844.333
fit0@Optim$iterations
fit0@Optim$ng
max(fit0@Optim$rg)
sqrt(sum(fit0@Optim$rg*-fit0@Optim$dir))
sqrt(sum(fit0@Optim$rg*fit0@Optim$rg))

penalties <- list(
  beta  = list(alpha = 0, lambda = 1, power = 1),
  class = list(alpha = 1),
  prob  = list(alpha = 1),
  sd    = list(alpha = 1)
)
Y <- as.matrix(empathy[, 9:10]) # Covariates
set.seed(2026)
fit <- lca(data = data, X = cbind(X, Y), model = fit0,
           item = rep("gaussian", ncol(data)),
           nclasses = 4L, penalties = F,
           control = list(opt = "newton",
                          step_maxit = 100,
                          tcg_maxit = 100),
           do.fit = TRUE)
fit@timing
fit@loglik # -1747.135
fit@penalized_loglik # -1750.669
fit@Optim$iterations
fit@Optim$convergence
fit@Optim$ng
max(fit@Optim$rg)
sqrt(sum(fit@Optim$rg*-fit@Optim$dir))
sqrt(sum(fit@Optim$rg*fit@Optim$rg))
fit@transformed_pars$beta

all.equal(fit0@parameters$items, fit@parameters$items)

beta <- fit@parameters$beta
vcov <- SE$vcov[1:9, 1:9]
matrix(sqrt(diag(vcov)), 3, 3)

new_se <- effects_coding(beta, vcov)
new_se$beta_new
matrix(new_se$se_new, 3, 4)
new_se <- move_intercept(beta, vcov)
new_se$beta_new
matrix(new_se$se_new, 3, 3)

# Plot model fit info:
fit

# Get fit indices:
getfit(fit)

# Inspect model objects:
latInspect(fit, what = "coefs", digits = 5)
latInspect(fit, what = "classes", digits = 5)
latInspect(fit, what = "profile", digits = 3)
latInspect(fit, what = "posterior", digits = 3)

predict(fit, new = rbind(c(2, 2, 2.428571, 2.142857),
                         c(1, 2, 3, 4)))
fitted(fit)

# Get standard errors:
SE <- se(fit, type = "standard", model = "model", digits = 4)
SE$table

# Get confidence intervals:
CI <- ci(fit, type = "standard", model = "model",
         confidence = 0.95, digits = 2)
CI$table

x <- plot(fit,
     type = "standard",
     what = "OR",
     effects = "coding",
     confidence = 0.95,
     show_est_ci = TRUE,
     est_ci_header_cex = 0.5,
     cex_y = 0.5,
     mfrow = c(2, 2),
     xlim = c(0, 5))

#### CFA ####

library(latent)
library(lavaan)

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9'

set.seed(2026)
fit <- lcfa(HolzingerSwineford1939, model = model,
            estimator = "ml",
            ordered = FALSE, std.lv = TRUE,
            mimic = "latent", do.fit = TRUE,
            control = list(opt = "newton",
                           step_maxit = 100,
                           tcg_maxit = 100))
fit@loss   # 0.283407
fit@loglik # -3427.131
fit@penalized_loglik # -3427.131
fit@Optim$iterations
fit@Optim$convergence
fit@timing
fit
SE <- se(fit, type = "standard", digits = 3)
SE$table_se

# With lavaan:
fit2 <- lavaan::cfa(model, data = HolzingerSwineford1939,
                    estimator = "ml",
                    # likelihood = "wishart",
                    std.lv = TRUE, std.ov = TRUE)
# Same loss value: OK
fit2@Fit@fx*2
fit@loss

fit@loglik
# WHAT'S THE DIFFERENCE?
fit2@loglik$loglik
fit2@h1$logl$loglik # SATURATED

lavaan::inspect(fit2, what = "se")

SE <- se(fit, digits = 3)
SE$table_se

fit@loss
fit@loglik # -0.283407
fit@penalized_loglik # -0.283407
fit@loss # 0.1574787
fit@Optim$iterations
fit@Optim$convergence
fit@timing

# Plot model fit info:
fit

# Get fit indices:
getfit(fit)

# Get a summary:
summary(fit)

# Inspect model objects:
latInspect(fit, what = "loadings", digits = 3)
latInspect(fit, what = "psi", digits = 3)
latInspect(fit, what = "theta", digits = 3)
latInspect(fit, what = "model", digits = 3)
latInspect(fit, what = "residuals", digits = 3)
latInspect(fit, what = "W", digits = 3)

library(lavaan)
fit2 <- cfa(model, data = HolzingerSwineford1939, estimator = "uls",
            std.lv = TRUE, std.ov = TRUE)
summary(fit2, fit.measures = FALSE)
fitMeasures(fit2, digits = 5)
inspect(fit2, what = "se")
inspect(fit2, what = "est")

SE <- se(fit, digits = 3)
SE$table_se

psi <- latInspect(fit, what = "psi", digits = 3)[[1]]
psi[lower.tri(psi, diag = TRUE)] %*%
         t(duplication(3, halflower = FALSE))

duplication(3, halflower = FALSE) %*% psi[lower.tri(psi, diag = TRUE)]

#### Multigroup CFA ####

library(latent)

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9'

fit <- lcfa(HolzingerSwineford1939, model = model,
            group = "school", estimator = "ml",
            ordered = FALSE, std.lv = TRUE,
            mimic = "latent", do.fit = TRUE)

fit@loglik # -3422.946 (ml)
fit@loss # 0.4055109 (uls) / 0.7677016 (ml)
fit@Optim$iterations
fit@Optim$convergence
fit@timing
latInspect(fit, what = "loglik")

# With lavaan:
fit2 <- lavaan::cfa(model, data = HolzingerSwineford1939,
                    group = NULL, estimator = "ml",
                    std.lv = TRUE, std.ov = TRUE)

fit2@loglik$loglik
fit2@Fit@fx*2
fit@loss

# Plot model fit info:
fit

# Get fit indices:
getfit(fit)

# Get a summary:
summary(fit)

# Inspect model objects:
latInspect(fit, what = "loadings", digits = 3)
latInspect(fit, what = "psi", digits = 3)
latInspect(fit, what = "uniquenesses", digits = 3)
latInspect(fit, what = "theta", digits = 3)
latInspect(fit, what = "model", digits = 3)
latInspect(fit, what = "residuals", digits = 3)

SE <- se(fit, type = "lavaan", model = "user", digits = 3)
SE$table_se
lavaan::inspect(fit2, what = "se")

library(lavaan)
fit2 <- cfa(model, data = HolzingerSwineford1939,
            estimator = "ml", group = "school",
            std.lv = TRUE, std.ov = TRUE)
fit2
summary(fit2, fit.measures = FALSE)
fitMeasures(fit2)
inspect(fit2, what = "est")
inspect(fit2, what = "std")

#### CFA (nonpositive definite) ####

library(latent)
library(lavaan)

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9
          x1 ~~ x5
          x1 ~~ x4
          x4 ~~ x5
          x4 ~~ x6'

# With lavaan:
fit2 <- cfa(model, data = HolzingerSwineford1939,
            estimator = "ml", std.lv = TRUE, std.ov = TRUE)
inspect(fit2, what = "est") # NEGATIVE VARIANCE
det(inspect(fit2, what = "est")$theta)
fit2@Fit@fx*2
fit2@loglik$loglik

# With latent:
set.seed(2026)
fit <- lcfa(data = HolzingerSwineford1939, model = model,
            estimator = "ml", positive = TRUE,
            penalties = list(logdet = list(w = 0.01)),
            ordered = FALSE, std.lv = TRUE,
            mimic = "latent", do.fit = TRUE,
            control = list(opt = "newton", maxit = 1000L,
                           cores = 20L, rstarts = 20L, eps = 1e-05,
                           tcg_maxit = 10))

fit@loglik # -3421.497 (ML)
fit@penalized_loglik # -3421.502 (ML)
fit@loss # 0.1419955 (ULS) / 0.2467385 (ML)
fit@Optim$iterations
fit@Optim$convergence
fit@Optim$ng
fit@timing

fit2@Fit@fx*2
fit@loss

# Plot model fit info:
fit

# Get fit indices:
getfit(fit)

# Inspect model objects:
latInspect(fit, what = "loadings", digits = 3)
latInspect(fit, what = "psi", digits = 3)
latInspect(fit, what = "uniquenesses", digits = 3)
det(latInspect(fit, what = "theta", digits = 3)[[1]])
round(latInspect(fit, what = "theta", digits = 3)[[1]], 3)
latInspect(fit, what = "model", digits = 3)

theta <- latInspect(fit, what = "theta", digits = 3)[[1]]
det(theta)
eigen(theta)$values
sum(eigen(theta)$values)
# solve(theta)

psi <- latInspect(fit, what = "psi", digits = 3)[[1]]
det(psi)
eigen(psi)$values
sum(eigen(psi)$values)
# solve(psi)

#### Multigroup CFA (nonpositive definite) ####

library(latent)
library(lavaan)

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9
          x1 ~~ x5
          x1 ~~ x4
          x4 ~~ x5
          x4 ~~ x6'

# With lavaan:
fit2 <- cfa(model, data = HolzingerSwineford1939,
            group = "school", estimator = "ml",
            std.lv = TRUE, std.ov = TRUE)
fit2
inspect(fit2, what = "est")
fitmeasures(fit2, fit.measures = c("cfi", "tli", "rmsea", "srmr"))

# With latent:
fit <- lcfa(data = HolzingerSwineford1939,
            model = model, group = "school",
            estimator = "ml",
            ordered = FALSE, std.lv = TRUE,
            mimic = "latent",
            positive = TRUE, penalties = TRUE,
            do.fit = TRUE,
            control = list(opt = "lbfgs", maxit = 284L,
                           cores = 1L, rstarts = 1L,
                           print = TRUE, print_interval = 30))

fit@loglik # -3415.092 (ML)
fit@penalized_loglik # -3415.095 (ML)
fit@loss # 0.1419955 (ULS) / 0.3327037 (ML)
fit@Optim$iterations
fit@Optim$convergence
fit@timing

# Plot model fit info:
fit

# Get fit indices:
getfit(fit)

# Inspect model objects:
latInspect(fit, what = "loadings", digits = 3)
latInspect(fit, what = "psi", digits = 3)
latInspect(fit, what = "uniquenesses", digits = 3)
latInspect(fit, what = "model", digits = 3)

ps <- fit@transformed_pars[[2]]$pj_psi
crossprod(ps)

theta <- latInspect(fit, what = "theta", digits = 3)[[1]]
det(theta)
eigen(theta)$values
sum(eigen(theta)$values)
# solve(theta)

psi <- latInspect(fit, what = "psi", digits = 3)[[1]]
det(psi)
eigen(psi)$values
sum(eigen(psi)$values)
# solve(psi)

#### Polychorics ####

library(latent)
samples <- unique(hexaco$sample) # industry mooc fire student dutch
Ns <- sapply(samples, FUN = function(x) sum(hexaco$sample == x))
names(Ns) <- samples

# Subset the items pertaining to the HEXACO-100
selection <- 5:104
selection <- 5:60
full <- hexaco[, selection]

mooc <- full[hexaco$sample == samples[2], ]
dim(mooc)
set.seed(2026)
fit <- lpoly(data = mooc, do.fit = TRUE, penalties = TRUE,
             control = list(opt = "lbfgs",
                            maxit = 500,
                            step_maxit = 50,
                            tcg_maxit = 30,
                            ss = 0.001,
                            eps = 1e-06,
                            print = TRUE,
                            print_interval = 10))
fit@loglik # -75695.53
fit@penalized_loglik # -75695.53
fit@Optim$iterations
fit@Optim$ng
fit@Optim$convergence
fit@Optim$f
max(fit@Optim$rg)
max(fit@Optim$g)
fit@timing

x <- Turbofuns:::PolychoricRM(as.matrix(mooc), estimate.acm = TRUE)
x$ACM
H[13:21, 13:21]

fit@modelInfo$poly_trans
fit@modelInfo$control_manifold
fit@modelInfo$control_estimator
fit@modelInfo$control$parameters

fit <- polyfast(as.matrix(mooc))
fit$iters
R <- fit$correlation
neg_value <- eigen(R)$values < 0

taus <- fit$thresholds
n <- fit$contingency_tables

#### CFA polychorics ####

samples <- unique(hexaco$sample) # industry mooc fire student dutch
Ns <- sapply(samples, FUN = function(x) sum(hexaco$sample == x))
names(Ns) <- samples

# Subset the items pertaining to the HEXACO-100
selection <- 5:104
full <- hexaco[, selection]

mooc <- full[hexaco$sample == samples[2], ]
dim(mooc)

model.EM <- "FEA =~ hexemfea146 + hexemfea170 + hexemfea74 + hexemfea2
             ANX =~ hexemanx128 + hexemanx8 + hexemanx80 + hexemanx176
             DEP =~ hexemdep62 + hexemdep182 + hexemdep134 + hexemdep158
             SEN =~ hexemsen44 + hexemsen164 + hexemsen20 + hexemsen68"

fit <- lcfa(model = model.EM, data = mooc,
            ordered = TRUE, estimator = "ml", do.fit = TRUE,
            control = NULL)
fit@loglik # -90154.77
fit@penalized_loglik # -90154.77
fit@loss # 90154.77
fit@Optim$iterations
fit@Optim$convergence
fit@timing

#### Check derivatives ####

control_manifold <- fit@modelInfo$control_manifold
control_transform <- fit@modelInfo$control_transform
control_estimator <- fit@modelInfo$control_estimator
control_optimizer <- fit@modelInfo$control
# control_optimizer$parameters[[1]] <- fit@Optim$parameters
# control_optimizer$transparameters[[1]] <- fit@Optim$transparameters
# control_optimizer$parameters[[1]] <- fit@modelInfo$control$parameters[[1]] +
#                                   rnorm(length(fit@modelInfo$control$parameters[[1]]), 0, 0.001)
# control_optimizer$transparameters[[1]] <- fit@modelInfo$control$transparameters[[1]]
x <- grad_comp(control_manifold, control_transform,
               control_estimator, control_optimizer,
               compute = "all",
               eps = 1e-07)
# x$f
# round(c(x$g) - c(x$numg), 5)
# max(abs(c(x$g) - c(x$numg)))
round(c(x$dg) - c(x$numdg), 5)
max(abs(c(x$dg) - c(x$numdg)))

x2 <- get_grad(control_manifold, control_transform,
               control_estimator, control_optimizer)
round(c(x2$g)-c(x$numg), 3)
max(abs(c(x$g) - c(x$numg)))

# x3 <- get_dgrad(control_manifold, control_transform,
#                 control_estimator, control_optimizer)
# round(c(x3$dg)-c(x$numdg), 3)
# round(c(x3$dg)/c(x$numdg), 3)

# Calculate the Hessian matrix using numerical approximations:
G <- function(parameters) {

  control_optimizer$parameters[[1]] <- parameters
  g <- get_grad(control_manifold = control_manifold,
                control_transform = control_transform,
                control_estimator = control_estimator,
                control_optimizer = control_optimizer)$g

  return(g)

}

H <- numDeriv::jacobian(func = G, x = control_optimizer$parameters[[1]])
H <- 0.5*(H + t(H)) # Force symmetry

x <- get_hess(control_manifold, control_transform,
              control_estimator, control_optimizer)
max(abs(H - x$h))
