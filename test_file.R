# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 09/04/2026

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
fit@parameters

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
latInspect(fit, what = "classes", digits = 3)

fit@loglik # -1841.336
fit@penalized_loglik # -1844.333
fit@Optim$iterations # 66
fit@Optim$ng
max(fit@Optim$rg)
fit@Optim$convergence # TRUE
fit@timing
fit@parameters

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
SE <- se(fit, type = "standard", model = "model", digits = 4)
SE$table

# Get confidence intervals:
CI <- ci(fit, type = "standard", model = "model",
         confidence = 0.95, digits = 2)
CI$table

#### Mixed LCA (multinomial and gaussian) ####

library(latent)
set.seed(2026)
fit <- lca(data = cancer[, 1:6], nclasses = 3L,
           item = c("gaussian", "gaussian",
                    "multinomial", "multinomial",
                    "gaussian", "gaussian"),
           control = list(opt = "lbfgs",
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
SE <- se(fit, type = "robust", model = "model", digits = 4)
SE$table

# Get confidence intervals:
CI <- ci(fit, type = "standard", model = "model",
         confidence = 0.95, digits = 2)
CI$table

# hypothesis(fit, "b1|2 - b1|3 = 0")

#### LCA with covariates (gaussian) ####

library(latent)
data <- empathy[, 1:6]
X <- as.matrix(empathy[, 7:8]) # Covariates

set.seed(2026)
fit1 <- lca(data = data, X = NULL, model = NULL,
            item = rep("gaussian", ncol(data)),
            nclasses = 4L, penalties = TRUE,
            control = list(opt = "lbfgs",
                           step_maxit = 100,
                           tcg_maxit = 100),
            do.fit = TRUE)
fit1@loglik           # -1841.336
fit1@penalized_loglik # -1844.333
fit1@Optim$iterations
fit1@Optim$ng
fit1@timing

set.seed(2026)
Y <- as.matrix(empathy[, 9:10]) # Covariates
fit2 <- lca(data = data, X = cbind(X, Y), model = fit1,
           item = rep("gaussian", ncol(data)),
           nclasses = 4L, penalties = TRUE,
           do.fit = TRUE,
           control = list(opt = "lbfgs",
                          step_maxit = 100,
                          tcg_maxit = 100))
fit2@loglik # -1747.135
fit2@penalized_loglik # -1750.566
fit2@Optim$iterations
fit2@Optim$convergence
fit2@Optim$ng
fit2@timing

# check that the measurement model was fixed:
all.equal(fit1@parameters[-1], fit2@parameters[-1])

# Standard errors:
SE1 <- se(fit1, type = "standard", model = "model", digits = 4)
SE1$se
SE2 <- se(fit2, type = "standard", model = "model", digits = 4)
SE2$se

# Effects-coding parameterization:
new_se <- effects_coding(fit2@parameters$beta, SE2$vcov)
new_se$beta
new_se$table_se
new_se$se

# Plot model fit2 info:
fit2

# Get fit2 indices:
getfit2(fit2)

# Inspect model objects:
latInspect(fit2, what = "coefs", digits = 5)
latInspect(fit2, what = "classes", digits = 5)
latInspect(fit2, what = "profile", digits = 3)
latInspect(fit2, what = "posterior", digits = 3)

predict(fit2, new = rbind(c(2, 2, 2.428571, 2.142857),
                         c(1, 2, 3, 4)))
fitted(fit2)

# Get confidence intervals:
CI <- ci(fit2, type = "standard", model = "model",
         confidence = 0.95, digits = 2)
CI$table

x <- plot.llca(fit2,
          type = "standard",
          what = "OR",
          effects = "coding",
          confidence = 0.95,
          show_est_ci = TRUE,
          est_ci_header_cex = 0.5,
          cex_y = 0.5,
          mfrow = c(2, 2),
          xlim = c(0, 5))

# new_se <- move_intercept(beta, vcov)
# new_se$beta_new
# matrix(new_se$se_new, 3, 3)

#### CFA ####

library(latent)
library(lavaan)

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9'

set.seed(2026)
fit <- lcfa(HolzingerSwineford1939, model = model,
            estimator = "ml", positive = FALSE,
            ordered = FALSE, acov = "standard",
            std.lv = FALSE, std.ov = FALSE,
            mimic = "latent", do.fit = TRUE,
            meanstructure = TRUE,
            # likelihood = "wishart",
            control = NULL)
fit@loss             # 0.283407
fit@loglik           # -3737.745
fit@penalized_loglik # -3737.745
fit@Optim$iterations
fit@Optim$convergence
fit@timing

# With lavaan:
fit2 <- lavaan::cfa(model, data = HolzingerSwineford1939,
                    estimator = "ml",
                    # likelihood = "wishart",
                    std.lv = FALSE, std.ov = FALSE,
                    meanstructure = FALSE)
# Same loss value: OK
fit2@Fit@fx*2      # 0.283407
fit@loss           # 0.283407
fit2@loglik$loglik # -3737.745
fit@loglik         # -3737.745

fit@Optim$SE$se
fit2@ParTable$se

# NACOV <- lavTech(fit2, "gamma")
# lapply(NACOV, FUN = dim)

#### Multigroup CFA ####

library(latent)

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9'

fit <- lcfa(HolzingerSwineford1939, model = model,
            group = "school", estimator = "ml",
            ordered = FALSE, std.lv = FALSE,
            std.ov = FALSE, likelihood = "normal",
            mimic = "latent", do.fit = TRUE)

fit@loss             # 0.3848882
fit@loglik           # -3682.198
fit@penalized_loglik # -3682.198
fit@Optim$iterations
fit@Optim$convergence
fit@timing

# With lavaan:
fit2 <- lavaan::cfa(model, data = HolzingerSwineford1939,
                    group = "school", estimator = "ml",
                    likelihood = "normal",
                    std.lv = FALSE, std.ov = FALSE)
fit2@loglik$loglik # -3682.198
fit@loglik         # -3682.198
fit2@Fit@fx*2      # 0.3848882
fit@loss           # 0.3848882

fit@Optim$SE$se
fit2@ParTable$se

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

# With latent:
set.seed(2026)
fit <- lcfa(data = HolzingerSwineford1939, model = model,
            estimator = "dwls", positive = TRUE,
            penalties = list(logdet = list(w = 0.01)),
            ordered = FALSE, std.lv = TRUE,
            mimic = "latent", do.fit = TRUE,
            control = NULL)

fit@loglik # -3422.761 (ML)
fit@penalized_loglik # -3422.766 (ML)
fit@loss # 0.1419955 (ULS) / 0.2467385 (ML)
fit@Optim$iterations
fit@Optim$convergence
fit@Optim$ng
fit@timing
fit@Optim$SE$se

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
            estimator = "ml", ordered = FALSE,
            std.lv = TRUE, mimic = "latent",
            positive = TRUE, penalties = TRUE,
            do.fit = TRUE, control = NULL)

fit@loglik # -3415.092 (ML)
fit@penalized_loglik # -3415.095 (ML)
fit@loss # 0.1419955 (ULS) / 0.3327037 (ML)
fit@Optim$iterations
fit@Optim$convergence
fit@timing

#### Polychorics ####

library(latent)
samples <- unique(hexaco$sample) # industry mooc fire student dutch
Ns <- sapply(samples, FUN = function(x) sum(hexaco$sample == x))
names(Ns) <- samples

# Subset the items pertaining to the HEXACO-100
selection <- 5:104
selection <- 5:10
full <- hexaco[, selection]

mooc <- full[hexaco$sample == samples[2], ]
dim(mooc)
set.seed(2026)
POLY <- polyfast(as.matrix(mooc))
taus <- lapply(POLY$thresholds, FUN = \(x) x[-c(1, length(x))])
fit <- lpoly(data = mooc,
             # model = list(taus = taus),
             method = "two-step",
             positive = TRUE,
             penalties = TRUE,
             do.fit = TRUE,
             control = list(opt = "grad",
                            maxit = 500,
                            step_maxit = 50,
                            tcg_maxit = 30,
                            ss = 0.001,
                            eps = 1e-06))
fit@loglik # -176520.8
fit@penalized_loglik # -176520.8
fit@Optim$iterations
fit@Optim$ng
fit@Optim$convergence
fit@Optim$f # 41.18545
max(fit@Optim$rg)
max(fit@Optim$g)
fit@timing

fit@modelInfo$control$parameters[[1]] <- fit@Optim$parameters
fit@modelInfo$control$transparameters[[1]] <- fit@Optim$transparameters
x <- get_hess(fit@modelInfo$control_manifold, fit@modelInfo$control_transform,
              fit@modelInfo$control_estimator, fit@modelInfo$control,
              cores = 32L)
ACOV <- solve(x$h)

Tur <- Turbofuns:::PolychoricRM(as.matrix(mooc), estimate.acm = TRUE)

#### CFA (polychorics) ####

library(latent)
library(lavaan)

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
            ordered = TRUE, estimator = "dwls",
            # meanstructure = TRUE,
            std.ov = FALSE, std.lv = FALSE,
            do.fit = TRUE, control = NULL)
fit@loglik           # -90154.77 (ml)
fit@penalized_loglik # -90154.77 (ml)
fit@loss             # 0.3817476 (dwls)
fit@Optim$iterations
fit@Optim$convergence
fit@timing

# With lavaan:
fit2 <- lavaan::cfa(model = model.EM, data = mooc,
                    ordered = TRUE,
                    estimator = "dwls",
                    # likelihood = "wishart",
                    std.ov = FALSE, std.lv = FALSE,
                    parameterization = "theta")
fit2@Fit@fx*2      # 0.4663271
fit@loss           # 0.3817476
fit2@loglik$loglik # -3737.745
fit@loglik         # -90154.77

fit@Optim$SE$se
fit2@ParTable$se

# thresholds only
TH <- lavInspect(fit2, "th")

#### Multigroup CFA (polychorics) ####

library(latent)
library(lavaan)

samples <- unique(hexaco$sample) # industry mooc fire student dutch
Ns <- sapply(samples, FUN = function(x) sum(hexaco$sample == x))
names(Ns) <- samples

# Subset the items pertaining to the HEXACO-100
selection <- 5:104
full <- hexaco[, selection]

mooc <- full[hexaco$sample == samples[2], ]
dim(mooc)
mooc$school <- rep(c("s1", "s2"), times = c(2000, 2286))

model.EM <- "FEA =~ hexemfea146 + hexemfea170 + hexemfea74 + hexemfea2
             ANX =~ hexemanx128 + hexemanx8 + hexemanx80 + hexemanx176
             DEP =~ hexemdep62 + hexemdep182 + hexemdep134 + hexemdep158
             SEN =~ hexemsen44 + hexemsen164 + hexemsen20 + hexemsen68"

fit <- lcfa(model = model.EM, data = mooc,
            ordered = TRUE, estimator = "dwls",
            group = "school",
            do.fit = TRUE, control = NULL)
fit@loglik # -45281.39
fit@penalized_loglik # -45281.39
fit@loss # 0.394728
fit@Optim$iterations
fit@Optim$convergence
fit@timing
fit@Optim$SE$se

# With lavaan:
fit2 <- lavaan::cfa(model = model.EM, data = mooc,
                    ordered = TRUE,
                    estimator = "dwls",
                    group = "school",
                    # likelihood = "wishart",
                    std.lv = TRUE, std.ov = TRUE,
                    parameterization = "theta")
# Same loss value: OK
fit2@Fit@fx*2
fit@loss

lavaan::inspect(fit2, what = "se")$lambda
round(fit@Optim$SE$table_se$lambda.group1, 3)
diag(lavaan::inspect(fit2, what = "se")$theta)
round(diag(fit@Optim$SE$table_se$theta.group1), 3)
lavaan::inspect(fit2, what = "se")$psi
round(fit@Optim$SE$table_se$psi.group1, 3)

W <- lavInspect(fit2, "wls.v")
W2 <- matrix(0, 16, 16)
W2[lower.tri(W2, diag = TRUE)] <- diag(W)
W2 <- (W2+t(W2))/2
round(W2, 2)
round(fit@data_list$correl[[1]]$W, 2)

#### CFA (FIML) ####

library(latent)
library(lavaan)

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9'

data_missing <- HolzingerSwineford1939
data_missing[1:19, "x1"] <- NA

set.seed(2026)
fit <- lcfa(data_missing,
            model = model,
            estimator = "ml",
            positive = FALSE,
            ordered = FALSE,
            std.lv = FALSE,
            # missing = "pairwise.complete.obs",
            missing = "fiml",
            do.fit = TRUE,
            std.ov = FALSE,
            # meanstructure = TRUE,
            control = NULL)

fit@loss   # 0.283407
fit@loglik # -3714.303
fit@penalized_loglik # -3714.303
fit@Optim$iterations
fit@Optim$convergence
fit@timing

# With lavaan:
fit2 <- cfa(model, data = data_missing,
            estimator = "ml",
            missing = "fiml",
            # likelihood = "wishart",
            # meanstructure = TRUE,
            std.lv = FALSE, std.ov = FALSE)

fit2@loglik$loglik
fit@loglik
fit2@Fit@fx
fit@loss

lavInspect(fit2, "se")
fit@Optim$SE$se

inspect(fit2, what = "est")$lambda
fit@transformed_pars$lambda.group

#### Check derivatives ####

control_manifold <- fit@modelInfo$control_manifold
control_transform <- fit@modelInfo$control_transform
control_estimator <- fit@modelInfo$control_estimator
control_optimizer <- fit@modelInfo$control_optimizer
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
round(c(x$g) - c(x$numg), 5)
max(abs(c(x$g) - c(x$numg)))
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

VCOV <- get_vcov(control_manifold, control_transform, control_estimator,
                 control_optimizer, x$h)
VCOV$vcov
diag(VCOV$vcov)

x <- get_jacob(control_manifold, control_transform,
          control_estimator, control_optimizer)
x

# saveRDS(list(fit@modelInfo$control_manifold,
#              fit@modelInfo$control_transform,
#              fit@modelInfo$control_estimator,
#              fit@modelInfo$control),
#         file = "C:/Users/marco/OneDrive/Documentos/deletethis.rds")
# X <- readRDS("C:/Users/marco/OneDrive/Documentos/deletethis.rds")
# X[[3]][[1]][3] <- NULL; X[[3]][[2]][3] <- NULL
# all.equal(fit@modelInfo$control_manifold, X[[1]])
# all.equal(fit@modelInfo$control_transform, X[[2]])
# all.equal(fit@modelInfo$control_estimator, X[[3]])
# all.equal(fit@modelInfo$control, X[[4]])
# fit@modelInfo$control$parameters
# X[[4]]$parameters
# fit@modelInfo$param
