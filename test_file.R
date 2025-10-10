# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 09/10/2025

#### Install latent ####

# devtools::install_github("marcosjnez/latent", force = TRUE)

#### LCA (multinomial) ####

library(latent)

# nmiss <- 30
# missrow <- sample(1:nrow(data), size = nmiss)
# misscol <- sample(1:ncol(data), size = nmiss, replace = TRUE)
# for(i in 1:nmiss) data[missrow[i], misscol[i]] <- NA
# control <- list(opt = "lbfgs", rstarts = 50L)

fit <- lca(data = gss82, nclasses = 3L,
           item = rep("multinomial", ncol(gss82)),
           penalties = TRUE, control = NULL, do.fit = TRUE)
fit@timing
fit@loglik # -2754.643
fit@penalized_loglik # -2759.507
fit@Optim$opt$iterations

# Plot model fit info:
fit

# Get fit indices:
getfit(fit)

# Get a summary:
summary(fit)

# Inspect model objects:
latInspect(fit, what = "coefs", digits = 3)
latInspect(fit, what = "classes", digits = 3)
latInspect(fit, what = "profile", digits = 3)
latInspect(fit, what = "posterior", digits = 3)
latInspect(fit, what = "table", digits = 3)
latInspect(fit, what = "pattern", digits = 3)

# Get standard errors:
SE <- se(fit, type = "standard", model = "user", digits = 4)
SE$table

# Get confidence intervals:
CI <- ci(fit, type = "standard", model = "user",
         confidence = 0.95, digits = 2)
CI$table

#### LCA (gaussian) ####

library(latent)

fit <- lca(data = empathy[, 1:6], nclasses = 4L,
           item = rep("gaussian", ncol(empathy[, 1:6])),
           penalties = TRUE, control = NULL, do.fit = TRUE)

fit@timing
fit@loglik # -1841.336
fit@penalized_loglik # -1844.333
fit@Optim$opt$iterations
fit@Optim$opt$convergence

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

# Control parameters for the optimizer:
# control <- list(opt = "em-lbfgs",
#                 rstarts_em = 50L, # Number of random starts with EM
#                 rstarts = 10L, # Number of random starts with LBFGS
#                 cores = 32L) # Number of cores
# The em-lbfgs optimizer runs EM 50 times, then pick up the 10 sets of
# parameter values with best loglik, and submits these sets to LBFGS to get
# the final convergence

fit <- lca(data = cancer[, 1:6], nclasses = 3L,
           item = c("gaussian", "gaussian",
                    "multinomial", "multinomial",
                    "gaussian", "gaussian"),
           penalties = list(class = list(alpha = 1), # Regularization for class probabilities
                            prob  = list(alpha = 1), # Regularization for conditional probabilities
                            sd    = list(alpha = 1)), # Regularization for standard deviations
           control = NULL, do.fit = TRUE)
fit@timing
fit@loglik # -5784.701
fit@penalized_loglik # -5795.573
fit@Optim$opt$iterations

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

#### LCA with covariates (gaussian) ####

# EM for covariates
# POBLQ in multigroup
# Automatic creation of structures

library(latent)

data <- empathy[, 1:6]
X <- as.matrix(empathy[, 7:8]) # Covariates

fit0 <- lca(data = data, X = NULL, model = NULL,
           item = rep("gaussian", ncol(data)),
           nclasses = 4L, penalties = TRUE,
           control = NULL, do.fit = TRUE)
fit0@timing
fit0@loglik # -1798.885
fit0@penalized_loglik # -1801.996
fit0@Optim$opt$iterations

fit <- lca(data = data, X = X, model = fit0,
           item = rep("gaussian", ncol(data)),
           nclasses = 4L, penalties = TRUE,
           control = NULL, do.fit = TRUE)
fit@timing
fit@loglik # -1798.885
fit@penalized_loglik # -1801.996
fit@Optim$opt$iterations

t(apply(fit@transformed_pars$beta, 1, FUN = \(x) x -sum(x)/4))

# Plot model fit info:
fit

# Get fit indices:
getfit(fit)

# Inspect model objects:
latInspect(fit, what = "coefs", digits = 5)
latInspect(fit, what = "classes", digits = 5)
latInspect(fit, what = "profile", digits = 3)
latInspect(fit, what = "posterior", digits = 3)

predict(fit, new = rbind(c(1, 1.571),
                         c(1, 2)))
fitted(fit)

# Get standard errors:
SE <- se(fit, type = "standard", model = "user", digits = 4)
SE$table

# Get confidence intervals:
CI <- ci(fit, type = "standard", model = "user",
         confidence = 0.95, digits = 2)
CI$table

#### CFA ####

library(latent)

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9'

fit <- cfast(HolzingerSwineford1939, model = model,
             estimator = "ml", cor = "pearson", do.fit = TRUE)
fit@loglik # -0.283407
fit@penalized_loglik # -0.283407
fit@Optim$opt$iterations

# control_manifold <- fit@Optim$control_manifold
# control_transform <- fit@Optim$control_transform
# control_estimator <- fit@Optim$control_estimator
# control <- fit@Optim$control
#
# x <- grad_comp(control_manifold, control_transform,
#                control_estimator, control, eps = 1e-04,
#                compute = "all")
# x$f
# round(c(x$g)-c(x$numg), 4)

fit@loglik # -0.283407 (ML)
fit@loss # 0.1574787 (ULS) / 0.283407 (ML)
fit@Optim$opt$iterations
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

library(lavaan)
fit2 <- cfa(model, data = HolzingerSwineford1939, estimator = "uls",
            std.lv = TRUE, std.ov = TRUE)
summary(fit2, fit.measures = FALSE)
fitMeasures(fit2, digits = 5)
inspect(fit2, what = "se")

SE <- se(fit, digits = 3)
SE$table_se

#### Multigroup CFA ####

library(latent)

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9'

# control <- list(opt = "lbfgs", maxit = 1000, rstarts = 1L, cores = 1L,
#                 eps = 1e-05)
fit <- cfast(HolzingerSwineford1939, model = model,
             group = "school", estimator = "ml", cor = "pearson",
             std.lv = TRUE, positive = FALSE,
             control = NULL, do.fit = TRUE)

fit@loglik # -0.1919254 (ml)
fit@loss # 0.2027554 (uls) / 0.3838508 (ml)
fit@Optim$opt$iterations

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
latInspect(fit, what = "rhat", digits = 3)

SE <- se(fit, type = "standard", model = "user", digits = 3)

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
fit2
inspect(fit2, what = "est")
fitmeasures(fit2, fit.measures = c("cfi", "tli", "rmsea", "srmr"))

# With latent:
fit <- cfast(data = HolzingerSwineford1939, model = model,
             estimator = "ml", cor = "pearson",
             positive = TRUE, penalties = TRUE,
             std.lv = TRUE, do.fit = TRUE,
             control = list(opt = "lbfgs", maxit = 1000L,
                            # logdetw = 1e-03,
                            cores = 10L, rstarts = 10L))

fit@loglik # -0.249485 (ML)
fit@penalized_loglik # -0.249485 (ML)
fit@loss # 0.1441886 (ULS) / 0.249485 (ML)
fit@Optim$opt$iterations
fit@Optim$opt$convergence

# Plot model fit info:
fit

# Get fit indices:
getfit(fit)

# Inspect model objects:
latInspect(fit, what = "loadings", digits = 3)
latInspect(fit, what = "psi", digits = 3)
latInspect(fit, what = "uniquenesses", digits = 3)
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
