# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 04/11/2025

#### Install latent ####

# devtools::install_github("marcosjnez/latent", force = TRUE)

#### LCA (multinomial) ####

library(latent)

fit <- lca(data = gss82, nclasses = 3,
           # multinomial = c("X1", "X2"),
           # poisson = ,
           # beta = ,
           # mimic = "LG",
           # model = formula(X1 ~ 1 + cluster1,
           #                 X2 ~ 1 + cluster2),
           item = rep("multinomial", ncol(gss82)),
           penalties = TRUE, control = NULL, do.fit = TRUE)
fit@timing
fit@loglik # -2754.643
fit@penalized_loglik # -2759.507
fit@Optim$opt$iterations
fit@Optim$opt$convergence
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
fit@Optim$opt$convergence
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

fit0 <- lca(data = data, X = X, model = NULL,
            item = rep("gaussian", ncol(data)),
            nclasses = 1:4L, penalties = TRUE,
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
fit@Optim$opt$convergence

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

set.seed(2025)
fit <- lcfa(HolzingerSwineford1939, model = model,
            estimator = "ml", do.fit = TRUE,
            cor = "pearson",
            # ordered = FALSE,
            # mimic = "lavaan",
            control = NULL)
fit@loss
fit@timing

# With lavaan:
fit2 <- lavaan::cfa(model, data = HolzingerSwineford1939,
                    estimator = "ml",
                    std.lv = TRUE, std.ov = TRUE,
                    likelihood = "wishart")
# Same loss value: OK
fit2@Fit@fx*2
fit@loss

fit@loglik
# WHAT'S THE DIFFERENCE?
fit2@loglik$loglik
fit2@h1$logl$loglik # SATURATED


control_manifold <- fit@modelInfo$control_manifold
control_transform <- fit@modelInfo$control_transform
control_estimator <- fit@modelInfo$control_estimator
control <- fit@modelInfo$control
control$parameters[[1]] <- fit@Optim$opt$parameters
control$transparameters[[1]] <- fit@Optim$opt$transparameters
x <- grad_comp(control_manifold, control_transform,
               control_estimator, control, eps = 1e-04,
               compute = "all")
x$f
round(c(x$g) - c(x$numg), 6)
round(c(x$g) / c(x$numg), 6)
round(c(x$dg) - c(x$numdg), 6)

fit@loss
fit@loglik # -0.283407
fit@penalized_loglik # -0.283407
fit@loss # 0.1574787
fit@Optim$opt$iterations
fit@Optim$opt$convergence
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

SE <- se(fit, digits = 3)
SE$table_se

#### Multigroup CFA ####

library(latent)

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9'

fit <- lcfa(HolzingerSwineford1939, model = model,
            group = "school", estimator = "uls", cor = "pearson",
            std.lv = TRUE, do.fit = TRUE, control = list(opt = "lbfgs"))

fit@timing
fit@Optim$opt$iterations

latInspect(fit, "weights")
latInspect(fit, "loglik")

fit@loglik # -3422.946 (ml)
fit@loss # 0.4055109 (uls) / 0.7677016 (ml)
fit@Optim$opt$iterations
fit@Optim$opt$convergence
fit@timing

# With lavaan:
fit2 <- lavaan::cfa(model, data = HolzingerSwineford1939,
                    group = "school", estimator = "ml",
                    std.lv = TRUE, std.ov = TRUE)

fit2@timing
fit2@Fit@fx*2
fit@loss

fit2@loglik$loglik

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
inspect(fit2, what = "est") # NEGATIVE VARIANCE
det(inspect(fit2, what = "est")$theta)
fit2@Fit@fx*2
fit2@loglik$loglik

# With latent:
fit <- lcfa(data = HolzingerSwineford1939, model = model,
            estimator = "ml", cor = "pearson",
            positive = TRUE,
            penalties = list(logdet = list(w = 0.001)),
            std.lv = TRUE, do.fit = TRUE,
            control = list(opt = "lbfgs", maxit = 100L,
                           cores = 20L, rstarts = 20L))

fit@loglik # -3421.613 (ML)
fit@penalized_loglik # -0.249485 (ML)
fit@loss # 0.1419955 (ULS) / 0.2467385 (ML)
fit@Optim$opt$iterations
fit@Optim$opt$convergence
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
            estimator = "ml", cor = "pearson",
            positive = TRUE, penalties = TRUE,
            std.lv = TRUE, control = NULL)

fit@loglik # -0.1678072 (ML)
fit@penalized_loglik # -0.1678072 (ML)
fit@loss # 0.1871894 (ULS) / 0.1678072 (ML)
fit@Optim$opt$iterations
fit@Optim$opt$convergence
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
             ANX ~= hexemanx128 + hexemanx8 + hexemanx80 + hexemanx176
             DEP ~= hexemdep62 + hexemdep182 + hexemdep134 + hexemdep158
             SEN ~= hexemsen44 + hexemsen164 + hexemsen20 + hexemsen68"

fit <- lcfa(model = model.EM, data = mooc,
            ordered = TRUE, estimator = "ml", do.fit = TRUE,
            control = NULL)
fit@loglik # -0.283407
fit@penalized_loglik # -0.283407
fit@loss # 0.1574787
fit@Optim$opt$iterations
fit@Optim$opt$convergence
fit@timing
