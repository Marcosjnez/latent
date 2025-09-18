#### Install latent ####

devtools::install_github("marcosjnez/latent", force = TRUE)

#### LCA (multinomial) ####

library(latent)

data <- gss82
item <- rep("multinomial", ncol(data))

nclasses <- 3L

nmiss <- 30
missrow <- sample(1:nrow(data), size = nmiss)
misscol <- sample(1:ncol(data), size = nmiss, replace = TRUE)
for(i in 1:nmiss) data[missrow[i], misscol[i]] <- NA
# control <- list(opt = "em", rstarts = 50L)

fit <- lca(data = data, item = item, nclasses = nclasses,
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
latInspect(fit, what = "classes", digits = 3)
latInspect(fit, what = "profile", digits = 3)
latInspect(fit, what = "posterior", digits = 3)
latInspect(fit, what = "table", digits = 3)
latInspect(fit, what = "pattern", digits = 3)

# Get confidence intervals:
CI <- ci(fit, type = "standard", model = "user",
         confidence = 0.95, digits = 2)
CI$table

#### LCA (gaussian) ####

library(latent)

nclasses <- 4L
data <- empathy[, 1:6]

item <- rep("gaussian", ncol(data))

fit <- lca(data = data, item = item, nclasses = nclasses,
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
latInspect(fit, what = "classes", digits = 3)
latInspect(fit, what = "profile", digits = 3)
latInspect(fit, what = "posterior", digits = 3)

# Get confidence intervals:
CI <- ci(fit, type = "standard", model = "user",
         confidence = 0.95, digits = 2)
CI$table

#### Mixed LCA (multinomial and gaussian) ####

library(latent)
data <- cancer[, 1:6]
names(data)

item <- c("gaussian", "gaussian",
          "multinomial", "multinomial",
          "gaussian", "gaussian")

nclasses <- 3L

# Regularization:
penalties <- list(
  class = list(alpha = 1), # Regularization for class probabilities
  prob  = list(alpha = 1), # Regularization for conditional probabilities
  sd    = list(alpha = 1) # Regularization for standard deviations
)

# Control parameters for the optimizer:
control <- list(opt = "em-lbfgs",
                rstarts_em = 50L, # Number of random starts with EM
                rstarts = 10L, # Number of random starts with LBFGS
                cores = 32L) # Number of cores
# The em-lbfgs optimizer runs EM 50 times, then pick up the 10 sets of
# parameter values with best loglik, and submits these sets to LBFGS to get
# the final convergence

fit <- lca(data = data, item = item, nclasses = nclasses,
           penalties = penalties, control = control,
           do.fit = TRUE)
fit@timing
fit@loglik # -5784.701
fit@penalized_loglik # -5795.573
fit@Optim$opt$iterations

# Plot model fit info:
fit

# Get fit indices:
getfit(fit)

# Inspect model objects:
latInspect(fit, what = "classes", digits = 3)
latInspect(fit, what = "profile", digits = 3)
latInspect(fit, what = "posterior", digits = 3)

# Get confidence intervals:
CI <- ci(fit, type = "standard", model = "user",
         confidence = 0.95, digits = 2)
CI$table

#### CFA ####

library(latent)

data <- HolzingerSwineford1939

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9'

estimator = "uls"
control <- list(opt = "lbfgs", maxit = 1000L, rstarts = 1L,
                cores = 1L, eps = 1e-05)
fit <- cfast(data, model = model,
             estimator = estimator, cor = "pearson",
             std.lv = TRUE, positive = FALSE,
             control = control, do.fit = TRUE)

fit@loglik # -0.283407 (ML)
fit@loss # 0.1574787 (ULS) / 0.283407 (ML)
fit@Optim$opt$iterations

# Plot model fit info:
fit

# Get fit indices:
getfit(fit)

# Get a summary:
summary(fit)

lapply(fit@transformed_pars, FUN = \(x) round(x$lambda, 3))
lapply(fit@transformed_pars, FUN = \(x) round(x$psi, 3))
lapply(fit@transformed_pars, FUN = \(x) round(x$theta, 3))

SE <- se(fit, type = "standard", model = "user", digits = 3)
round(cbind(SE$se, SE$z, SE$pval), 3)
getfit(fit)

library(lavaan)
fit2 <- cfa(model, data = HolzingerSwineford1939, estimator = "uls",
            std.lv = TRUE, std.ov = TRUE)
summary(fit2, fit.measures = FALSE)
fitMeasures(fit2)
inspect(fit2, what = "est")
inspect(fit2, what = "std")

#### Multigroup CFA ####

library(latent)

data <- lavaan::HolzingerSwineford1939

model <- 'visual  =~ c(a,a)*x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9'

estimator = "uls"
group <- "school"
control <- list(opt = "lbfgs", maxit = 1000, rstarts = 1L, cores = 1L,
                eps = 1e-05)
fit <- cfast(data, model = model, group = group,
             estimator = estimator, cor = "pearson",
             std.lv = TRUE, positive = FALSE,
             control = control, do.fit = TRUE)

fit@loglik # -0.7809653 (ml)
fit@loss # 0.4269743 (uls) / 0.7809653 (ml)
fit@Optim$opt$iterations

fit@transformed_pars
