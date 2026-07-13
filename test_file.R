# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 12/07/2026

#### Store a dataset ####

# usethis::use_data(gss82, overwrite = TRUE)

#### Build the package ####

# roxygen2::roxygenise()

#### Install latent ####

# devtools::install_github("marcosjnez/latent", force = TRUE)

#### LCA (multinomial) ####

library(latent)
set.seed(2026)

gss82$EDUCR <- as.integer(gss82$EDUCR)-1L
fit <- lca(data = gss82,
           nclasses = 3L,
           multinomial = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"),
           covariates = c("RACE", "SEX", "EDUCR", "AGE"),
           # outcomes = "MARITAL",
           # adjustment = "bk",
           # classification = "modal",
           # model = list("UNDERSTA ~~ COOPERAT
           #               PURPOSE ~~ COOPERAT"),
           penalties = list(class = list(alpha=1),
                            prob  = list(alpha=0)),
           do.fit = TRUE)
latInspect(fit, what = "loglik")
# loglik: -3891.252 # penalized_loglik: -3892.478
# loglik: -3879.167 # penalized_loglik: -3880.371 ("UNDERSTA ~~ COOPERAT
#                                                   PURPOSE ~~ COOPERAT")

# Print model fit info:
fit

# Get fit indices:
getfit(fit)

# Get a summary:
summary(fit)

# Bivariate residuals:
lbvr(fit)

# Plot:
plot(fit)
plot(fit, type = "coefficients", what = "OR")

# Predictions if there are covariates:
df <- data.frame(RACE = c("BLACK", "WHITE"),
                 SEX = c("MALE", "MALE"),
                 EDUCR = c(3, 3),
                 AGE = c(35, 35))
predict(fit, new = df)

# Inspect model objects:
latInspect(fit, what = "convergence")
latInspect(fit, what = "profile")
latInspect(fit, what = "coefs")
latInspect(fit, what = "classification")
latInspect(fit, what = "pattern")
latInspect(fit, what = "table")

# Get standard errors:
SE <- se(fit, type = "standard", digits = 4)
SE$table

# Get confidence intervals:
CI <- ci(fit, type = "standard", confidence = 0.95, digits = 2)
CI$table

#### LCA (gaussian) ####

library(latent)
set.seed(2026)

fit <- lca(data = empathy,
           nclasses = 4L,
           gaussian = c("ec1", "ec2", "ec3", "ec4", "ec5", "ec6"),
           # covariates = c("pt1", "pt2", "pt3", "pt4"),
           # outcomes = c("pt1", "pt2"),
           # model = list("ec2 ~~ ec3"),
           penalties = TRUE,
           control = list(rstarts = 50L, cores = 32L),
           do.fit = TRUE)
latInspect(fit, what = "loglik")
# loglik: -1841.336 # penalized_loglik: -1844.333
# loglik: -1814.483 # penalized_loglik: -1817.526 ("ec2 ~~ ec3")

# Print model fit info:
fit

# Get fit indices:
getfit(fit)

# Bivariate residuals:
lbvr(fit)

# Inspect model objects:
latInspect(fit, what = "convergence")
latInspect(fit, what = "profile")
latInspect(fit, what = "coefs")
latInspect(fit, what = "posterior")

# Plot:
plot(fit)

# Get standard errors:
SE <- se(fit, type = "standard", digits = 4)
SE$table

# Get confidence intervals:
CI <- ci(fit, type = "standard", confidence = 0.95, digits = 2) # FIX THIS
CI$table

#### Mixed LCA (multinomial and gaussian) ####

library(latent)
set.seed(2026)
penalties <- list(
  # beta  = list(alpha = 0, lambda = 0, power = 0),
  # beta  = list(alpha = 0),
  class = list(alpha = 1),
  prob  = list(alpha = 1),
  var   = list(alpha = 1),
  Sigma = list(alpha = 1)
) # FIX defaults in penalties
fit <- lca(data = cancer,
           nclasses = 3L,
           gaussian = c("Age", "WeightIndex", "SystolicBloodPressure",
                        "DiastolicBloodPressure"),
           multinomial = c("PerformanceRating", "CardiovascularDiseaseHistory"),
           penalties = penalties,
           do.fit = TRUE)
latInspect(fit, what = "loglik")
# loglik: -5784.701 # penalized_loglik: -5795.573

# Print model fit info:
fit

# Get fit indices:
getfit(fit)

# Inspect model objects:
latInspect(fit, what = "fit.matrix")
latInspect(fit, what = "coefs")
latInspect(fit, what = "classes")
latInspect(fit, what = "profile")
latInspect(fit, what = "posterior")

# Get standard errors:
SE <- se(fit, type = "standard", digits = 4) # FIX stdv
SE$table

# Get confidence intervals:
CI <- ci(fit, type = "standard", confidence = 0.95, digits = 2)
CI$table

# hypothesis(fit, "b1|2 - b1|3 = 0")

plot(fit)

#### LCA with covariates (gaussian) ####

fit <- lca(data = empathy,
            nclasses = 4L,
            gaussian = c("ec1", "ec2", "ec3", "ec4", "ec5", "ec6"),
            covariates = c("pt1", "pt2", "pt3", "pt4"),
            # outcomes = list(gaussian = c("pt5")),
            adjustment = "bk",
            # classification = "modal",
            # model = list("ec2 ~~ ec3"),
            penalties = TRUE,
            do.fit = TRUE)
latInspect(fit$structural, what = "loglik")
# loglik: -1747.135 # penalized_loglik: -1750.566
# loglik: -2049.840 # penalized_loglik: -2053.322 (outcomes = list(gaussian = c("pt5")))
# SE <- se(fit, type = "standard", digits = 4)
# SE$se

# # Effects-coding parameterization:
# new_se <- effects_coding(fit$structural@parameters$beta, SE$vcov)
# new_se$beta
# new_se$table_se
# new_se$se

# new_se <- move_intercept(beta, vcov)
# new_se$beta_new
# matrix(new_se$se_new, 3, 3)

# Print model fit:
fit

# Get fit indices:
getfit(fit)

plot.llcalist(fit)
plot.llcalist(fit, type = "coefficients", what = "OR", xlim = c(0, 7))

# Inspect model objects:
latInspect(fit, what = "loglik")
# loglik: -2049.840 # penalized_loglik: -2053.322
latInspect(fit, what = "coefs")
latInspect(fit, what = "classes")
latInspect(fit, what = "profile")
latInspect(fit, what = "posterior")

predict(fit, new = rbind(c(2, 2, 2.428571, 2.142857),
                         c(1, 2, 3, 4)))
fitted(fit)

# Get confidence intervals:
CI <- ci(fit, type = "standard", confidence = 0.95, digits = 2)
CI$table

#### CFA ####

library(latent)
library(lavaan)

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9'

set.seed(2026)
estimator <- "ml"
std.ov <- FALSE
std.lv <- FALSE
meanstructure <- TRUE
likelihood <- "normal"
acov <- "standard"

fit <- lcfa(HolzingerSwineford1939,
            model = model,
            estimator = estimator,
            acov = acov,
            std.ov = std.ov,
            std.lv = std.lv,
            meanstructure = meanstructure,
            likelihood = likelihood,
            se = TRUE,
            control = NULL,
            do.fit = TRUE)

latInspect(fit, "est")
latInspect(fit, "loglik") # loglik           -3737.745
                          # penalized_loglik -3737.745
                          # loglik_base      -4211.418
                          # loglik_sat       -3695.092
getfit(fit)

# With lavaan:
fit2 <- lavaan::cfa(data = HolzingerSwineford1939,
                    model = model,
                    estimator = estimator,
                    std.lv = std.lv,
                    std.ov = std.ov,
                    meanstructure = meanstructure,
                    likelihood = likelihood,
                    do.fit = TRUE)
# Same loss value: OK
fit2@Fit@fx*2      # 0.283407
fit2@loglik$loglik # -3737.745
latInspect(fit, "loss")
latInspect(fit, "loglik")

inspect(fit2, "est")
fit@parameters

fitmeasures(fit2)
getfit(fit)
latInspect(fit, "loss")
fitMeasures(fit2, "unrestricted.logl")
fitMeasures(fit2, "baseline.chisq")
fitmeasures(fit2, "rmr")

# NACOV <- lavTech(fit2, "gamma")
# lapply(NACOV, FUN = dim)

#### Multigroup CFA ####

library(latent)

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9'

estimator <- "ml"
std.ov <- FALSE
std.lv <- FALSE
meanstructure <- TRUE
likelihood <- "normal"
acov <- "standard"

fit <- lcfa(HolzingerSwineford1939,
            model = model,
            group = "school",
            estimator = estimator,
            acov = acov,
            std.ov = std.ov,
            std.lv = std.lv,
            meanstructure = meanstructure,
            likelihood = likelihood,
            se = TRUE,
            control = NULL,
            do.fit = TRUE)

latInspect(fit, "est")
latInspect(fit, "loglik") # loglik           -3682.198
                          # penalized_loglik -3682.198
                          # loglik_base      -4150.500
                          # loglik_sat       -3624.272

latInspect(fit, what = "loss")

# With lavaan:
fit2 <- lavaan::cfa(data = HolzingerSwineford1939,
                    model = model,
                    group = "school",
                    estimator = estimator,
                    std.ov = std.ov,
                    std.lv = std.lv,
                    meanstructure = meanstructure,
                    likelihood = likelihood,
                    do.fit = TRUE)
fit2@loglik$loglik # -3682.198
fit2@Fit@fx*2      # 0.3848882
fitMeasures(fit2, "unrestricted.logl")

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

# With latent:
estimator <- "ml"
std.ov <- FALSE
std.lv <- FALSE
meanstructure <- TRUE
likelihood <- "normal"
acov <- "standard"
fit <- lcfa(data = HolzingerSwineford1939, model = model,
            estimator = "ml", ordered = FALSE,
            positive = TRUE,
            # penalties = TRUE,
            penalties = list(logdet = list(w = 0.001)),
            std.lv = FALSE, std.ov = FALSE,
            likelihood = "normal",
            do.fit = TRUE, control = NULL)

latInspect(fit, "loglik") # loglik           -3732.196
                          # penalized_loglik -3732.198
                          # loglik_base      -
                          # loglik_sat       -

# fit@Optim$SE$se
det(fit@transformed_pars$theta)
det(fit@transformed_pars$psi)

# With lavaan:
fit2 <- cfa(model, data = HolzingerSwineford1939,
            estimator = "ml", std.lv = FALSE, std.ov = FALSE)
inspect(fit2, what = "est") # NEGATIVE VARIANCE
det(inspect(fit2, what = "est")$theta)
fit2@Fit@fx*2

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

fit@loglik # -3674.3 (ML)
fit@penalized_loglik # -3674.298 (ML)
fit@loss # 0.1419955 (ULS) / 0.3343969 (ML)
fit@penalized_loss # 0.3363794 (ML)
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
             method = "one-step",
             positive = F,
             penalties = F,
             do.fit = TRUE,
             control = list(opt = "grad",
                            # subfix = ".group1",
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
fit@Optim$f # 41.18626
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

fit <- lcfa(data = mooc,
            model = model.EM,
            ordered = TRUE,
            estimator = "dwls",
            std.ov = FALSE,
            std.lv = FALSE,
            meanstructure = FALSE,
            se = TRUE,
            control = NULL,
            do.fit = TRUE)
fit@loglik           # -90154.77 (ml)
fit@penalized_loglik # -90154.77 (ml)
fit@loss             # 0.3845136 (dwls)
fit@Optim$iterations
fit@Optim$convergence
fit@timing

# With lavaan:
fit2 <- lavaan::cfa(data = mooc,
                    model = model.EM,
                    ordered = TRUE,
                    estimator = "dwls",
                    meanstructure = FALSE,
                    std.ov = FALSE,
                    std.lv = FALSE,
                    # parameterization = "theta",
                    parameterization = "delta",
                    do.fit = TRUE)
fit2@Fit@fx*2      # 0.4663271
fit@loss           # 0.3817476
fit2@loglik$loglik # -3737.745
fit@loglik         # -90154.77

lavaan::inspect(fit2, "est")$delta
round(fit@parameters$delta, 3)

diag(lavaan::inspect(fit2, "est")$theta)
diag(round(fit@parameters$theta., 3))

lavaan::inspect(fit2, "est")$lambda
round(fit@parameters$lambda., 3)

lavaan::inspect(fit2, "est")$psi
round(fit@parameters$psi., 3)

# fit@Optim$SE$se
# fit2@ParTable$se

#### CFA (Yule correlation) ####

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
            ordered = "yule", estimator = "dwls",
            positive = FALSE, penalties = FALSE,
            std.ov = FALSE,
            std.lv = FALSE,
            meanstructure = FALSE,
            se = TRUE,
            do.fit = TRUE,
            control = NULL)
fit@loglik           # -90154.77 (ml)
fit@penalized_loglik # -90154.77 (ml)
fit@loss             # 497.0507 (dwls)
fit@Optim$iterations
fit@Optim$convergence
fit@timing

fit@Optim$SE$se

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
data_missing[1:3, "x1"] <- NA

set.seed(2026)
fit <- lcfa(data_missing,
            model = model,
            estimator = "ml",
            positive = FALSE,
            ordered = FALSE,
            std.lv = FALSE,
            # missing = "pairwise.complete.obs",
            missing = "fiml",
            std.ov = FALSE,
            se = TRUE,
            do.fit = TRUE,
            control = NULL)

fit@loss   # 0.9767859
fit@loglik # -3733.769
fit@penalized_loglik # -3733.769
fit@Optim$iterations
fit@Optim$convergence
fit@timing

latInspect(fit, what = "loglik")

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
fit2@ParTable$se
fit@Optim$SE$se

#### lpearson ####

library(latent)

data <- HolzingerSwineford1939[, paste("x", 1:9, sep = "")]
fit <- lpearson(data = data, std.ov = FALSE,
                acov = "standard", likelihood = "normal",
                missing = "pairwise.complete.obs",
                do.fit = TRUE)
fit@parameters

#### Yule correlations ####

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

fit <- lyule(data = as.matrix(mooc),
             do.fit = TRUE)

#### Check derivatives ####

control_manifold <- fit@modelInfo$control_manifold
control_transform <- fit@modelInfo$control_transform
control_estimator <- fit@modelInfo$control_estimator
control_optimizer <- fit@modelInfo$control_optimizer

x <- grad_comp(control_manifold, control_transform,
               control_estimator, control_optimizer,
               compute = "all",
               eps = 1e-07)
x$f # 4362.65 # 13800.13

round(c(x$g) - c(x$numg), 5)
max(abs(c(x$g) - c(x$numg)))
round(c(x$dg) - c(x$numdg), 5)
max(abs(c(x$dg) - c(x$numdg)))

Optim <- optimizer(control_manifold, control_transform,
                   control_estimator, control_optimizer)
Optim$f

x2 <- get_grad(control_manifold, control_transform,
               control_estimator, control_optimizer)
round(c(x2$g)-c(x$numg), 3)
max(abs(c(x$g) - c(x$numg)))

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

#### To-do ####
# Fix class ordering by size
# Create a class and methods for measurement + structural objects
# loglik from BK¿?
# SE for ML-modal-prop
