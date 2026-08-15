# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 15/08/2026

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
set.seed(20216)
fit <- lca(data = gss82,
           nclasses = 3L,
           multinomial = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"),
           covariates = c("RACE", "SEX", "EDUCR", "AGE"),
           outcomes = "MARITAL",
           # adjustment = "bk",
           # classification = "modal",
           # model = list("UNDERSTA ~~ COOPERAT
           #               PURPOSE ~~ COOPERAT"),
           # start = start,
           penalties = list(class = list(alpha=1),
                            prob  = list(alpha=1)),
           # control = list(opt = "em", rstarts = 30, cores = 30,
           #                maxit = 50L, eps = 1e-05, step_maxit = 30L,
           #                mopt = "grad", mstep_maxit = 20L, mstep_eps = 1e-05),
           do.fit = TRUE)
latInspect(fit, what = "loglik")
latInspect(fit, what = "convergence")
latInspect(fit, what = "elapsed")
# loglik: -3891.252 # penalized_loglik: -3892.478
# loglik: -3879.167 # penalized_loglik: -3880.371 ("UNDERSTA ~~ COOPERAT
#                                                   PURPOSE ~~ COOPERAT")

dconstraints <- constraints_derivs(fit)
labs <- c(fit@modelInfo$trans$PURPOSE)
dconstraints[labs, ]
colnames(dconstraints)
dim(dconstraints)

jacob <- jacobian(fit)
labs1 <- c(fit@modelInfo$trans$PURPOSE[1, ])
labs2 <- c(fit@modelInfo$trans$PURPOSE[4, ])
jacob[labs1, labs2]

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
SE <- se(fit, type = "information", digits = 4)
SE$table

# Get confidence intervals:
CI <- ci(fit, type = "information", confidence = 0.95, digits = 2)
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
SE <- se(fit, type = "information", digits = 4)
SE$table

# Get confidence intervals:
CI <- ci(fit, type = "information", confidence = 0.95, digits = 2) # FIX THIS
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
           # control = list(opt = "em", rstarts = 30, cores = 30,
           #                maxit = 50L, eps = 1e-05, step_maxit = 30L,
           #                mopt = "grad", mstep_maxit = 20L, mstep_eps = 1e-05),
           do.fit = TRUE)
latInspect(fit, what = "loglik")
# loglik: -5784.701 # penalized_loglik: -5795.573
latInspect(fit, what = "convergence")
fit@Optim$elapsed

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
SE <- se(fit, type = "information", digits = 4)
SE$table

# Get confidence intervals:
CI <- ci(fit, type = "information", confidence = 0.95, digits = 2)
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
plot.llcalist(fit, type = "coefficients", what = "OR")

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

# Get standard errors:
SE <- se(fit, type = "standard", digits = 4)
SE$table

# Get confidence intervals:
CI <- ci(fit, type = "standard", confidence = 0.95, digits = 2)
CI$table

#### CFA ####

library(latent)
library(lavaan)

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9'
S <- cov(HolzingerSwineford1939[, paste("x", 1:9, sep = "")])
means <- colMeans(HolzingerSwineford1939[, paste("x", 1:9, sep = "")])

set.seed(2026)
estimator <- "ml"
std.ov <- FALSE
std.lv <- FALSE
meanstructure <- TRUE
likelihood <- "normal"
acov <- "standard"

fit <- lcfa(model = model,
            data = HolzingerSwineford1939,
            # sample.cov = S,
            # sample.mean = means,
            # sample.nobs = 301,
            estimator = estimator,
            acov = acov,
            std.ov = std.ov,
            std.lv = std.lv,
            meanstructure = meanstructure,
            likelihood = likelihood,
            se = TRUE,
            control = NULL,
            do.fit = TRUE)

# With lavaan:
fit2 <- lavaan::cfa(model = model,
                    data = HolzingerSwineford1939,
                    # sample.cov = S,
                    # sample.mean = means,
                    # sample.nobs = 301,
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
latInspect(fit, "loglik") # loglik           -3737.745
                          # penalized_loglik -3737.745
                          # loglik_base      -4211.418
                          # loglik_sat       -3695.092

inspect(fit2, "est")
latInspect(fit, "est")
fit@Optim$SE$se
fit2@ParTable$se

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

# latInspect(fit, "est")
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
lavaan::fitMeasures(fit2, "unrestricted.logl")

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
            estimator = estimator, ordered = FALSE,
            positive = TRUE,
            # penalties = TRUE,
            penalties = list(logdet = list(w = 0.001)),
            std.lv = std.lv, std.ov = std.ov,
            likelihood = likelihood,
            do.fit = TRUE, control = NULL)

latInspect(fit, "loglik") # loglik           -3732.196
                          # penalized_loglik -3732.197
                          # loglik_base      -4211.418
                          # loglik_sat       -3695.092

latInspect(fit, what = "est")

# fit@Optim$SE$se
det(fit@transformed_pars$theta)
det(fit@transformed_pars$psi)

# With lavaan:
fit2 <- cfa(model, data = HolzingerSwineford1939,
            estimator = estimator, std.lv = std.lv, std.ov = std.ov,
            likelihood = likelihood)
fitMeasures(fit2, "logl")
fitMeasures(fit2, "unrestricted.logl")
fitMeasures(fit2, "baseline.chisq")

inspect(fit2, what = "est") # NEGATIVE VARIANCE
det(inspect(fit2, what = "est")$theta)

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

# With latent:
estimator <- "ml"
std.ov <- FALSE
std.lv <- FALSE
meanstructure <- TRUE
likelihood <- "normal"
acov <- "standard"
fit <- lcfa(data = HolzingerSwineford1939,
            model = model, group = "school",
            estimator = estimator, ordered = FALSE,
            std.ov = std.ov, std.lv = std.lv,
            positive = TRUE, penalties = TRUE,
            meanstructure = meanstructure,
            likelihood = likelihood,
            acov = acov,
            do.fit = TRUE, control = NULL)
latInspect(fit, "loglik") # loglik           -3674.300
                          # penalized_loglik -3674.302
                          # loglik_base      -4150.500
                          # loglik_sat       -3624.272

fit@Optim$iterations
fit@Optim$convergence
fit@timing

# With lavaan:
fit2 <- cfa(model, data = HolzingerSwineford1939,
            group = "school", estimator = estimator,
            std.lv = std.lv, std.ov = std.ov,
            meanstructure = meanstructure,
            likelihood = likelihood,)
fitmeasures(fit2, fit.measures = c("cfi", "tli", "rmsea", "srmr"))
fitMeasures(fit2, "logl")
fitMeasures(fit2, "unrestricted.logl")
fitMeasures(fit2, "baseline.chisq")

inspect(fit2, what = "est")$Pasteur$theta # NEGATIVE VARIANCE
det(inspect(fit2, what = "est")$Pasteur$theta)

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

estimator <- "dwls"
std.ov <- FALSE
std.lv <- FALSE
meanstructure <- TRUE
likelihood <- "normal"
acov <- "standard"
fit <- lcfa(data = mooc,
            model = model.EM,
            ordered = TRUE,
            estimator = estimator,
            std.ov = std.ov,
            std.lv = std.lv,
            meanstructure = meanstructure,
            acov = acov,
            se = TRUE,
            control = NULL,
            do.fit = TRUE)
latInspect(fit, "loss") # loss           0.3845136
                        # penalized_loss 0.3845136
                        # loss_base      6.4196570
                        # loss_sat       0.0000000
fit@Optim$iterations
fit@Optim$convergence
fit@timing
latInspect(fit, what = "lambda")

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

lavaan::inspect(fit2, "est")$delta
round(fit@parameters$delta, 3)

diag(lavaan::inspect(fit2, "est")$theta)
diag(round(fit@parameters$theta, 3))

lavaan::inspect(fit2, "est")$lambda
round(fit@parameters$lambda, 3)

lavaan::inspect(fit2, "est")$psi
round(fit@parameters$psi, 3)

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

estimator <- "dwls"
std.ov <- FALSE
std.lv <- FALSE
meanstructure <- TRUE
likelihood <- "normal"
acov <- "standard"
fit <- lcfa(model = model.EM, data = mooc,
            ordered = "yule", estimator = estimator,
            positive = FALSE, penalties = FALSE,
            std.ov = std.ov,
            std.lv = std.lv,
            meanstructure = meanstructure,
            se = TRUE,
            do.fit = TRUE,
            control = NULL)
latInspect(fit, "loss") # loss           0.1159708
                        # penalized_loss 0.1159708
                        # loss_base      2.6771636
                        # loss_sat       0.0000000
fit@Optim$iterations
fit@Optim$convergence
fit@timing

latInspect(fit, what = "lambda")
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

#### lmean ####

# Estimate means of variables
library(latent)

df <- HolzingerSwineford1939[, paste("x", 1:9, sep = "")]
set.seed(2026)
std.ov <- FALSE

fit <- lmean(data = df, model = NULL, std.ov = std.ov)
fit@parameters$means
fit@Optim$SE$se
# colMeans(df)
# sqrt(apply(df, MARGIN = 2, var)/nrow(df))

#### lpearson ####

# Pearson correlations

library(latent)

df <- HolzingerSwineford1939[, paste("x", 1:9, sep = "")]
std.ov <- FALSE
likelihood <- "normal"
acov <- "standard"
fit <- lpearson(data = df, std.ov = std.ov,
                acov = acov, likelihood = likelihood,
                missing = "pairwise.complete.obs",
                do.fit = TRUE)
fit@parameters
fit@Optim$SE$se

#### lpoly ####

# Polychoric correlations

library(latent)
samples <- unique(hexaco$sample) # industry mooc fire student dutch
Ns <- sapply(samples, FUN = function(x) sum(hexaco$sample == x))
names(Ns) <- samples

# Subset the items pertaining to the HEXACO-100
selection <- 5:104
selection <- 5:10
full <- hexaco[, selection]

mooc <- full[hexaco$sample == samples[2], ]
# dim(mooc)
# POLY <- polyfast(as.matrix(mooc))
# taus <- lapply(POLY$thresholds, FUN = \(x) matrix(x[-c(1, length(x))],
#                                                   ncol = 1))
#
# names(taus) <- paste("taus", colnames(mooc), sep = "")
# for(j in 1:length(taus)) {
#   colnames(taus[[j]]) <- colnames(mooc)[j]
#   rownames(taus[[j]]) <- 1:length(taus[[j]])
# }

set.seed(2026)
fit <- lpoly(data = mooc,
             method = "one-step",
             positive = FALSE,
             penalties = FALSE,
             # model = fit@parameters[-7],
             # start = fit@parameters,
             do.fit = TRUE,
             control = list(opt = "grad",
                            # subfix = ".group1",
                            maxit = 500,
                            step_maxit = 50,
                            tcg_maxit = 30,
                            ss = 0.001,
                            eps = 1e-06))
# fit@loglik # -176520.8
# fit@penalized_loglik # -176520.8
fit@Optim$iterations
fit@Optim$ng
fit@Optim$convergence
fit@Optim$f # 41.18545 # 41.18626
max(fit@Optim$rg)
max(fit@Optim$g)
fit@timing
fit@modelInfo$param
fit@modelInfo$parameters_labels

fit@parameters
fit@Optim$SE$se
VCOV <- vcov(fit)
VCOV$se
# FIX: get standard errors for taus

# Tur <- Turbofuns:::PolychoricRM(as.matrix(mooc), estimate.acm = TRUE)

#### lyule ####

# Yule correlations

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

fit <- lyule(data = as.matrix(mooc), do.fit = TRUE)
fit@parameters
fit@Optim$SE$se

# Compare with lpoly
fit <- lpoly(data = mooc,
             method = "one-step",
             positive = FALSE,
             penalties = FALSE,
             do.fit = TRUE)
fit@parameters$S
fit@Optim$SE$se

#### lrotate ####

# Factor rotation

library(latent)

set.seed(2026)

# Simulate data:
nfactors <- 3L
nitems <- 12
sim <- simfactor(nfactors = 3, nitems = nitems/nfactors,
                 correlations = 0.40, crossloadings = 0.30)
scores <- MASS::mvrnorm(1e3, rep(0, nrow(sim$R)), Sigma = sim$R)
s <- cor(scores)

estimator <- "uls"
rotation <- "oblimin"
projection <- "poblq"
# Fit efa with bifactor:
fit <- bifactor::efast(s, nfactors = nfactors, estimator = estimator,
                       rotation = rotation, projection = projection,
                       oblq_factors = c(2),
                       gamma = 0, random_starts = 10L, cores = 1L)
fit$rotation$f
rownames(fit$efa$lambda) <- paste("X", 1:nitems, sep = "")
colnames(fit$efa$lambda) <- paste("F", 1:nfactors, sep = "")
lambda <- list(fit$efa$lambda)

set.seed(2028)
p <- nitems
q <- nfactors
target <- diag(nfactors) %x% rep(1, nitems/nfactors)
weight <- 1-target
psitarget <- matrix(1, q, q)
psitarget[1:2, 1:2] <- 0
psiweight <- 1-psitarget; diag(weight) <- 0
constraints <- matrix(0, nfactors, nfactors)
constraints[1:2, 1:2] <- 1; diag(constraints) <- 0
fit <- lrotate(lambda, rotation = rotation, projection = projection,
               target = target, weight = weight,
               psitarget = psitarget, psiweight = psiweight, w = 1,
               gamma = 0.00, epsilon = 0.01, k = 0, a = 1, b = 0.5,
               constraints = constraints,
               do.fit = TRUE, control = list(opt = "newton", rstarts = 100L))
fit@Optim$f
fit@Optim$iterations
fit@Optim$convergence
fit@Optim$ng
fit@Optim$elapsed

#### lefa ####

# Exploratory factor analysis

library(latent)

set.seed(2026)

# Simulate data:
nfactors <- 3L
nitems <- 12
sim <- simfactor(nfactors = 3, nitems = nitems/nfactors,
                 correlations = 0.40, crossloadings = 0.30)
scores <- MASS::mvrnorm(1e3, rep(0, nrow(sim$R)), Sigma = sim$R)
s <- cor(scores)

estimator <- "uls"
rotation <- "oblimin"
projection <- "poblq"
# Fit efa with bifactor:
fit <- bifactor::efast(s, nfactors = nfactors, estimator = estimator,
                       rotation = rotation, projection = projection,
                       oblq_factors = c(2),
                       gamma = 0, random_starts = 10L, cores = 1L)
fit$rotation$f

fit <- lefa(data = scores, #sample.cov = s,
            nfactors = nfactors, estimator = estimator,
            rotation = rotation, projection = projection)



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
# SE for ML-modal-prop
