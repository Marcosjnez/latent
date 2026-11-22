# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 19/11/2026

#### Store a dataset ####

# usethis::use_data(dataset, overwrite = TRUE)

#### Build the package ####

# roxygen2::roxygenise()

#### Install latent ####

# pak::pak("marcosjnez/latent")
# devtools::install_github("marcosjnez/latent", force = TRUE)

# Install a previous version:
# remotes::install_github("marcosjnez/latent", ref = "6e3603e")

#### LCA (multinomial) ####

library(latent)
gss <- latent::gss82
gss$EDUCR <- as.integer(gss$EDUCR)-1L
set.seed(20216)
fit <- lca(data = gss,
           nclasses = 3L,
           multinomial = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"),
           # covariates = c("RACE", "SEX", "EDUCR", "AGE"),
           # outcomes = "MARITAL",
           # adjustment = "bk",
           # classification = "modal",
           # model = list("UNDERSTA ~~ COOPERAT
           #               PURPOSE ~~ COOPERAT"),
           # start = start,
           penalties = list(class = list(alpha=1),
                            prob  = list(alpha=0)),
           # control = list(opt = "lbfgs", step = "armijo",
           #                rstarts = 30, cores = 30),
           control = list(opt = "em", maxit = 1000L, eps = 1e-10,
                          mopt = "lbfgs", mstep_maxit = 50L, mstep_eps = 1e-05,
                          step = "armijo", step_maxit = 30L,
                          rstarts = 30, cores = 30),
           do.fit = TRUE)
latInspect(fit, what = "loglik")
latInspect(fit, what = "convergence")
latInspect(fit, what = "elapsed")
# loglik: -3891.252 # penalized_loglik: -3892.478
# With both residual dependencies and direct L-BFGS:
# loglik: -3879.167 # penalized_loglik: -3880.371

jacob <- jacobian(fit, parameters = fit@modelInfo$trans[names(fit@modelInfo$param)])
dim(jacob)
labs1 <- c(fit@modelInfo$trans$PURPOSE)
labs2 <- c(fit@modelInfo$trans$PURPOSE)
jacob[labs1, labs2]

dconstr <- constraints_derivs(fit, parameters = fit@modelInfo$trans[names(fit@modelInfo$param)])
dim(dconstr$dconstr)
labs1 <- c(fit@modelInfo$trans$ACCURACY)
dconstr$dconstr[labs1, , drop = FALSE]

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
plot(fit, type = "coefficients", what = "OR", se_type = "information")

# With the covariates above, predict class membership for two respondents:
if(length(fit@dataList$covariates_names)) {
  df <- data.frame(RACE = c("BLACK", "WHITE"),
                   SEX = c("MALE", "MALE"),
                   EDUCR = c(3, 3),
                   AGE = c(35, 35))
  predict(fit, new = df)
} else {
  head(predict(fit)) # Common class probabilities in the intercept-only model
}

# Inspect model objects:
latInspect(fit, what = "convergence")
latInspect(fit, what = "profile")
latInspect(fit, what = "coefs")
latInspect(fit, what = "classification")
latInspect(fit, what = "pattern")
latInspect(fit, what = "table")

# Get standard errors:
SE <- se(fit, type = "information", digits = 4,
         parameters = fit@modelInfo$trans[names(fit@modelInfo$param)])
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
           # adjustment = "bk",
           penalties = TRUE,
           control = list(rstarts = 30L, cores = 30L),
           # control = list(opt = "em", maxit = 50L, eps = 1e-05,
           #                mopt = "grad", mstep_maxit = 1L, mstep_eps = 1e-05,
           #                step = "armijo", step_maxit = 30L,
           #                rstarts = 30, cores = 30),
           do.fit = TRUE)
latInspect(fit, what = "loglik")
latInspect(fit, what = "convergence")
latInspect(fit, what = "elapsed")
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
CI <- ci(fit, type = "information", confidence = 0.95, digits = 2)
CI$table
stopifnot(isTRUE(all.equal(CI$upper-CI$lower, 2*qnorm(0.975)*SE$se)))

#### Mixed LCA (multinomial and gaussian) ####

library(latent)
set.seed(2026)
penalties <- list(
  beta  = list(alpha = 0),
  class = list(alpha = 1),
  prob  = list(alpha = 1),
  var   = list(alpha = 1),
  Sigma = list(alpha = 1)
) # The same defaults are used with penalties = TRUE.
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

# Wald test for equality of the Class 2 and Class 3 intercepts:
labels <- fit@modelInfo$trans$beta[1, 2:3]
contrast <- c(1, -1)
delta <- sum(contrast*fit@Optim$transparameters[labels])
wald <- delta^2/drop(contrast %*% SE$vcov[labels, labels] %*% contrast)
pchisq(wald, df = 1, lower.tail = FALSE)

plot(fit)

#### LCA with covariates (gaussian) ####

library(latent)
set.seed(2026)

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
latInspect(fit, what = "loglik")
# loglik: -1747.135 # penalized_loglik: -1750.566
# loglik: -2049.840 # penalized_loglik: -2053.322 (outcomes = list(gaussian = c("pt5")))
latInspect(fit, what = "convergence")

# Print model fit:
fit

# Get fit indices:
getfit(fit)

plot(fit)
plot(fit, type = "coefficients", what = "OR", se_type = "information")

# Inspect model objects:
latInspect(fit, what = "loglik")
latInspect(fit, what = "coefs")
latInspect(fit, what = "classes")
latInspect(fit, what = "profile")
latInspect(fit, what = "posterior")

predict(fit, new = rbind(c(2, 2, 2.428571, 2.142857),
                         c(1, 2, 3, 4)))
fitted(fit)

# Get standard errors:
SE <- se(fit, type = "information", digits = 4)
SE$table

# Effects coding and changing the reference to Class 2:
beta <- fit@transformed_pars$beta
labels <- fit@modelInfo$trans$beta
effects_coding(beta, SE$vcov[c(labels), c(labels)])$table_se
free <- c(labels[, -1L, drop = FALSE])
move_intercept(beta, SE$vcov[free, free], column = 2L)$table_se

# Get confidence intervals:
CI <- ci(fit, type = "information", confidence = 0.95, digits = 2)
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
estimator <- "uls"
std.ov <- FALSE
std.lv <- FALSE
meanstructure <- TRUE
likelihood <- "normal"

fit <- lcfa(model = model,
            data = HolzingerSwineford1939,
            # sample.cov = S,
            # sample.mean = means,
            # sample.nobs = 301,
            estimator = estimator,
            std.ov = std.ov,
            std.lv = std.lv,
            meanstructure = meanstructure,
            likelihood = likelihood,
            se = "standard",
            # control = list(opt = "lbfgs", step = "armijo"),
            do.fit = TRUE)
latInspect(fit, what = "convergence")
fit@Optim$elapsed
fit@modelInfo$param

# fit@modelInfo$param

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
# The default ULS criteria use different covariance divisors and weights.
# latent uses S * (n-1)/n and counts both triangles in its squared residuals.
latent_loss <- latInspect(fit, "loss")["loss", "overall"]
c(latent = latent_loss, lavaan_twice_fmin = 2 * unname(fitMeasures(fit2, "fmin")))

# Supply the same sample moments and weights for a direct comparison.
# lavaan stacks the means followed by the unique entries of the covariance.
n <- nrow(HolzingerSwineford1939)
W <- matrix(1, nrow(S), ncol(S))
diag(W) <- 0.5
WLS.V <- diag(c(rep(0.5, length(means)), W[lower.tri(W, diag = TRUE)]))
fit2_matched <- lavaan::cfa(model, sample.cov = S * (n-1)/n, sample.mean = means,
                            sample.nobs = n, sample.cov.rescale = FALSE,
                            ordered = FALSE, estimator = "DWLS", WLS.V = WLS.V,
                            std.lv = std.lv, std.ov = std.ov,
                            meanstructure = meanstructure, se = "none")
# lavaan's least-squares test statistic uses n-1 observations.
lavaan_loss <- lavaan::lavInspect(fit2_matched, "test")$standard$stat/(n-1)
c(latent = latent_loss, lavaan_matched = lavaan_loss)
stopifnot(abs(latent_loss - lavaan_loss) < 1e-7)

lavaan::lavInspect(fit2, "est")
latInspect(fit, "est", sort = FALSE)
fit@Optim$SE$se
lavaan::lavInspect(fit2, "se")
fitMeasures(fit2, c("baseline.chisq", "rmr"))
getfit(fit)

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

fit <- lcfa(HolzingerSwineford1939,
            model = model,
            group = "school",
            estimator = estimator,
            std.ov = std.ov,
            std.lv = std.lv,
            meanstructure = meanstructure,
            likelihood = likelihood,
            se = "standard",
            control = NULL,
            do.fit = TRUE)

latInspect(fit, "est")
# fit@modelInfo$param

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
loss <- latInspect(fit, "loss")["loss", "overall"]
ll <- latInspect(fit, "loglik")
c(latent = loss, lavaan = 2 * unname(fitMeasures(fit2, "fmin")))
c(latent = ll["loglik", "overall"], lavaan = as.numeric(logLik(fit2)))
stopifnot(abs(loss - 2 * fitMeasures(fit2, "fmin")) < 1e-7,
          abs(ll["loglik", "overall"] - as.numeric(logLik(fit2))) < 1e-5,
          abs(ll["loglik_sat", "overall"] - fitMeasures(fit2, "unrestricted.logl")) < 1e-5)

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
std.lv <- TRUE
meanstructure <- TRUE
likelihood <- "normal"

set.seed(2026)
fit <- lcfa(data = HolzingerSwineford1939,
            model = model,
            estimator = estimator,
            meanstructure = meanstructure,
            ordered = FALSE,
            positive = TRUE,
            std.lv = std.lv,
            std.ov = std.ov,
            likelihood = likelihood,
            # penalties = TRUE,
            penalties = list(logdet = list(w = 0.001)),
            se = "standard",
            do.fit = TRUE,
            control = list(rstarts = 10L, se_method = "KKT", opt = "newton"))

# The positivity constraint and log-determinant penalty change the optimum.
latInspect(fit, "loglik")
# Different Gram factors give the same covariance matrix. The information
# calculation can warn about variances of their individual coordinates.
# Inspect the reported model-parameter SEs separately from those coordinates.
X <- se(fit, parameters = fit@modelInfo$trans[c("lambda", "psi", "theta")])
lapply(X$table_se, FUN = round, digits = 3L)

latInspect(fit, what = "est")

# fit@Optim$SE$se
det(fit@transformed_pars$theta)
det(fit@transformed_pars$psi)

# With lavaan:
fit2 <- cfa(model, data = HolzingerSwineford1939,
            estimator = estimator, std.lv = std.lv, std.ov = std.ov,
            meanstructure = meanstructure, likelihood = likelihood)
fitMeasures(fit2, "logl")
fitMeasures(fit2, "unrestricted.logl")
fitMeasures(fit2, "baseline.chisq")

diag(lavaan::lavInspect(fit2, "est")$theta) # The unconstrained fit has a negative variance.
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

# Obtain starting values from the unrestricted covariance fit. lavaan
# warns about a negative residual variance in this example.
fit2 <- lavaan::cfa(model, data = HolzingerSwineford1939,
                    group = "school", estimator = estimator,
                    std.lv = std.lv, std.ov = std.ov,
                    meanstructure = meanstructure, likelihood = likelihood)
est <- lavaan::lavInspect(fit2, "est")
start <- list()
for(g in names(est)) {
  theta <- est[[g]]$theta
  smallest <- min(eigen(theta, symmetric = TRUE, only.values = TRUE)$values)
  diag(theta) <- diag(theta) + max(0, 0.1 - smallest)
  start[[paste0("xtheta.", g)]] <- chol(theta)
  start[[paste0("xpsi.", g)]] <- chol(est[[g]]$psi)
  start[[paste0("lambda.", g)]] <- est[[g]]$lambda
}

set.seed(2026)
fit <- lcfa(data = HolzingerSwineford1939,
            model = model,
            group = "school",
            estimator = estimator,
            ordered = FALSE,
            std.ov = std.ov,
            std.lv = std.lv,
            positive = TRUE,
            penalties = TRUE,
            meanstructure = meanstructure,
            likelihood = likelihood,
            se = "standard",
            do.fit = TRUE,
            control = list(opt = "newton", step = "trust", tcg_maxit = 100L,
                           start = start, rstarts = 1L, se_method = "RH"))
stopifnot(fit@Optim$convergence)
# Expected log likelihood is about -3674.300. Redundant Gram coordinates
# can produce variance warnings even at a converged fit.
latInspect(fit, "loglik")

fit@Optim$iterations
fit@Optim$convergence
fit@timing

# Compare the unconstrained fit:
fitmeasures(fit2, fit.measures = c("cfi", "tli", "rmsea", "srmr"))
fitMeasures(fit2, "logl")
fitMeasures(fit2, "unrestricted.logl")
fitMeasures(fit2, "baseline.chisq")

diag(lavaan::lavInspect(fit2, "est")$Pasteur$theta) # Unconstrained residual variances
det(inspect(fit2, what = "est")$Pasteur$theta)

#### CFA (polychorics) ####

library(lavaan)
library(latent)

samples <- unique(hexaco$sample) # industry mooc fire student dutch
Ns <- sapply(samples, FUN = function(x) sum(hexaco$sample == x))
names(Ns) <- samples

# Subset the items pertaining to the HEXACO-100
selection <- 5:104
full <- hexaco[, selection]
mooc <- full[hexaco$sample == "mooc", ]
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

fit <- lcfa(data = mooc,
            model = model.EM,
            ordered = TRUE,
            estimator = estimator,
            std.ov = std.ov,
            std.lv = std.lv,
            meanstructure = meanstructure,
            se = "standard",
            control = NULL,
            do.fit = TRUE)

# With lavaan:
fit2 <- lavaan::cfa(data = mooc,
                    model = model.EM,
                    ordered = TRUE,
                    estimator = estimator,
                    meanstructure = meanstructure,
                    std.ov = std.ov,
                    std.lv = std.lv,
                    parameterization = "delta",
                    do.fit = TRUE)

# The two-step polychoric weight estimates differ between the packages.
latent_loss <- latInspect(fit, "loss")["loss", "overall"]
n <- nrow(mooc)
c(latent = latent_loss, lavaan_default = lavaan::lavInspect(fit2, "test")$standard$stat/(n-1))

# latent sums both halves of the correlation matrix and half-weights the
# threshold residuals. Use those weights in lavaan's unique-statistic order.
est <- fit@modelInfo$control_estimator[[1L]]
WLS.V <- diag(c(est$w_thresholds/2, est$W[lower.tri(est$W)]))
fit2_matched <- lavaan::cfa(model.EM, data = mooc, ordered = TRUE,
                            estimator = estimator, meanstructure = meanstructure,
                            std.lv = std.lv, std.ov = std.ov,
                            parameterization = "delta", WLS.V = WLS.V)
lavaan_loss <- lavaan::lavInspect(fit2_matched, "test")$standard$stat/(n-1)
c(latent = latent_loss, lavaan_matched = lavaan_loss)
stopifnot(abs(latent_loss - lavaan_loss) < 1e-6)

diag(lavaan::lavInspect(fit2_matched, "est")$theta)
diag(latInspect(fit, "est", sort = FALSE)$theta)

#### CFA (Yule correlation) ####

library(latent)
library(lavaan)

samples <- unique(hexaco$sample) # industry mooc fire student dutch
Ns <- sapply(samples, FUN = function(x) sum(hexaco$sample == x))
names(Ns) <- samples

# Subset the items pertaining to the HEXACO-100
selection <- 5:104
full <- hexaco[, selection]

mooc <- full[hexaco$sample == "mooc", ]
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

fit <- lcfa(model = model.EM,
            data = mooc,
            ordered = "yule",
            estimator = estimator,
            positive = FALSE,
            penalties = FALSE,
            std.ov = std.ov,
            std.lv = std.lv,
            meanstructure = meanstructure,
            se = "standard",
            do.fit = TRUE,
            control = NULL)
latInspect(fit, "loss")
fit@Optim$iterations
fit@Optim$convergence
fit@timing

latInspect(fit, what = "lambda")
# This Yule example uses an approximate inverse for the first-stage Hessian
# when propagating the sampling uncertainty into the factor parameters.
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

mooc <- full[hexaco$sample == "mooc", ]
dim(mooc)
mooc$school <- rep(c("s1", "s2"), times = c(2000L, nrow(mooc) - 2000L))

model.EM <- "FEA =~ hexemfea146 + hexemfea170 + hexemfea74 + hexemfea2
             ANX =~ hexemanx128 + hexemanx8 + hexemanx80 + hexemanx176
             DEP =~ hexemdep62 + hexemdep182 + hexemdep134 + hexemdep158
             SEN =~ hexemsen44 + hexemsen164 + hexemsen20 + hexemsen68"

fit <- lcfa(model = model.EM,
            data = mooc,
            ordered = TRUE,
            estimator = "dwls",
            group = "school",
            se = "standard",
            do.fit = TRUE,
            control = NULL)
loss <- latInspect(fit, "loss")["loss", ]

# With lavaan, use the same marker indicators and delta parameterization.
fit2 <- lavaan::cfa(model = model.EM, data = mooc, ordered = TRUE,
                    estimator = "dwls", group = "school", std.lv = FALSE,
                    meanstructure = TRUE, parameterization = "delta")
n <- as.numeric(lavaan::lavInspect(fit2, "nobs"))
default_loss <- lavaan::lavInspect(fit2, "test")$standard$stat.group/(n-1)
c(latent = unname(loss["overall"]), lavaan_default = weighted.mean(default_loss, n))

WLS.V <- lapply(fit@modelInfo$control_estimator, function(est) {
  diag(c(est$w_thresholds/2, est$W[lower.tri(est$W)]))
})
fit2_matched <- lavaan::cfa(model.EM, data = mooc, ordered = TRUE,
                            estimator = "dwls", group = "school", std.lv = FALSE,
                            meanstructure = TRUE, parameterization = "delta",
                            WLS.V = WLS.V)
group_loss <- lavaan::lavInspect(fit2_matched, "test")$standard$stat.group/(n-1)
# latent reports each group's contribution weighted by its sample fraction.
comparison <- cbind(latent = loss[-length(loss)], lavaan_matched = group_loss*n/sum(n))
comparison
stopifnot(max(abs(comparison[, 1L] - comparison[, 2L])) < 1e-6)

# Inspect sampling uncertainty using the group-specific parameter structures.
blocks <- with(fit@dataList$data_param, c(lambda_group, theta_group, psi_group))
latent_se <- se(fit, parameters = fit@modelInfo$trans[blocks], digits = NULL)$table_se
lavaan_se <- lavaan::lavInspect(fit2, "se")
latent_se[c("lambda.s1", "theta.s1", "psi.s1")]
lavaan_se$s1[c("lambda", "theta", "psi")]

# Compare correlation weights separately from thresholds, in the same pair order.
W <- lavaan::lavInspect(fit2, "wls.v")
lapply(seq_along(W), function(g) {
  native <- fit@modelInfo$control_estimator[[g]]$W
  weights <- diag(W[[g]])
  pairs <- grepl("~~", names(weights), fixed = TRUE)
  cbind(latent = native[lower.tri(native)], lavaan = weights[pairs])
})

#### CFA (FIML) ####

library(latent)
library(lavaan)

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9'

data_missing <- HolzingerSwineford1939
data_missing[1:5, "x1"] <- NA
data_missing[9:20, "x5"] <- NA
data_missing[167:189, "x6"] <- NA
data_missing[244:251, "x9"] <- NA
data_missing[280, "x8"] <- NA
# rvars <- sample(paste("x", c(1, 4, 6), sep = ""), size = nmiss, replace = TRUE)
# rsubjs <- sample(1:1000, size = nmiss)
# for(j in 1:nmiss) data_missing[rsubjs[j], rvars[j]] <- NA

estimator <- "ml"
std.ov <- FALSE
std.lv <- FALSE
meanstructure <- TRUE
likelihood <- "normal"

set.seed(2026)
fit <- lcfa(data = data_missing,
            model = model,
            estimator = estimator,
            missing = "ml",
            group = "school",
            meanstructure = meanstructure,
            std.lv = std.lv,
            std.ov = std.ov,
            likelihood = likelihood,
            # missing = "pairwise.complete.obs",
            positive = FALSE,
            se = "standard",
            do.fit = TRUE)

# fit@modelInfo$param

# With lavaan:
fit2 <- cfa(model = model,
            data = data_missing,
            estimator = estimator,
            missing = "fiml",
            group = "school",
            meanstructure = meanstructure,
            std.lv = std.lv,
            std.ov = std.ov,
            likelihood = likelihood)

# The direct pattern likelihood has a different additive constant in its
# optimization loss. Compare the observed-data log likelihood on the same cases.
ll <- latInspect(fit, "loglik")["loglik", "overall"]
c(latent = ll, lavaan = as.numeric(logLik(fit2)))
stopifnot(abs(ll - as.numeric(logLik(fit2))) < 1e-5)

latInspect(fit, "est", sort = FALSE)
lavInspect(fit2, "est")

lavInspect(fit2, "se")
fit2@ParTable$se
blocks <- with(fit@dataList$data_param, c(lambda_group, theta_group, psi_group))
x <- se(fit, parameters = fit@modelInfo$trans[blocks])
x$table_se

#### lmean ####

# Estimate means of variables
library(latent)

df <- HolzingerSwineford1939[, paste("x", 1:9, sep = "")]
set.seed(2026)
std.ov <- FALSE

fit <- lmean(data = df, model = NULL, std.ov = std.ov)
fit@parameters$means
fit@Optim$SE$se
stopifnot(all.equal(c(fit@parameters$means), unname(colMeans(df))),
          all.equal(unname(fit@Optim$SE$se), sqrt(unname(diag(cov(df)))/nrow(df))))

#### lpearson ####

# Pearson covariances (normal likelihood uses denominator N)

library(latent)

df <- HolzingerSwineford1939[, paste("x", 1:9, sep = "")]
std.ov <- FALSE
likelihood <- "normal"
VCOV_type <- "standard"
fit <- lpearson(data = df, std.ov = std.ov,
                VCOV = VCOV_type, likelihood = likelihood,
                missing = "pairwise.complete.obs",
                do.fit = TRUE)
fit@parameters
fit@Optim$SE$se
stopifnot(all.equal(unname(fit@parameters$S), unname(cov(df))*(nrow(df)-1L)/nrow(df)))

#### lpoly ####

# Polychoric correlations

library(latent)
samples <- unique(hexaco$sample) # industry mooc fire student dutch
Ns <- sapply(samples, FUN = function(x) sum(hexaco$sample == x))
names(Ns) <- samples

# First six HEXACO-100 items in the MOOC sample
selection <- 5:10
full <- hexaco[, selection]

mooc <- full[hexaco$sample == "mooc", ]
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
             # start = taus, # Optional thresholds computed above
             do.fit = TRUE,
             control = list(opt = "grad",
                            # subfix = ".group1",
                            maxit = 500,
                            step_maxit = 50,
                            tcg_maxit = 30,
                            ss = 0.001,
                            eps = 1e-06))
# The composite log-likelihood is summed over cases; the objective is per case.
loglik <- fit@Optim$outputs$estimators$doubles[[1L]][2L] # -176520.8
stopifnot(all.equal(fit@Optim$f, -loglik/nrow(mooc)))
fit@Optim$iterations
fit@Optim$ng
fit@Optim$convergence
fit@Optim$f # 41.18545
max(abs(fit@Optim$rg))
max(abs(fit@Optim$g))
fit@timing
fit@modelInfo$param
fit@modelInfo$parameters_labels

fit@parameters
fit@Optim$SE$se
taus <- fit@modelInfo$trans[paste0("taus", names(mooc))]
# lpoly() stores the two-step ACOV, including the thresholds.
VCOV <- vcov(fit, v = fit@Optim$SE$VCOV, parameters = taus)
VCOV$se
stopifnot(all.equal(VCOV$se, fit@Optim$SE$se[unlist(taus, use.names = FALSE)]))

# Tur <- Turbofuns::PolychoricRM(as.matrix(mooc), estimate.acm = TRUE)

#### lyule ####

# Yule correlations

library(latent)

samples <- unique(hexaco$sample) # industry mooc fire student dutch
Ns <- sapply(samples, FUN = function(x) sum(hexaco$sample == x))
names(Ns) <- samples

# First six HEXACO-100 items in the MOOC sample
selection <- 5:10
full <- hexaco[, selection]

mooc <- full[hexaco$sample == "mooc", ]
dim(mooc)

fit_yule <- lyule(data = as.matrix(mooc), do.fit = TRUE)
fit_yule@parameters
fit_yule@Optim$SE$se

# Compare with lpoly
fit <- lpoly(data = mooc,
             method = "one-step",
             positive = FALSE,
             penalties = FALSE,
             do.fit = TRUE)
fit@parameters$S
fit@Optim$SE$se

#### lrotate ####

library(latent)

# Compare rotated loadings and factor correlations after matching columns and signs.
rotation_differences <- function(fit, lambda, phi) {
  L <- fit@transformed_pars$lambda_rotated
  P <- fit@transformed_pars$psi_rotated
  congruence <- psych::factor.congruence(L, lambda, digits = 10)
  index <- max.col(abs(congruence), ties.method = "first")
  stopifnot(!anyDuplicated(index))
  signs <- sign(congruence[cbind(seq_along(index), index)])
  errors <- c(loadings = max(abs(L - sweep(lambda[, index, drop = FALSE], 2, signs, "*"))),
              correlations = max(abs(P - phi[index, index, drop = FALSE]*outer(signs, signs))))
  stopifnot(max(errors) < 2e-4)
  errors
}

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9'

estimator <- "ml"
std.ov <- TRUE
std.lv <- TRUE
meanstructure <- FALSE
likelihood <- "normal"

fit_cfa <- lcfa(data = HolzingerSwineford1939,
                model = model,
                estimator = estimator,
                std.ov = std.ov,
                std.lv = std.lv,
                meanstructure = meanstructure,
                likelihood = likelihood,
                orthogonal = TRUE)
fit_rotation <- lrotate(fit = fit_cfa,
                        projection = "oblq",
                        rotation = "oblimin",
                        control = list(se_method = "KKT"),
                        se = TRUE)
fit_rotation@Optim$f
fit_rotation@parameters
fit_rotation@Optim$SE
fit_rotation@Optim$SE$table_se
# fit_rotation@modelInfo$param

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
fit_b <- bifactor::efast(s, nfactors = nfactors, estimator = estimator,
                         rotation = rotation, projection = projection,
                         oblq_factors = 2L,
                         gamma = 0, random_starts = 10L, cores = 1L)
lambda <- fit_b$efa$lambda
rownames(lambda) <- paste0("X", seq_len(nitems))
colnames(lambda) <- paste0("F", seq_len(nfactors))

set.seed(2028)
fit <- lrotate(lambda = lambda, rotation = rotation, projection = projection,
               oblique = 2L, gamma = 0,
               do.fit = TRUE, control = list(opt = "newton", rstarts = 100L))
c(latent = fit@Optim$f, bifactor = fit_b$rotation$f)
rotation_differences(fit, fit_b$rotation$lambda, fit_b$rotation$phi)
fit@Optim$iterations
fit@Optim$convergence
fit@Optim$ng
fit@Optim$elapsed

#### lefa ####

library(latent)

df <- HolzingerSwineford1939[, paste("x", 1:9, sep = "")]

nfactors <- 3L
estimator <- "ml"
std.ov <- TRUE
std.lv <- TRUE
meanstructure <- FALSE
likelihood <- "normal"
rotation <- "oblimin"
projection <- "oblq"
constraints <- matrix(0, nfactors, nfactors)
constraints[1:2, 1:2] <- 1
diag(constraints) <- 0

set.seed(10)
fit <- lefa(data = df,
            nfactors = nfactors,
            estimator = estimator,
            std.ov = std.ov,
            std.lv = std.lv,
            meanstructure = meanstructure,
            likelihood = likelihood,
            orthogonal = TRUE,
            rotation = rotation,
            projection = projection,
            # constraints = constraints,
            # oblique = list(2:3),
            control.efa = list(rstarts = 10L, se_method = "KKT",
                               orth.lambda = TRUE),
            control.rotation = list(rstarts = 10L, se_method = "KKT",
                                    opt = "newton", step = "trust"),
            se = TRUE)
# fit@Optim$SE$table_se
# fit@modelInfo$param$lambda
latInspect(fit, what = "est")
latInspect(fit, what = "convergence")
fit@Optim$iterations
fit@Optim$elapsed

# Use the same ML covariance matrix. efast() would standardize raw data to unit variance.
S <- fit@extra$efa@dataList$sample.cov[[1L]]
fit_b <- bifactor::efast(S, nfactors = nfactors, nobs = nrow(df), estimator = estimator,
                         rotation = rotation, projection = projection,
                         gamma = 0, random_starts = 10L, cores = 1L)
fit_b$elapsed/1e9
fit@extra$efa@Optim$f
fit_b$efa$f

fit@Optim$f
fit_b$rotation$f

rotation_differences(fit, fit_b$rotation$lambda, fit_b$rotation$phi)

# EFA constraints derivatives:
parameters <- fit@extra$efa@modelInfo$param$lambda
dconstr <- constraints_derivs(fit@extra$efa, parameters = parameters)
dconstr$dconstr
dim(dconstr$dconstr)
dim(dconstr$d2constr)

# Constraints in rotated parameters:
gradient <- get_grad(fit)$grad
names(gradient) <- fit@modelInfo$transparameters_labels
L <- fit@transformed_pars$lambda_rotated
P <- fit@transformed_pars$psi_rotated
L_labels <- fit@modelInfo$trans$lambda_rotated
P_labels <- fit@modelInfo$trans$psi_rotated
GL <- matrix(gradient[L_labels], nrow = nrow(L), ncol = ncol(L))
GP <- matrix(gradient[P_labels], nrow = nrow(P), ncol = ncol(P))
stationarity <- 2*GP - solve(P, crossprod(GL, L))
stopifnot(max(abs(stationarity[row(P) != col(P)])) < 1e-6)

# Rotation constraints derivatives:
parameters <- fit@modelInfo$trans$X
dconstr <- constraints_derivs(fit, parameters = parameters)
dconstr$dconstr
dim(dconstr$dconstr)
dim(dconstr$d2constr)

#### Exploratory bifactor ####

library(latent)
library(GPArotation)

set.seed(2026)
ulambda <- rorth(9, 4)

# Bi-geomin:

rot <- bigeominQ(ulambda, delta = 0.01, randomStarts = 10L)
fit <- lrotate(lambda = ulambda, projection = "oblq",
               rotation = list(geomin = list(factors = 2:4, epsilon = 0.01)),
               control = list(rstarts = 10L))
fit@Optim$iterations
fit@Optim$f
rotation_differences(fit, rot$loadings, rot$Phi)

# Bi-quartimin:

rot <- bifactorQ(ulambda, randomStarts = 10L)
fit <- lrotate(lambda = ulambda, projection = "oblq",
               rotation = list(oblimin = list(factors = 2:4)),
               control = list(rstarts = 10L))
fit@Optim$iterations
fit@Optim$f
rotation_differences(fit, rot$loadings, rot$Phi)

#### SEM ####

library(latent)
library(lavaan)

model <- '
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + a*y6 + b*y7 + c*y8

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60

  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'
set.seed(2026)
estimator <- "ml"
std.ov <- FALSE
std.lv <- FALSE
meanstructure <- TRUE
likelihood <- "normal"

fit <- lcfa(model = model,
            data = PoliticalDemocracy,
            estimator = estimator,
            std.ov = std.ov,
            std.lv = std.lv,
            meanstructure = meanstructure,
            likelihood = likelihood,
            se = FALSE,
            do.fit = TRUE)
fit@modelInfo$param
latInspect(fit, what = "convergence")
fit@Optim$elapsed
fit@transformed_pars$B
fit@modelInfo$control_optimizer$deltaparam
fit@parameters

fit2 <- sem(model, data = PoliticalDemocracy,
            estimator = estimator,
            std.ov = std.ov,
            std.lv = std.lv,
            meanstructure = meanstructure,
            likelihood = likelihood)
# Compare the ML discrepancy and the log likelihoods of the same three models.
reference <- lavaan::fitMeasures(fit2, c("fmin", "logl", "unrestricted.logl",
                                         "baseline.chisq"))
ll <- latent::latInspect(fit, "loglik")[, "overall"]
comparison <- rbind(
  latent = c(latent::latInspect(fit, "loss")["loss", "overall"],
             ll["loglik"], ll["loglik_base"], ll["loglik_sat"]),
  lavaan = c(2*reference["fmin"], reference["logl"],
             reference["unrestricted.logl"] - reference["baseline.chisq"]/2,
             reference["unrestricted.logl"]))
colnames(comparison) <- c("loss", "loglik", "loglik_base", "loglik_sat")
comparison
stopifnot(max(abs(comparison[1L, ] - comparison[2L, ])) < 1e-5)

lavaan::lavInspect(fit2, "est")$beta
round(latent::latInspect(fit, "est", sort = FALSE)$B, 3)
fit2

#### latscores ####

library(latent)
fit <- lcfa(HolzingerSwineford1939[, paste("x", 1:9, sep = "")],
            model = 'visual =~ x1 + x2 + x3
                     textual =~ x4 + x5 + x6
                     speed =~ x7 + x8 + x9')
fs <- latscores(fit, method = "Green")
fs@parameters$scores
fs@Optim$SE$table_se$scores
latInspect(fs, "scores")               # Display order and signs.
se(fs)$table_se                        # Native order and signs.
vcov(fs, person = 1)                   # Small error covariance for person 1.
fs_new <- latscores(fit,
                    data = HolzingerSwineford1939[, paste("x", 1:9, sep = "")],
                    method = "regression")

#### Hidden Markov Model ####

library(latent)

# Data: 120 consecutive observations from one participant
data("speed", package = "depmixS4")
dat <- head(speed, 120)
dat$id <- "A"
dat$time <- seq_len(nrow(dat))

# 1. latent
set.seed(123)
fit_latent <- lca(
  data = dat,
  nclasses = 2,
  gaussian = "rt",
  time = "time",
  penalties = FALSE,
  control = list(opt = "lbfgs", rstarts = 10)
)

latInspect(fit_latent, what = "loglik")
latInspect(fit_latent, what = "convergence")
latInspect(fit_latent, what = "elapsed")

initial <- latInspect(fit_latent, "initial", sort = FALSE)
forward <- latInspect(fit_latent, "forward", sort = FALSE)
backward <- latInspect(fit_latent, "backward", sort = FALSE)
smoothed <- latInspect(fit_latent, "smoothed", sort = FALSE)
path <- latInspect(fit_latent, "viterbi")
# Columns are previous states and rows are next states:
posterior_transition <- latInspect(fit_latent, "posterior_transition")
stopifnot(is(fit_latent, "llca"),
          max(abs(forward * backward - smoothed)) < 1e-10,
          max(abs(rowSums(smoothed) - 1)) < 1e-10,
          length(path) == nrow(dat))

set.seed(123)
fit_em <- lca(
  data = dat,
  nclasses = 2,
  gaussian = "rt",
  id = "id",
  time = "time",
  penalties = FALSE,
  control = list(opt = "em", mopt = "lbfgs", rstarts = 10,
                 eps = 1e-10, mstep_eps = 1e-7)
)

# 2. depmixS4
mod <- depmixS4::depmix(
  rt ~ 1,
  data = dat,
  nstates = 2,
  family = gaussian()
)

set.seed(123)
fit_depmix <- depmixS4::fit(mod, verbose = FALSE)

# 3. Compare log-likelihoods
ll_latent <- latInspect(fit_latent, "loglik")["loglik", "overall"]
ll_em <- latInspect(fit_em, "loglik")["loglik", "overall"]
ll_depmix <- as.numeric(depmixS4::logLik(fit_depmix))

c(
  latent = ll_latent,
  latent_em = ll_em,
  depmixS4 = ll_depmix,
  difference = ll_latent - ll_depmix
)

# 4. Order states by mean response time. Each column is a previous state.
latent_order <- order(fit_latent@transformed_pars$rt["mean", ])
em_order <- order(fit_em@transformed_pars$rt["mean", ])
states <- seq_along(latent_order)
depmix_order <- order(vapply(states, function(k) {
  depmixS4::getpars(depmixS4::getmodel(fit_depmix, "response", k, 1L))[1L]
}, numeric(1L)))
depmix_transition <- vapply(states, function(k) {
  depmixS4::getpars(depmixS4::getmodel(fit_depmix, "transition", k))
}, numeric(length(states)))
A_latent <- latInspect(fit_latent, "transition", sort = FALSE)[latent_order, latent_order]
A_em <- latInspect(fit_em, "transition", sort = FALSE)[em_order, em_order]
A_depmix <- depmix_transition[depmix_order, depmix_order]
A_latent-A_depmix
stopifnot(abs(ll_latent-ll_depmix) < 1e-5, max(abs(A_latent-A_depmix)) < 1e-4,
          abs(ll_em-ll_depmix) < 1e-5, max(abs(A_em-A_depmix)) < 1e-4)

# 5. Compare state probabilities and the decoded path after aligning the labels
filter_depmix <- depmixS4::posterior(fit_depmix, type = "filtering")
smooth_depmix <- depmixS4::posterior(fit_depmix, type = "smoothing")
path_depmix <- depmixS4::posterior(fit_depmix, type = "global")
stopifnot(max(abs(forward[, latent_order] - filter_depmix[, depmix_order])) < 5e-4,
          max(abs(smoothed[, latent_order] - smooth_depmix[, depmix_order])) < 5e-4,
          all(match(path, latent_order) == match(path_depmix, depmix_order)))

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
