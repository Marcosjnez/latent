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
