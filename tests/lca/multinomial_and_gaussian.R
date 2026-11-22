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
