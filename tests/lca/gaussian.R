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
