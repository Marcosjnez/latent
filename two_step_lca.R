library(latent)

data <- empathy[, 1:6]
# Covariates:
X <- empathy[, 7:8]
# Covariates may also include factors

# Step one, fit the measurement model without the covariates:
fit0 <- lca(data = data, X = NULL, model = NULL,
            item = rep("gaussian", ncol(data)),
            nclasses = 4L, penalties = TRUE,
            control = list(maxit = 10000))

# Step two, fit the model with covariates fixing the measurement part:
fit <- lca(data = data, X = X, model = fit0,
           item = rep("gaussian", ncol(data)),
           nclasses = 4L, penalties = TRUE,
           control = list(maxit = 10000))
fit@timing
fit@loglik # -1798.885
fit@penalized_loglik # -1801.996
fit@Optim$opt$iterations
fit@Optim$opt$convergence

# Plot model fit info:
fit

# Get fit indices:
getfit(fit)

# Inspect model objects:
latInspect(fit, what = "coefs", digits = 5)
latInspect(fit, what = "classes", digits = 5)
latInspect(fit, what = "profile", digits = 3)
latInspect(fit, what = "posterior", digits = 3)
