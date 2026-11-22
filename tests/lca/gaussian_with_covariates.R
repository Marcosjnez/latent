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
