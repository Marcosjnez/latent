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
