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
