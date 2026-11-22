library(latent)

# Marginal thresholds use every observed response. Each correlation uses
# the jointly observed responses and is invariant to ordinal score spacing.
x <- c(-3, -3, 0.5, 0.5, 9, 9, -3, 0.5, 9, NA, 9, NA)
y <- c(1, 7, 1, 7, 7, 1, 1, NA, 7, 7, NA, 1)
ordinal <- cbind(x, y)
fit <- polyfast(ordinal)
expected <- as.numeric(stats::qnorm(cumsum(table(x))[1:2]/sum(!is.na(x))))
stopifnot(isTRUE(all.equal(fit$thresholds[[1]][2:3], expected, tolerance = 1e-7)),
          identical(unname(do.call(rbind, fit$contingency_tables[[1]])),
                    matrix(as.integer(table(x, y)), 3, 2)))
recoded <- apply(ordinal, 2L, function(z) match(z, sort(unique(z))))
stopifnot(isTRUE(all.equal(fit$correlation, polyfast(recoded)$correlation)),
          polyfast(cbind(x, -x))$correlation[2L, 1L] < -0.99)

yule <- lyule(data.frame(ordinal))
complete <- as.data.frame(ordinal[complete.cases(ordinal), ])
complete[] <- lapply(complete, factor)
reference <- lyule(complete)
stopifnot(isTRUE(all.equal(yule@parameters$S, reference@parameters$S)),
          isTRUE(all.equal(yule@Optim$SE$se, reference@Optim$SE$se)))

# Gaussian covariance and correlation uncertainty have closed forms.
S <- matrix(c(2, 1, 1, 4), 2L)
stopifnot(asymptotic_normal(S)[2L, 2L] == S[1L, 1L]*S[2L, 2L] + S[2L, 1L]^2,
          isTRUE(all.equal(as.numeric(asymptotic_normal(S, FALSE, FALSE)),
                           (1 - stats::cov2cor(S)[2L, 1L]^2)^2)))

# Unequal numbers of observed pairs enter point estimates and standard errors.
data <- data.frame(x, y)
fit <- lpearson(data)
counts <- crossprod(!is.na(ordinal))
S <- stats::cov(ordinal, use = "pairwise.complete.obs")*(counts-1L)/counts
stopifnot(isTRUE(all.equal(fit@parameters$S, S)),
          isTRUE(all.equal(unname(fit@Optim$SE$VCOV[2L, 2L]),
                           (S[1L, 1L]*S[2L, 2L] + S[2L, 1L]^2)/counts[2L, 1L])),
          lpearson(data, missing = "complete.obs")@dataList$nobs == sum(complete.cases(data)))
