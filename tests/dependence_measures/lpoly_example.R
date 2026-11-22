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
# The one-step sandwich covariance includes the thresholds.
VCOV <- vcov(fit, v = fit@Optim$SE$VCOV, parameters = taus)
VCOV$se
stopifnot(all.equal(VCOV$se, fit@Optim$SE$se[unlist(taus, use.names = FALSE)]))

# Tur <- Turbofuns::PolychoricRM(as.matrix(mooc), estimate.acm = TRUE)
