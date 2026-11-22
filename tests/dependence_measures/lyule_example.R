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
