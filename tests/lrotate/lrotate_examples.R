library(latent)

# Compare rotated loadings and factor correlations after matching columns and signs.
rotation_differences <- function(fit, lambda, phi) {
  L <- fit@transformed_pars$lambda_rotated
  P <- fit@transformed_pars$psi_rotated
  congruence <- psych::factor.congruence(L, lambda, digits = 10)
  index <- max.col(abs(congruence), ties.method = "first")
  stopifnot(!anyDuplicated(index))
  signs <- sign(congruence[cbind(seq_along(index), index)])
  errors <- c(loadings = max(abs(L - sweep(lambda[, index, drop = FALSE], 2, signs, "*"))),
              correlations = max(abs(P - phi[index, index, drop = FALSE]*outer(signs, signs))))
  stopifnot(max(errors) < 2e-4)
  errors
}

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9'

estimator <- "ml"
std.ov <- TRUE
std.lv <- TRUE
meanstructure <- FALSE
likelihood <- "normal"

fit_cfa <- lcfa(data = HolzingerSwineford1939,
                model = model,
                estimator = estimator,
                std.ov = std.ov,
                std.lv = std.lv,
                meanstructure = meanstructure,
                likelihood = likelihood,
                orthogonal = TRUE)
fit_rotation <- lrotate(fit = fit_cfa,
                        projection = "oblq",
                        rotation = "oblimin",
                        control = list(se_method = "KKT"),
                        se = TRUE)
fit_rotation@Optim$f
fit_rotation@parameters
fit_rotation@Optim$SE
fit_rotation@Optim$SE$table_se
# fit_rotation@modelInfo$param

# Factor rotation

library(latent)

set.seed(2026)

# Simulate data:
nfactors <- 3L
nitems <- 12
sim <- simfactor(nfactors = 3, nitems = nitems/nfactors,
                 correlations = 0.40, crossloadings = 0.30)
scores <- MASS::mvrnorm(1e3, rep(0, nrow(sim$R)), Sigma = sim$R)
s <- cor(scores)

estimator <- "uls"
rotation <- "oblimin"
projection <- "poblq"
# Fit efa with bifactor:
fit_b <- bifactor::efast(s, nfactors = nfactors, estimator = estimator,
                         rotation = rotation, projection = projection,
                         oblq_factors = 2L,
                         gamma = 0, random_starts = 10L, cores = 1L)
lambda <- fit_b$efa$lambda
rownames(lambda) <- paste0("X", seq_len(nitems))
colnames(lambda) <- paste0("F", seq_len(nfactors))

set.seed(2028)
fit <- lrotate(lambda = lambda, rotation = rotation, projection = projection,
               oblique = 2L, gamma = 0,
               do.fit = TRUE, control = list(opt = "newton", rstarts = 100L))
c(latent = fit@Optim$f, bifactor = fit_b$rotation$f)
rotation_differences(fit, fit_b$rotation$lambda, fit_b$rotation$phi)
fit@Optim$iterations
fit@Optim$convergence
fit@Optim$ng
fit@Optim$elapsed
