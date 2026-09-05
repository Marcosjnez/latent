#### EFA ####

library(latent)

# Simulate data:
set.seed(2026)
sim <- simfactor(nfactors = 3, nitems = 3, correlations = 0.40,
                 estimator = "ml", fit = "rmsr", misfit = 0)

nfactors <- 3L
estimator <- "uls"
std.ov <- TRUE
std.lv <- TRUE
meanstructure <- FALSE
likelihood <- "normal"
rotation <- "target"
projection <- "oblq"
target <- ifelse(sim$lambda > 0, 1, 0)

N <- 300L
scores <- MASS::mvrnorm(N, rep(0, nrow(sim$R_error)), Sigma = sim$R_error)
colnames(scores) <- paste("x", 1:9, sep = "")
scores <- as.data.frame(scores)
# scores$school <- rep(c("s1", "s2"), each = 500)

fit <- lefa(data = scores,
            nfactors = nfactors,
            orthogonal = TRUE,
            estimator = estimator,
            std.ov = std.ov,
            std.lv = std.lv,
            meanstructure = meanstructure,
            likelihood = likelihood,
            rotation = rotation,
            projection = projection,
            target = target,
            weight = 1-target,
            # positive = TRUE,
            # penalties = list(logdet = list(w = 0.001)),
            # do.fit = FALSE,
            control = list(rstarts = 10L, se_method = "KKT"),
            rotation.control = list(rstarts = 10L, se_method = "KKT"),
            se = TRUE)

#### KKT vs Riemannian Hessian ####

labels <- fit@modelInfo$parameters_labels

# control_optimizer <- fit@modelInfo$control_optimizer
# control_optimizer$parameters[[1L]] <- fit@Optim$parameters
# control_optimizer$transparameters[[1L]] <- fit@Optim$transparameters

# RH <- latent:::get_rhess(
#   control_manifold = fit@modelInfo$control_manifold,
#   control_transform = fit@modelInfo$control_transform,
#   control_estimator = fit@modelInfo$control_estimator,
#   control_optimizer = control_optimizer
# )
RH <- get_rhess(fit)

P_Riemannian <- RH$P

P_KKT <- latent:::invert_hessian_latent(
  fit = fit,
  H = hessian(fit),
  labels = labels,
  object_name = "rotation Hessian"
)

max(abs(P_KKT-P_Riemannian))
