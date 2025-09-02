#### LCA (multinomial) ####

library(latent)

data <- gss82
item <- rep("multinomial", ncol(data))
nclasses <- 3L
control <- list(opt = "lbfgs", maxit = 1000, rstarts = 32L,
                cores = 32L, eps = 1e-05)
penalties <- list(
  class = list(alpha = 1),
  prob  = list(alpha = 1),
  sd    = list(alpha = 0)
)
fit <- lca(data = data, item = item, nclasses = nclasses,
           penalties = penalties, control = control, do.fit = TRUE)
fit@timing
fit@loglik # -2754.643 / -2754.639
fit@penalized_loglik # -2759.507 / -2758.241
fit@Optim$opt$iterations
fit@transformed_pars
fit@parameters

control_manifold <- fit@Optim$control_manifold
control_transform <- fit@Optim$control_transform
control_estimator <- fit@Optim$control_estimator
control_optimizer <- fit@Optim$control
control_optimizer$parameters[[1]] <- fit@Optim$opt$parameters
control_optimizer$transparameters[[1]] <- fit@Optim$opt$transparameters
computations <- grad_comp(control_manifold = control_manifold,
                          control_transform = control_transform,
                          control_estimator = control_estimator,
                          control_optimizer = control_optimizer,
                          compute = "all", eps = 1e-04)
computations$f
fit@penalized_loglik
fit@loglik

SE <- se(fit, type = "standard", model = "user", digits = 4)
SE$table
SE$table_se

CI <- ci(fit, type = "standard", model = "user",
         confidence = 0.95, digits = 2)
CI$table
CI$lower_table
CI$upper_table

#### LCA (gaussian) ####

library(latent)

nclasses <- 4L
data <- empathy[, 1:6]
item <- rep("gaussian", ncol(data))
control <- list(opt = "em", maxit = 1000L, rstarts = 32L, cores = 32L,
                eps = 1e-05)
penalties <- list(
  class = list(alpha = 1),
  prob  = list(alpha = 0),
  sd    = list(alpha = 1)
)
fit <- lca(data = data, item = item, nclasses = nclasses,
           penalties = penalties, control = control, do.fit = TRUE)
fit@timing
fit@loglik # -1841.336
fit@penalized_loglik # -1844.333
fit@Optim$opt$iterations
fit@transformed_pars

SE <- se(fit, type = "standard", model = "user", digits = 3)
SE$table
SE$table_se

CI <- ci(fit, type = "standard", model = "user", digits = 2,
         confidence = 0.95)
CI$table
CI$lower_table
CI$upper_table

#### Mixed LCA (multinomial and gaussian) ####

library(latent)
data <- cancer[, 1:6]
names(data)
item <- c("gaussian", "gaussian",
          "multinomial", "multinomial",
          "gaussian", "gaussian")
nclasses <- 3L
control <- list(opt = "lbfgs", maxit = 1000L,
                rstarts = 32L, cores = 32L,
                eps = 1e-05, step_maxit = 30)
penalties <- list(
  class = list(alpha = 1),
  prob  = list(alpha = 1),
  sd    = list(alpha = 1)
)
fit <- lca(data = data, item = item, nclasses = 3,
           penalties = penalties, control = control, do.fit = TRUE)
fit@timing
fit@loglik # -5784.701
fit@penalized_loglik # -5795.573
fit@Optim$opt$iterations
fit@transformed_pars

SE <- se(fit, type = "standard", model = "user", digits = 4)
SE$table
SE$table_se

CI <- ci(fit, type = "standard", model = "user",
         confidence = 0.95, digits = 2)
CI$table
CI$lower_table
CI$upper_table
