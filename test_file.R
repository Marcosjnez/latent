#### LCA (multinomial) ####

library(latent)

data <- gss82
item <- rep("multinomial", ncol(data))
nclasses <- 3L
control <- list(opt = "lbfgs", maxit = 1000, rstarts = 16L, cores = 1L,
                eps = 1e-05)
penalties <- list(
  class = list(alpha = 1),
  prob  = list(alpha = 1),
  sd    = list(alpha = 0)
)
# nmiss <- 30
# missrow <- sample(1:nrow(data), size = nmiss)
# misscol <- sample(1:ncol(data), size = nmiss, replace = TRUE)
# for(i in 1:nmiss) data[missrow[i], misscol[i]] <- NA
fit <- lca(data = data, item = item, nclasses = nclasses,
           penalties = penalties, control = NULL, do.fit = TRUE)
fit@timing
fit@loglik # -2754.643
fit@penalized_loglik # -2759.507
fit@Optim$opt$iterations
fit@transformed_pars

SE <- se(fit, type = "robust", model = "user", digits = 4)
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
control <- list(opt = "lbfgs", maxit = 1000L, rstarts = 32L, cores = 32L,
                eps = 1e-05)
penalties <- list(
  class = list(alpha = 1),
  prob  = list(alpha = 0),
  sd    = list(alpha = 1)
)
# nmiss <- 30
# missrow <- sample(1:nrow(data), size = nmiss)
# misscol <- sample(1:ncol(data), size = nmiss, replace = TRUE)
# for(i in 1:nmiss) data[missrow[i], misscol[i]] <- NA
fit <- lca(data = data, item = item, nclasses = nclasses,
           penalties = penalties, control = control, do.fit = TRUE)
fit@timing
fit@loglik # -1841.336
fit@penalized_loglik # -1844.333
fit@Optim$opt$iterations
fit@transformed_pars

SE <- se(fit, type = "standard", model = "user", digits = 2)
SE$table
SE$table_se

CI <- ci(fit, type = "standard", model = "user",
         confidence = 0.95, digits = 2)
CI$table
CI$lower_table
CI$upper_table

#### Mixed LCA (multinomial and gaussian) ####

library(latent)
data <- cancer[, 1:6]
names(data)
item <- c("gaussian", "gaussian", "multinomial", "multinomial",
          "gaussian", "gaussian")
nclasses <- 3L
control <- list(opt = "lbfgs", maxit = 1000L, rstarts = 32L, cores = 32L,
                eps = 1e-05, step_maxit = 30)
penalties <- list(
  class = list(alpha = 1),
  prob  = list(alpha = 1),
  sd    = list(alpha = 1)
)
# nmiss <- 30
# missrow <- sample(1:nrow(data), size = nmiss)
# misscol <- sample(1:ncol(data), size = nmiss, replace = TRUE)
# for(i in 1:nmiss) data[missrow[i], misscol[i]] <- NA
fit <- lca(data = data, item = item, nclasses = nclasses,
           penalties = penalties, control = control, do.fit = TRUE)
fit@timing
fit@loglik # -5784.701
fit@penalized_loglik # -5795.573
fit@Optim$opt$iterations
fit@transformed_pars

SE <- se(fit, type = "robust", model = "user", digits = 4)
SE$table
SE$table_se

CI <- ci(fit, type = "standard", model = "user",
         confidence = 0.95, digits = 2)
CI$table
CI$lower_table
CI$upper_table
