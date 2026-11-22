library(latent)

df <- HolzingerSwineford1939[, paste("x", 1:9, sep = "")]

nfactors <- 3L
estimator <- "ml"
std.ov <- TRUE
std.lv <- TRUE
meanstructure <- FALSE
likelihood <- "normal"
rotation <- "oblimin"
projection <- "oblq"
constraints <- matrix(0, nfactors, nfactors)
constraints[1:2, 1:2] <- 1
diag(constraints) <- 0

set.seed(10)
fit <- lefa(data = df,
            nfactors = nfactors,
            estimator = estimator,
            std.ov = std.ov,
            std.lv = std.lv,
            meanstructure = meanstructure,
            likelihood = likelihood,
            orthogonal = TRUE,
            rotation = rotation,
            projection = projection,
            # constraints = constraints,
            # oblique = list(2:3),
            control.efa = list(rstarts = 10L, se_method = "KKT",
                               orth.lambda = TRUE),
            control.rotation = list(rstarts = 10L, se_method = "KKT",
                                    opt = "newton", step = "trust"),
            se = TRUE)
# fit@Optim$SE$table_se
# fit@modelInfo$param$lambda
latInspect(fit, what = "est")
latInspect(fit, what = "convergence")
fit@Optim$iterations
fit@Optim$elapsed

# Use the same ML covariance matrix. efast() would standardize raw data to unit variance.
S <- fit@extra$efa@dataList$sample.cov[[1L]]
fit_b <- bifactor::efast(S, nfactors = nfactors, nobs = nrow(df), estimator = estimator,
                         rotation = rotation, projection = projection,
                         gamma = 0, random_starts = 10L, cores = 1L)
fit_b$elapsed/1e9
fit@extra$efa@Optim$f
fit_b$efa$f

fit@Optim$f
fit_b$rotation$f

rotation_differences(fit, fit_b$rotation$lambda, fit_b$rotation$phi)

# EFA constraints derivatives:
parameters <- fit@extra$efa@modelInfo$param$lambda
dconstr <- constraints_derivs(fit@extra$efa, parameters = parameters)
dconstr$dconstr
dim(dconstr$dconstr)
dim(dconstr$d2constr)

# Constraints in rotated parameters:
gradient <- get_grad(fit)$grad
names(gradient) <- fit@modelInfo$transparameters_labels
L <- fit@transformed_pars$lambda_rotated
P <- fit@transformed_pars$psi_rotated
L_labels <- fit@modelInfo$trans$lambda_rotated
P_labels <- fit@modelInfo$trans$psi_rotated
GL <- matrix(gradient[L_labels], nrow = nrow(L), ncol = ncol(L))
GP <- matrix(gradient[P_labels], nrow = nrow(P), ncol = ncol(P))
stationarity <- 2*GP - solve(P, crossprod(GL, L))
stopifnot(max(abs(stationarity[row(P) != col(P)])) < 1e-6)

# Rotation constraints derivatives:
parameters <- fit@modelInfo$trans$X
dconstr <- constraints_derivs(fit, parameters = parameters)
dconstr$dconstr
dim(dconstr$dconstr)
dim(dconstr$d2constr)
