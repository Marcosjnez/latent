library(latent)

# Independent Poisson regressions with a common design matrix.
set.seed(295)
X <- data.frame(x = rnorm(80), group = gl(2, 40))
Y <- cbind(y = rpois(80, exp(.3 + .2*X$x + .3*as.numeric(X$group))),
           z = rpois(80, 3))
fit <- pois_reg(Y, X)
models <- lapply(as.data.frame(Y), function(y) glm(y ~ ., X, family = poisson()))
stopifnot(max(abs(fit@parameters$beta - sapply(models, coef))) < 1e-6,
          abs(fit@Optim$f + sum(vapply(models, logLik, numeric(1)))) < 1e-8,
          max(abs(unname(se(fit)$se) -
                  unlist(lapply(models, function(m) sqrt(diag(vcov(m)))), use.names = FALSE))) < 1e-4,
          length(pois_reg(Y, do.fit = FALSE)@Optim) == 0L)

# Check the generated correlation matrix itself against the target discrepancy.
for(method in c("cudeck", "yuan")) {
  for(estimator in c("uls", "ml")) {
    set.seed(101)
    sim <- simfactor(error_method = method, estimator = estimator, fit = "cfi", misfit = .97)
    S <- sim$R_error
    R <- sim$R
    if(estimator == "ml") {
      null <- -as.numeric(determinant(S)$modulus)
      loss <- as.numeric(determinant(R)$modulus) + null + sum(S*solve(R)) - nrow(R)
    } else {
      null <- sum(S[lower.tri(S)]^2)
      loss <- sum((R-S)^2)/2
    }
    stopifnot(abs(loss-sim$delta) < 1e-5, abs(1-loss/null-.97) < 1e-5,
              all(diag(S) == 1), all(diag(chol(S)) > 0))
    residual <- S-R
    if(estimator == "ml") residual <- solve(R) %*% residual %*% solve(R)
    score <- residual %*% sim$lambda %*% sim$Phi
    stopifnot(max(abs(score[sim$lambda != 0])) < 1e-10)
  }
}
stopifnot(ncol(simfactor(lambda = matrix(.6, 4, 1), misfit = "zero")$Phi) == 1L)
sim <- simfactor(ngenerals = 2, nfactors = 3, pure = TRUE, crossloadings = .2, misfit = "zero")
stopifnot(sum(rowSums(sim$lambda != 0) == 1) == 6L)
set.seed(42)
original <- simfactor(misfit = 0, correlations = .25)
set.seed(42)
cross <- simfactor(misfit = 0, correlations = .25, crossloadings = .3)
stopifnot(max(abs(original$uniquenesses-cross$uniquenesses)) < 1e-12)
k <- grid_search(1e-7, diag(c(1, -1e-12)))
stopifnot(abs(k-log1p(k) - 1e-7) < 1e-14)

# Named, nested and empty blocks retain the requested label order.
stopifnot(identical(create_parameters(list()), list()),
          identical(extract_unique_values(list(a = c("v", "v"), b = list(s = "w")),
                                          list(b = list(s = 3), a = c(NA, 5))), c(v = 5, w = 3)))
labels <- matrix(letters[1:4], 2, 2)
stopifnot(identical(create_parameters(list(list(name = "A", type = "matrix",
                                                dim = c(2, 2), labels = labels)))$A, labels))

# Each rotation component reports its contribution to the fitted objective.
rotation <- lrotate(lambda = matrix(rnorm(24), 8, 3), rotation = "varimax",
                    projection = "orth", se = FALSE)
loadings <- rotation@transformed_pars$lambda_rotated
centered <- sweep(loadings^2, 2L, colMeans(loadings^2))
loss <- -sum(centered^2)/4
stopifnot(abs(rotation@Optim$f-loss) < 1e-9,
          abs(rotation@Optim$outputs$estimators$doubles[[1L]][1L]-loss) < 1e-9)

# Missing responses select the matching submatrix of the fitted covariance.
X <- rbind(c(1, NA, 2), c(NA, 3, 1), c(2, NA, 0), c(NA, NA, NA))
L <- matrix(c(.6, .5, .7), 3, 1)
S <- tcrossprod(L) + diag(3)
dataList <- list(ngroups = 1L, centered = list(X), method = "regression", missing = "available",
                 models = list(list(lambda = L, latent_cov = matrix(1), model = S, latent_means = .3)),
                 scores_group = "F", group_label = "group1")
modelInfo <- list(control_optimizer = list(tol = 1e-10), parameters_labels = letters[1:4])
scores <- fit_latscores(dataList, modelInfo)
expected <- apply(X, 1, function(x) {
  obs <- !is.na(x)
  if(!any(obs)) return(.3)
  .3 + drop(x[obs] %*% solve(S[obs, obs], L[obs]))
})
stopifnot(max(abs(scores$scores$F - expected)) < 1e-12,
          identical(scores$Optim$case_pattern[[1]], c(1L, 2L, 1L, 3L)))
