library(latent)

set.seed(216)
control <- list(rstarts = 1L, cores = 1L, maxit = 2L)

# Partially observed dependent indicators contribute their observed marginal.
# These levels also produce two distinct joint cells with the same pasted label.
items <- expand.grid(a = c("A.B", "A"), b = c("C", "B.C"), c = c("x", "y"))
d <- items[rep(seq_len(nrow(items)), 3L), ]
d$a[1:3] <- NA
d$b[2:4] <- NA
d$c[3] <- NA
fit <- lca(d, 2L, multinomial = names(d), model = list("a ~~ b", "b ~~ c"),
           penalties = FALSE, control = control, verbose = FALSE)
joint <- fit@transformed_pars$a.b.c[seq_len(nrow(items)), , drop = FALSE]
priors <- latInspect(fit, "classes", sort = FALSE)
expected <- apply(d, 1L, function(response) {
  keep <- rep(TRUE, nrow(items))
  for(j in which(!is.na(response))) keep <- keep & items[[j]] == response[j]
  log(sum(colSums(joint[keep, , drop = FALSE]) * priors))
})
stopifnot(max(abs(latInspect(fit, "loglik_case") - expected)) < 1e-9,
          is.na(getfit(fit)["L2"]))

# Weights and conditioning strata agree between model fit and inspection.
d <- expand.grid(a = 0:1, b = 0:1, group = 1:3)
d <- d[rep(seq_len(nrow(d)), 2L), ]
w <- rep(c(1, 3, 2), length.out = nrow(d))
calls <- 0L
fit <- lca({calls <- calls + 1L; d}, 2L, multinomial = c("a", "b"),
           covariates = "group", weights = w, adjustment = "none",
           penalties = FALSE, control = control, verbose = FALSE)
tab <- latInspect(fit, "summary")
indices <- getfit(fit, digits = NULL)
expected_totals <- tapply(w, d$group, sum)
expected_L2 <- 2*sum(tab$Observed * log(tab$Observed/tab$Estimated))
post <- latInspect(fit, "posterior", sort = FALSE)
modal <- max.col(post, ties.method = "first")
assigned <- post[cbind(seq_len(nrow(post)), modal)]
diagnosis <- lclass_diag(fit, type = "all", digits = 10)
stopifnot(calls == 1L,
          max(abs(tapply(tab$Estimated, tab$group, sum) - expected_totals)) < 1e-9,
          abs(indices["L2"] - expected_L2) < 1e-9,
          abs(indices["ICL"] + indices["BIC"] - 2*sum(w*log(assigned))) < 1e-9,
          max(abs(diagnosis$Sum.Posterior$count - colSums(post*w))) < 1e-8)

adjusted <- lca(d, 2L, multinomial = c("a", "b"), covariates = "group",
                weights = w, adjustment = "ml", classification = "prop",
                penalties = FALSE, control = control, verbose = FALSE)
stopifnot(abs(sum(adjusted@dataList$pattern_weights) - sum(w)) < 1e-9)

# Literal NA labels and labels containing spaces remain distinct patterns.
d <- data.frame(a = c("a b", "a", "NA", NA, "a b"), b = c("c", "b c", "x", "x", "c"))
p <- latent:::make_patterns(d, seq_len(nrow(d)), any_continuous = FALSE)
stopifnot(identical(p$short2full, c(1L, 2L, 3L, 4L, 1L)),
          identical(p$pattern_weights, c(6, 2, 3, 4)))

# The sandwich score covariance uses case scores evaluated at the penalized fit.
d <- data.frame(y = rnorm(30L))
fit <- lca(d, 1L, gaussian = "y", penalties = TRUE, control = control, verbose = FALSE)
mu <- fit@transformed_pars$y["mean", 1L]
variance <- fit@transformed_pars$y["var", 1L]
scores <- cbind((mu-d$y)/variance, 0.5*(1-(d$y-mu)^2/variance))
expected <- crossprod(scores)*nrow(d)/(nrow(d)-1L)
stopifnot(max(abs(robust(fit)$B - expected)) < 1e-8)
limits <- ci(fit, confidence = 1 - .Machine$double.eps/2)
stopifnot(all(is.finite(limits$upper)), all(is.finite(limits$lower)))

# Compare HMM posteriors and likelihoods to all state paths of each trajectory.
d <- data.frame(id = c(2, 1, 1, 2, 2), time = c(2, 2, 1, 1, 3), a = c(0, 1, 1, 0, 1))
fit <- lca(d, 2L, multinomial = "a", id = "id", time = "time",
           penalties = FALSE, control = control, verbose = FALSE)
initial <- latInspect(fit, "initial", sort = FALSE)
A <- latInspect(fit, "transition", sort = FALSE)
B <- latInspect(fit, "item", sort = FALSE)$a
post <- latInspect(fit, "posterior", sort = FALSE)
likelihood <- latInspect(fit, "loglik_case", sort = FALSE)
for(id in unique(d$id)) {
  rows <- which(d$id == id)
  rows <- rows[order(d$time[rows])]
  times <- seq_along(rows)
  paths <- expand.grid(rep(list(seq_len(ncol(A))), length(times)))
  prob <- apply(paths, 1L, function(state) {
    initial[as.character(id), state[1L]] * prod(A[cbind(state[-1L], head(state, -1L))]) *
      prod(B[cbind(d$a[rows] + 1L, state)])
  })
  expected <- vapply(paths, function(state) sum(prob[state == 1L])/sum(prob), numeric(1L))
  stopifnot(abs(likelihood[as.character(id)] - log(sum(prob))) < 1e-9,
            max(abs(post[rows, 1L] - expected)) < 1e-9)
}

# Prediction uses the contrasts from estimation, even if the R options change.
d <- expand.grid(a = 0:1, b = 0:1, group = LETTERS[1:3])
original <- options(contrasts = c("contr.sum", "contr.poly"))
fit <- lca(d, 2L, multinomial = c("a", "b"), covariates = "group", adjustment = "none",
           penalties = FALSE, control = control, verbose = FALSE)
options(contrasts = c("contr.treatment", "contr.poly"))
rows <- c(2L, 5L)
pred <- predict(fit, newdata = d[rows, "group", drop = FALSE])
stopifnot(max(abs(as.matrix(pred[, c("Class1", "Class2")]) - fitted(fit)[rows, ])) < 1e-10)
options(original)

# Coefficient plots support a single odds ratio and honor the selected penalty.
fit <- lca(gss82, 2L, multinomial = c("PURPOSE", "ACCURACY", "COOPERAT"),
           control = list(rstarts = 1L, cores = 1L, maxit = 50L), verbose = FALSE)
grDevices::pdf(NULL)
odds <- plot(fit, type = "coefficients", effects = "dummy")
fits <- structure(list(fit, fit), class = "llcalist")
criteria <- plot(getfit(fits))
grDevices::dev.off()
stopifnot(length(odds$parameters) == 1L, all(is.finite(odds$se)),
          identical(names(criteria), c("nclasses", "AICp", "BICp")))

# Multistage derivatives reuse the fitted estimates without drawing new starts.
fit <- lca(gss82, 2L, multinomial = c("PURPOSE", "ACCURACY", "COOPERAT"),
           covariates = "SEX", adjustment = "bk",
           control = list(rstarts = 2L, cores = 1L, maxit = 100L), verbose = FALSE)
seed <- .Random.seed
beta_se <- se(fit, parameters = fit@modelInfo$trans["beta"])$se
stopifnot(identical(seed, .Random.seed), all(is.finite(beta_se)))
