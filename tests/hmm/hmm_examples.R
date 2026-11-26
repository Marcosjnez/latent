library(latent)

# Data: 120 consecutive observations from one participant
data("speed", package = "depmixS4")
dat <- head(speed, 120)
dat$id <- "A"
dat$time <- seq_len(nrow(dat))

# 1. latent
set.seed(123)
fit_latent <- lca(
  data = dat,
  nclasses = 2,
  gaussian = "rt",
  time = "time", id = "id",
  penalties = FALSE,
  control = list(opt = "lbfgs", rstarts = 10)
)

latInspect(fit_latent, what = "loglik")
latInspect(fit_latent, what = "convergence")
latInspect(fit_latent, what = "elapsed")

initial <- latInspect(fit_latent, "initial", sort = FALSE)
forward <- latInspect(fit_latent, "forward", sort = FALSE)
backward <- latInspect(fit_latent, "backward", sort = FALSE)
smoothed <- latInspect(fit_latent, "smoothed", sort = FALSE)
path <- latInspect(fit_latent, "viterbi")
# Columns are previous states and rows are next states:
posterior_transition <- latInspect(fit_latent, "posterior_transition", sort = FALSE)
transition <- latInspect(fit_latent, "transition", sort = FALSE)
stopifnot(is(fit_latent, "llca"),
          max(abs(forward[["A"]] * backward[["A"]] - smoothed[["A"]])) < 1e-10,
          max(abs(rowSums(smoothed[["A"]]) - 1)) < 1e-10,
          length(path[["A"]]) == nrow(dat))

# 2. depmixS4
mod <- depmixS4::depmix(
  rt ~ 1,
  data = dat,
  nstates = 2,
  family = gaussian()
)

set.seed(123)
fit_depmix <- depmixS4::fit(mod, verbose = FALSE)

# 3. Compare log-likelihoods
ll_latent <- latInspect(fit_latent, "loglik")["loglik", "overall"]
ll_depmix <- as.numeric(depmixS4::logLik(fit_depmix))

c(
  latent = ll_latent,
  depmixS4 = ll_depmix,
  difference = ll_latent - ll_depmix
)

# 4. Order states by mean response time. Each column is a previous state.
latent_order <- order(fit_latent@transformed_pars$rt["mean", ])
states <- seq_along(latent_order)
depmix_order <- order(vapply(states, function(k) {
  depmixS4::getpars(depmixS4::getmodel(fit_depmix, "response", k, 1L))[1L]
}, numeric(1L)))
depmix_transition <- vapply(states, function(k) {
  depmixS4::getpars(depmixS4::getmodel(fit_depmix, "transition", k))
}, numeric(length(states)))
A_latent <- transition[latent_order, latent_order]
A_depmix <- depmix_transition[depmix_order, depmix_order]
A_latent-A_depmix
stopifnot(abs(ll_latent-ll_depmix) < 1e-5, max(abs(A_latent-A_depmix)) < 1e-4)

# 5. Compare state probabilities and the decoded path after aligning the labels
filter_depmix <- depmixS4::posterior(fit_depmix, type = "filtering")
smooth_depmix <- depmixS4::posterior(fit_depmix, type = "smoothing")
path_depmix <- depmixS4::posterior(fit_depmix, type = "global")
stopifnot(max(abs(forward[["A"]][, latent_order] - filter_depmix[, depmix_order])) < 5e-4,
          max(abs(smoothed[["A"]][, latent_order] - smooth_depmix[, depmix_order])) < 5e-4,
          all(match(path[["A"]], latent_order) == match(path_depmix, depmix_order)))

fit_hmm <- lca(
  data = datasets::Indometh,
  nclasses = 2L,
  gaussian = "conc",
  time = "time", id = "Subject",
  control = list(opt = "lbfgs")
)

initial <- latInspect(fit_hmm, "initial", sort = FALSE)
forward <- latInspect(fit_hmm, "forward", sort = FALSE)
backward <- latInspect(fit_hmm, "backward", sort = FALSE)
smoothed <- latInspect(fit_hmm, "smoothed", sort = FALSE)
path <- latInspect(fit_hmm, "viterbi")
# Columns are previous states and rows are next states:
posterior_transition <- latInspect(fit_hmm, "posterior_transition",
                                   sort = FALSE, polish = TRUE)
transition <- latInspect(fit_hmm, "transition", sort = FALSE,
                         polish = TRUE)

# Subject 2: state probabilities at each observation, and the joint
# probabilities for the transition from 0.25 to 0.5 hours.
smoothed[["2"]]
move <- posterior_transition[["2"]][["0.5"]]
stopifnot(max(abs(rowSums(move) - smoothed[["2"]]["0.5", ])) < 1e-10,
          max(abs(colSums(move) - smoothed[["2"]]["0.25", ])) < 1e-10)
