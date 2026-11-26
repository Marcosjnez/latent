library(latent)

# Input order interleaves IDs and times. Subject C has no transition.
d <- data.frame(id = c("B", "C", "A", "B", "A", "B"),
                time = c(5, 2, 8, 1, 3, 4),
                y = c(0.7, 0.2, -1, -0.9, 1.2, 0.9))
set.seed(25)
fit <- lca(d, nclasses = 2L, gaussian = "y", time = "time", id = "id",
           control = list(rstarts = 1L, cores = 1L, maxit = 10L))
inspect <- function(what, sort = FALSE) latInspect(fit, what, sort = sort)
initial <- inspect("initial")
forward <- inspect("forward")
backward <- inspect("backward")
smoothed <- inspect("smoothed")
path <- inspect("viterbi")
modal <- inspect("states")
moves <- inspect("posterior_transition")
A <- inspect("transition")
stopifnot(is.matrix(A), identical(forward, inspect("filtered")),
          identical(initial, inspect("fullclasses")), identical(smoothed, inspect("posterior")))

for(result in list(initial, forward, backward, smoothed, path, modal, moves)) {
  stopifnot(identical(names(result), unique(d$id)))
}
sorted <- inspect("posterior_transition", TRUE)
stopifnot(identical(moves, latInspect(fit, "posterior_transition", sort = FALSE, polish = FALSE)),
          identical(A, latInspect(fit, "transition", sort = FALSE, polish = FALSE)))
states <- colnames(A)
for(id in unique(d$id)) {
  obs <- d[d$id == id, ]
  obs <- obs[order(obs$time), ]
  times <- as.character(obs$time)
  post <- smoothed[[id]]
  stopifnot(identical(dim(initial[[id]]), c(1L, 2L)),
            identical(rownames(initial[[id]]), id),
            identical(rownames(forward[[id]]), times), identical(rownames(backward[[id]]), times),
            identical(rownames(post), times), identical(names(path[[id]]), times),
            identical(names(modal[[id]]), times), all(modal[[id]] == max.col(post, ties.method = "first")),
            max(abs(forward[[id]] * backward[[id]] - post)) < 1e-10,
            length(moves[[id]]) == nrow(obs)-1L)

  # Sum the probabilities of all possible paths to verify each move and state.
  n <- nrow(obs)
  paths <- as.matrix(expand.grid(rep(list(1:2), n)))
  emission <- exp(fit@transformed_pars$loglik[fit@dataList$data$id == id, , drop = FALSE])
  prob <- apply(paths, 1L, function(z) {
    initial[[id]][1L, z[1L]] * prod(emission[cbind(seq_len(n), z)]) *
      prod(A[cbind(z[-1L], z[-n])])
  })
  prob <- prob/sum(prob)
  stopifnot(all(path[[id]] == paths[which.max(prob), ]),
            max(abs(post[, 1L] - colSums((paths == 1L) * prob))) < 1e-10)
  if(n == 1L) next
  stopifnot(identical(names(moves[[id]]), times[-1L]))
  for(t in 2:n) {
    move <- moves[[id]][[times[t]]]
    expected <- outer(1:2, 1:2, Vectorize(function(i, j) {
      sum(prob[paths[, t] == i & paths[, t-1L] == j])
    }))
    stopifnot(identical(rownames(move), states), identical(colnames(move), states),
              max(abs(move - expected)) < 1e-10,
              max(abs(sorted[[id]][[times[t]]][states, states] - move)) < 1e-10)
  }
}

# Omitting id creates one trajectory.
single <- d[d$id == "B", ]
fit_one <- lca(single, nclasses = 2L, gaussian = "y", time = "time",
               control = list(rstarts = 1L, cores = 1L, maxit = 1L))
stopifnot(identical(names(latInspect(fit_one, "forward")), "1"))

# Polishing adds a printed explanation while preserving ID/time indexing.
polished <- latInspect(fit, "posterior_transition", polish = TRUE)
stopifnot(identical(names(polished), names(sorted)),
          identical(polished[["B"]][["4"]], sorted[["B"]][["4"]]),
          identical(polished[["C"]], sorted[["C"]]))
display <- capture.output(print(polished))
for(label in c("ID = B", "Joint state probabilities from time 1 to time 4",
               "Joint state probabilities from time 4 to time 5",
               "One observation. No transitions.", "To state", "From state")) {
  stopifnot(any(grepl(label, display, fixed = TRUE)))
}

shared <- latInspect(fit, "transition", sort = FALSE, polish = TRUE)
stopifnot(is.matrix(shared), identical(as.numeric(shared), as.numeric(A)),
          identical(dimnames(shared), dimnames(A)),
          any(grepl("Transition probabilities shared by all trajectories",
                    capture.output(print(shared)), fixed = TRUE)),
          identical(smoothed, latInspect(fit, "smoothed", sort = FALSE, polish = TRUE)))

fits <- structure(list(two_states = fit), class = "llcalist")
collection <- latInspect(fits, "posterior_transition", polish = TRUE)
stopifnot(identical(collection[["two_states"]], polished))
