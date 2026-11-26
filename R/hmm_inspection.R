# HMM inspection in trajectory order and in the fitted state order.

.print_hmm_transitions <- function(x, ...) {
  print_matrix <- function(m) {
    names(dimnames(m)) <- c("To state", "From state")
    print.default(m, ...)
  }
  if(is.matrix(x)) {
    cat("Transition probabilities shared by all trajectories\n")
    print_matrix(unclass(x))
  } else {
    first_time <- attr(x, "first_time")
    for(q in seq_along(x)) {
      if(q > 1L) cat("\n")
      cat("ID = ", names(x)[q], "\n", sep = "")
      moves <- x[[q]]
      if(!length(moves)) cat("One observation. No transitions.\n")
      from <- first_time[q]
      for(j in seq_along(moves)) {
        to <- names(moves)[j]
        cat("Joint state probabilities from time ", from, " to time ", to, "\n", sep = "")
        print_matrix(moves[[j]])
        from <- to
      }
    }
  }
  invisible(x)
}

.hmm_by_id <- function(value, d) {
  times <- as.character(d$data[[d$time]])
  rows <- split(seq_along(times), rep.int(seq_along(d$lengths), d$lengths))
  names(rows) <- as.character(d$trajectory_ids)
  if(is.matrix(value)) {
    rownames(value) <- times
    return(lapply(rows, function(i) value[i, , drop = FALSE]))
  }
  names(value) <- times
  lapply(rows, function(i) value[i])
}

.hmm_log_row_sums <- function(x) {
  # log(rowSums(exp(x))), allowing an impossible transition or emission.
  shift <- x[cbind(seq_len(nrow(x)), max.col(x, ties.method = "first"))]
  shift[shift == -Inf] <- 0
  shift + log(rowSums(exp(x - shift)))
}

.hmm_viterbi <- function(log_emission, log_transition, initial, first, last) {
  nclasses <- ncol(log_emission)
  states <- seq_len(nclasses)
  path <- integer(nrow(log_emission))
  predecessor <- matrix(0L, max(last - first + 1L), nclasses)

  for(q in seq_along(first)) {
    start <- first[q]
    end <- last[q]
    score <- log(initial[q, ]) + log_emission[start, ]
    if(end > start) {
      for(s in seq.int(start + 1L, end)) {
        # Rows are next states; columns are previous states.
        joint <- sweep(log_transition, 2L, score - max(score), "+")
        previous <- max.col(joint, ties.method = "first")
        predecessor[s - start + 1L, ] <- previous
        score <- joint[cbind(states, previous)] + log_emission[s, ]
      }
    }
    path[end] <- which.max(score)
    if(end > start) {
      for(s in seq.int(end, start + 1L)) {
        path[s - 1L] <- predecessor[s - start + 1L, path[s]]
      }
    }
  }
  path
}

.hmm_quantities <- function(fit, what, with_smoothed = FALSE) {
  parameters <- fit@transformed_pars
  log_emission <- parameters$loglik
  if(what != "viterbi") {
    filtered <- fit@Optim$outputs$estimators$matrices[[1L]][[1L]]
    dim(filtered) <- dim(log_emission)
    if(what == "forward") return(filtered)
  }
  # Common emission factors cancel from posterior probabilities and paths.
  peak <- max.col(log_emission, ties.method = "first")
  log_emission <- log_emission - log_emission[cbind(seq_len(nrow(log_emission)), peak)]
  log_transition <- log(parameters$transition)
  nclasses <- ncol(log_emission)
  first <- fit@dataList$initial_rows
  last <- cumsum(fit@dataList$lengths)
  if(what == "viterbi") {
    return(.hmm_viterbi(log_emission, log_transition, parameters$class, first, last))
  }

  by_previous <- t(log_transition)
  need_smoothed <- with_smoothed || what == "smoothed"
  if(need_smoothed) smoothed <- filtered # At trajectory ends the two coincide.
  if(what == "posterior_transition") {
    result <- setNames(vector("list", length(first)), as.character(fit@dataList$trajectory_ids))
    times <- as.character(fit@dataList$data[[fit@dataList$time]])
  } else if(what == "backward") {
    result <- matrix(1, nrow(log_emission), nclasses)
  }

  for(q in seq_along(first)) {
    start <- first[q]
    end <- last[q]
    if(what == "posterior_transition") {
      destinations <- start + seq_len(end - start)
      result[[q]] <- setNames(vector("list", length(destinations)), times[destinations])
    }
    if(end == start) next

    log_filtered <- log(filtered[start:end, , drop = FALSE])
    # Recover log probabilities if conversion to the stored probability scale
    # lost precision. Later observations can still favor such a state.
    if(any(log_filtered < log(.Machine$double.xmin))) {
      log_prob <- log(parameters$class[q, ]) + log_emission[start, ]
      for(s in seq.int(start, end)) {
        if(s > start) {
          joint <- sweep(log_transition, 2L, log_prob, "+")
          log_prob <- .hmm_log_row_sums(joint) + log_emission[s, ]
        }
        log_prob <- log_prob - max(log_prob)
        log_prob <- log_prob - log(sum(exp(log_prob)))
        log_filtered[s - start + 1L, ] <- log_prob
      }
      if(need_smoothed) smoothed[end, ] <- exp(log_prob)
    }

    log_backward <- numeric(nclasses)
    for(s in seq.int(end, start + 1L)) {
      future <- log_emission[s, ] + log_backward
      joint <- sweep(by_previous, 2L, future - max(future), "+")
      log_backward <- .hmm_log_row_sums(joint)
      log_forward <- log_filtered[s - start, ]
      log_state <- log_forward + log_backward
      shift <- max(log_state)
      normalizer <- shift + log(sum(exp(log_state - shift)))

      if(need_smoothed) smoothed[s - 1L, ] <- exp(log_state - normalizer)
      if(what == "backward") {
        # P(future observations | current state) / P(future | past and present).
        # The final value is one; multiplying by the filter gives the smoother.
        result[s - 1L, ] <- exp(log_backward - normalizer)
      } else if(what == "posterior_transition") {
        result[[q]][[s - start]] <- t(exp(joint + log_forward - normalizer))
      }
    }
  }
  if(what == "smoothed") return(smoothed)
  if(with_smoothed) return(list(value = result, smoothed = smoothed))
  result
}
