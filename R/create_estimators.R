# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 20/03/2026

get_estimators <- function(estimators, structures) {

  # Collect the unique parameter labels that are not fixed values:
  vector_structures <- unname(unique(c(unlist(structures))))

  # Handle the estimators:
  estimators_and_labels <- estimators

  # Check that a list was provided:
  if(!is.list(estimators_and_labels)) {
    stop("estimators should be a list")
  }
  nestimators <- length(estimators_and_labels) # Number of estimators
  # Pick the estimators:
  estimators <- lapply(estimators_and_labels, FUN = \(x) x$estimator)
  # Pick the parameter labels going to each estimator:
  inputs <- lapply(estimators_and_labels, FUN = \(x) x$parameters)
  # Pick the extra objects going to each estimator:
  dots <- lapply(estimators_and_labels, FUN = \(x) x$extra)

  # For each estimator:
  result <- vector("list", length = nestimators)
  k <- 1L
  for(i in 1:nestimators) {

    estimator <- estimators[[i]]

    # Choose which extra objects from 'dots' should be kept for this estimator:
    est_objects <- switch(estimator,
                          beta_loglik        = c("X"),
                          binomial_loglik    = c("X", "Ntrials"),
                          exponential_loglik = c("X"),
                          gamma_loglik       = c("X"),
                          gaussian_loglik    = c("alpha", "means", "N", "sds"),
                          laplace_loglik     = c("X"),
                          poisson_loglik     = c("X"),
                          t_loglik           = c("X"),
                          weibull_loglik     = c("X"),
                          cf                 = c("p", "q", "k"),
                          geomin             = c("p", "q", "epsilon"),
                          lclf               = c("p", "q", "epsilon"),
                          oblimin            = c("p", "q", "gamma"),
                          target             = c("target", "weight"),
                          varimax            = c("p", "q"),
                          varimin            = c("p", "q"),
                          xtarget            = c("target", "weight", "w",
                                                 "psitarget", "psiweight"),
                          ridge              = c("lambda", "power", "N"),
                          lreg               = c("y", "X"),
                          polycor            = c("p", "n", "N"),
                          lca                = c("S", "J", "I", "weights"),
                          bayesconst1        = c("K", "alpha", "N", "U"),
                          bayesconst2        = c("K", "alpha", "N", "pihat"),
                          bayesconst3        = c("K", "alpha", "N", "varshat"),
                          logdetmat          = c("lower_indices", "logdetw", "p"),
                          cfa_dwls           = c("q", "w", "R", "W"),
                          cfa_ml             = c("R", "w", "n"),
                          cfa_fml            = c("R", "w", "n"),
                          cfa_dwls_error     = c("w", "p", "W"),
                          cfa_fml_error      = c("w", "p", "n"),
                          cfa_ml_error       = c("w", "p", "n"),
                          stop("Unknown estimator: ", estimator)
    )

    # Pick only those objects from 'dots':
    extra <- dots[[i]][est_objects]

    # Ensure the required extras are present and named:
    if (length(est_objects) > 0L) {
      missing <- setdiff(est_objects, names(extra))
      if (length(missing) > 0L) {
        stop("Missing required object(s) for estimator '", estimator, "': ",
             paste(missing, collapse = ", "))
      }
    }

    ninputs <- length(inputs[[i]]) # Number of estimators
    indices <- labels <- vector("list", length = ninputs)
    for(j in 1:ninputs) {

      # Collect the unique subset of parameter labels that are not fixed values:
      if(is.list(inputs[[i]][j])) {
        labels_vector <- unname(c(unlist(inputs[[i]][[j]])))
      } else {
        labels_vector <- unname(c(unlist(structures[inputs[[i]][[j]]])))
      }

      # Get the indices of the vector_structures that are in the labels_vector:
      m <- match(labels_vector, vector_structures)
      if (anyNA(m)) { # Check for wrong parameter labels
        stop("Some parameters were not found in structures: ",
             paste(labels_vector[is.na(m)], collapse = ", "))
      }

      indices[[j]] <- m-1L # C++ indexing starts at 0
      labels[[j]] <- labels_vector

    }

    result[[i]] <- c(
      list(estimator = estimator,
           indices  = indices,
           labels = labels),
      extra
    )

  }

  return(result)

}
