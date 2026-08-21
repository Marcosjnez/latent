# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 21/08/2026

create_estimators <- function(estimators, structures) {

  # Collect the unique parameter labels that are not fixed values:
  vector_structures <- unname(unique(c(unlist(structures))))

  # Handle the estimators:
  estimators_and_labels <- estimators

  # Check that a list was provided:
  if(!is.list(estimators_and_labels)) {
    stop("estimators should be a list")
  }

  nestimators <- length(estimators_and_labels)
  estimators <- lapply(estimators_and_labels, FUN = \(x) x$estimator)
  inputs <- lapply(estimators_and_labels, FUN = \(x) x$parameters)
  dots <- lapply(estimators_and_labels, FUN = \(x) x$extra)

  result <- vector("list", length = nestimators)

  for(i in seq_len(nestimators)) {

    estimator <- estimators[[i]]

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
                          lca                = c("S", "I", "weights"),
                          lcaEM              = c("S", "I", "weights"),
                          bayesconst1        = c("K", "alpha", "N", "U"),
                          bayesconst2        = c("K", "alpha", "N", "pihat"),
                          bayesconst3        = c("K", "alpha", "N", "varshat"),
                          bayesconst4        = c("K", "alpha", "J", "D"),
                          logdetmat          = c("lower_indices", "logdetw", "p"),
                          logdetR            = c("lower_indices", "logdetw", "p"),
                          cfa_dwls           = c("q", "w", "p", "W"),
                          cfa_means_dwls     = c("q", "w", "p", "W", "w_means"),
                          cfa_ml             = c("p", "w", "n"),
                          cfa_fml            = c("p", "w", "n"),
                          cfa_means_fml      = c("p", "w", "n"),
                          cfa_means_ml       = c("p", "w", "n"),
                          stop("Unknown estimator: ", estimator))

    extra <- dots[[i]][c(est_objects, "double_names", "matrix_names")]

    if(length(est_objects) > 0L) {
      missing <- setdiff(est_objects, names(extra))
      if(length(missing) > 0L) {
        stop("Missing required object(s) for estimator '", estimator, "': ",
             paste(missing, collapse = ", "))
      }
    }

    ninputs <- length(inputs[[i]])
    indices <- vector("list", length = ninputs)

    for(j in seq_len(ninputs)) {

      if(is.list(inputs[[i]][j])) {
        labels_vector <- unname(c(unlist(inputs[[i]][[j]])))
      } else {
        labels_vector <- unname(c(unlist(structures[inputs[[i]][[j]]])))
      }

      m <- match(labels_vector, vector_structures)

      if(anyNA(m)) {
        stop("Some parameters were not found in structures: ",
             paste(labels_vector[is.na(m)], collapse = ", "))
      }

      indices[[j]] <- m-1L

    }

    result[[i]] <- c(list(estimator = estimator,
                          indices = indices),
                     extra)

  }

  #### Result ####

  return(result)

}
