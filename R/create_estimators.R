# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 21/08/2026

create_estimators <- function(estimators, structures) {

  #### Check inputs ####

  if(!is.list(estimators)) {
    stop("estimators should be a list")
  }

  vector_structures <- unname(unique(c(unlist(structures))))
  nestimators <- length(estimators)

  estimator_names <- lapply(estimators,
                            FUN = \(x) x$estimator)
  inputs <- lapply(estimators,
                   FUN = \(x) x$parameters)
  dots <- lapply(estimators,
                 FUN = \(x) x$extra)

  metadata_names <- c("group_index", "group_label",
                      "pattern_index", "role",
                      "covariance_name", "means_name")

  experimental_estimators <- c(
    "beta_loglik",
    "binomial_loglik",
    "exponential_loglik",
    "gamma_loglik",
    "laplace_loglik",
    "t_loglik",
    "weibull_loglik"
  )

  #### Create estimator control objects ####

  result <- vector("list", length = nestimators)

  for(i in seq_len(nestimators)) {

    estimator <- estimator_names[[i]]

    if(estimator %in% experimental_estimators) {
      stop("Estimator '", estimator, "' is present only as experimental ",
           "source code and is not registered because its derivative ",
           "interface is incomplete.")
    }

    est_objects <- switch(
      estimator,
      gaussian_loglik    = c("alpha", "means", "N", "sds"),
      poisson_loglik     = c("X"),
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
      stop("Unknown estimator: ", estimator)
    )

    extra_i <- dots[[i]]

    if(is.null(extra_i)) {
      extra_i <- list()
    }

    extra <- extra_i[c(est_objects,
                       "double_names",
                       "matrix_names")]

    if(length(est_objects) > 0L) {

      missing_objects <- setdiff(est_objects, names(extra))

      if(length(missing_objects) > 0L) {
        stop("Missing required object(s) for estimator '",
             estimator, "': ",
             paste(missing_objects, collapse = ", "))
      }

    }

    metadata <- extra_i[
      intersect(metadata_names, names(extra_i))
    ]

    #### Parameter indices ####

    ninputs <- length(inputs[[i]])
    indices <- vector("list", length = ninputs)

    for(j in seq_len(ninputs)) {

      if(is.list(inputs[[i]][j])) {
        labels_vector <- unname(c(unlist(inputs[[i]][[j]])))
      } else {
        labels_vector <- unname(c(unlist(
          structures[inputs[[i]][[j]]]
        )))
      }

      match_indices <- match(labels_vector,
                             vector_structures)

      if(anyNA(match_indices)) {
        stop("Some parameters were not found in structures: ",
             paste(labels_vector[is.na(match_indices)],
                   collapse = ", "))
      }

      indices[[j]] <- match_indices-1L

    }

    result[[i]] <- c(
      list(estimator = estimator,
           indices = indices),
      extra,
      metadata
    )

  }

  #### Result ####

  return(result)

}
