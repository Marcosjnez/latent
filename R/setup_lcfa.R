# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 08/04/2026

create_cfa_datalist <- function(data, model = NULL, cor = "pearson",
                                estimator = "ml", ordered = FALSE,
                                group = NULL, sample.cov = NULL, nobs = NULL,
                                positive = FALSE, penalties = TRUE,
                                missing = "pairwise.complete.obs",
                                std.lv = TRUE, std.ov = FALSE,
                                acov = "standard", message = FALSE,
                                likelihood = NULL, meanstructure = TRUE,
                                args = NULL,
                                ...) {

  cor <- tolower(cor)
  estimator <- tolower(estimator)
  acov <- tolower(acov)
  missing <- tolower(missing)

  if (is.null(group)) {
    ngroups <- 1L
    group_label <- "group"
  } else {
    group_label <- unique(data[[group]])
    ngroups <- length(group_label)
  }

  item_names <- extract_item_names_lavaan(model, ngroups = ngroups)

  get_group_data <- function(i) {
    if (ngroups > 1L) {
      data[data[[group]] == group_label[i], item_names[[i]], drop = FALSE]
    } else {
      data[, item_names[[i]], drop = FALSE]
    }
  }

  unwrap_single <- function(x) {
    if (length(x) == 0L || all(vapply(x, is.null, logical(1)))) {
      return(NULL)
    }
    if (ngroups == 1L) x[[1L]] else x
  }

  normalize_to_list <- function(x, ngroups) {
    if (is.null(x)) {
      return(vector("list", ngroups))
    }
    if (ngroups == 1L && !is.list(x)) {
      return(list(x))
    }
    x
  }

  sample.cov <- normalize_to_list(sample.cov, ngroups)

  if (!is.null(nobs) && ngroups == 1L && !is.list(nobs)) {
    nobs_in <- list(nobs)
  } else {
    nobs_in <- nobs
  }

  X <- vector("list", ngroups)
  correl <- vector("list", ngroups)
  nobs_list <- vector("list", ngroups)
  NACOV <- vector("list", ngroups)
  WLS.V <- vector("list", ngroups)
  thresholds <- vector("list", ngroups)

  if (!all(vapply(sample.cov, is.null, logical(1)))) {

    for (i in seq_len(ngroups)) {
      X[[i]] <- get_group_data(i)

      if (!is.null(nobs_in)) {
        nobs_list[[i]] <- nobs_in[[i]]
      } else {
        nobs_list[[i]] <- nrow(X[[i]])
      }

      rownames(sample.cov[[i]]) <- colnames(sample.cov[[i]]) <- item_names[[i]]

      correl[[i]] <- list(
        S = sample.cov[[i]],
        W = {
          p <- nrow(sample.cov[[i]])
          W <- matrix(1, nrow = p, ncol = p)
          diag(W) <- 1
          W
        },
        nobs = nobs_list[[i]],
        item_names = item_names[[i]]
      )
    }

  } else {

    sample.cov <- vector("list", ngroups)

    for (i in seq_len(ngroups)) {
      X[[i]] <- get_group_data(i)
      nobs_list[[i]] <- nrow(X[[i]])

      correl[[i]] <- lcov(
        data = X[[i]],
        item_names = item_names[[i]],
        cor = cor,
        estimator = estimator,
        acov = acov,
        nobs = nobs_list[[i]],
        missing = missing,
        std.ov = std.ov,
        likelihood = likelihood,
        meanstructure = meanstructure
      )

      correl[[i]]$nobs <- nobs_list[[i]]
      correl[[i]]$item_names <- item_names[[i]]

      sample.cov[[i]] <- correl[[i]][[1]]$S
      rownames(sample.cov[[i]]) <- colnames(sample.cov[[i]]) <-
        correl[[i]][[1]]$item_names

      NACOV[[i]] <- correl[[i]][[1]]$ACOV

      if (!is.null(correl[[i]][[1]]$thresholds)) {
        thresholds[[i]] <- correl[[i]][[1]]$thresholds
      }

      WLS.V[[i]] <- diag(c(correl[[i]][[1]]$W))
    }
  }

  LAV <- lavaan::cfa(
    model = model,
    sample.cov = unwrap_single(sample.cov),
    sample.nobs = unwrap_single(nobs_list),
    group = group,
    NACOV = unwrap_single(NACOV),
    WLS.V = unwrap_single(WLS.V),
    ordered = ordered,
    std.lv = std.lv,
    std.ov = std.ov,
    meanstructure = TRUE,
    do.fit = FALSE,
    warn = FALSE,
    ...
  )

  LAV@Options$positive <- positive

  lavmodel <- LAV@Model
  ngroups <- LAV@Model@ngroups

  item_label <- LAV@Data@ov.names
  nobs_list <- LAV@Data@nobs
  # group_label <- LAV@Data@group.label
  factor_label <- replicate(ngroups, list(LAV@Model@dimNames[[1]][[2]]))

  model_out <- getmodel_fromlavaan(LAV)
  if (ngroups == 1L) {
    model_out <- list(model_out)
  }

  nitems <- as.list(lavmodel@nvar)
  npatterns <- lapply(nitems, function(p) 0.5 * p * (p + 1))
  nfactors <- lapply(model_out, function(x) ncol(x$lambda))

  if (is.null(args)) {
    args <- as.list(match.call(expand.dots = TRUE))[-1]
  }

  data_list <- list()
  data_list$ngroups <- ngroups
  data_list$data <- data
  data_list$data_per_group <- X
  data_list$nobs <- nobs_list
  data_list$nitems <- nitems
  data_list$npatterns <- npatterns
  data_list$nfactors <- nfactors
  data_list$correl <- correl
  data_list$positive <- positive
  data_list$estimator <- estimator
  data_list$cor <- cor
  data_list$group_label <- group_label
  data_list$item_label <- item_label
  data_list$factor_label <- factor_label
  data_list$LAV <- LAV
  data_list$args <- args
  data_list$model <- model_out
  data_list$sample.cov <- sample.cov
  data_list$NACOV <- NACOV
  data_list$WLS.V <- WLS.V
  data_list$thresholds <- thresholds

  return(data_list)

}

create_cfa_model <- function(data_list, model, control) {

  # Generate the model syntax and initial parameter values

  list2env(data_list, envir = environment())

  # Initialize the objects to store the initial parameters:
  param <- trans <- vector("list")
  fixed <- nonfixed <- fixed_values_list <- vector("list")

  #### Model for the transformed parameters ####

  # Initialize the target matrices for positive-definite constraints:
  target_psi <- target_theta <- targets <- vector("list", length = ngroups)
  rest <- 0L

  lambda_group <- paste("lambda.", data_list$group_label, sep = "")
  psi_group <- paste("psi.", data_list$group_label, sep = "")
  theta_group <- paste("theta.", data_list$group_label, sep = "")
  xpsi_group <- paste("xpsi.", data_list$group_label, sep = "")
  xtheta_group <- paste("xtheta.", data_list$group_label, sep = "")
  model_group <- paste("model.", data_list$group_label, sep = "")
  nu_group <- paste("nu.", data_list$group_label, sep = "")
  tau_group <- paste("tau.", data_list$group_label, sep = "")
  S_group <- M_group <- vector("list", length = ngroups)
  for(i in 1:ngroups) {
      for(j in 1:correl[[i]]$npatterns) {
        S_group[[i]][[j]] <- paste("S.", data_list$group_label[i],
                                   ".pattern", j, sep = "")
        M_group[[i]][[j]] <- paste("means.", data_list$group_label[i],
                                       ".pattern", j, sep = "")
    }
  }

  # Transformed parameters:
  list_struct <- vector("list")
  k <- 1L
  for(i in 1:ngroups) {

    # Lambda:
    list_struct[[k]] <- list(name = lambda_group[i],
                             type = "matrix",
                             dim = c(nitems[[i]], nfactors[[i]]),
                             rownames = item_label[[i]],
                             colnames = factor_label[[i]])
    k <- k+1L

    # Create additional parameters if there are positive-definite constraints:
    if(positive) {

      # Theta:
      list_struct[[k]] <- list(name = xtheta_group[i],
                               type = "matrix",
                               dim = c(nitems[[i]], nitems[[i]]),
                               rownames = item_label[[i]],
                               colnames = item_label[[i]])
      k <- k+1L

      # Psi:
      list_struct[[k]] <- list(name = xpsi_group[i],
                               type = "matrix",
                               dim = c(nfactors[[i]], nfactors[[i]]),
                               rownames = factor_label[[i]],
                               colnames = factor_label[[i]])
      k <- k+1L

    }

    # Theta:
    list_struct[[k]] <- list(name = theta_group[i],
                             type = "matrix",
                             dim = c(nitems[[i]], nitems[[i]]),
                             rownames = item_label[[i]],
                             colnames = item_label[[i]],
                             symmetric = TRUE)
    k <- k+1L

    # Psi:
    list_struct[[k]] <- list(name = psi_group[i],
                             type = "matrix",
                             dim = c(nfactors[[i]], nfactors[[i]]),
                             rownames = factor_label[[i]],
                             colnames = factor_label[[i]],
                             symmetric = TRUE)
    k <- k+1L

    # Model matrix:
    list_struct[[k]] <- list(name = model_group[i],
                             type = "matrix",
                             dim = c(nitems[[i]], nitems[[i]]),
                             rownames = item_label[[i]],
                             colnames = item_label[[i]],
                             symmetric = TRUE)
    k <- k+1L

    if(control$meanstructure) {
      # Model means vector:
      list_struct[[k]] <- list(name = nu_group[i],
                               type = "matrix",
                               dim = c(nitems[[i]], 1),
                               rownames = item_label[[i]],
                               colnames = "intrcp")
      k <- k+1L
    }

    # S and M:

      for(j in 1:correl[[i]]$npatterns) {

        nitems_ij <- nrow(correl[[i]][[j]]$S)
        item_label_ij <- correl[[i]][[j]]$item_names
        list_struct[[k]] <- list(name = S_group[[i]][[j]],
                                 type = "matrix",
                                 dim = c(nitems_ij, nitems_ij),
                                 rownames = item_label_ij,
                                 colnames = item_label_ij,
                                 symmetric = TRUE)
        k <- k+1L

        if(control$meanstructure) {
          nitems_ij <- length(correl[[i]][[j]]$item_names)
          list_struct[[k]] <- list(name = M_group[[i]][[j]],
                                   type = "matrix",
                                   dim = c(nitems_ij, 1),
                                   rownames = item_label_ij,
                                   colnames = "intrcp")
          k <- k+1L
        }

    }

  }

  trans <- create_parameters(list_struct)

  #### Replace latent labels by lavaan labels ####

  for(i in 1:ngroups) {

    # Get the positions of parameters and fixed values:

    # Ensure the same label ordering than in model:
    group_i <- c(lambda_group[i], theta_group[i], psi_group[i], nu_group[i])

    nonfixed[group_i] <- lapply(model[[i]][1:4], FUN = \(x) {
      which(is.na(suppressWarnings(as.numeric(x))))
    })

    fixed[group_i] <- lapply(model[[i]][1:4], FUN = \(x) {
      which(!is.na(suppressWarnings(as.numeric(x))))
    })

    fixed_values_list[group_i] <- lapply(model[[i]][1:4], FUN = \(x) {
      numerals <- suppressWarnings(as.numeric(x))
      inds <- which(!is.na(numerals))
      return(numerals[inds])
    })

    trans[[lambda_group[i]]][nonfixed[[lambda_group[i]]]] <- model[[i]]$lambda[nonfixed[[lambda_group[i]]]]
    trans[[theta_group[i]]][nonfixed[[theta_group[i]]]] <- model[[i]]$theta[nonfixed[[theta_group[i]]]]
    trans[[psi_group[i]]][nonfixed[[psi_group[i]]]] <- model[[i]]$psi[nonfixed[[psi_group[i]]]]
    if(control$meanstructure) {
      trans[[nu_group[i]]][nonfixed[[nu_group[i]]]] <- model[[i]]$nu[nonfixed[[nu_group[i]]]]
    }

  }

  #### Model for the parameters ####

  for(i in 1:ngroups) {

    param[[lambda_group[i]]] <- trans[[lambda_group[i]]]
    # Insert fixed values in the model:
    param[[lambda_group[i]]][fixed[[lambda_group[i]]]] <- model[[i]]$lambda[fixed[[lambda_group[i]]]]

    if(positive) {

      # Theta:
      param[[xtheta_group[i]]] <- trans[[xtheta_group[i]]]
      # Psi:
      param[[xpsi_group[i]]] <- trans[[xpsi_group[i]]]

    } else {

      # Theta:
      param[[theta_group[i]]] <- trans[[theta_group[i]]]
      # Insert fixed values in the model:
      param[[theta_group[i]]][fixed[[theta_group[i]]]] <- model[[i]]$theta[fixed[[theta_group[i]]]]
      if(control$deltaparam) {
        diag(param[[theta_group[i]]]) <- "1"
      }

      # Psi:
      param[[psi_group[i]]] <- trans[[psi_group[i]]]
      # Insert fixed values in the model:
      param[[psi_group[i]]][fixed[[psi_group[i]]]] <- model[[i]]$psi[fixed[[psi_group[i]]]]

    }

    # Fix the sample covariance matrix:
    if(control$free_S) {
      param[unlist(S_group[[i]])] <- trans[unlist(S_group[[i]])]
    } else {
        for(j in 1:correl[[i]]$npatterns) {
          param[[S_group[[i]][[j]]]] <- correl[[i]][[j]]$S
      }
    }

    if(cor %in% c("poly", "polys", "polychoric", "polychorics")) {
      fix_diag <- TRUE
    } else {
      fix_diag <- FALSE
    }

    # Fix the diagonal of the sample covariance matrix:
    if(fix_diag) {
      for(j in 1:correl[[i]]$npatterns) {
        diag(param[[S_group[[i]][[j]]]]) <- 1
      }
    }

    if(control$meanstructure) {

      if(fix_diag) {
        for(j in 1:correl[[i]]$npatterns) {
          param[[M_group[[i]][[j]]]] <- matrix(rep(0, length(correl[[i]]$item_names)),
                                               ncol = 1L)
        }
        param[[nu_group[[i]]]] <- matrix(rep(0, nitems[[i]]), ncol = 1L)
      } else {
        for(j in 1:correl[[i]]$npatterns) {
          param[[M_group[[i]][[j]]]] <- matrix(correl[[i]][[j]]$means, ncol = 1L)
        }
        param[[nu_group[i]]] <- trans[[nu_group[i]]]
      }

      if(control$free_M) {
        param[unlist(M_group[[i]])] <- trans[unlist(M_group[[i]])]
      }
      # # Insert fixed values in the model:
      # param[[nu_group[i]]][fixed[[nu_group[i]]]] <- model[[i]]$nu[fixed[[nu_group[i]]]]

    }

    # Create the target matrices for positive-definite constraints:
    if(positive) {

      target_theta[[i]] <- matrix(0, nrow = nitems[[i]], ncol = nitems[[i]])
      target_theta[[i]][nonfixed[[theta_group[i]]]] <- 1

      target_psi[[i]] <- matrix(0, nrow = nfactors[[i]], ncol = nfactors[[i]])
      target_psi[[i]][nonfixed[[psi_group[i]]]] <- 1

      # UGLY FIX THIS TO COUNT THE DEGREES OF FREEDOM AUTOMATICALLY
      q <- nfactors[[i]]
      p <- nitems[[i]]
      lower_theta <- lower.tri(diag(p), diag = TRUE)
      lower_psi <- lower.tri(diag(q), diag = TRUE)
      nconstraints <- sum(unlist(c(target_theta[[i]][lower_theta],
                                   target_psi[[i]][lower_psi])) == 0)
      rest <- rest + 0.5*q*(q-1) + 0.5*p*(p-1) + nconstraints

    }
  }

  #### Create the initial values for the parameters ####

  # Collect the unique nontransformed parameters and the unique transformed parameters:

  init_param <- vector("list", length = control$rstarts)

  for(rs in 1:control$rstarts) {

    init_param[[rs]] <- vector("list")

    for(i in 1:ngroups) {

      init_param[[rs]][[lambda_group[i]]] <- rorth(nitems[[i]], nfactors[[i]])
      init_param[[rs]][[lambda_group[i]]][fixed[[lambda_group[i]]]] <- fixed_values_list[[lambda_group[i]]]

      if(positive) {

        init_param[[rs]][[xtheta_group[i]]] <- rpoblq(nitems[[i]], nitems[[i]], constraints = target_theta[[i]])
        init_param[[rs]][[xpsi_group[i]]] <- rpoblq(nfactors[[i]], nfactors[[i]], constraints = target_psi[[i]])

        init_param[[rs]][[theta_group[i]]] <- crossprod(init_param[[rs]][[xtheta_group[i]]])
        init_param[[rs]][[psi_group[i]]] <- crossprod(init_param[[rs]][[xpsi_group[i]]])

      } else {

        init_param[[rs]][[theta_group[i]]] <- diag(runif(nitems[[i]]))
        init_param[[rs]][[theta_group[i]]][fixed[[theta_group[i]]]] <- fixed_values_list[[theta_group[i]]]

        init_param[[rs]][[psi_group[i]]] <- diag(nfactors[[i]])
        init_param[[rs]][[psi_group[i]]][fixed[[psi_group[i]]]] <- fixed_values_list[[psi_group[i]]]

      }

      Lambda <- init_param[[rs]][[lambda_group[i]]]
      Theta <- init_param[[rs]][[theta_group[i]]]
      Psi <- init_param[[rs]][[psi_group[i]]]
      init_param[[rs]][[model_group[i]]] <- Lambda %*% Psi %*% t(Lambda) + Theta

      if(control$meanstructure) {
        init_param[[rs]][[nu_group[i]]] <- matrix(colMeans(data_list$data_per_group[[i]],
                                                           na.rm = TRUE), ncol = 1L)
      }

      for(j in 1:correl[[i]]$npatterns) {
        init_param[[rs]][[S_group[[i]][[j]]]] <- correl[[i]][[j]]$S
        if(control$meanstructure) {
          init_param[[rs]][[M_group[[i]][[j]]]] <- matrix(correl[[i]][[j]]$means,
                                                          ncol = 1L)
        }
      }

    }

  }

  #### Custom initial values ####

  # # Replace initial starting values by custom starting values:
  #
  # if(!is.null(control$start)) {
  #
  #   nm <- names(control$start)
  #   nm <- nm[!vapply(control$start, is.null, logical(1))]
  #
  #   for (i in seq_len(control$rstarts)) {
  #     common_nm <- intersect(nm, names(init_param[[i]]))
  #     for (j in common_nm) {
  #       init_param[[i]][[j]] <- insert_object(init_param[[i]][[j]],
  #                                             control$start[[j]])
  #     }
  #   }
  #
  # }

  #### Return ####

  result <- list(param = param,
                 trans = trans,
                 init_param = init_param,
                 target_psi = target_psi,
                 target_theta = target_theta,
                 rest = rest)

  return(result)

}

create_cfa_modelInfo <- function(data_list, full_model, control) {

  # Generate control_manifold, control_transform, and control_estimator

  list2env(data_list, envir = environment())
  list2env(full_model, envir = environment())

  lambda_group <- paste("lambda.", data_list$group_label, sep = "")
  psi_group <- paste("psi.", data_list$group_label, sep = "")
  theta_group <- paste("theta.", data_list$group_label, sep = "")
  xpsi_group <- paste("xpsi.", data_list$group_label, sep = "")
  xtheta_group <- paste("xtheta.", data_list$group_label, sep = "")
  model_group <- paste("model.", data_list$group_label, sep = "")
  nu_group <- paste("nu.", data_list$group_label, sep = "")
  S_group <- M_group <- vector("list", length = ngroups)
  for(i in 1:ngroups) {
    for(j in 1:correl[[i]]$npatterns) {
      S_group[[i]][[j]] <- paste("S.", data_list$group_label[i],
                                 ".pattern", j, sep = "")
      M_group[[i]][[j]] <- paste("means.", data_list$group_label[i],
                                 ".pattern", j, sep = "")
    }
  }

  #### Manifolds ####

  manifolds <- list()
  k <- 1L

  for(i in 1:ngroups) {

    manifolds[[k]] <- list(manifold = "euclidean",
                           parameters = lambda_group[i])
    k <- k+1L

    manifolds[[k]] <- list(manifold = "euclidean",
                           parameters = nu_group[i])
    k <- k+1L


    manifolds[[k]] <- list(manifold = "euclidean",
                           parameters = unlist(M_group[[i]]))
    k <- k+1L

    manifolds[[k]] <- list(manifold = "euclidean",
                           parameters = unlist(S_group[[i]]))
    k <- k+1L

    if(positive) {

      manifolds[[k]] <- list(manifold = "poblq",
                             parameters = xpsi_group[i],
                             extra = list(
                               p = nfactors[[i]],
                               q = nfactors[[i]],
                               constraints = target_psi[[i]]))
      k <- k+1L

      manifolds[[k]] <- list(manifold = "poblq",
                             parameters = xtheta_group[i],
                             extra = list(
                               p = nitems[[i]],
                               q = nitems[[i]],
                               constraints = target_theta[[i]]))
      k <- k+1L

    } else {

      manifolds[[k]] <- list(manifold = "euclidean",
                             parameters = psi_group[i])
      k <- k+1L

      manifolds[[k]] <- list(manifold = "euclidean",
                             parameters = theta_group[i])
      k <- k+1L

    }

  }

  control_manifold <- create_manifolds(manifolds = manifolds,
                                       structures = param)

  #### Transformations ####

  transforms <- list()
  dots <- list()
  k <- 1L

  for(i in 1:ngroups) {

    if(positive) {

      lower_psi <- lower.tri(trans[[psi_group[i]]], diag = TRUE)
      lower_theta <- lower.tri(trans[[theta_group[i]]], diag = TRUE)

      dots$p <- nrow(trans[[psi_group[i]]])
      transforms[[k]] <- list(transform = "crossprod",
                             parameters_in = xpsi_group[i],
                             parameters_out = psi_group[i],
                             extra = dots)
      k <- k+1L

      dots$p <- nrow(trans[[theta_group[i]]])
      transforms[[k]] <- list(transform = "crossprod",
                              parameters_in = xtheta_group[i],
                              parameters_out = theta_group[i],
                              extra = dots)
      k <- k+1L

    }

    if(control$deltaparam) {

      dots$p <- nrow(trans[[theta_group[i]]])
      dots$q <- nrow(trans[[psi_group[i]]])
      transforms[[k]] <- list(transform = "deltaparam",
                              parameters_in = c(lambda_group[i],
                                                psi_group[i]),
                              parameters_out = list(diag(trans[[theta_group[i]]])),
                              extra = dots)
      k <- k+1L

    }

    # Model matrix correlation:
    lower_psi <- lower.tri(trans[[psi_group[i]]], diag = TRUE)
    lower_theta <- lower.tri(trans[[theta_group[i]]], diag = TRUE)
    lower_diag <- lower.tri(trans[[model_group[i]]], diag = TRUE)

    dots$p <- nitems[[i]]
    dots$q <- nfactors[[i]]
    transforms[[k]] <- list(transform = "factor_cor",
                            parameters_in = c(lambda_group[i],
                                              psi_group[i],
                                              theta_group[i]),
                            parameters_out = model_group[i],
                            extra = dots)
    k <- k+1L

  }

  control_transform <- create_transforms(transforms = transforms,
                                         structures = trans)

  #### Estimators ####

  estimators <- list()
  k <- 1L

  for(i in 1:ngroups) {

    estimator <- tolower(estimator)
    cfa_estimator <- switch(estimator,
                            uls = "cfa_dwls",
                            dwls  = "cfa_dwls",
                            ml = "cfa_fml",
                            fml  = "cfa_fml",
                            means_fml = "cfa_means_fml",
                            means_dwls = "cfa_means_dwls",
                            means_uls = "cfa_means_dwls",
                            stop("Unknown estimator: ", estimator)
    )

    for(j in 1:correl[[i]]$npatterns) {

      pick <- correl[[i]][[j]]$item_names
      S_group_ij <- S_group[[i]][[j]]
      M_group_ij <- M_group[[i]][[j]]
      estimators[[k]] <- list(estimator = cfa_estimator,
                              parameters = list(c(trans[[model_group[i]]][pick, pick]),
                                                c(trans[[S_group_ij]]),
                                                c(trans[[nu_group[i]]][pick, ]),
                                                c(trans[[M_group_ij]])),
                              extra = list(W = correl[[i]][[j]]$W,
                                           w = correl[[i]][[j]]$nobs /
                                             sum(unlist(nobs)),
                                           w_means = correl[[i]][[j]]$w_means,
                                           q = nrow(trans[[psi_group[i]]]),
                                           p = nrow(correl[[i]][[j]]$S),
                                           n = correl[[i]][[j]]$nobs))

      k <- k + 1L

    }

  }

  if(positive & control$reg) {

    for(i in 1:ngroups) {

      # For the psi matrix:

      lower_indices <- which(lower.tri(trans[[psi_group[i]]], diag = TRUE))
      estimators[[k]] <- list(estimator = "logdetmat",
                              parameters = psi_group[i],
                              extra = list(lower_indices = lower_indices-1L,
                                           p = nrow(trans[[psi_group[i]]]),
                                           logdetw = control$penalties$logdet$w))
      k <- k+1L

      # For the theta matrix:

      lower_indices <- which(lower.tri(trans[[theta_group[i]]], diag = TRUE))
      estimators[[k]] <- list(estimator = "logdetmat",
                              parameters = theta_group[i],
                              extra = list(lower_indices = lower_indices-1L,
                                           p = nrow(trans[[theta_group[i]]]),
                                           logdetw = control$penalties$logdet$w))
      k <- k+1L

    }

  }

  control_estimator <- create_estimators(estimators = estimators,
                                         structures = trans)

  #### Pass the initial values to vectors ####

  idx_transformed <- unlist(lapply(control_transform,
                                   FUN = \(x) unlist(x$indices_out)+1L))
  inits <- create_init(trans, param, init_param,
                       idx_transformed = idx_transformed, control)

  parameters <- inits$parameters
  parameters_labels <- names(parameters[[1]])
  nparam <- length(parameters_labels)

  transparameters <- inits$transparameters
  transparameters_labels <- names(transparameters[[1]])
  ntrans <- length(transparameters_labels)

  trans2param <- match(parameters_labels, transparameters_labels)

  #### Set up the optimizer ####

  # Create defaults for the control of the optimizer:
  control_optimizer <- control
  control_optimizer$parameters <- parameters
  control_optimizer$transparameters <- transparameters
  control_optimizer$init_param <- init_param
  control_optimizer$transparam2param <- trans2param-1L

  #### Collect all the model information ####

  modelInfo <- list(param = param,
                    trans = trans,
                    nparam = nparam - rest,
                    ntrans = ntrans,
                    parameters_labels = parameters_labels,
                    transparameters_labels = transparameters_labels,
                    dof = sum(unlist(npatterns)) - nparam + rest,
                    control_manifold = control_manifold,
                    control_transform = control_transform,
                    control_estimator = control_estimator,
                    control_optimizer = control_optimizer)

  #### Return ####

  return(modelInfo)

}
