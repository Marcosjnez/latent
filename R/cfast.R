#' @title
#' Confirmatory factor analysis.
#' @export
cfast <- function(data, lambda, phi = NULL, psi = NULL, cor = "pearson",
                  estimator = "uls", rotate = NULL, missing = "pairwise.complete.cases",
                  nobs = NULL, group = NULL, invariance = "metric",
                  control = NULL, positive = FALSE, random_starts = 1L, cores = 1L) {

  # cor = "pearson"; estimator = "uls";
  # missing = "pairwise.complete.cases"
  # nobs = 500; control = NULL;
  # data <- sim$R_error
  # positive = TRUE
  # group = NULL; invariance = "metric"
  # lambda <- target
  # phi <- targetphi
  # psi <- targetpsi

  # S is a list of correlation matrices (one for each group)
  # lambda is a list of matrices for the loadings indicating the parameters
  # (as characters) and fixed values (numeric)
  # phi is a list of matrices for the factor correlations indicating the
  # parameters (as characters) and fixed values (numeric)
  # psi is a list of matrices for the uniquenesses and correlated errors
  # indicating the parameters (as characters) and fixed values (numeric)

  # The targets must be fully specified, that is, every entry in a target
  # must contain either a character or a numeric value

  # Transform single matrices into a list with a matrix:
  if(is.matrix(data)) {

    n <- 1
    data1 <- data
    data <- list()
    data[[1]] <- data1

  } else if(is.data.frame(data) & !is.null(group)) {

    groups <- unique(data$group)
    n <- length(groups)
    if(n == 1L) stop("Group invariance is not possible with only one group")
    data1 <- data
    data <- list()
    remove <- which(colnames(data1) == group)

    # Repeat the same model for each group:
    for(i in 1:n) {
      data[[i]] <- as.matrix(data1[data1$group == groups[i], -remove])
      nobs <- vector(length = n)
      nobs[i] <- nrows(data[[i]])
    }

    if(invariance == "metric") {
      for(i in 1:n) {
        lambda[[i]] = lambda[[1]]
      }
    } else if(invariance == "scalar") {
      stop("scalar invariance not available yet")
      for(i in 1:n) {
        lambda[[i]] = lambda[[1]]
      }
    } else if(invariance == "residual") {
      for(i in 1:n) {
        lambda[[i]] = lambda[[1]]
        psi[[i]] = psi[[1]]
      }
    } else {
      stop("Unkown invariance")
    }

  }

  if(is.matrix(lambda)) {
    for(i in 1:n) {
      lambda1 <- lambda
      lambda <- list()
      lambda[[i]] <- lambda1
    }
  }

  if(is.null(phi)) {
    for(i in 1:n) {
      q <- ncol(lambda[[i]])
      phi <- list()
      phi[[i]] <- diag(q) # Identity by default
    }
  } else if(is.matrix(phi)) {
    phi1 <- phi
    phi <- list()
    phi[[1]] <- phi1
  }

  # By default, estimate the uniquenesses:
  if(is.null(psi)) {
    for(i in 1:n) {
      p <- nrow(lambda[[i]])
      psi[[i]] <- matrix(0, p, p)
      diag(psi[[i]]) <- paste(i, "psi", 1:p, sep = "")
    }
  } else if(is.matrix(psi)) {
    psi1 <- psi
    psi <- list()
    psi[[1]] <- psi1
  }

  n <- length(data) # Number of groups
  p <- q <- vector(length = n) # Number of variables (p) and factors (q) in each group
  Lambda <- Phi <- Psi <-
    indexes_lambda <- indexes_phi <- indexes_psi <-
    indexes_target <- indexes_targetphi <- indexes_targetpsi <-
    lambda_hat <- phi_hat <- psi_hat <- uniquenesses_hat <-
    R <- Rhat <- residuals <- fill_Phi_Target <- fill_Psi_Target <-
    free_indices_phi <- free_indices_psi <- list()
  indexes_factorvars <- indexes_uniquenesses <- c()
  original_targetphi <- phi
  original_targetpsi <- psi

  if(positive) {
    for(i in 1:n) {
      fill_Phi_Target[[i]] <- which(is.na(suppressWarnings(as.numeric(phi[[i]]))))
      length_phi <- prod(dim(phi[[i]]))
      phi[[i]] <- matrix(paste(i, "proj", 1:length_phi, sep = ""),
                               nrow(phi[[i]]), ncol(phi[[i]]))
      fill_Psi_Target[[i]] <- which(is.na(suppressWarnings(as.numeric(psi[[i]]))))
      length_psi <- prod(dim(psi[[i]]))
      psi[[i]] <- matrix(paste(i, "psi_proj", 1:length_psi, sep = ""),
                         nrow(psi[[i]]), ncol(psi[[i]]))
    }
  }

  # Find the unique elements in all the targets:
  uniques <- unique(unlist(c(lambda, phi, psi)))
  # Get the parameters (those elements that are not digits):
  z <- suppressWarnings(as.numeric(uniques))
  parameter_vector <- uniques[which(is.na(z))]
  # Get the fixed values (those elements that are digits):
  fixed_vector <- uniques[which(!is.na(z))]
  # Initialize the vector of initial parameter estimates:
  init <- vector(length = length(parameter_vector))

  for(i in 1:n) { # For each group...

    p[i] <- nrow(lambda[[i]]) # Number of variables in group i
    q[i] <- ncol(lambda[[i]]) # Number of factors in group i
    free_indices_phi[[i]] <- free_indices_psi[[i]] <- integer(0)

    # Find which elements in the targets correspond to a parameter:
    indexes_target[[i]] <- which(lambda[[i]] %in% parameter_vector) # Which lambdas are estimated in group i
    if(positive) {
      indexes_targetphi[[i]] <- which(phi[[i]] %in% parameter_vector) # Which phis are estimated in group i
      free_indices_phi[[i]] <- which(is.na(suppressWarnings(as.numeric(diag(original_targetphi[[i]])))))
      indexes_targetpsi[[i]] <- which(psi[[i]] %in% parameter_vector) # Which psis are estimated in group i
      free_indices_psi[[i]] <- which(is.na(suppressWarnings(as.numeric(diag(original_targetpsi[[i]])))))
    } else {
      indexes_targetphi[[i]] <- which(phi[[i]] %in% parameter_vector & lower.tri(phi[[i]], diag = TRUE)) # Which phis are estimated in group i
      indexes_targetpsi[[i]] <- which(psi[[i]] %in% parameter_vector & lower.tri(psi[[i]], diag = TRUE)) # Which psis are estimated in group i
    }

    # Get the indexes for the factor variances:
    indexes_factorvars <- c(indexes_factorvars,
                            which(parameter_vector %in% diag(phi[[i]])))

    # Get the indexes for the uniquenesses:
    indexes_uniquenesses <- c(indexes_uniquenesses,
                              which(parameter_vector %in% diag(psi[[i]])))

    # Relate the parameters in the targets to the parameters in the parameter vector:
    if(length(indexes_target[[i]]) == 0) {
      indexes_lambda[[i]] <- logical(0)
    } else {
      indexes_lambda[[i]] <- match(lambda[[i]][indexes_target[[i]]], parameter_vector) # Which lambdas are estimated in group i
    }
    if(length(indexes_targetphi[[i]]) == 0) {
      indexes_phi[[i]] <- logical(0)
    } else {
      indexes_phi[[i]] <- match(phi[[i]][indexes_targetphi[[i]]], parameter_vector) # Which phis are estimated in group i
    }
    if(length(indexes_targetpsi[[i]]) == 0) {
      indexes_psi[[i]] <- logical(0)
    } else {
      indexes_psi[[i]] <- match(psi[[i]][indexes_targetpsi[[i]]], parameter_vector) # Which psis are estimated in group i
    }

    # Find which elements in the targets for correspond to a fixed value:
    indexes_fixtarget <- which(lambda[[i]] %in% fixed_vector) # Which lambdas are fixed in group i
    indexes_fixtargetphi <- which(phi[[i]] %in% fixed_vector) # Which phis are fixed in group i
    indexes_fixtargetpsi <- which(psi[[i]] %in% fixed_vector) # Which psis are fixed in group i

    # Relate the elements in the targets to the fixed values in the fixed vector:
    if(length(indexes_fixtarget) == 0) {
      indexes_fixlambda <- logical(0)
    } else {
      indexes_fixlambda <- match(lambda[[i]][indexes_fixtarget], fixed_vector) # Which lambdas are fixed in group i
    }
    if(length(indexes_fixtargetphi) == 0) {
      indexes_fixphi <- logical(0)
    } else {
      indexes_fixphi <- match(phi[[i]][indexes_fixtargetphi], fixed_vector) # Which phis are fixed in group i
    }
    if(length(indexes_fixtargetpsi) == 0) {
      indexes_fixpsi <- logical(0)
    } else {
      indexes_fixpsi <- match(psi[[i]][indexes_fixtargetpsi], fixed_vector) # Which psis are fixed in group i
    }

    # non-specified elements in lambda are 0:
    lambda_hat[[i]] <- matrix(0, p[i], q[i])
    lambda_hat[[i]][indexes_fixtarget] <- fixed_vector[indexes_fixlambda]
    # non-specified elements in phi are zero if off-diagonal and 1 if diagonal:
    phi_hat[[i]] <- matrix(0, q[i], q[i]); #diag(phi_hat[[i]]) <- 1
    psi_hat[[i]] <- matrix(0, p[i], p[i]); #diag(psi_hat[[i]]) <- 1
    if(positive) {
      phi_hat[[i]][fill_Phi_Target[[i]]] <- 1
      psi_hat[[i]][fill_Psi_Target[[i]]] <- 1
    }
    # non-specified elements in phi and psi are zero if off-diagonal and estimated if diagonal:
    phi_hat[[i]][indexes_fixtargetphi] <- fixed_vector[indexes_fixphi]
    psi_hat[[i]][indexes_fixtargetpsi] <- fixed_vector[indexes_fixpsi]

    class(lambda_hat[[i]]) <- "numeric"
    class(phi_hat[[i]]) <- "numeric"
    class(psi_hat[[i]]) <- "numeric"

    # Initial lambda and uniqueness values based on the eigendecomposition of
    # the reduced correlation matrix:
    if(nrow(data[[i]]) == ncol(data[[i]])) {
      S <- data[[i]]
    } else {
      S <- cor(data[[i]])
    }
    u <- 1/diag(solve(S))
    diag(S) <- u
    e <- eigen(S)
    D <- matrix(0, q[i], q[i])
    diag(D) <- sqrt(e$values[1:q[i]])
    V <- e$vectors[, 1:q[i]]
    VD <- V %*% D
    VV <- VD %*% t(VD)
    init[indexes_lambda[[i]]] <- VV[indexes_target[[i]]]
    init[indexes_phi[[i]]] <- diag(q[i])[indexes_targetphi[[i]]]
    init[indexes_psi[[i]]] <- diag(u)[indexes_targetpsi[[i]]]

  }

  indexes_factorvars <- unique(indexes_factorvars)
  indexes_uniquenesses <- unique(indexes_uniquenesses)
  lambda_p <- length(unique(unlist(indexes_lambda)))
  phi_p <- length(unique(unlist(indexes_phi)))
  psi_p <- length(unique(unlist(indexes_psi)))

  lower_psi <- rep(-0.995, psi_p) # Lower bounds for correlated residuals
  upper_psi <- rep(0.995, psi_p) # Upper bounds for correlated residuals
  lower <- c(rep(-10, lambda_p), rep(-0.995, phi_p), lower_psi)
  upper <- c(rep(10, lambda_p), rep(0.995, phi_p), upper_psi)
  lower[indexes_uniquenesses] <- 0.001
  upper[indexes_uniquenesses] <- 0.999
  lower[indexes_factorvars] <- 0.001
  upper[indexes_factorvars] <- 10
  control$lower <- lower
  control$upper <- upper
  control$target_positive <- c(indexes_factorvars, indexes_uniquenesses)

  # x <- init
  if(!is.null(control$init)) {
    x <- init
  } else {
    x <- init #c(stats::runif(lambda_p), rep(0.5, phi_p), rep(0.5, psi_p))
  }

  if(positive) {
    projection = "positive"
  } else {
    projection <- "id"
  }

  if(length(cor) == 1L) cor <- rep(cor, n)
  if(length(estimator) == 1L) estimator <- rep(estimator, n)
  if(length(projection) == 1L) projection <- rep(projection, n)
  if(length(missing) == 1L) missing <- rep(missing, n)
  if(length(nobs) == 1L) nobs <- rep(nobs, n)

  fit <- cfa(parameters = x, X = data, nfactors = q, nobs = nobs,
             lambda = lambda_hat, phi = phi_hat, psi = psi_hat,
             lambda_indexes = indexes_lambda,
             phi_indexes = indexes_phi,
             psi_indexes = indexes_psi,
             target_indexes = indexes_target,
             targetphi_indexes = indexes_targetphi,
             targetpsi_indexes = indexes_targetpsi,
             free_indices_phi = free_indices_phi,
             free_indices_psi = free_indices_psi,
             cor = cor, estimator = estimator,
             projection = projection,
             missing = rep(missing, n),
             se = "robust",
             control = control)
  # c(fit$cfa$parameters)[indexes_phi[[1]]]
  # fit$cfa$phi
  # fit$cfa$f
  # fit$cfa$iterations

  # Arrange the parameter estimates in the lambda, phi, and psi matrices:
  # for(i in 1:n) {
  #
  #   # Arrange lambda parameter estimates:
  #   lambda_hat[[i]] <- matrix(fit$cfa$lambda[[i]], p[i], q[i])
  #
  #   # Arrange phi parameter estimates:
  #   phi_hat[[i]] <- matrix(fit$cfa$phi[[i]], q[i], q[i])
  #
  #   # Arrange psi parameter estimates:
  #   psi_hat[[i]] <- matrix(fit$cfa$psi[[i]], p[i], p[i])
  #   uniquenesses_hat[[i]] <- diag(psi_hat[[i]])
  #
  #   # Model matrix:
  #   R[[i]] <- matrix(fit$cfa$R[[i]], p[i], p[i])
  #   Rhat[[i]] <- matrix(fit$cfa$Rhat[[i]], p[i], p[i])
  #   residuals[[i]] <- R[[i]] - Rhat[[i]]
  #
  # }

  # if(!is.null(rotate)) {
  #   rotation <- rotate$rotation
  #   projection <- rotate$projection
  #   gamma <- rotate$gamma
  #   epsilon <- rotate$gamma
  #   if(!is.null(rotate$k)) k <- rotate$k
  #   w <- rotate$w
  #   Target <- rotate$Target
  #   Weight <- rotate$Weight
  #   PhiTarget <- rotate$PhiTarget
  #   PhiWeight <- rotate$PhiWeight
  #   blocks <- rotate$blocks
  #   block_weights <- rotate$block_weights
  #   oblq_factors <- rotate$oblq_factors
  #   normalization <- rotate$normalization
  #   rot_control <- rotate$rot_control
  #   random_starts <- random_starts
  #   cores <- cores
  #
  #   for(i in 1:n) {
  #     lambda <- fit$cfa$lambda[[i]]
  #     rot <- rotate(lambda = lambda, rotation = rotation, projection = projection,
  #                   gamma = gamma, epsilon = epsilon, k = k, w = w, Target = Target,
  #                   Weight = Weight, PhiTarget = PhiTarget, PhiWeight = PhiWeight,
  #                   blocks = blocks, block_weights = block_weights,
  #                   oblq_factors = oblq_factors, normalization = normalization,
  #                   rot_control = rot_control, random_starts = random_starts, cores = cores)
  #     fit$rotation[[i]] <- rot
  #   }
  # }

  return(fit)

}

