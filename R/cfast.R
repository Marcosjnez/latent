#' @title
#' Confirmatory factor analysis.
#' @export
cfast <- function(data, target, targetphi = NULL, targetpsi = NULL, cor = "pearson",
                  estimator = "uls", missing = "pairwise.complete.cases",
                  nobs = NULL, group = NULL, invariance = "metric",
                  control = NULL, positive = FALSE) {

  # cor = "pearson"; estimator = "uls";
  # missing = "pairwise.complete.cases"
  # nobs = 500; control = NULL;
  # X <- S
  # positive = TRUE
  # group = NULL; invariance = "metric"

  # S is a list of correlation matrices (one for each group)
  # target is a list of matrices for the loadings indicating the parameters
  # (as characters) and fixed values (numeric)
  # targetphi is a list of matrices for the factor correlations indicating the
  # parameters (as characters) and fixed values (numeric)
  # targetpsi is a list of matrices for the uniquenesses and correlated errors
  # indicating the parameters (as characters) and fixed values (numeric)

  # The targets must be fully specified, that is, every entry in a target
  # must contain either a character or a numeric value

  if(is.matrix(target)) {
    target1 <- target
    target <- list()
    target[[1]] <- target1
  }

  if(is.null(targetphi)) {
    for(i in 1:n) {
      q <- ncol(target[[i]])
      targetphi[[i]] <- diag(q)
    }
  } else if(is.matrix(targetphi)) {
    targetphi1 <- targetphi
    targetphi <- list()
    targetphi[[1]] <- targetphi1
  }

  if(is.null(targetpsi)) {
    for(i in 1:n) {
      p <- nrow(target[[i]])
      targetpsi[[i]] <- matrix(0, p, p)
      diag(targetpsi[[i]]) <- paste(i, "psi", 1:p, sep = "")
    }
  } else if(is.matrix(targetpsi)) {
    targetpsi1 <- targetpsi
    targetpsi <- list()
    targetpsi[[1]] <- targetpsi1
  }

  if(is.data.frame(data) & !is.null(group)) {

    groups <- unique(data$group)
    n <- length(groups)
    if(n == 1L) stop("Group invariance is not possible with only one group")
    data1 <- data
    data <- list()
    remove <- which(colnames(data1) == "group")

    for(i in 1:n) {
      data[[i]] <- as.matrix(data1[data1$group == groups[i], -remove])
      nobs <- vector(length = n)
      nobs[i] <- nrows(data[[i]])
    }

    if(invariance == "metric") {
      for(i in 1:n) {
        target[[i]] = target[[1]]
      }
    } else if(invariance == "scalar") {
      for(i in 1:n) {
        target[[i]] = target[[1]]
      }
    } else if(invariance == "residual") {
      for(i in 1:n) {
        target[[i]] = target[[1]]
        targetpsi[[i]] = targetpsi[[1]]
      }
    } else {
      stop("Unkown invariance")
    }

  }

  n <- length(data) # Number of groups
  p <- q <- vector(length = n) # Number of variables (p) and factors (q) in each group
  Lambda <- Phi <- Psi <-
    indexes_lambda <- indexes_phi <- indexes_psi <-
    indexes_target <- indexes_targetphi <- indexes_targetpsi <-
    lambda_hat <- phi_hat <- psi_hat <- uniquenesses_hat <-
    R <- Rhat <- residuals <- fill_Phi_Target <- free_indices_phi <- list()
  indexes_factorvars <- indexes_uniquenesses <- c()
  original_targetphi <- targetphi

  if(positive) {
    for(i in 1:n) {
      fill_Phi_Target[[i]] <- which(is.na(suppressWarnings(as.numeric(targetphi[[i]]))))
      length_phi <- prod(dim(targetphi[[i]]))
      targetphi[[i]] <- matrix(paste(i, "proj", 1:length_phi, sep = ""),
                               nrow(targetphi[[i]]), ncol(targetphi[[i]]))
    }
  }

  # Find the unique elements in all the targets:
  uniques <- unique(unlist(c(target, targetphi, targetpsi)))
  # Get the parameters (those elements that are not digits):
  z <- suppressWarnings(as.numeric(uniques))
  parameter_vector <- uniques[which(is.na(z))]
  # Get the fixed values (those elements that are digits):
  fixed_vector <- uniques[which(!is.na(z))]
  # Initialize the vector of initial parameter estimates:
  init <- vector(length = length(parameter_vector))

  for(i in 1:n) { # For each group...

    p[i] <- nrow(target[[i]]) # Number of variables in group i
    q[i] <- ncol(target[[i]]) # Number of factors in group i
    free_indices_phi[[i]] <- integer(0)

    # Find which elements in the targets correspond to a parameter:
    indexes_target[[i]] <- which(target[[i]] %in% parameter_vector) # Which lambdas are estimated in group i
    if(positive) {
      indexes_targetphi[[i]] <- which(targetphi[[i]] %in% parameter_vector) # Which phis are estimated in group i
      free_indices_phi[[i]] <- which(is.na(suppressWarnings(as.numeric(diag(original_targetphi[[i]])))))
    } else {
      indexes_targetphi[[i]] <- which(targetphi[[i]] %in% parameter_vector & lower.tri(targetphi[[i]], diag = TRUE)) # Which phis are estimated in group i
    }
    indexes_targetpsi[[i]] <- which(targetpsi[[i]] %in% parameter_vector & lower.tri(targetpsi[[i]], diag = TRUE)) # Which psis are estimated in group i

    # Get the indexes for the factor variances:
    indexes_factorvars <- c(indexes_factorvars,
                            which(parameter_vector %in% diag(targetphi[[i]])))

    # Get the indexes for the uniquenesses:
    indexes_uniquenesses <- c(indexes_uniquenesses,
                              which(parameter_vector %in% diag(targetpsi[[i]])))

    # Relate the parameters in the targets to the parameters in the parameter vector:
    if(length(indexes_target[[i]]) == 0) {
      indexes_lambda[[i]] <- logical(0)
    } else {
      indexes_lambda[[i]] <- match(target[[i]][indexes_target[[i]]], parameter_vector) # Which lambdas are estimated in group i
    }
    if(length(indexes_targetphi[[i]]) == 0) {
      indexes_phi[[i]] <- logical(0)
    } else {
      indexes_phi[[i]] <- match(targetphi[[i]][indexes_targetphi[[i]]], parameter_vector) # Which phis are estimated in group i
    }
    if(length(indexes_targetpsi[[i]]) == 0) {
      indexes_psi[[i]] <- logical(0)
    } else {
      indexes_psi[[i]] <- match(targetpsi[[i]][indexes_targetpsi[[i]]], parameter_vector) # Which psis are estimated in group i
    }

    # Find which elements in the targets for correspond to a fixed value:
    indexes_fixtarget <- which(target[[i]] %in% fixed_vector) # Which lambdas are fixed in group i
    indexes_fixtargetphi <- which(targetphi[[i]] %in% fixed_vector) # Which phis are fixed in group i
    indexes_fixtargetpsi <- which(targetpsi[[i]] %in% fixed_vector) # Which psis are fixed in group i

    # Relate the elements in the targets to the fixed values in the fixed vector:
    if(length(indexes_fixtarget) == 0) {
      indexes_fixlambda <- logical(0)
    } else {
      indexes_fixlambda <- match(target[[i]][indexes_fixtarget], fixed_vector) # Which lambdas are fixed in group i
    }
    if(length(indexes_fixtargetphi) == 0) {
      indexes_fixphi <- logical(0)
    } else {
      indexes_fixphi <- match(targetphi[[i]][indexes_fixtargetphi], fixed_vector) # Which phis are fixed in group i
    }
    if(length(indexes_fixtargetpsi) == 0) {
      indexes_fixpsi <- logical(0)
    } else {
      indexes_fixpsi <- match(targetpsi[[i]][indexes_fixtargetpsi], fixed_vector) # Which psis are fixed in group i
    }

    # non-specified elements in lambda are 0:
    lambda_hat[[i]] <- matrix(0, p[i], q[i])
    lambda_hat[[i]][indexes_fixtarget] <- fixed_vector[indexes_fixlambda]
    # non-specified elements in phi are zero if off-diagonal and 1 if diagonal:
    phi_hat[[i]] <- matrix(0, q[i], q[i]); #diag(phi_hat[[i]]) <- 1
    if(positive) phi_hat[[i]][fill_Phi_Target[[i]]] <- 1
    phi_hat[[i]][indexes_fixtargetphi] <- fixed_vector[indexes_fixphi]
    # non-specified elements in psi are zero if off-diagonal and estimated if diagonal:
    psi_hat[[i]] <- matrix(0, p[i], p[i])
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
             cor = cor, estimator = estimator,
             projection = projection,
             missing = rep(missing, n),
             control = control)
  # c(fit$cfa$parameters)[indexes_phi[[1]]]
  # fit$cfa$phi
  # fit$cfa$f
  # fit$cfa$iterations

  # Arrange the parameter estimates in the lambda, phi, and psi matrices:
  for(i in 1:n) {

    # Arrange lambda parameter estimates:
    lambda_hat[[i]] <- matrix(fit$cfa$lambda[[i]], p[i], q[i])

    # Arrange phi parameter estimates:
    phi_hat[[i]] <- matrix(fit$cfa$phi[[i]], q[i], q[i])

    # Arrange psi parameter estimates:
    psi_hat[[i]] <- matrix(fit$cfa$psi[[i]], p[i], p[i])
    uniquenesses_hat[[i]] <- diag(psi_hat[[i]])

    # Model matrix:
    R[[i]] <- matrix(fit$cfa$R[[i]], p[i], p[i])
    Rhat[[i]] <- matrix(fit$cfa$Rhat[[i]], p[i], p[i])
    residuals[[i]] <- R[[i]] - Rhat[[i]]

  }

  return(fit)

}

