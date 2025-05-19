#' @title
#' Confirmatory factor analysis.
#' @export
cfast <- function(data, model = NULL, lambda = NULL, phi = NULL, psi = NULL, cor = "pearson",
                  estimator = "uls", rotate = NULL, missing = "pairwise.complete.cases",
                  nobs = NULL, group = NULL, invariance = "none",
                  control = NULL, std.lv = FALSE, positive = FALSE) {

  # cor = "pearson";
  # group = NULL;
  # invariance = "none";
  # missing = "pairwise.complete.cases"
  # data <- sim$R_error
  # positive = FALSE
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

  control <- lca_control(control)

  extract_fit <- lavaan::cfa(model = model, data = data, std.lv = std.lv,
                             do.fit = FALSE, group = group)
  item_names <- unique(extract_fit@ParTable$rhs[extract_fit@ParTable$op == "=~"])
  dataset <- data[, c(item_names, group)]
  matrices <- getmodel_cfa(extract_fit)

  if(is.null(group)) {

    lambda <- matrices$lambda
    phi <- matrices$psi
    psi <- matrices$theta
    ngroups <- 1 # Number of groups
    data1 <- as.matrix(dataset)
    dataset <- list()
    dataset[[1]] <- data1
    nobs <- nrow(data1)

  } else {

    lambda <- lapply(matrices, FUN = \(x) x$lambda)
    phi <- lapply(matrices, FUN = \(x) x$psi)
    psi <- lapply(matrices, FUN = \(x) x$theta)

    groups <- unique(data[, group]) # Number of groups
    ngroups <- length(groups)
    if(ngroups == 1L) stop("Group invariance is not possible with only one group")
    data1 <- dataset
    dataset <- list()

    # Repeat the same model for each group:
    nobs <- vector(length = ngroups)
    for(i in 1:ngroups) {
      dataset[[i]] <- as.matrix(data1[data1[, group] == groups[i], item_names])
      nobs[i] <- nrow(dataset[[i]])
    }

    if(invariance == "none") {
    } else if(invariance == "metric") {
      for(i in 1:ngroups) {
        lambda[[i]] = lambda[[1]]
      }
    } else if(invariance == "scalar") {
      stop("scalar invariance not available yet")
      for(i in 1:ngroups) {
        lambda[[i]] = lambda[[1]]
      }
    } else if(invariance == "residual") {
      for(i in 1:ngroups) {
        lambda[[i]] = lambda[[1]]
        psi[[i]] = psi[[1]]
      }
    } else {
      stop("Unkown type of invariance")
    }

  }

  if(length(estimator) == 1L) estimator <- rep(estimator, ngroups)

  # Arrange the matrices of parameters:

  if(is.null(lambda)) {
    stop("Provide a model, please")
  } else if(is.matrix(lambda)) {
    for(i in 1:ngroups) {
      lambda1 <- lambda
      lambda <- list()
      lambda[[i]] <- lambda1
    }
  }

  # By default, everything is orthogonal:
  if(is.null(phi)) {
    for(i in 1:ngroups) {
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
    for(i in 1:ngroups) {
      p <- nrow(lambda[[i]])
      psi[[i]] <- matrix(0, p, p)
      diag(psi[[i]]) <- paste(i, "psi", 1:p, sep = "")
    }
  } else if(is.matrix(psi)) {
    psi1 <- psi
    psi <- list()
    psi[[1]] <- psi1
  }

  positive <- FALSE
  # if(positive) {
  #   for(i in 1:ngroups) {
  #     fill_Phi_Target[[i]] <- which(is.na(suppressWarnings(as.numeric(phi[[i]]))))
  #     length_phi <- prod(dim(phi[[i]]))
  #     phi[[i]] <- matrix(paste(i, "proj", 1:length_phi, sep = ""),
  #                        nrow(phi[[i]]), ncol(phi[[i]]))
  #     fill_Psi_Target[[i]] <- which(is.na(suppressWarnings(as.numeric(psi[[i]]))))
  #     length_psi <- prod(dim(psi[[i]]))
  #     psi[[i]] <- matrix(paste(i, "psi_proj", 1:length_psi, sep = ""),
  #                        nrow(psi[[i]]), ncol(psi[[i]]))
  #   }
  # }

  # Find the unique elements in all the targets:
  uniques <- unique(unlist(matrices))
  # Get the parameters (those elements that are not digits):
  z <- suppressWarnings(as.numeric(uniques))
  parameter_vector <- uniques[which(is.na(z))]
  # Get the fixed values (those elements that are digits):
  fixed_vector <- uniques[which(!is.na(z))]

  estimator_setup <- list()
  for(i in 1:ngroups) { # For each group, get the setup

    if(estimator[i] == "uls" || estimator[i] == "dwls") {
      p <- nrow(lambda[[i]])
      W <- matrix(1, nrow = p, ncol = p)
      estimator_setup[[i]] <- getargs_cfa_dwls(parameter_vector, fixed_vector,
                                               dataset[[i]], lambda[[i]], phi[[i]], psi[[i]],
                                               W = W, positive = FALSE)
    } else if(estimator[i] == "ml") {
      estimator_setup[[i]] <- getargs_cfa_ml(parameter_vector, fixed_vector,
                                             dataset[[i]], lambda[[i]], phi[[i]], psi[[i]],
                                             positive = FALSE)
    } else {
      stop("Not yet there")
    }

  }

  # if(positive) {
  #   projection = "positive"
  # } else {
  #   projection <- "id"
  # }

  indices <- unique(unlist(lapply(estimator_setup, function(x) x$indices)))
  manifolds <- c("euclidean")
  arguments <- list()
  arguments[[1]] <- list(indices = indices)
  manifold_setup <- setup_all_manifolds(manifolds, arguments)
  # manifold_setup <- latent2:::setup_all_manifolds(manifolds, arguments)

  nparam <- length(indices)
  transform_setup <- list()
  transform_setup[[1]] <- list(transform = "identity",
                               indices = 1:nparam - 1,
                               target_indices = 1:nparam - 1)

  # Collect the initital values of all the entries of lambda, phi and psi:
  if(is.null(control$init)) { # Extract the vector of parameter matrices
    init <- unlist(lapply(estimator_setup, FUN = \(x) unlist(x$matrices)))
  }

  # Get the indices of the first nonduplicated elements that are estimated:
  vmatrices <- unlist(matrices)
  nondupli <- which(!duplicated(vmatrices) & is.na(suppressWarnings(as.numeric(vmatrices))))

  # Relate these first nonduplicated elements to the values in init:
  parameters <- vector("list", length = control$rstarts)
  for(i in 1:control$rstarts) {
    parameters[[i]] <- init[nondupli]
  }

  if(control$maxit < 1) {

    x <- vector("list")
    x$model <- matrices
    x$transform_setup <- transform_setup
    x$estimator_setup <- estimator_setup
    x$manifold_setup <- manifold_setup
    x$control <- control
    x$parameters <- parameters

    return(x)

  }

  control$parameters <- parameters

  x <- optimizer(control_transform = transform_setup,
                 control_estimator = estimator_setup,
                 control_manifold = manifold_setup,
                 control_optimizer = control)

  x$init_parameters <- parameters
  x$transform_setup <- transform_setup
  x$estimator_setup <- estimator_setup
  x$manifold_setup <- manifold_setup
  x$control <- control

  outputs <- vector("list", length = ngroups)

  for(i in 1:ngroups) {

    p <- nrow(x$estimator_setup[[i]]$lambda)
    q <- ncol(x$estimator_setup[[i]]$lambda)

    # Arrange lambda parameter estimates:
    outputs[[i]]$lambda <- matrix(x$matrices[[i]][[1]], p, q)

    # Arrange phi parameter estimates:
    outputs[[i]]$phi <- matrix(x$matrices[[i]][[2]], q, q)

    # Arrange psi parameter estimates:
    outputs[[i]]$psi <- matrix(x$matrices[[i]][[3]], p, p)
    # uniquenesses_hat[[i]] <- diag(psi_hat[[i]])

    # Model matrix:
    outputs[[i]]$model <- matrix(x$matrices[[i]][[4]], p, p)
    outputs[[i]]$residuals <- matrix(x$matrices[[i]][[5]], p, p)

  }

  x$model <- matrices
  x$outputs <- outputs

  # fit <- cfa(parameters = x, X = data, nfactors = q, nobs = nobs,
  #            lambda = lambda_hat, phi = phi_hat, psi = psi_hat,
  #            lambda_indexes = indexes_lambda,
  #            phi_indexes = indexes_phi,
  #            psi_indexes = indexes_psi,
  #            target_indexes = indexes_target,
  #            targetphi_indexes = indexes_targetphi,
  #            targetpsi_indexes = indexes_targetpsi,
  #            free_indices_phi = free_indices_phi,
  #            free_indices_psi = free_indices_psi,
  #            cor = cor,
  #            estimator = estimator,
  #            projection = projection,
  #            missing = rep(missing, ngroups),
  #            se = "robust",
  #            control = control)
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
  #   rstarts <- control$rstarts
  #   cores <- control$cores
  #
  #   for(i in 1:n) {
  #     lambda <- fit$cfa$lambda[[i]]
  #     rot <- rotate(lambda = lambda, rotation = rotation, projection = projection,
  #                   gamma = gamma, epsilon = epsilon, k = k, w = w, Target = Target,
  #                   Weight = Weight, PhiTarget = PhiTarget, PhiWeight = PhiWeight,
  #                   blocks = blocks, block_weights = block_weights,
  #                   oblq_factors = oblq_factors, normalization = normalization,
  #                   rot_control = rot_control, random_starts = rstarts, cores = cores)
  #     fit$rotation[[i]] <- rot
  #   }
  # }

  return(x)

}

