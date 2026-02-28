# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 22/02/2026

pois_reg <- function(Y, X = NULL, penalties = FALSE, do.fit = TRUE, control = NULL) {

  ## store original call
  mc  <- match.call()

  control$penalties <- penalties
  control_optimizer <- lcfa_control(control)

  #### Initial input checks ####

  nobs <- nrow(Y)
  q <- ncol(Y)

  # Process the covariates:
  if(is.null(X)) {

    X <- matrix(1, nrow = nobs, ncol = 1L)
    colnames(X) <- "(Intercept)"

  } else {

    # Check that X is either a data.frame or a matrix:
    if(!is.data.frame(X) & !is.matrix(X)) {
      stop("data must be a matrix or data.frame")
    }

    if(nrow(X) != nobs) {
      stop("Number of cases in the data and covariates do not match")
    }

    # Transform characters into factors:
    X_df <- as.data.frame(X)
    X_df[] <- lapply(X_df, function(x) if (is.character(x)) factor(x) else x)

    # Create the design matrix:
    X <- model.matrix(~ . + 1, X_df)
    # Center the variables:
    # X[, -1] <- apply(X[, -1], MARGIN = 2, FUN = \(x) x-mean(x))

    # Put an underscore between the variable names and their level names:
    for (v in names(X_df)[sapply(X_df, is.factor)]) {
      i <- which(startsWith(colnames(X), v))
      colnames(X)[i] <- paste0(v, "_", make.names(levels(X_df[[v]]))[seq_along(i)])
    }

  }

  p <- ncol(X)

  # Combine the data and the covariates:
  dt <- cbind(Y, X)
  pcov <- ncol(X) # Number of covariates

  # Convert data to a data.table object:
  dt <- data.table::as.data.table(dt)

  ## Collect some information from the data ##
  counts_dt <- dt[, .(index = .I[1], count = .N), by = names(dt)]
  # Data matrix with the unique response patterns:
  patterns <- as.matrix(counts_dt[, names(dt), with = FALSE])[, -(1:pcov) , drop = FALSE]
  # Covariates with unique response patterns:
  cov_patterns2 <- as.matrix(counts_dt[, names(dt), with = FALSE])[, 1:pcov , drop = FALSE]
  # Number of unique response patterns:
  npatterns <- nrow(counts_dt)
  # Counts of each response pattern:
  weights <- counts_dt$count
  # Indices to map the original data to the matrix of unique patterns:
  full2short <- counts_dt$index
  # Indices to map the matrix of unique patterns to the original data:
  short2full <- match(do.call(paste, dt),
                      do.call(paste, counts_dt[, -c("index", "count"), with = FALSE]))

  ## Collect some information from the covariates ##
  cov_dt <- data.table::as.data.table(X)
  counts_cov_dt <- cov_dt[, .(index = .I[1], count = .N), by = names(cov_dt)]
  # Covariates with unique patterns:
  cov_patterns <- as.matrix(counts_cov_dt[, names(cov_dt), with = FALSE])
  # Number of unique covariate patterns:
  ncov_patterns <- nrow(cov_patterns)
  # Indices to map the original data to the matrix of unique patterns:
  cov_full2short <- counts_cov_dt$index
  # Indices to map the matrix of unique patterns to the original data:
  cov_short2full <- match(do.call(paste, cov_dt),
                          do.call(paste, counts_cov_dt[, -c("index", "count"),
                                                       with = FALSE]))

  #### Data list ####

  # Put in a list the objects generated form the data:
  data_list <- vector("list")
  # data_list$data <- data
  data_list$Y <- Y
  data_list$X <- X
  data_list$nobs <- nobs
  data_list$patterns <- patterns
  data_list$cov_patterns2 <- cov_patterns2
  data_list$npatterns <- npatterns
  data_list$weights <- weights
  data_list$full2short <- full2short
  data_list$short2full <- short2full
  data_list$cov_patterns <- cov_patterns
  data_list$ncov_patterns <- ncov_patterns
  data_list$cov_full2short <- cov_full2short
  data_list$cov_short2full <- cov_short2full

  #### Create the model ####

  param <- trans <- list()

  labels_beta <- paste("beta", rep(1:p, times = q), "|",
                  rep(1:q, each = p), sep = "")
  param$beta <- matrix(labels_beta, nrow = p, ncol = q)
  rownames(param$beta) <- colnames(X)
  colnames(param$beta) <- paste("Y", 1:q, sep = "")
  trans$beta <- param$beta

  labels_linpred <- paste("linpred", rep(1:nobs, times = q), "|",
                          rep(1:q, each = nobs), sep = "")
  trans$linpred <- matrix(labels_linpred, nrow = nobs, ncol = q)
  colnames(trans$linpred) <- paste("Y", 1:q, sep = "")

  labels_lambda <- paste("lambda", rep(1:nobs, times = q), "|",
                       rep(1:q, each = nobs), sep = "")
  trans$lambda <- matrix(labels_lambda, nrow = nobs, ncol = q)
  colnames(trans$lambda) <- paste("Y", 1:q, sep = "")

  #### Arrange labels ####

  # Arrange parameter labels:
  vector_param <- unname(unique(unlist(param)))

  # Select the unique, nonnumeric labels:
  nonfixed_pars <- which(is.na(suppressWarnings(as.numeric(vector_param))))
  parameters_labels <- vector_param[nonfixed_pars]
  nparam <- length(parameters_labels)

  # Arrange transparameter labels:
  vector_trans <- unname(unlist(trans))
  transparameters_labels <- unique(vector_trans)
  ntrans <- length(transparameters_labels)

  #### Relate the transformed parameters to the parameters ####

  param2trans <- match(transparameters_labels, parameters_labels)
  param2trans <- param2trans[!is.na(param2trans)]
  # Relate the parameters to the transformed parameters:
  trans2param <- match(parameters_labels, transparameters_labels)

  #### Create the initial values for the parameters ####

  init_trans <- vector("list")
  init_trans[[1]] <- trans
  init_trans[[1]]$beta <- matrix(1, nrow = p, ncol = q)
  init_trans[[1]]$linpred <- X %*% init_trans[[1]]$beta
  init_trans[[1]]$lambda <- exp(init_trans[[1]]$linpred)

  #### Create the vectors of parameters and transformed parameters ####

  parameters <- transparameters <- vector("list",
                                          length = control_optimizer$rstarts)
  # Indices of the unique transparameters in init_trans:
  trans_inds <- match(transparameters_labels, vector_trans)
  init_inds <- match(parameters_labels, vector_trans)

  for(rs in 1:control_optimizer$rstarts) {

    transparameters[[rs]] <- unlist(init_trans[[rs]])[trans_inds]
    names(transparameters[[rs]]) <- transparameters_labels
    parameters[[rs]] <- unlist(init_trans[[rs]])[init_inds]
    names(parameters[[rs]]) <- parameters_labels

  }

  #### Set up the optimizer ####

  # Create defaults for the control of the optimizer:
  control_optimizer$parameters <- parameters
  control_optimizer$transparameters <- transparameters
  control_optimizer$param2transparam <- param2trans-1L
  control_optimizer$transparam2param <- trans2param-1L

  #### Manifolds ####

  mani_and_labs <- list(
    list("euclidean", param)
  )
  control_manifold <- create_manifolds(manifolds_and_labels = mani_and_labs,
                                       param_structures = param)

  #### Transformations ####

  trans_and_labs <- list()
  dots <- list()

  # betas to linpreds:
  dots$X <- X
  trans_and_labs[[1]] <- extra_transforms(transform = "column_space",
                                          labels_in = list(trans$beta),
                                          labels_out = list(trans$linpred),
                                          dots)

  # linpred to lambdas:
  trans_and_labs[[2]] <- extra_transforms(transform = "exponential",
                                          labels_in = list(trans$linpred),
                                          labels_out = list(trans$lambda),
                                          dots)

  control_transform <- create_transforms(transforms_and_labels = trans_and_labs,
                                         param_structures = trans)

  #### Estimators ####

  control_estimator <- list()

  indices <- list(match(trans$lambda, transparameters_labels)-1L)
  control_estimator[[1]] <- list(estimator = "poisson_loglik",
                                 indices = indices,
                                 X = Y)

  #### Collect all the model information ####

  # Model information:
  modelInfo <- list(
    nobs = nobs,
    nparam = nparam,
    npatterns = npatterns,
    dof = sum(unlist(npatterns)) - nparam,
    ntrans = ntrans,
    parameters_labels = parameters_labels,
    transparameters_labels = transparameters_labels,
    param = param,
    trans = trans,
    control_manifold = control_manifold,
    control_transform = control_transform,
    control_estimator = control_estimator,
    control = control_optimizer)

  #### Fit ####

  control_optimizer$cores <- min(control_optimizer$rstarts,
                                 control_optimizer$cores)
  # Fit the model:
  Optim <- optimizer(control_manifold = control_manifold,
                     control_transform = control_transform,
                     control_estimator = control_estimator,
                     control_optimizer = control_optimizer)

  #### Estimated model structures ####

  # Create the structures of untransformed parameters:
  indices_pars <- match(modelInfo$parameters_labels,
                        unlist(modelInfo$param))

  vv <- rep(0, times = length(unlist(modelInfo$param)))
  vv[indices_pars] <- Optim$parameters
  parameters <- fill_list_with_vector(modelInfo$param, vv)
  parameters <- allnumeric(parameters)
  # FIXED PARAMETERS?

  # Create the structures of transformed parameters:
  indices_trans <- match(modelInfo$transparameters_labels,
                         unlist(modelInfo$trans))
  vv <- rep(0, times = length(unlist(modelInfo$trans)))
  vv[indices_trans] <- Optim$transparameters
  transformed_pars <- fill_list_with_vector(modelInfo$trans, vv)
  transformed_pars <- allnumeric(transformed_pars)

  #### Process the fit information ####

  loss <- -Optim$outputs$estimators$doubles[[1]][[1]]
  penalized_loss <- -Optim$outputs$estimators$doubles[[1]][[1]]
  loglik <- Optim$outputs$estimators$doubles[[1]][[1]]
  penalized_loglik <- Optim$outputs$estimators$doubles[[1]][[1]]

  #### Result ####

  elapsed <- Optim$elapsed
  result <- new("lcfa",
                version            = as.character( packageVersion('latent') ),
                call               = mc,
                timing             = elapsed,
                data_list          = data_list,
                modelInfo          = modelInfo,
                Optim              = Optim,
                parameters         = parameters,
                transformed_pars   = transformed_pars,
                loglik             = loglik,
                penalized_loglik   = penalized_loglik,
                loss               = loss,
                penalized_loss     = penalized_loss
  )

  return(result)

}

