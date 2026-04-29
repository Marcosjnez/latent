# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 29/04/2026

lyule <- function(data,
                  model = NULL,
                  do.fit = TRUE,
                  control = NULL,
                  ...) {

  ## store original call
  mc  <- match.call()

  # Check the arguments to control_optimizer and create defaults:
  control <- lyule_control(control)

  #### Create the data_list ####

  data_list <- create_lyule_datalist(data, control)
  list2env(data_list, envir = environment())

  #### Create the model ####

  full_model <- create_lyule_model(data_list = data_list,
                                   model = model,
                                   control = control)
  list2env(full_model, envir = environment())

  #### Create the modelInfo ####

  modelInfo <- create_lyule_modelInfo(data_list = data_list,
                                      full_model = full_model,
                                      control = control)

  #### Fit the model ####

  if(!do.fit) {

    lcfa_list <- new("lcfa",
                     version            = as.character(packageVersion('latent')),
                     call               = mc, # matched call
                     timing             = numeric(), # timing information
                     data_list          = data_list,
                     modelInfo          = modelInfo,
                     Optim              = list(),
                     parameters         = list(),
                     transformed_pars   = list(),
                     loglik             = numeric(), # loglik values
                     penalized_loglik   = numeric(),
                     loss               = numeric(),
                     penalized_loss     = numeric()
    )

    return(lcfa_list)

  }

  est_and_se <- yule_cor_full_rcpp(data)
  S <- est_and_se$cor
  rownames(S) <- colnames(S) <- item_names

  Optim <- list()
  Optim$f <- 0

  Optim$parameters <- S[lower.tri(S, diag = FALSE)]
  Optim$transparameters <- S[lower.tri(S, diag = TRUE)]
  names(Optim$parameters) <- modelInfo$parameters_labels
  names(Optim$transparameters) <- modelInfo$transparameters_labels

  #### Standard errors ####

  lower_idx <- lower.tri(est_and_se$cor, diag = FALSE)
  Optim$SE$ACOV <- diag(est_and_se$se[lower_idx]^2)*data_list$nobs

  rownames(Optim$SE$ACOV) <- colnames(Optim$SE$ACOV) <-
    modelInfo$parameters_labels
  Optim$SE$se <- sqrt(diag(Optim$SE$ACOV) / data_list$nobs)

  # Collect all the information about the optimization:

  elapsed <- 0

  #### Estimated model structures ####

  # Create the structures of transformed parameters:
  transformed_pars <- fill_in(modelInfo$trans, Optim$transparameters)

  # Create the structures of untransformed parameters:
  parameters <- transformed_pars[names(modelInfo$param)]

  #### Process the fit information ####

  loss <- Optim$f
  penalized_loss <- loss
  loglik <- sum(unlist(lapply(Optim$outputs$estimators$doubles,
                              FUN = \(x) x[[2]])))
  penalized_loglik <- loglik

  #### latent object ####

  result <- new("latent",
                version            = as.character(packageVersion('latent')),
                call               = mc,
                timing             = elapsed,
                dataList           = data_list,
                modelInfo          = modelInfo,
                Optim              = Optim,
                parameters         = parameters,
                transformed_pars   = transformed_pars,
                loglik             = loglik,
                penalized_loglik   = penalized_loglik,
                loss               = loss,
                penalized_loss     = penalized_loss,
                extra              = list()
  )

  #### Return ####

  return(result)

}

create_lyule_datalist <- function(data, control) {

  data_list <- vector("list")
  data_list$data <- data
  data_list$nobs <- nrow(data)
  data_list$nitems <- ncol(data)
  data_list$npatterns <- nrow(data)
  data_list$item_names <- colnames(data)

  return(data_list)

}

create_lyule_model <- function(data_list, model, control) {

  # Generate the model syntax and initial parameter values

  list2env(data_list, envir = environment())
  list2env(control, envir = environment())

  # Initialize the objects to store the initial parameters:
  param <- trans <- vector("list")

  S_matrix <- paste("S.", control$subfix, sep = "")

  #### Model for the transformed parameters ####

  # Transformed parameters:
  list_struct <- vector("list")
  k <- 1L

  # Covariance matrix:
  list_struct[[k]] <- list(name = S_matrix,
                           type = "matrix",
                           dim = c(nitems, nitems),
                           rownames = item_names,
                           colnames = item_names,
                           symmetric = TRUE)
  k <- k+1L

  trans <- create_parameters(list_struct)

  #### Model for the parameters ####

  param <- trans
  diag(param[[S_matrix]]) <- 1

  #### Create the initial values for the parameters ####

  init_param <- vector("list", length = control$rstarts)
  for(rs in 1:control$rstarts) {

    init_param[[rs]] <- vector("list")
    init_param[[rs]][[S_matrix]] <- diag(1, nrow = data_list$nitems)

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
                 init_param = init_param)

  return(result)

}

create_lyule_modelInfo <- function(data_list, full_model, control) {

  # Generate control_manifold, control_transform, and control_estimator

  list2env(data_list, envir = environment())
  list2env(full_model, envir = environment())
  list2env(control, envir = environment())

  S_matrix <- paste("S.", control$subfix, sep = "")

  #### Manifolds ####

  manifolds <- list(
    list(manifold = "euclidean",
         parameters = S_matrix)
  )

  control_manifold <- create_manifolds(manifolds = manifolds,
                                       structures = param)

  #### Transformations ####

  transforms <- list()

  control_transform <- create_transforms(transforms = transforms,
                                         structures = trans)

  #### Estimators ####

  estimators <- list()

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
                    nparam = nparam,
                    ntrans = ntrans,
                    parameters_labels = parameters_labels,
                    transparameters_labels = transparameters_labels,
                    dof = npatterns - nparam,
                    control_manifold = control_manifold,
                    control_transform = control_transform,
                    control_estimator = control_estimator,
                    control_optimizer = control_optimizer)

  #### Return ####

  return(modelInfo)

}

lyule_control <- function(control) {

  if(is.null(control$step_maxit)) {
    control$step_maxit <- 30L
  }

  if(is.null(control$c1)) {
    control$c1 <- 0.5
  }

  if(is.null(control$c2)) {
    control$c2 <- 0.5
  }

  if(is.null(control$step_eps)) {
    control$step_eps <- 1e-09
  }

  if(is.null(control$df_eps)) {
    control$df_eps <- 1e-09
  }

  if(is.null(control$M)) {
    control$M <- 100L
  }

  if(is.null(control$eps)) {
    control$eps <- 1e-04
  }

  if(is.null(control$ss_fac)) {
    control$ss_fac <- 2
  }

  if(is.null(control$maxit)) {
    control$maxit <- 1000L
  }

  if(is.null(control$rstarts)) {
    control$rstarts <- 1L
  }

  if(is.null(control$cores)) {
    control$cores <- 1L
  }

  if(is.null(control$tcg_maxit)) {
    control$tcg_maxit <- 10L
  }

  if(is.null(control$ss)) {
    # Step sizes should be small so taus are not far from sensible bounds
    control$ss <- 0.001
  }

  # Fixed rstarts and cores:
  control$rstarts <- 1L
  control$cores <- 1L

  return(control)

}
