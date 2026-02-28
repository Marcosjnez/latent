# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 31/01/2026
#'
#' @title
#' Rotate the lambda matrix of an orthogonal factor model.
#'
#' @usage
#'
#' lrotate(lambda, projection = "oblq", rotation = "oblimin",
#'  group = NULL, positive = FALSE, penalties = TRUE,
#'  do.fit = TRUE, control = NULL, ...)
#'
#' @param lambda List, loading matrices for each group.
#' @param projection String. Can be "orth", "oblq", or "poblq".
#' @param rotation String. Name of the variable that splits the data in different groups.
#' @param group String. Name of the variable that splits the data in different groups.
#' @param positive Force a positive-definite solution. Defaults to FALSE.
#' @param penalties list of penalty terms for the parameters.
#' @param do.fit TRUE to fit the model and FALSE to return only the model setup. Defaults to TRUE.
#' @param control List of control parameters for the optimization algorithm. See 'details' for more information.
#' @param ... Additional arguments.
#'
#' @details \code{lrotate} estimates confirmatory factor models.
#'
#' @return List with the following objects:
#' \item{version}{Version number of 'latent' when the model was estimated.}
#' \item{call}{Code used to estimate the model.}
#' \item{ModelInfo}{Model information.}
#' \item{Optim}{Output of the optimizer.}
#' \item{parameters}{Structure with all model parameters.}
#' \item{transparameters}{Structure with all transformed model parameters.}
#' \item{loglik}{Logarithm likelihood of the model.}
#' \item{penalized_loglik}{Logarithm likelihood + logarithm priors of the model.}
#'
#' @examples
#'
#' \dontrun{
#'
#' fit <- lrotate(lambda = , projection = "oblq", rotation = "oblimin")
#' summary(fit, digits = 3L)
#'}
#'
#' @export
lrotate <- function(lambda, projection = "oblq", rotation = "oblimin",
                    group = NULL, positive = FALSE, penalties = TRUE,
                    do.fit = TRUE, control = NULL, ...) {

  # Capture everything in ... as a named list
  dots <- list(...)

  #### Arrange the data ####

  # Check orthogonality
  if(projection == "oblq" || projection == "poblq") {
    orthogonal <- FALSE
  } else if(projection == "orth") {
    orthogonal <- TRUE
  } else {
    stop("Unknown projection")
  }

  # Check the arguments to control_optimizer and create defaults:
  estimator <- rotation
  control$penalties <- penalties
  control$estimator <- tolower(rotation)
  control <- lcfa_control(control)
  if(control$opt == "lbfgs") control$opt <- "newton"

  # Data and structure information:
  ngroups        <- length(lambda)
  group_label    <- names(lambda)
  item_label   <- lapply(lambda, FUN = \(x) rownames(x))
  factor_label   <- lapply(lambda, FUN = \(x) colnames(x))
  nitems <- lapply(lambda, FUN = nrow)
  nfactors <- lapply(lambda, FUN = ncol)

  data_list <- vector("list")
  data_list$ngroups <- ngroups
  data_list$data <- lambda
  data_list$nitems <- nitems
  data_list$nfactors <- nfactors
  data_list$positive <- positive
  data_list$estimator <- estimator
  data_list$group_label <- group_label
  data_list$item_label <- item_label
  data_list$factor_label <- factor_label
  data_list$orthogonal <- orthogonal

  ## store original call
  mc  <- match.call()

  #### Create the model ####

  # # Generate the model syntax and initial parameter values
  #
  # list2env(data_list, envir = environment())

  # Initialize the objects to store the initial parameters:
  rot_param <- rot_trans <- model <- vector("list", length = ngroups)
  fixed <- fixed_values <- nonfixed <- vector("list", length = ngroups)

  for(i in 1:ngroups) {

    model[[i]]$ulambda <- lambda[[i]]
    X_labels <- paste("g", i, ".X[", rep(1:nfactors[[i]], times = nfactors[[i]]),
                      ",", rep(1:nfactors[[i]], each = nfactors[[i]]), "]", sep = "")
    model[[i]]$X <- matrix(X_labels, nrow = nfactors[[i]], ncol = nfactors[[i]])

    # Transformed parameters:

    # Get the positions of parameters and fixed values:

    nonfixed[[i]] <- lapply(model[[i]], FUN = \(x) {
      which(is.na(suppressWarnings(as.numeric(x))))
    })

    fixed[[i]] <- lapply(model[[i]], FUN = \(x) {
      which(!is.na(suppressWarnings(as.numeric(x))))
    })

    fixed_values[[i]] <- lapply(model[[i]], FUN = \(x) {
      numerals <- suppressWarnings(as.numeric(x))
      inds <- which(!is.na(numerals))
      return(numerals[inds])
    })

    # Unrotated lambda:
    ulambda_labels <- paste("g", i, ".ulambda[", rep(1:nitems[[i]], times = nfactors[[i]]),
                            ",", rep(1:nfactors[[i]], each = nitems[[i]]), "]", sep = "")
    rot_trans[[i]]$ulambda <- matrix(ulambda_labels, nrow = nitems[[i]], ncol = nfactors[[i]])

    # X:
    rot_trans[[i]]$X <- model[[i]]$X

    if(!orthogonal) {
      # Xinv:
      Xinv_labels <- paste("g", i, ".Xinv[", rep(1:nfactors[[i]], times = nfactors[[i]]),
                           ",", rep(1:nfactors[[i]], each = nfactors[[i]]), "]", sep = "")
      rot_trans[[i]]$Xinv <- matrix(Xinv_labels, nrow = nfactors[[i]], ncol = nfactors[[i]])
    }

    # Rotated lambda:
    lambda_labels <- paste("g", i, ".lambda[", rep(1:nitems[[i]], times = nfactors[[i]]),
                           ",", rep(1:nfactors[[i]], each = nitems[[i]]), "]", sep = "")
    rot_trans[[i]]$lambda <- matrix(lambda_labels, nrow = nitems[[i]], ncol = nfactors[[i]])

    # Latent correlations:
    psi_labels <- paste("g", i, ".psi[", rep(1:nfactors[[i]], times = nfactors[[i]]),
                           ",", rep(1:nfactors[[i]], each = nfactors[[i]]), "]", sep = "")
    rot_trans[[i]]$psi <- matrix(psi_labels, nrow = nfactors[[i]], ncol = nfactors[[i]])

    # Untransformed parameters:

    rot_param[[i]]$ulambda <- model[[i]]$ulambda
    rot_param[[i]]$X <- rot_trans[[i]]$X

  }

  #### Arrange labels ####

  # Arrange parameter labels:
  vector_param <- unname(unique(unlist(rot_param)))

  # Select the unique, nonnumeric labels:
  nonfixed_pars <- which(is.na(suppressWarnings(as.numeric(vector_param))))
  parameters_labels <- vector_param[nonfixed_pars]
  nparam <- length(parameters_labels)

  # Arrange transparameter labels:
  vector_trans <- unname(unlist(rot_trans))
  transparameters_labels <- unique(vector_trans)
  ntrans <- length(transparameters_labels)

  #### Relate the transformed parameters to the parameters ####

  param2trans <- match(transparameters_labels, parameters_labels)
  param2trans <- param2trans[!is.na(param2trans)]
  # Relate the parameters to the transformed parameters:
  trans2param <- match(parameters_labels, transparameters_labels)

  #### Create the initial values for the parameters ####

  # Collect the unique nontransformed parameters and the unique transformed parameters:

  init_trans <- vector("list", length = control$rstarts)

  for(rs in 1:control$rstarts) {

    init_trans[[rs]] <- vector("list", length = ngroups)

    for(i in 1:ngroups) {

      init_trans[[rs]][[i]]$ulambda <- lambda[[i]]
      init_trans[[rs]][[i]]$X <- rorth(nfactors[[i]], nfactors[[i]])
      # init_trans[[rs]][[i]]$X <- roblq(nfactors[[i]], nfactors[[i]])

      if(orthogonal) {
        init_trans[[rs]][[i]]$lambda <- lambda[[i]] %*%
          init_trans[[rs]][[i]]$X
      } else {
        init_trans[[rs]][[i]]$Xinv <- solve(init_trans[[rs]][[i]]$X)
        init_trans[[rs]][[i]]$lambda <- lambda[[i]] %*%
          t(init_trans[[rs]][[i]]$Xinv)
      }

      init_trans[[rs]][[i]]$psi <- crossprod(init_trans[[rs]][[i]]$X)

    }

  }

  #### Create the vectors of parameters and transformed parameters ####

  parameters <- transparameters <- vector("list", length = control$rstarts)
  # Indices of the unique transparameters in init_trans:
  trans_inds <- match(transparameters_labels, vector_trans)
  init_inds <- match(parameters_labels, vector_trans)

  for(rs in 1:control$rstarts) {

    transparameters[[rs]] <- unlist(init_trans[[rs]])[trans_inds]
    names(transparameters[[rs]]) <- transparameters_labels
    parameters[[rs]] <- unlist(init_trans[[rs]])[init_inds]
    names(parameters[[rs]]) <- parameters_labels

  }

  #### Set up the optimizer ####

  # Create defaults for the control of the optimizer:
  control$parameters <- parameters
  control$transparameters <- transparameters
  control$param2transparam <- param2trans-1L
  control$transparam2param <- trans2param-1L

  full_model <- list(parameters_labels = parameters_labels,
                     nparam = nparam,
                     transparameters_labels = transparameters_labels,
                     ntrans = ntrans,
                     rot_param = rot_param,
                     rot_trans = rot_trans,
                     fixed = fixed,
                     nonfixed = nonfixed,
                     init_trans = init_trans,
                     control = control)

  #### Create the structures ####

  #### Manifolds ####

  mani_and_labs <- list()
  k <- 1L

  for(i in 1:ngroups) {

    dots$p <- nfactors[[i]]
    dots$q <- nfactors[[i]]
    # Get the extra objects required for the manifold:
    mani_and_labs[[k]] <- extra_manifolds(projection,
                                          rot_param[[i]]$X,
                                          dots)

    k <- k+1L

  }

  control_manifold <- create_manifolds(manifolds_and_labels = mani_and_labs,
                                       param_structures = rot_param)

  #### Transformations ####

  trans_and_labs <- list()
  k <- 1L

  for(i in 1:ngroups) {

    if(orthogonal) {

      # Rotated lambda:

      # Get the extra objects required for the transformation:
      dots$p <- nitems[[i]]
      dots$q <- nfactors[[i]]
      trans_and_labs[[k]] <- extra_transforms(transform = "XY",
                                              labels_in = list(rot_trans[[i]]$ulambda,
                                                               rot_trans[[i]]$X),
                                              labels_out = list(rot_trans[[i]]$lambda),
                                              dots)
      k <- k+1L

    } else {

      # Inverse of X:

      # Get the extra objects required for the transformation:
      dots$p <- nfactors[[i]]
      trans_and_labs[[k]] <- extra_transforms(transform = "matrix_inverse",
                                              labels_in = list(rot_trans[[i]]$X),
                                              labels_out = list(rot_trans[[i]]$Xinv),
                                              dots)
      k <- k+1L

      # Rotated lambda:

      # Get the extra objects required for the transformation:
      dots$p <- nitems[[i]]
      dots$q <- nfactors[[i]]
      trans_and_labs[[k]] <- extra_transforms("XYt",
                                              labels_in = list(rot_trans[[i]]$ulambda,
                                                               rot_trans[[i]]$Xinv),
                                              labels_out = list(rot_trans[[i]]$lambda),
                                              dots)
      k <- k+1L

    }

    # Latent covariances:

    # Get the extra objects required for the transformation:
    dots$p <- nfactors[[i]]
    trans_and_labs[[k]] <- extra_transforms("crossprod",
                                            labels_in = list(rot_trans[[i]]$X),
                                            labels_out = list(rot_trans[[i]]$psi[lower.tri(rot_trans[[i]]$psi, diag = TRUE)]), # LOWER DIAG
                                            dots)
    k <- k+1L

  }

  control_transform <- create_transforms(transforms_and_labels = trans_and_labs,
                                         param_structures = rot_trans)

  #### Estimators ####

  control_estimator <- list()
  k <- 1L

  for(i in 1:ngroups) {

    # Rotation criteria:
    p <- nitems[[i]]
    q <- nfactors[[i]]
    lower_psi <- lower.tri(diag(q, q), diag = TRUE)

    lambda_labels <- c(rot_trans[[i]]$lambda)
    psi_labels <- rot_trans[[i]]$psi[lower_psi]
    indices <- list(match(lambda_labels, transparameters_labels)-1L,
                    match(psi_labels, transparameters_labels)-1L)
    control_estimator[[k]] <- list(estimator = rotation,
                                   indices = indices,
                                   p = p,
                                   q = q,
                                   ...)
    k <- k+1L

  }

  #### Structures ####

  structures <- list(control_manifold = control_manifold,
                     control_transform = control_transform,
                     control_estimator = control_estimator)

  #### Collect all the model information ####

  # Model information:
  modelInfo <- list(
                    # nobs = nobs,
                    # nparam = nparam - rest,
                    # npatterns = npatterns,
                    # dof = sum(unlist(npatterns)) - nparam + rest,
                    ntrans = ntrans,
                    parameters_labels = parameters_labels,
                    transparameters_labels = transparameters_labels,
                    rot_param = rot_param,
                    rot_trans = rot_trans,
                    rotation = rotation,
                    control_manifold = control_manifold,
                    control_transform = control_transform,
                    control_estimator = control_estimator,
                    control = control)

  #### Fit the model ####

  if(!do.fit) {

    lcfa_list <- new("lcfa",
                     version            = as.character( packageVersion('latent') ),
                     call               = mc,
                     timing             = numeric(),
                     data_list          = data_list,
                     modelInfo          = modelInfo,
                     Optim              = list(),
                     parameters         = list(),
                     transformed_pars   = list(),
                     loglik             = numeric(),
                     penalized_loglik   = numeric(),
                     loss               = numeric(),
                     penalized_loss     = numeric()
    )

    return(lcfa_list)

  }

  control$cores <- min(control$rstarts, control$cores)
  # Fit the model:
  x <- optimizer(control_manifold = control_manifold,
                 control_transform = control_transform,
                 control_estimator = control_estimator,
                 control_optimizer = control)

  # Collect all the information about the optimization:

  Optim <- x
  elapsed <- x$elapsed

  #### Estimated model structures ####

  # Create the structures of untransformed parameters:
  indices_pars <- match(modelInfo$parameters_labels,
                        unlist(modelInfo$rot_param))

  vv <- rep(0, times = length(unlist(modelInfo$rot_param)))
  vv[indices_pars] <- Optim$parameters
  parameters <- fill_list_with_vector(modelInfo$rot_param, vv)
  parameters <- allnumeric(parameters)
  # FIXED PARAMETERS?

  # Create the structures of transformed parameters:
  indices_trans <- match(modelInfo$transparameters_labels,
                         unlist(modelInfo$rot_trans))
  vv <- rep(0, times = length(unlist(modelInfo$rot_trans)))
  vv[indices_trans] <- Optim$transparameters
  transformed_pars <- fill_list_with_vector(modelInfo$rot_trans, vv)
  transformed_pars <- allnumeric(transformed_pars)

  #### Process the fit information ####

  # Initialize the objects to be returned:
  loss <- penalized_loss <- loglik <- penalized_loglik <- penalty <-
    vector("list", length = ngroups)

  # For each group, extract the loss, penalized loss, loglik and penalized loglik
  for(i in 1:ngroups) {

    loss[[i]] <- c(x$outputs$estimators$doubles[[1]][[1]])
    loglik[[i]] <- c(x$outputs$estimators$doubles[[1]][[2]])
    penalized_loss[[i]] <- loss[[i]]
    penalized_loglik[[i]] <- loglik[[i]]

  }

  loss <- sum(unlist(loss))
  penalized_loss <- sum(unlist(penalized_loss))
  loglik <- sum(unlist(loglik))
  penalized_loglik <- sum(unlist(penalized_loglik))

  #### Result ####

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

