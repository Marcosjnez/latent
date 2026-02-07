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

  # Check orthogonality
  if(projection == "oblq") {
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

  # # Generate control_manifold, control_transform, and control_estimator
  #
  # list2env(data_list, envir = environment())
  # list2env(full_model, envir = environment())

  #### Create the structures ####

  #### Manifolds ####

  control_manifold <- list()

  # # Euclidean for unrotated lambdas:
  # ulambdas <- unlist(lapply(rot_param, FUN = \(x) x$ulambda))
  # indices_ulambda <- match(unique(ulambdas), parameters_labels)
  # indices_ulambda <- indices_ulambda[!is.na(indices_ulambda)]
  # labels <- parameters_labels[indices_ulambda]
  # indices <- list(indices_ulambda-1L)
  # control_manifold[[1]] <- list(manifold = "euclidean",
  #                               parameters = labels,
  #                               indices = indices)
  # k <- 2L

  k <- 1L

  for(i in 1:ngroups) {

    if(orthogonal) {
      manifoldX <- "orth"
    } else {
      manifoldX <- "oblq"
    }

    # Orthogonal/Oblique X in each group:
    Xs <- c(rot_param[[i]]$X)
    indices_Xs <- match(unique(Xs), parameters_labels)
    indices_Xs <- indices_Xs[!is.na(indices_Xs)]
    labels <- parameters_labels[indices_Xs]
    indices <- list(indices_Xs-1L)
    control_manifold[[k]] <- list(manifold = manifoldX,
                                  parameters = labels,
                                  indices = indices,
                                  q = nfactors[[i]])
    k <- k+1L

  }

  #### Transformations ####

  control_transform <- list()
  k <- 1L

  for(i in 1:ngroups) {

    if(orthogonal) {

      # Rotated lambda:
      labels_ulambda <- c(rot_trans[[i]]$ulambda)
      labels_X <- c(rot_trans[[i]]$X)
      labels_in <- list(labels_ulambda, labels_X)
      labels_out <- c(rot_trans[[i]]$lambda)
      indices_in_ulambda <- match(labels_ulambda, transparameters_labels)-1L
      indices_in_X <- match(labels_X, transparameters_labels)-1L
      indices_in <- list(indices_in_ulambda, indices_in_X)
      indices_out <- list(match(labels_out, transparameters_labels)-1L)
      control_transform[[k]] <- list(transform = "XY",
                                     labels_in = labels_in,
                                     indices_in = indices_in,
                                     labels_out = labels_out,
                                     indices_out = indices_out,
                                     p = nitems[[i]],
                                     q = nfactors[[i]])
      k <- k+1L

    } else {

      # Inverse of X:
      labels_in <- c(rot_trans[[i]]$X)
      labels_out <- c(rot_trans[[i]]$Xinv)
      indices_in <- list(match(labels_in, transparameters_labels)-1L)
      indices_out <- list(match(labels_out, transparameters_labels)-1L)
      control_transform[[k]] <- list(transform = "matrix_inverse",
                                     labels_in = labels_in,
                                     indices_in = indices_in,
                                     labels_out = labels_out,
                                     indices_out = indices_out,
                                     p = nfactors[[i]])
      k <- k+1L

      # Rotated lambda:
      labels_ulambda <- c(rot_trans[[i]]$ulambda)
      labels_Xinv <- c(rot_trans[[i]]$Xinv)
      labels_in <- list(labels_ulambda, labels_Xinv)
      labels_out <- c(rot_trans[[i]]$lambda)
      indices_in_ulambda <- match(labels_ulambda, transparameters_labels)-1L
      indices_in_Xinv <- match(labels_Xinv, transparameters_labels)-1L
      indices_in <- list(indices_in_ulambda, indices_in_Xinv)
      indices_out <- list(match(labels_out, transparameters_labels)-1L)
      control_transform[[k]] <- list(transform = "XYt",
                                     labels_in = labels_in,
                                     indices_in = indices_in,
                                     labels_out = labels_out,
                                     indices_out = indices_out,
                                     p = nitems[[i]],
                                     q = nfactors[[i]])
      k <- k+1L

    }

    # # Latent covariances:
    lower_psi <- lower.tri(rot_trans[[i]]$psi, diag = TRUE)
    labels_in <- c(rot_trans[[i]]$X)
    indices_in <- list(match(labels_in, transparameters_labels)-1L)
    labels_out <- c(rot_trans[[i]]$psi[lower_psi])
    indices_out <- list(match(labels_out, transparameters_labels)-1L)
    control_transform[[k]] <- list(transform = "crossprod",
                                   labels_in = labels_in,
                                   indices_in = indices_in,
                                   labels_out = labels_out,
                                   indices_out = indices_out,
                                   p = nfactors[[i]],
                                   q = nfactors[[i]])
    k <- k+1L

  }

  #### Estimators ####

  control_estimator <- list()
  k <- 1L

  for(i in 1:ngroups) {

    # Rotation criteria:
    p <- nitems[[i]]
    q <- nfactors[[i]]
    lower_psi <- lower.tri(diag(q, q), diag = TRUE)

    # if(is.null(gamma)) {
    #   gamma <- 0.00
    # }
    #
    # if(is.null(epsilon)) {
    #   epsilon <- 0.01
    # }

    # if(rotation == "oblimin") {

      lambda_labels <- c(rot_trans[[i]]$lambda)
      psi_labels <- rot_trans[[i]]$psi[lower_psi]
      indices <- list(match(lambda_labels, transparameters_labels) - 1L,
                      match(psi_labels, transparameters_labels) - 1L)
      control_estimator[[k]] <- list(estimator = rotation,
                                     labels = labels,
                                     indices = indices,
                                     p = p,
                                     q = q,
                                     ...)
      k <- k+1L

    # } else if(rotation == "geomin") {
    #
    #   labels <- c(rot_trans[[i]]$lambda)
    #   indices <- list(match(labels, transparameters_labels) - 1L)
    #   control_estimator[[k]] <- list(estimator = "geomin",
    #                                  labels = labels,
    #                                  indices = indices,
    #                                  epsilon = 0.01,
    #                                  p = nitems[[i]],
    #                                  q = nfactors[[i]])
    #   k <- k+1L
    #
    # } else if(rotation == "target") {
    #
    #   labels <- c(rot_trans[[i]]$lambda)
    #   indices <- list(match(labels, transparameters_labels) - 1L)
    #   control_estimator[[k]] <- list(estimator = "target",
    #                                  labels = labels,
    #                                  indices = indices,
    #                                  target = target,
    #                                  weight = weight)
    #   k <- k+1L
    #
    # } else if(rotation == "xtarget") {
    #
    #   labels <- c(rot_trans[[i]]$lambda)
    #   indices <- list(match(labels, transparameters_labels) - 1L)
    #   control_estimator[[k]] <- list(estimator = "xtarget",
    #                                  labels = labels,
    #                                  indices = indices,
    #                                  target = target,
    #                                  weight = weight,
    #                                  psitarget = psitarget,
    #                                  psiweight = psiweight,
    #                                  w = w)
    #   k <- k+1L
    #
    # } else if(rotation == "varimax") {
    #
    #   labels <- c(rot_trans[[i]]$lambda)
    #   indices <- list(match(labels, transparameters_labels) - 1L)
    #   control_estimator[[k]] <- list(estimator = "varimax",
    #                                  labels = labels,
    #                                  indices = indices,
    #                                  p = p,
    #                                  q = q)
    #   k <- k+1L
    #
    # } else if(rotation == "varimin") {
    #
    #   labels <- c(rot_trans[[i]]$lambda)
    #   indices <- list(match(labels, transparameters_labels) - 1L)
    #   control_estimator[[k]] <- list(estimator = "varimin",
    #                                  labels = labels,
    #                                  indices = indices,
    #                                  p = p,
    #                                  q = q)
    #   k <- k+1L
    #
    # } else {
    #   stop("Unknown rotation. Available rotations: 'oblimin', 'geomin',
    #        'target', 'xtarget', 'varimax', 'varimin'")
    # }

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
                        unlist(modelInfo$rot_param)) # FIX THIS

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
                call               = mc, # matched call
                timing             = elapsed, # timing information
                data_list          = data_list,
                modelInfo          = modelInfo,
                Optim              = Optim,
                parameters         = parameters,
                transformed_pars   = transformed_pars,
                loglik             = loglik, # loglik values
                penalized_loglik   = penalized_loglik,
                loss               = loss,
                penalized_loss     = penalized_loss
  )

  return(result)

}

