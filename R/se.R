#' @title
#' Standard Errors
#' @description
#'
#' Compute standard errors.
#'
#' @usage
#'
#' se(fit)
#'
#' @param fit model.
#'
#' @details Compute standard errors.
#'
#' @return List with the following objects:
#' \item{vcov}{Variance-covariance matrix between the parameters.}
#' \item{se}{Standard errors.}
#'
#' @references
#'
#' None yet.
#'
#' @export
se <- function(fit) {

  data <- fit$opt$data
  model <- fit$modelInfo$model_vector # ALLOW CUSTOM MODELS
  nclasses <- length(fit$parameters$classes)
  control_manifold <- fit$opt$manifold_setup
  control_estimator <- fit$opt$estimator_setup
  control_optimizer <- fit$opt$control
  # control <- list(opt = "em", maxit = 1L, rstarts = 1L, cores = 1L)
  # control$force.combi <- TRUE # Use latentloglik_combination.h even when opt = "em"
  # empty_fit <- lca(data = data, model = model, nclasses = nclasses,
  #                  control = control, do.fit = FALSE)
  # control_manifold <- empty_fit$opt$manifold_setup
  # control_estimator <- empty_fit$opt$estimator_setup
  # control_optimizer <- empty_fit$opt$control

  gauss <- "gaussian" %in% model
  multin <- "multinomial" %in% model

  # parameters <- unlist(fit$parameters)
  parameters <- c(fit$opt$parameters)
  H <- numDeriv::jacobian(func = latent2::grad_comp,
                          x = parameters,
                          control_manifold = control_manifold,
                          control_estimator = control_estimator,
                          control_optimizer = control_optimizer)
  non_alnum_indices <- grep("^(?!-?\\d+(\\.\\d+)?$)",
                            unlist(fit$modelInfo$model), perl = TRUE)
  rownames(H) <- colnames(H) <- unlist(fit$modelInfo$model)[non_alnum_indices]

  # Create the constraints
  extract_column_constraints <- function(x) {
    # Recursive function to process each element
    if (is.list(x)) {
      # If it's a list, recurse and flatten
      unlist(lapply(x, extract_column_constraints), recursive = FALSE)
    } else if (is.matrix(x)) {
      # If it's a matrix, extract each column as a vector
      lapply(seq_len(ncol(x)), function(j) rep(1, times = length(x[, j])))
    } else if (is.vector(x)) {
      # If it's a vector, return as single-item list
      list(rep(1, times = length(x)))
    } else {
      stop("Unsupported element type in structure")
    }
  }

  if(multin & gauss) {
    multinom_indices <- which(model == "multinomial")
    vector_constraints <- extract_column_constraints(fit$modelInfo$model)
  } else if(multin) {
    vector_constraints <- extract_column_constraints(fit$modelInfo$model)
  } else if(gauss) {
    gauss_indices <- which(model == "gaussian")
    vector_constraints <- extract_column_constraints(fit$modelInfo$model$classes)
  }

  # Total number of rows is the number of parameters
  nparam <- length(unlist(fit$parameters))
  n_rows <- nparam
  n_cols <- length(vector_constraints)

  # Initialize matrix with zeros
  constraints <- matrix(0, nrow = n_rows, ncol = n_cols)

  # Fill each column block
  start_row <- 1
  for (i in seq_along(vector_constraints)) {
    v <- vector_constraints[[i]]
    len <- length(v)
    constraints[start_row:(start_row + len - 1), i] <- v
    start_row <- start_row + len
  }

  # REMOVE PROBABILITIES CLOSE TO ZERO in H and constraints:

  # remove <- which(unlist(empty_fit$modelInfo$model) == "-1")
  # remove <- which(abs(unlist(fit$parameters)) < 1e-05)
  # constraints <- constraints[-remove, ]
  # H <- H[-remove, ]
  # H <- H[, -remove]

  constraints <- constraints[non_alnum_indices, ]
  K <- t(constraints)
  D <- H + t(K) %*% K
  D_inv <- solve(D)
  C <- D_inv - D_inv %*% t(K) %*% solve(K %*% D_inv %*% t(K)) %*% K %*% D_inv
  se <- sqrt(diag(C))

  result <- list()
  result$se <- se
  result$vcov <- C

  return(result)

}
