# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 05/03/2026

pois_reg <- function(Y, X = NULL, penalties = FALSE, do.fit = TRUE, control = NULL) {

  mc <- match.call()
  Y <- as.matrix(Y)
  if(!is.numeric(Y) || !nrow(Y) || !ncol(Y) ||
     any(!is.finite(Y) | Y < 0 | Y != floor(Y))) {
    stop("Y must contain non-negative integer counts.")
  }
  if(!is.null(X) && ((!is.data.frame(X) && !is.matrix(X)) || nrow(X) != nrow(Y))) {
    stop("X must be a matrix or data.frame with one row per observation.")
  }
  if(!isFALSE(penalties)) stop("Poisson regression currently supports penalties = FALSE.")
  if(!is.null(control) && !is.list(control)) stop("control must be a list.")
  control$positive <- FALSE
  control$penalties <- FALSE
  control <- lcfa_control(control)

  if(is.null(X)) {
    X <- matrix(1, nrow(Y), 1L, dimnames = list(NULL, "(Intercept)"))
  } else {
    X <- as.data.frame(X)
    X <- make_design_matrix(names(X), X)
    if(any(!is.finite(X))) stop("Covariates must be observed and finite.")
  }
  nobs <- nrow(Y)
  p <- ncol(X)
  q <- ncol(Y)
  response_names <- colnames(Y)
  if(is.null(response_names)) response_names <- paste0("Y", seq_len(q))

  # Keep the observed response/covariate patterns and their row mappings.
  joint <- make_patterns(as.data.frame(cbind(unname(Y), unname(X))), any_continuous = FALSE)
  cov <- make_patterns(as.data.frame(X), any_continuous = FALSE)
  dataList <- list(Y = Y, X = X, nobs = nobs, npatterns = joint$npatterns,
                   patterns = Y[joint$full2short, , drop = FALSE],
                   cov_patterns2 = X[joint$full2short, , drop = FALSE],
                   weights = joint$pattern_weights, full2short = joint$full2short,
                   short2full = joint$short2full, cov_patterns = X[cov$full2short, , drop = FALSE],
                   ncov_patterns = cov$npatterns, cov_full2short = cov$full2short,
                   cov_short2full = cov$short2full)

  #### Set up the model ####

  trans <- create_parameters(list(
    list(name = "beta", type = "matrix", dim = c(p, q),
         rownames = colnames(X), colnames = response_names),
    list(name = "linpred", type = "matrix", dim = c(nobs, q), colnames = response_names),
    list(name = "lambda", type = "matrix", dim = c(nobs, q), colnames = response_names)
  ))
  param <- trans["beta"]
  manifolds <- create_manifolds(list(list(manifold = "euclidean", parameters = "beta")), param)
  transforms <- create_transforms(list(
    list(transform = "column_space", parameters_in = "beta", parameters_out = "linpred",
         extra = list(X = X)),
    list(transform = "exponential", parameters_in = "linpred", parameters_out = "lambda")
  ), trans)
  estimators <- create_estimators(list(list(estimator = "poisson_loglik", parameters = "lambda",
                                           extra = list(X = Y))), trans)
  beta <- matrix(0, p, q)
  beta[1L, ] <- log(pmax(colMeans(Y), 0.5))
  starts <- rep(list(list(beta = beta)), control$rstarts)
  initial <- create_init(trans, param, starts, transforms, control)
  control$parameters <- initial$parameters
  control$transparameters <- initial$transparameters
  control$transparam2param <- initial$trans2param - 1L
  control$cores <- min(control$cores, control$rstarts)
  modelInfo <- list(param = param, trans = trans, nobs = nobs, nparam = initial$nparam,
                    ntrans = initial$ntrans, npatterns = joint$npatterns,
                    dof = nobs*q - initial$nparam, parameters_labels = initial$parameters_labels,
                    transparameters_labels = initial$transparameters_labels,
                    control_manifold = manifolds, control_transform = transforms,
                    control_estimator = estimators, control_optimizer = control)

  result <- new("latent", version = as.character(packageVersion("latent")),
                call = mc, dataList = dataList, modelInfo = modelInfo)
  if(do.fit) {
    Optim <- optimizer(manifolds, transforms, estimators, control)
    names(Optim$parameters) <- initial$parameters_labels
    names(Optim$transparameters) <- initial$transparameters_labels
    result@Optim <- Optim
    result@timing <- Optim$elapsed
    result@transformed_pars <- fill_in(trans, Optim$transparameters)
    result@parameters <- result@transformed_pars["beta"]
  }
  result
}
