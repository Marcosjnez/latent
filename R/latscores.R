# Factor scores for latent: closed-form estimates and prediction errors.
# Works both in R/ and via source("latscores.R") after library(latent).

if (!methods::isClass("latent")) {
  if (!requireNamespace("latent", quietly = TRUE)) {
    stop("Install latent before sourcing latscores.R.")
  }
}
if (!methods::isClass("latscores")) {
  methods::setClass("latscores", contains = "latent")
}

#' Factor Scores and Their Prediction Errors
#'
#' Compute continuous factor scores from a fitted lcfa or lefa model, using
#' closed-form weights. No numerical optimization or score Hessian is needed.
#'
#' @param fit A fitted lcfa or lefa object.
#' @param data Raw observations; NULL uses the original fitting data. A matrix
#'   or data frame is accepted. For multiple groups, use the fitting group
#'   column, or a named list of matrices/data frames (one per fitted group).
#'   Columns are matched by name. An unnamed matrix must be in indicator order.
#' @param method One of "regression", "tenBerge", "Bartlett", "Harman", or
#'   "Green", ignoring case. "determinant", "Krijnen", and "McDonald" are
#'   aliases for tenBerge; "OLS" is an alias for Harman.
#' @param se Logical. Store model-based prediction errors. These treat fitted
#'   measurement-model parameters as fixed; they do not include uncertainty
#'   from estimating loadings, covariances, means, or standardization constants.
#' @param missing "patterns" uses each person's observed indicators;
#'   "complete" returns NA scores for incomplete rows without deleting them.
#'   Regression permits any pattern, including no observed indicators (the
#'   fitted latent mean). Other methods require full column rank within every
#'   scored pattern and stop if that requirement is not met.
#' @param center,scale Optional raw-data centering and scaling constants, as
#'   vectors in indicator order (or named vectors), or lists by fitted group.
#'   Scalars are recycled. Defaults reproduce the fitted model's scale and
#'   fitted observed means, using training data, never the new scoring sample.
#'   Explicit values override these defaults: scores are computed from
#'   sweep(sweep(data, 2, center, "-"), 2, scale, "/"). For a covariance-only
#'   fit without recorded means, supply center (0 for already centered data).
#' @param do.fit Logical. FALSE prepares the dataList and modelInfo only.
#' @param message Logical. Print a short progress message.
#' @param control Optional list; tol (default 1e-10) controls relative matrix
#'   rank/positive-definiteness checks. No optimizer settings are needed.
#'
#' @details
#' Write Sigma = Lambda Phi Lambda' + Theta, where Phi is latent_cov, including
#'   structural regressions when present. Scores are alpha + B'(y - mu).
#'   The full residual covariance Theta is retained.
#'
#' Regression minimizes mean squared prediction error without a covariance
#' constraint. Green minimizes its trace subject to B' Sigma B = Phi.
#' tenBerge is the determinant-optimal covariance-preserving version given in
#' the SAM workshop. Both preserve the fitted covariance in the population
#' model; the covariance of a finite sample of scores generally differs.
#' For orthonormal factors these methods coincide with Anderson-Rubin scores.
#' Bartlett is conditionally unbiased generalized least squares; Harman is
#' conditionally unbiased ordinary least squares.
#'
#' For every method, the reported prediction-error covariance is
#'   E = (B' Lambda - I) Phi (B' Lambda - I)' + B' Theta B.
#' The score SEs are sqrt(diag(E)), with no division by sample size. Under
#' Gaussian regression scoring E is also the posterior factor covariance.
#' For Green and tenBerge it is marginal prediction MSE, not an individual
#' posterior covariance or a parameter-estimation standard error.
#'
#' Only one small covariance matrix is stored per group/missingness pattern.
#' No covariance matrix of dimension (number of people times factors)^2 is
#' constructed. se() and vcov() use dedicated methods rather than the generic
#' multistep/multistage parameter-estimation machinery.
#'
#' Native factor order and signs are retained in all slots, as in lcfa/lefa.
#' latInspect(..., sort = TRUE) applies the source model's display ordering
#' and signs. se() and vcov() retain native coordinates.
#'
#' The formulas require continuous responses. Ordinal/polychoric/Yule models
#' are rejected because scoring numeric category codes is not equivalent to
#' scoring the model's underlying continuous responses.
#'
#' @return An S4 latscores object inheriting from latent, with dataList,
#'   modelInfo, Optim, parameters, transformed_pars, and extra slots.
#'   parameters$scores is the score matrix for a single group; multiple-group
#'   blocks are named scores.<group>. Optim$weights contains weight matrices
#'   by group and observed-data pattern, in model units. Optim$SE$table_se has
#'   score-shaped SE matrices, and Optim$SE$se is their labeled vector.
#'   Optim$SE$VCOV is a nested list of small prediction-error covariance blocks;
#'   Optim$SE$case_pattern identifies the block belonging to each person.
#'   extra$model retains fit. With se = FALSE, se() can compute SEs later.
#'
#' @references
#' Green (1969). Best linear composites with a specified structure.
#' https://doi.org/10.1007/BF02289359
#'
#' Krijnen, Wansbeek, and ten Berge (1996). Best linear predictors for factor
#' scores. https://doi.org/10.1080/03610929608831883
#'
#' ten Berge, Krijnen, Wansbeek, and Shapiro (1999). Some new results on
#' correlation-preserving factor scores prediction methods.
#' https://doi.org/10.1016/S0024-3795(97)10007-6
#'
#' @examples
#' \dontrun{
#' library(latent)
#' fit <- lcfa(HolzingerSwineford1939,
#'             model = 'visual =~ x1 + x2 + x3
#'                      textual =~ x4 + x5 + x6
#'                      speed =~ x7 + x8 + x9')
#' fs <- latscores(fit, method = "Green")
#' fs@parameters$scores
#' fs@Optim$SE$table_se$scores
#' latInspect(fs, "scores")               # Display order and signs.
#' se(fs)$table_se                        # Native order and signs.
#' vcov(fs, person = 1)                   # Small error covariance for person 1.
#' fs_new <- latscores(fit, data = HolzingerSwineford1939[1:5, ],
#'                     method = "regression")
#' }
#' @export
latscores <- function(fit, data = NULL, method = "regression", se = TRUE,
                      missing = "patterns", center = NULL, scale = NULL,
                      do.fit = TRUE, message = FALSE, control = NULL) {

  if (!inherits(fit, "lcfa") && !inherits(fit, "lefa")) {
    stop("fit must be a fitted lcfa or lefa object.")
  }
  if (!length(fit@Optim$transparameters)) stop("The input model is unfitted.")
  if (identical(fit@Optim$convergence, FALSE)) {
    warning("The input model did not converge; inspect it before using scores.")
  }
  for (x in list(se, do.fit, message)) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop("se, do.fit, and message must be TRUE or FALSE.")
    }
  }
  if (!is.character(method) || length(method) != 1L || is.na(method)) {
    stop("method must be a single character string.")
  }
  method <- switch(tolower(method), regression = "regression",
                   tenberge =, determinant =, krijnen =, mcdonald = "tenBerge",
                   bartlett = "Bartlett", harman =, ols = "Harman",
                   green = "Green", stop("Unknown factor-score method."))
  missing <- match.arg(missing, c("patterns", "complete"))
  if (is.null(control)) control <- list()
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || anyDuplicated(names(control)) ||
       any(!names(control) %in% "tol")))) {
    stop("control must be a list containing only tol.")
  }
  if (is.null(control$tol)) control$tol <- 1e-10
  if (!is.numeric(control$tol) || length(control$tol) != 1L ||
      !is.finite(control$tol) || control$tol <= 0 || control$tol >= 1) {
    stop("control$tol must be a number strictly between 0 and 1.")
  }
  control <- c(control, list(method = method, missing = missing, se = se,
                             opt = "analytic", rstarts = 1L, cores = 1L))
  .latscores_register()
  mc <- match.call()
  start <- proc.time()[["elapsed"]]

  dataList <- create_latscores_dataList(fit, data, center, scale, control)
  modelInfo <- create_latscores_modelInfo(dataList, control)
  result <- methods::new("latscores", version = fit@version, call = mc,
                         timing = numeric(), dataList = dataList,
                         modelInfo = modelInfo, Optim = list(),
                         parameters = list(), transformed_pars = list(),
                         extra = list(model = fit))
  if (!do.fit) return(result)
  if (message) cat("Computing", method, "factor scores.\n")

  output <- fit_latscores(dataList, modelInfo)
  result@parameters <- result@transformed_pars <- output$scores
  result@Optim <- output$Optim
  result@modelInfo$control_optimizer$parameters <- list(output$Optim$parameters)
  result@modelInfo$control_optimizer$transparameters <- list(output$Optim$transparameters)
  result@modelInfo$control_optimizer$init_param <- list(output$scores)
  result@modelInfo$control_optimizer$transparam2param <- seq_along(output$Optim$parameters) - 1L
  result@dataList$npatterns <- lapply(output$Optim$patterns, length)
  if (se) result@Optim$SE <- compute_se_latscores(result)
  result@timing <- result@Optim$elapsed <- proc.time()[["elapsed"]] - start
  result
}

# Data preprocessing is kept separate from the score formulas.
create_latscores_dataList <- function(fit, data, center, scale, control) {
  source <- if (inherits(fit, "lefa")) latent::lefa_source_lcfa(fit) else fit
  d <- source@dataList
  if (!identical(d$cor, "pearson")) {
    stop("These closed-form scores require a continuous (Pearson) model.")
  }
  factor_map <- if (inherits(fit, "lefa")) fit@modelInfo$data_param else d$data_param
  source_map <- d$data_param
  ng <- d$ngroups
  groups <- if (ng == 1L) "1" else as.character(d$group_label)
  blocks <- if (ng == 1L) "scores" else paste0("scores.", groups)
  std <- isTRUE(source@modelInfo$control_optimizer$std.ov)
  original <- d$args$data
  has_original <- is.data.frame(original) || is.matrix(original)
  stored_std <- std && d$missing %in% c("ml", "fiml") && !has_original
  default_data <- is.null(data)
  if (is.null(data)) data <- if (has_original) original else d$data
  if (is.null(data)) stop("No raw data were retained; supply data.")
  if (stored_std && !default_data && (is.null(center) || is.null(scale))) {
    stop("Original standardization constants were not retained; supply center and scale.")
  }

  # Accept one raw table or a list explicitly matched to the fitted groups.
  split_data <- function(x) {
    if (is.list(x) && !is.data.frame(x)) {
      if (length(x) != ng) stop("Supply one data matrix per fitted group.")
      if (ng > 1L) {
        if (is.null(names(x)) || anyDuplicated(names(x)) ||
            !setequal(names(x), groups)) stop("Data-list names must match fitted group labels.")
        x <- x[groups]
      }
      return(x)
    }
    if (!is.data.frame(x) && !is.matrix(x)) stop("data must be a matrix, data frame, or group list.")
    if (ng == 1L) return(list(x))
    if (is.null(d$group) || !d$group %in% colnames(x)) {
      stop("Supply the fitting group column, or a named data list by group.")
    }
    g <- as.character(x[, d$group])
    if (anyNA(g) || any(!g %in% groups)) stop("Unknown or missing group labels in data.")
    lapply(groups, function(a) x[g == a, , drop = FALSE])
  }
  raw <- split_data(data)
  training <- if (has_original) split_data(original) else d$data_per_group
  if (is.null(training)) training <- rep(list(NULL), ng)

  block <- function(pars, map, key, g, fallback = NULL) {
    for (field in c(key, fallback)) {
      name <- map[[field]]
      if (length(name) >= g && !is.na(name[g]) && !is.null(pars[[name[g]]])) {
        return(pars[[name[g]]])
      }
    }
    NULL
  }
  option <- function(x, g, items) {
    if (is.null(x)) return(NULL)
    if (is.list(x)) {
      if (length(x) != ng) stop("Centering/scaling lists must have one entry per group.")
      if (!is.null(names(x)) && ng > 1L) {
        if (anyDuplicated(names(x)) || !setequal(names(x), groups)) {
          stop("Centering/scaling list names must match fitted groups.")
        }
        x <- x[groups]
      }
      x <- x[[g]]
    }
    if (!is.null(names(x)) && length(x) != 1L) {
      if (anyDuplicated(names(x)) || !all(items %in% names(x))) {
        stop("Named centering/scaling vectors must include every indicator.")
      }
      x <- x[items]
    }
    if (!is.numeric(x) || !length(x) || any(!is.finite(x)) ||
        !length(x) %in% c(1L, length(items))) stop("Invalid centering/scaling constants.")
    rep_len(as.numeric(x), length(items))
  }
  select <- function(x, items) {
    if (!is.matrix(x) && !is.data.frame(x)) stop("Each data group must be a matrix or data frame.")
    if (is.null(colnames(x))) {
      if (ncol(x) != length(items)) stop("An unnamed matrix must have one column per indicator.")
      colnames(x) <- items
    }
    if (anyNA(colnames(x)) || anyDuplicated(colnames(x)) ||
        !all(items %in% colnames(x))) stop("Data columns must include the model's unique indicator names.")
    x <- x[, items, drop = FALSE]
    if (is.data.frame(x) && !all(vapply(x, is.numeric, logical(1L)))) {
      stop("All scored indicators must be numeric continuous variables.")
    }
    x <- as.matrix(x)
    if (!length(x)) storage.mode(x) <- "double"
    if (!is.numeric(x) || any(is.infinite(x))) stop("Indicators must be numeric without infinite values.")
    x
  }

  models <- centered <- centers <- scales <- display <- stats::setNames(vector("list", ng), groups)
  sorted <- latent::sort_latent(fit)
  orderings <- attr(sorted, "order", exact = TRUE)
  signs <- attr(sorted, "signs", exact = TRUE)
  for (g in seq_len(ng)) {
    L <- block(fit@transformed_pars, factor_map, "lambda_group", g)
    # Disturbance moments are a valid fallback only without structural paths.
    P <- block(fit@transformed_pars, factor_map, "latent_cov_group", g,
               if (!isTRUE(d$structural_model)) "psi_group")
    a <- block(fit@transformed_pars, factor_map, "latent_means_group", g,
               if (!isTRUE(d$structural_model)) "alpha_group")
    T <- block(source@transformed_pars, source_map, "theta_group", g)
    mu <- block(source@transformed_pars, source_map, "meanshat_group", g)
    if (!is.matrix(L) || !is.numeric(L) || any(!is.finite(L))) stop("Invalid fitted loadings.")
    p <- nrow(L); q <- ncol(L)
    items <- rownames(L); factors <- colnames(L)
    if (is.null(items) || is.null(factors) || p < 1L || q < 1L) stop("The fitted model lacks indicator/factor names.")
    if (!identical(dim(P), c(q, q)) || !identical(dim(T), c(p, p))) stop("Incompatible fitted covariance dimensions.")
    P <- .latscores_covariance(P, control$tol, "latent covariance")
    T <- .latscores_covariance(T, control$tol, "residual covariance")
    S <- .latscores_covariance(L %*% P %*% t(L) + T, control$tol, "model covariance", positive = TRUE)
    if (is.null(a)) {
      if (isTRUE(d$meanstructure)) stop("The fitted latent means are unavailable; refit the model.")
      a <- rep(0, q)
    }
    if (length(a) != q || any(!is.finite(a))) stop("Invalid fitted latent means.")
    if (is.null(mu)) {
      nu <- block(source@transformed_pars, source_map, "nu_group", g)
      if (is.null(nu)) {
        if (isTRUE(d$meanstructure)) stop("The fitted observed means/intercepts are unavailable.")
        nu <- rep(0, p)
      }
      mu <- c(nu) + c(L %*% a)
    }
    if (length(mu) != p || any(!is.finite(mu))) stop("Invalid fitted observed means.")
    X <- select(raw[[g]], items)
    train <- if (is.null(training[[g]])) NULL else select(training[[g]], items)
    c0 <- option(center, g, items)
    s0 <- option(scale, g, items)

    # Obtain the original training moments only when they are actually needed.
    train_mean <- if (!is.null(train) && nrow(train)) colMeans(train, na.rm = TRUE) else NULL
    if (is.null(train_mean)) train_mean <- option(d$args$sample.mean, g, items)
    if (is.null(s0)) {
      s0 <- rep(1, p)
      if (std && !stored_std) {
        s0 <- if (!is.null(train) && nrow(train)) {
          apply(train, 2L, stats::sd, na.rm = TRUE)
        } else sqrt(diag(d$sample.cov.input[[g]][items, items, drop = FALSE]))
      }
    }
    if (is.null(c0)) {
      if (stored_std) {
        c0 <- c(mu)
      } else if (std) {
        if (is.null(train_mean)) stop("Original observed means are unavailable; supply center.")
        c0 <- train_mean + s0 * c(mu)
      } else if (isTRUE(d$meanstructure)) {
        c0 <- c(mu)
      } else {
        if (is.null(train_mean)) stop("No mean structure or training means; supply center (0 for centered data).")
        c0 <- train_mean
      }
    }
    if (length(c0) != p || length(s0) != p || any(!is.finite(c0)) ||
        any(!is.finite(s0)) || any(s0 <= 0)) stop("Invalid training centering/scaling constants.")
    centers[[g]] <- stats::setNames(c0, items)
    scales[[g]] <- stats::setNames(s0, items)
    centered[[g]] <- sweep(sweep(X, 2L, c0, "-"), 2L, s0, "/")
    raw[[g]] <- X
    models[[g]] <- list(lambda = L, latent_cov = P, theta = T,
                         model = S, latent_means = stats::setNames(c(a), factors))
    index <- if (is.null(orderings)) seq_len(q) else orderings[[g]]
    direction <- signs[[factor_map$lambda_group[g]]]
    direction <- if (is.null(direction)) rep(1, q) else direction[1L, ]
    display[[g]] <- list(index = index, sign = direction)
  }
  names(raw) <- groups
  list(data = data, data_per_group = raw, centered = centered, models = models,
       center = centers, scale = scales, display = display, ngroups = ng,
       group_label = groups, scores_group = blocks, method = control$method,
       missing = control$missing, nobs = lapply(raw, nrow),
       nitems = lapply(raw, ncol), nfactors = lapply(models, function(x) ncol(x$lambda)),
       factor_label = lapply(models, function(x) colnames(x$lambda)))
}

create_latscores_modelInfo <- function(dataList, control) {
  specification <- lapply(seq_len(dataList$ngroups), function(g) {
    spec <- list(name = dataList$scores_group[g], type = "matrix",
         dim = c(dataList$nobs[[g]], dataList$nfactors[[g]]),
         rownames = rownames(dataList$data_per_group[[g]]),
         colnames = dataList$factor_label[[g]])
    if (dataList$nobs[[g]] == 0L) {
      spec$labels <- matrix(character(0L), 0L, dataList$nfactors[[g]])
    }
    spec
  })
  param <- latent::create_parameters(specification)
  labels <- unname(unlist(param))
  manifolds <- lapply(names(param), function(x) list(manifold = "euclidean", parameters = x))
  list(param = param, trans = param, nparam = length(labels), ntrans = length(labels),
       parameters_labels = labels, transparameters_labels = labels,
       data_param = list(scores_group = names(param)), propagate_uncertainty = FALSE,
       dof = NA_integer_, control_manifold = latent::create_manifolds(manifolds, param),
       control_transform = list(), control_estimator = list(), control_optimizer = control)
}

# Validate covariance matrices; clip only negative roundoff, never a Heywood case.
.latscores_covariance <- function(x, tol, name, positive = FALSE) {
  if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x))) stop("Invalid ", name, ".")
  size <- max(abs(x), .Machine$double.xmin)
  if (max(abs(x - t(x))) > tol * size * 10) stop(name, " must be symmetric.")
  x <- (x + t(x))/2
  e <- eigen(x, symmetric = TRUE)
  cutoff <- tol * max(abs(e$values), .Machine$double.xmin)
  if (min(e$values) < -cutoff || (positive && min(e$values) <= cutoff)) {
    stop(name, " must be ", if (positive) "positive definite" else "positive semidefinite", ".")
  }
  if (any(e$values < 0)) {
    fixed <- tcrossprod(sweep(e$vectors, 2L, sqrt(pmax(e$values, 0)), "*"))
    dimnames(fixed) <- dimnames(x)
    x <- fixed
  }
  x
}

.latscores_power <- function(x, power, tol, name) {
  x <- .latscores_covariance(x, tol, name, positive = TRUE)
  e <- eigen(x, symmetric = TRUE)
  sweep(e$vectors, 2L, e$values^power, "*") %*% t(e$vectors)
}

# Weights have p rows and q columns: centered_data %*% weights.
latscores_weights <- function(lambda, phi, sigma, method, tol = 1e-10) {
  if (method == "Harman") {
    K <- .latscores_covariance(crossprod(lambda), tol, "loading crossproduct", TRUE)
    return(t(solve(K, t(lambda))))
  }
  R <- chol(sigma)
  Z <- forwardsolve(t(R), lambda)
  G <- backsolve(R, Z)
  if (method == "regression") return(G %*% phi)
  if (method == "Bartlett") {
    K <- .latscores_covariance(crossprod(Z), tol, "loading information matrix", TRUE)
    return(t(solve(K, t(G))))
  }
  H <- .latscores_power(phi, 0.5, tol, "latent covariance")
  D <- if (method == "Green") phi %*% H else H
  # The polar factor U V' equals C (C'C)^(-1/2), C = Sigma^(-1/2) Lambda D.
  # Thin SVD avoids squaring the condition number in the inverse square root.
  sv <- svd(Z %*% D, nu = ncol(lambda), nv = ncol(lambda))
  if (length(sv$d) < ncol(lambda) || min(sv$d) <= tol * max(sv$d)) {
    stop("Covariance-preserving scores require full loading rank in this pattern.")
  }
  backsolve(R, sv$u %*% t(sv$v)) %*% H
}

fit_latscores <- function(dataList, modelInfo) {
  scores <- weights <- patterns <- memberships <- vector("list", dataList$ngroups)
  tol <- modelInfo$control_optimizer$tol
  for (g in seq_len(dataList$ngroups)) {
    X <- dataList$centered[[g]]
    m <- dataList$models[[g]]
    q <- ncol(m$lambda)
    F <- matrix(NA_real_, nrow(X), q, dimnames = list(rownames(X), colnames(m$lambda)))
    key <- if (nrow(X)) apply(!is.na(X), 1L, paste0, collapse = "") else character(0L)
    keys <- unique(key)
    pats <- wts <- vector("list", length(keys))
    membership <- match(key, keys)
    for (k in seq_along(keys)) {
      rows <- which(membership == k)
      obs <- which(!is.na(X[rows[1L], ]))
      pats[[k]] <- list(rows = rows, observed = obs)
      if (dataList$missing == "complete" && length(obs) < ncol(X)) next
      if (!length(obs)) {
        if (dataList$method != "regression") stop("Only regression scores support a completely missing row.")
        B <- matrix(0, 0L, q)
      } else {
        B <- tryCatch(latscores_weights(m$lambda[obs, , drop = FALSE], m$latent_cov,
                                        m$model[obs, obs, drop = FALSE], dataList$method, tol),
                      error = function(e) stop("Group ", dataList$group_label[g],
                        ", observed pattern ", k, ": ", conditionMessage(e), call. = FALSE))
      }
      dimnames(B) <- list(rownames(m$lambda)[obs], colnames(m$lambda))
      wts[[k]] <- B
      F[rows, ] <- sweep(X[rows, obs, drop = FALSE] %*% B, 2L, m$latent_means, "+")
    }
    names(pats) <- names(wts) <- if (length(keys)) paste0("pattern", seq_along(keys)) else character(0L)
    names(membership) <- rownames(X)
    scores[[g]] <- F; weights[[g]] <- wts
    patterns[[g]] <- pats; memberships[[g]] <- membership
  }
  names(scores) <- dataList$scores_group
  names(weights) <- names(patterns) <- names(memberships) <- dataList$group_label
  values <- stats::setNames(unname(unlist(scores)), modelInfo$parameters_labels)
  list(scores = scores, Optim = list(f = 0, convergence = TRUE, iterations = 0L,
       elapsed = 0, parameters = values, transparameters = values, SE = list(),
       weights = weights, patterns = patterns, case_pattern = memberships))
}

compute_se_latscores <- function(fit) {
  if (!length(fit@parameters)) stop("The score object has not been fitted.")
  tables <- fit@parameters
  covariances <- fit@Optim$weights
  for (g in seq_along(tables)) {
    tables[[g]][] <- NA_real_
    m <- fit@dataList$models[[g]]
    q <- ncol(m$lambda)
    for (k in seq_along(covariances[[g]])) {
      B <- fit@Optim$weights[[g]][[k]]
      if (is.null(B)) next
      pattern <- fit@Optim$patterns[[g]][[k]]
      obs <- pattern$observed
      D <- crossprod(B, m$lambda[obs, , drop = FALSE]) - diag(q)
      E <- D %*% m$latent_cov %*% t(D) + crossprod(B, m$theta[obs, obs, drop = FALSE] %*% B)
      dimnames(E) <- list(colnames(m$lambda), colnames(m$lambda))
      E <- .latscores_covariance(E, fit@modelInfo$control_optimizer$tol, "prediction-error covariance")
      covariances[[g]][[k]] <- E
      tables[[g]][pattern$rows, ] <- matrix(sqrt(pmax(diag(E), 0)),
                                            length(pattern$rows), q, byrow = TRUE)
    }
  }
  values <- stats::setNames(unname(unlist(tables)), fit@modelInfo$parameters_labels)
  list(table_se = tables, se = values, VCOV = covariances, vcov = covariances,
       case_pattern = fit@Optim$case_pattern, type = "prediction",
       parameter_uncertainty = FALSE)
}

#' @rdname latscores
#' @param type Only "prediction" is supported for score uncertainty.
#' @param digits Optional rounding of the displayed SE matrices, not the
#'   stored SE vector or covariance blocks. NULL returns full precision.
#' @param ... Reserved; unsupported arguments produce an error.
#' @method se latscores
#' @export
se.latscores <- function(fit, type = "prediction", digits = NULL, ...) {
  if (length(list(...))) stop("Unsupported arguments to se.latscores().")
  if (!identical(type, "prediction")) stop("Use type = 'prediction'; model-parameter uncertainty is not included.")
  out <- if (length(fit@Optim$SE)) fit@Optim$SE else compute_se_latscores(fit)
  if (!is.null(digits)) {
    if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
        digits < 0 || digits != as.integer(digits)) stop("digits must be NULL or a non-negative integer.")
    out$table_se <- lapply(out$table_se, round, digits = digits)
  }
  out
}

#' @rdname latscores
#' @param object A latscores object.
#' @param person Optional row number within group; return its q by q
#'   prediction-error covariance. Without person, return the compact nested
#'   list of covariance blocks; case-to-pattern indices are in Optim$case_pattern.
#' @param group Fitted group label or index; required for person in multigroup fits.
#' @method vcov latscores
#' @export
vcov.latscores <- function(object, person = NULL, group = NULL, ...) {
  if (length(list(...))) stop("Unsupported arguments to vcov.latscores().")
  out <- se.latscores(object)
  if (is.null(person) && is.null(group)) return(out$VCOV)
  if (is.null(group)) {
    if (object@dataList$ngroups > 1L) stop("Specify group for a multigroup score object.")
    group <- 1L
  }
  if (is.character(group)) group <- match(group, object@dataList$group_label)
  if (length(group) != 1L || !is.numeric(group) || !is.finite(group) ||
      group != as.integer(group) || group < 1L || group > object@dataList$ngroups) stop("Unknown group.")
  if (is.null(person)) return(out$VCOV[[group]])
  if (!is.numeric(person) || length(person) != 1L || !is.finite(person) ||
      person != as.integer(person) || person < 1L || person > object@dataList$nobs[[group]]) stop("Invalid person row number.")
  block <- out$case_pattern[[group]][person]
  E <- out$VCOV[[group]][[block]]
  if (is.null(E)) stop("This person has no scores under missing = 'complete'.")
  E
}

#' @rdname latscores
#' @param what One of "est", "scores", "se", "vcov", "weights",
#'   "case.pattern", "structures", "param", "trans", or "model".
#' @param sort Apply the input model's display ordering and signs.
#' @method latInspect latscores
#' @export
latInspect.latscores <- function(fit, what = "est", sort = TRUE) {
  if (!is.logical(sort) || length(sort) != 1L || is.na(sort)) stop("sort must be TRUE or FALSE.")
  if (!is.character(what) || length(what) != 1L || is.na(what)) stop("what must be one character string.")
  what <- tolower(what)
  order_columns <- function(x, signed = TRUE) {
    if (sort) {
      for (g in seq_along(x)) {
        view <- fit@dataList$display[[g]]
        x[[g]] <- x[[g]][, view$index, drop = FALSE]
        if (signed) x[[g]] <- sweep(x[[g]], 2L, view$sign, "*")
      }
    }
    x
  }
  if (what == "model") return(fit@extra$model)
  if (what %in% c("case.pattern", "case_pattern")) return(fit@Optim$case_pattern)
  if (what == "structures") return(list(parameters = order_columns(fit@parameters),
    transformed_pars = order_columns(fit@transformed_pars),
    param = order_columns(fit@modelInfo$param, FALSE),
    trans = order_columns(fit@modelInfo$trans, FALSE)))
  if (what %in% c("param", "trans", "labels")) return(order_columns(fit@modelInfo$trans, FALSE))
  if (!length(fit@parameters)) stop("The score object has not been fitted.")
  kind <- switch(what, est =, estimates =, parameters =, scores =, fscores = "scores",
                  se = "se", vcov =, covariance = "vcov", weights = "weights",
                  stop("Unknown score inspection request: ", what))
  out <- switch(kind, scores = fit@parameters, se = se.latscores(fit)$table_se,
                 vcov = vcov.latscores(fit), weights = fit@Optim$weights)
  if (kind %in% c("scores", "se")) {
    out <- order_columns(out, signed = kind == "scores")
  } else if (sort) {
    for (g in seq_along(out)) {
      index <- fit@dataList$display[[g]]$index
      sign <- fit@dataList$display[[g]]$sign
      out[[g]] <- lapply(out[[g]], function(x) {
        if (is.null(x)) return(NULL)
        if (kind == "vcov") return(x[index, index, drop = FALSE] * outer(sign, sign))
        sweep(x[, index, drop = FALSE], 2L, sign, "*")
      })
    }
  }
  if (what %in% c("scores", "fscores", "se") && length(out) == 1L) return(out[[1L]])
  out
}

# Keep Hessian-based generic inference from allocating enormous score matrices.
.latscores_no_derivatives <- function(x, ...) {
  stop("Factor-score uncertainty uses se() or vcov(); model-parameter derivatives are not defined for latscores objects.")
}

.latscores_register <- function() {
  # Explicit registration also supports a standalone sourced file without a
  # NAMESPACE edit. Pass function objects so dispatch works in local environments.
  registerS3method("se", "latscores", se.latscores, envir = asNamespace("latent"))
  registerS3method("latInspect", "latscores", latInspect.latscores, envir = asNamespace("latent"))
  registerS3method("vcov", "latscores", vcov.latscores, envir = asNamespace("stats"))
  for (generic in c("hessian", "information", "robust", "jacobian", "constraints_derivs")) {
    registerS3method(generic, "latscores", .latscores_no_derivatives, envir = asNamespace("latent"))
  }
  invisible(NULL)
}

methods::setMethod("show", "latscores", function(object) {
  cat(if (length(object@parameters)) "Factor scores" else "Unfitted factor-score specification", "\n")
  cat("  Method:", object@dataList$method, "\n")
  cat("  Observations:", sum(unlist(object@dataList$nobs)), "\n")
  cat("  Groups:", object@dataList$ngroups, "\n")
  cat("  Prediction errors:", if (length(object@Optim$SE)) "computed" else "not stored", "\n")
  invisible(object)
})

# latent currently has an .onAttach hook but no .onLoad hook. Register on load
# so an object restored with readRDS() works before any latscores() call.
if (isNamespace(environment())) {
  .onLoad <- function(libname, pkgname) .latscores_register()
} else {
  .latscores_register()
}
