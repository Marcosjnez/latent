# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 11/07/2026
#'
#' Predict latent class membership probabilities
#'
#' Compute latent class membership probabilities from a fitted latent class
#' model, using the class-membership regression coefficients and, when present,
#' observed covariates.
#'
#' @description
#' `predict.llca()` returns prior latent class membership probabilities implied
#' by the fitted class-membership model. Predictions depend only on the
#' covariates and the estimated `beta` coefficients. Measurement indicators and
#' distal outcomes are not used.
#'
#' `predict.llcalist()` applies the prediction method to collections of latent
#' class models. When the object contains named measurement and structural fits,
#' only predictions from the structural fit are returned.
#'
#' @usage
#' \method{predict}{llca}(model, new = NULL, ...)
#'
#' \method{predict}{llcalist}(model, new = NULL, ...)
#'
#' @param model A fitted object of class `"llca"` or `"llcalist"`.
#' @param new Optional `data.frame` or matrix containing new values for all
#'   covariates used in the fitted model. Column names must match the original
#'   covariate names. When `NULL`, predictions are computed for the data used to
#'   fit the model. For a model without covariates, `new` only determines the
#'   number and row names of the returned predictions.
#' @param ... Additional arguments. Currently not used.
#'
#' @details
#' For models with covariates, the method reconstructs the design matrix using
#' the variable types and factor levels from the fitted data. Its columns are
#' then matched to the design matrix used during estimation before calculating
#' the linear predictors and applying the softmax transformation.
#'
#' For models without covariates, the fitted latent class probabilities are
#' unconditional. The same probability vector is returned for every original
#' observation, or for every row of `new` when new data are supplied.
#'
#' The returned predictions are prior class-membership probabilities based only
#' on the class-membership model. They are different from posterior class
#' probabilities, which also condition on the observed indicators and can be
#' obtained with `latInspect(model, what = "posterior")`.
#'
#' For an `"llcalist"` containing named measurement and structural fits, the
#' measurement fit is skipped. For other collections containing both models with
#' and without covariates, only models with covariates are retained. For
#' class-enumeration objects, predictions are returned for every relevant fitted
#' model.
#'
#' @return
#' For an `"llca"` object, a `data.frame` containing the non-intercept columns
#' of the prediction design matrix, followed by one latent class probability
#' column per class. If the model has no covariates, only the class probability
#' columns are returned.
#'
#' For an `"llcalist"`, a single `data.frame` is returned when only one relevant
#' structural model is present. Otherwise, a list of prediction `data.frame`s is
#' returned with class `"predict.llcalist"`.
#'
#' @references
#' None yet.
#'
#' @examples
#' \dontrun{
#' fit <- lca(data = gss82, nclasses = 3L,
#'            multinomial = c("PURPOSE", "ACCURACY"))
#'
#' predict(fit)
#'
#' fit_cov <- lca(data = empathy, nclasses = 3L,
#'                gaussian = c("ec1", "ec2", "ec3"),
#'                covariates = c("sex", "pt1"), adjustment = "none")
#'
#' predict(fit_cov)
#' predict(fit_cov, new = data.frame(sex = c("female", "male"),
#'                                   pt1 = c(2.1, 3.4)))
#' }
#'
#' @method predict llca
#' @export
predict.llca <- function(model, new = NULL, ...) {

  if(!inherits(model, "llca")) {
    stop("model must inherit from class 'llca'.")
  }

  if(length(model@transformed_pars) == 0L ||
     is.null(model@transformed_pars$beta)) {
    stop("model must be fitted before predictions can be computed.")
  }

  covariates_names <- model@dataList$covariates_names
  if(is.null(covariates_names)) covariates_names <- character(0L)

  if(!is.null(new) && !is.data.frame(new) && !is.matrix(new)) {
    stop("new must be NULL, a data.frame, or a matrix.")
  }

  if(length(covariates_names) == 0L) {

    if(is.null(new)) {
      nnew <- nrow(model@dataList$data)
      prediction_names <- rownames(model@dataList$data)
    } else {
      nnew <- nrow(new)
      prediction_names <- rownames(new)
    }

    covariates <- matrix(1, nrow = nnew, ncol = 1L,
                         dimnames = list(prediction_names, "(Intercept)"))

  } else if(is.null(new)) {

    covariates <- model@dataList$design
    rownames(covariates) <- rownames(model@dataList$data)

  } else {

    new_names <- colnames(new)
    X_df <- as.data.frame(new)

    if(is.null(new_names) && ncol(X_df) == length(covariates_names)) {
      names(X_df) <- covariates_names
    }

    missing_covariates <- setdiff(covariates_names, names(X_df))
    if(length(missing_covariates) > 0L) {
      stop("The following covariates are missing in new: ",
           paste(missing_covariates, collapse = ", "), ".")
    }

    X_df <- X_df[, covariates_names, drop = FALSE]
    fitted_data <- model@dataList$data[, covariates_names, drop = FALSE]

    for(name in covariates_names) {

      original <- fitted_data[[name]]
      values <- X_df[[name]]

      if(is.factor(original) || is.character(original)) {

        original_factor <- if(is.factor(original)) droplevels(original) else factor(original)
        original_levels <- levels(original_factor)
        unknown_levels <- setdiff(unique(as.character(values)), original_levels)
        unknown_levels <- unknown_levels[!is.na(unknown_levels)]

        if(length(unknown_levels) > 0L) {
          stop("Unknown level(s) in covariate '", name, "': ",
               paste(unknown_levels, collapse = ", "), ".")
        }

        X_df[[name]] <- factor(values, levels = original_levels,
                               ordered = is.ordered(original_factor))

      } else if(is.numeric(original) && !is.numeric(values)) {
        stop("Covariate '", name, "' must be numeric.")
      }

    }

    if(anyNA(X_df)) {
      stop("new cannot contain missing values in the model covariates.")
    }

    mf <- model.frame(~ . + 1, data = X_df, na.action = na.pass)
    covariates <- model.matrix(~ . + 1, data = mf)

    fac_cols <- names(X_df)[vapply(X_df, is.factor, logical(1L))]
    terms_X <- terms(~ . + 1, data = X_df)

    for(v in fac_cols) {
      term_id <- match(v, attr(terms_X, "term.labels"))
      i <- which(attr(covariates, "assign") == term_id)

      if(length(i) > 0L) {
        old_names <- colnames(covariates)[i]
        level_names <- sub(paste0("^", v), "", old_names)
        colnames(covariates)[i] <- paste0(v, "_", level_names)
      }
    }

    expected_names <- colnames(model@dataList$design)
    extra_names <- setdiff(colnames(covariates), expected_names)

    if(length(extra_names) > 0L) {
      stop("The new data produce unexpected design-matrix column(s): ",
           paste(extra_names, collapse = ", "), ".")
    }

    missing_names <- setdiff(expected_names, colnames(covariates))
    if(length(missing_names) > 0L) {
      missing_matrix <- matrix(0, nrow = nrow(covariates),
                               ncol = length(missing_names),
                               dimnames = list(rownames(covariates), missing_names))
      covariates <- cbind(covariates, missing_matrix)
    }

    covariates <- covariates[, expected_names, drop = FALSE]

  }

  Z0 <- covariates %*% model@transformed_pars$beta
  row_max <- apply(Z0, MARGIN = 1L, FUN = max)
  Z0 <- sweep(Z0, MARGIN = 1L, STATS = row_max, FUN = "-")
  P0 <- exp(Z0)
  P0 <- P0 / rowSums(P0)

  class_names <- colnames(model@transformed_pars$class)
  if(is.null(class_names)) class_names <- paste("Class", seq_len(ncol(P0)), sep = "")
  colnames(P0) <- class_names
  rownames(P0) <- rownames(covariates)

  keep <- colnames(covariates) != "(Intercept)"
  if(any(keep)) {
    result <- data.frame(covariates[, keep, drop = FALSE], P0,
                         check.names = FALSE)
  } else {
    result <- data.frame(P0, check.names = FALSE)
  }

  #### Result ####

  return(result)

}

#' @rdname predict.llca
#' @method predict llcalist
#' @export
predict.llcalist <- function(model, new = NULL, ...) {

  if(!inherits(model, "llcalist")) {
    stop("model must inherit from class 'llcalist'.")
  }

  if(length(model) == 0L) {
    stop("model cannot be an empty llcalist.")
  }

  if("structural" %in% names(model) && inherits(model[["structural"]], "llca")) {
    result <- predict.llca(model[["structural"]], new = new, ...)

    #### Result ####

    return(result)
  }

  is_llca <- vapply(model, inherits, logical(1L), what = "llca")
  has_covariates <- logical(length(model))
  for(i in seq_along(model)) {
    if(is_llca[i]) {
      covariates_names <- model[[i]]@dataList$covariates_names
      has_covariates[i] <- length(covariates_names) > 0L
    }
  }

  if(any(has_covariates)) {
    keep <- which(has_covariates)
  } else {
    keep <- seq_along(model)
  }

  if(length(keep) == 1L && is_llca[keep]) {
    result <- predict.llca(model[[keep]], new = new, ...)

    #### Result ####

    return(result)
  }

  result <- vector("list", length = length(keep))
  for(i in seq_along(keep)) {
    x <- model[[keep[i]]]

    if(inherits(x, "llca")) {
      result[[i]] <- predict.llca(x, new = new, ...)
    } else if(inherits(x, "llcalist")) {
      result[[i]] <- predict.llcalist(x, new = new, ...)
    } else {
      stop("All selected elements of model must inherit from class 'llca' or 'llcalist'.")
    }
  }

  if(!is.null(names(model))) names(result) <- names(model)[keep]
  class(result) <- "predict.llcalist"

  #### Result ####

  return(result)

}
