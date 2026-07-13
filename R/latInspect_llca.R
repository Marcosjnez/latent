# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 13/07/2026
#'
#' Inspect fitted latent class models
#'
#' Extract model-implied class probabilities, posterior probabilities,
#' class-conditional response profiles, fit statistics, response-pattern
#' summaries, classification-error matrices, and optimization diagnostics from
#' fitted \code{"llca"} objects.
#'
#' @param fit An object of class \code{"llca"} returned by \code{\link{lca}}.
#' @param what Character string indicating the object to extract. Available
#'   values and aliases include:
#'   \describe{
#'     \item{\code{"profile"}}{Posterior class sizes and class-conditional
#'       indicator profiles. If distal outcomes are present, their
#'       class-conditional profiles are included in an \code{outcomes} element.}
#'     \item{\code{"classconditional"}, \code{"item"}, \code{"items"}}{
#'       Class-conditional indicator parameters. Gaussian indicators contain
#'       means and standard deviations; multinomial indicators contain
#'       probabilities only.}
#'     \item{\code{"outcome"}, \code{"outcomes"}}{Class-conditional distal
#'       outcome parameters.}
#'     \item{\code{"class"}, \code{"classes"}}{Average model-implied prior class
#'       probabilities, weighted by response-pattern weights.}
#'     \item{\code{"fullclasses"}}{Pattern-specific prior class probabilities.}
#'     \item{\code{"beta"}, \code{"coef"}, \code{"coefs"}}{Class-membership
#'       regression coefficients.}
#'     \item{\code{"posterior"}}{Case-specific posterior class probabilities.}
#'     \item{\code{"state"}, \code{"states"}}{Modal posterior class assignments.}
#'     \item{\code{"respconditional"}}{For fully multinomial measurement models,
#'       probabilities of latent-class membership conditional on each response
#'       category.}
#'     \item{\code{"probcat"}}{For fully multinomial measurement models,
#'       marginal response-category probabilities.}
#'     \item{\code{"data"}}{Processed data stored in the fitted object.}
#'     \item{\code{"measurement"}}{Processed indicator data.}
#'     \item{\code{"pattern"}, \code{"patterns"}}{Unique model-variable patterns
#'       and their observed weights.}
#'     \item{\code{"table"}, \code{"summary"}}{Response-pattern summary table.}
#'     \item{\code{"loglik_case"}}{Case-specific log-likelihood contributions.}
#'     \item{\code{"loglik_pattern"}}{Pattern-specific weighted
#'       log-likelihood contributions.}
#'     \item{\code{"loss"}, \code{"loglik"}, \code{"fit_matrix"}}{Optimizer loss,
#'       log-likelihood, or the complete matrix of fit components.}
#'     \item{\code{"classification"}}{Modal and proportional
#'       classification-error matrices. Rows represent latent classes and
#'       columns represent assigned classes.}
#'     \item{\code{"convergence"}}{Optimizer convergence information.}
#'     \item{\code{"gradient"}}{Gradient diagnostics.}
#'     \item{\code{"timing"}, \code{"elapsed"}}{Elapsed optimization time.}
#'   }
#'   Matching is case-insensitive.
#' @param digits Non-negative integer retained for compatibility with other
#'   inspection methods. Numeric results are returned without rounding so they
#'   can safely be used in subsequent computations.
#'
#' @details
#' Posterior probabilities and log-likelihood contributions are stored
#' internally by unique model-variable pattern and expanded back to the retained
#' observations when case-level results are requested.
#'
#' Pattern weights are used when calculating average class probabilities,
#' posterior class sizes, and classification-error matrices. Consequently,
#' frequency-weighted models produce weighted inspection results.
#'
#' Gaussian item matrices are read by row name, using \code{"mean"} and
#' \code{"stdv"}. Multinomial item matrices are assumed to contain probabilities
#' in their first rows and log-probability parameters in their following rows;
#' only the probability rows are returned in response profiles.
#'
#' @return The object requested through \code{what}. See the description of
#'   \code{what} for the possible return types.
#'
#' @references
#' None yet.
#'
#' @method latInspect llca
#' @export
latInspect.llca <- function(fit,
                            what = "profile",
                            digits = 4L) {

  #### Check inputs ####

  if(!inherits(fit, "llca")) {
    stop("fit must inherit from class 'llca'.")
  }

  if(!is.character(what) || length(what) != 1L || is.na(what)) {
    stop("what must be a single character string.")
  }

  if(length(digits) != 1L ||
     is.na(digits) ||
     digits < 0L ||
     digits != as.integer(digits)) {
    stop("digits must be a non-negative integer.")
  }

  if(length(fit@Optim) == 0L ||
     is.null(fit@Optim$outputs) ||
     is.null(fit@Optim$outputs$estimators)) {
    stop("The llca object does not contain fitted optimization results.")
  }

  what <- tolower(what)

  list2env(fit@dataList, envir = environment())
  list2env(fit@modelInfo, envir = environment())

  # Older and outcome-free llca objects may not contain outcomes_names.
  # Always create it explicitly before extracting conditional parameters.
  outcomes_names <- fit@dataList$outcomes_names
  if(is.null(outcomes_names)) outcomes_names <- character(0L)

  nclasses <- ncol(trans$class)
  total_weight <- sum(pattern_weights)

  #### Pattern-level likelihood and posterior quantities ####

  # Log-likelihood contribution of each response pattern:
  loglik_case_patterns <- fit@Optim$outputs$estimators$vectors[[1]][[1]]
  names(loglik_case_patterns) <- pattern_names

  # Weighted log-likelihood contribution of each response pattern:
  loglik_patterns <- pattern_weights * loglik_case_patterns
  names(loglik_patterns) <- pattern_names

  # Model-implied frequency on the total-weight scale:
  estimated <- exp(loglik_case_patterns) * total_weight

  # Posterior class probabilities by response pattern:
  posterior_patterns <- exp(matrix(fit@Optim$outputs$estimators$matrices[[1]][[1]],
                                   nrow = npatterns, ncol = nclasses))
  rownames(posterior_patterns) <- pattern_names
  colnames(posterior_patterns) <- paste0("P(Class", seq_len(nclasses), "|data)")

  # Modal posterior class by response pattern:
  state_patterns <- max.col(posterior_patterns, ties.method = "first")
  names(state_patterns) <- pattern_names

  #### Response-pattern summary ####

  patterns <- data[full2short, variables, drop = FALSE]

  summary_table <- cbind(patterns, Observed = pattern_weights,
                         Estimated = estimated, State = state_patterns,
                         Posterior = posterior_patterns,
                         loglik_case = loglik_case_patterns,
                         loglik_patterns = loglik_patterns)
  summary_table <- as.data.frame(summary_table)

  if(ncol(patterns) > 0L) {
    pattern_order <- do.call(order, patterns)
    summary_table <- summary_table[pattern_order, , drop = FALSE]
  }

  rownames(summary_table) <- paste0("pattern", seq_len(nrow(summary_table)))

  #### Expand pattern-level quantities to retained observations ####

  loglik_case <- loglik_case_patterns[short2full]
  posterior <- posterior_patterns[short2full, , drop = FALSE]
  state <- state_patterns[short2full]

  rownames(posterior) <- rownames(data)
  names(state) <- rownames(data)
  names(loglik_case) <- rownames(data)

  #### Average class probabilities ####

  weighted_classes <- sweep(fit@transformed_pars$class, MARGIN = 1L,
                            STATS = pattern_weights, FUN = "*")
  classes <- colSums(weighted_classes) / total_weight

  weighted_posterior <- sweep(posterior_patterns, MARGIN = 1L,
                              STATS = pattern_weights, FUN = "*")
  posterior_classes <- colSums(weighted_posterior) / total_weight

  #### Class-conditional parameters ####

  extract_conditional <- function(variable_names) {

    variable_names <- intersect(variable_names, names(fit@transformed_pars))
    result <- fit@transformed_pars[variable_names]

    gaussian_names_all <- c(gaussian$gaussian_names, mvgaussian$mvgaussian_names)

    for(name in names(result)) {

      x <- result[[name]]

      if(name %in% gaussian_names_all) {

        required_rows <- c("mean", "stdv")

        if(is.null(rownames(x)) ||
           !all(required_rows %in% rownames(x))) {
          stop(
            "Gaussian parameter matrix '", name,
            "' does not contain rows named 'mean' and 'stdv'."
          )
        }

        result[[name]] <- x[required_rows, , drop = FALSE]

      } else if(name %in% multinomial$multinomial_names) {

        j <- match(name, multinomial$multinomial_names)

        if(!is.null(multinomial$multinomial_prob_rows)) {
          prob_rows <- multinomial$multinomial_prob_rows[[j]]
        } else {
          prob_rows <- seq_len(nrow(x) / 2L)
        }

        result[[name]] <- x[prob_rows, , drop = FALSE]

      }

    }

    return(result)

  }

  ClassConditional <- extract_conditional(indicators_names)
  OutcomeConditional <- extract_conditional(outcomes_names)

  #### Additional outputs for fully multinomial measurement models ####

  RespConditional <- list()
  probCat <- list()

  if(length(variable_type) > 0L &&
     all(variable_type == "multinomial")) {

    # Marginal response-category probabilities:
    probCat <- lapply(ClassConditional, FUN = function(mat) {

      jointp <- sweep(t(mat), MARGIN = 1L, STATS = classes, FUN = "*")
      colSums(jointp)

    })

    # Class probabilities conditional on each response category:
    RespConditional <- lapply(ClassConditional, FUN = function(mat) {

      jointp <- sweep(t(mat), MARGIN = 1L, STATS = classes, FUN = "*")
      category_prob <- colSums(jointp)

      result <- sweep(t(jointp), MARGIN = 1L, STATS = category_prob, FUN = "/")
      result[!is.finite(result)] <- NA_real_

      return(result)

    })

  }

  profile <- list(class_size = posterior_classes,
                  indicators = ClassConditional)

  if(length(OutcomeConditional) > 0L) profile$outcomes <- OutcomeConditional

  #### Fit components ####

  doubles <- fit@Optim$outputs$estimators$doubles
  names(doubles) <- sapply(fit@modelInfo$control_estimator, FUN = \(x) x$estimator)
  names_doubles <- fit@Optim$outputs$estimators$names_doubles
  doubles <- Map(setNames, doubles, names_doubles)

  # Get fit information form the estimator structures. When doing so, set to
  # zero any fit measure that is missing in the estimators:

  get_zero <- function(x, name) {
    # Function to set to zero missing fit measures

    if(is.null(names(x)) || !name %in% names(x)) {
      value <- 0
    } else {
      value <- x[[name]]
    }

    return(value)

  }

  # Get all the fits:
  fit_matrix <- vapply(
    doubles,
    FUN = function(x) {
      loss        <- get_zero(x, "loss")
      loglik      <- get_zero(x, "loglik")
      penalty     <- get_zero(x, "penalty")
      loss_base   <- get_zero(x, "loss_baseline")
      loss_sat    <- get_zero(x, "loss_saturated")
      loglik_base <- get_zero(x, "loglik_baseline")
      loglik_sat  <- get_zero(x, "loglik_saturated")

      result <- c(
        loss             = loss,
        loss_base        = loss_base,
        loss_sat         = loss_sat,
        loglik           = loglik,
        loglik_base      = loglik_base,
        loglik_sat       = loglik_sat,
        penalty          = penalty,
        penalized_loss   = loss + penalty,
        penalized_loglik = loglik - penalty
      )

      #### Result ####
      return(result)
    },
    FUN.VALUE = numeric(9L)
  )

  estimator_names <- unlist(lapply(control_estimator, FUN = \(x) x$double_names),
                            use.names = FALSE)

  if(length(estimator_names) != ncol(fit_matrix)) {
    estimator_names <- paste0("estimator", seq_len(ncol(fit_matrix)))
  }

  colnames(fit_matrix) <- estimator_names
  fit_matrix <- cbind(fit_matrix, overall = rowSums(fit_matrix))

  #### Select requested result ####

  if(what == "profile") {

    result <- profile

  } else if(what %in% c("classconditional", "items", "item")) {

    result <- ClassConditional

  } else if(what %in% c("outcome", "outcomes")) {

    result <- OutcomeConditional

  } else if(what %in% c("class", "classes", "cluster", "clusters")) {

    result <- classes

  } else if(what == "fullclasses") {

    result <- fit@transformed_pars$class

  } else if(what %in% c("beta", "betas", "coef", "coefs",
                        "coefficient", "coefficients")) {

    result <- fit@transformed_pars$beta

  } else if(what == "respconditional") {

    result <- RespConditional

  } else if(what == "convergence") {

    result <- data.frame(Iterations = fit@Optim$iterations,
                         convergence = fit@Optim$convergence,
                         grad.norm = fit@Optim$ng)

  } else if(what == "gradient") {

    gradient_info <- data.frame(name = parameters_labels,
                                gradient = fit@Optim$g,
                                rgradient = fit@Optim$rg,
                                dir = fit@Optim$dir)

    result <- list(grad.norm = fit@Optim$ng, gradient.info = gradient_info)

  } else if(what == "data") {

    result <- data

  } else if(what == "measurement") {

    result <- data[, indicators_names, drop = FALSE]

  } else if(what %in% c("pattern", "patterns")) {

    result <- data.frame(patterns, Observed = pattern_weights)

    if(ncol(patterns) > 0L) {
      pattern_order <- do.call(order, patterns)
      result <- result[pattern_order, , drop = FALSE]
    }

    rownames(result) <- paste0("pattern", seq_len(nrow(result)))

  } else if(what %in% c("table", "summary")) {

    result <- summary_table

  } else if(what == "posterior") {

    result <- posterior

  } else if(what %in% c("state", "states")) {

    result <- state

  } else if(what %in% c("loglik_case", "loglik.case", "case")) {

    result <- loglik_case

  } else if(what %in% c("loglik_pattern", "loglik.pattern")) {

    result <- loglik_patterns

  } else if(what %in% c("loss", "losses")) {

    result <- fit_matrix[c("loss", "penalized_loss"), "overall", drop = FALSE]

  } else if(what %in% c("loglik", "ll")) {

    result <- fit_matrix[c("loglik", "penalized_loglik"), "overall", drop = FALSE]

  } else if(what %in% c("fit_matrix", "fit.matrix")) {

    result <- fit_matrix

  } else if(what == "probcat") {

    result <- probCat

  } else if(what == "classification") {

    denominator <- colSums(weighted_posterior)

    class_error_prop <- crossprod(posterior_patterns, weighted_posterior)
    class_error_prop <- sweep(class_error_prop, MARGIN = 1L,
                              STATS = denominator, FUN = "/")

    assigned_modal <- diag(nclasses)[state_patterns, , drop = FALSE]
    weighted_modal <- sweep(assigned_modal, MARGIN = 1L,
                            STATS = pattern_weights, FUN = "*")

    class_error_modal <- crossprod(posterior_patterns, weighted_modal)
    class_error_modal <- sweep(class_error_modal, MARGIN = 1L,
                               STATS = denominator, FUN = "/")

    class_error_prop[!is.finite(class_error_prop)] <- NA_real_
    class_error_modal[!is.finite(class_error_modal)] <- NA_real_

    class_names <- paste0("Class", seq_len(nclasses))
    assigned_names <- paste0("Assigned", seq_len(nclasses))

    rownames(class_error_prop) <- class_names
    colnames(class_error_prop) <- assigned_names
    rownames(class_error_modal) <- class_names
    colnames(class_error_modal) <- assigned_names

    result <- list(class_error_modal = class_error_modal,
                   class_error_prop = class_error_prop)

  } else if(what == "diagnosis") {

    result <- lclass_diag(fit, type = "all", digits = digits)

  } else if(what %in% c("timing", "elapsed")) {

    result <- fit@Optim$elapsed

  } else {

    stop(
      "Unknown value of what: '", what,
      "'. See ?latInspect.llca for the supported values."
    )

  }

  #### Return ####

  return(result)

}

#' @rdname latInspect.llca
#' @param model An object of class \code{"llcalist"} containing fitted
#'   \code{"llca"} objects.
#' @method latInspect llcalist
#' @export
latInspect.llcalist <- function(model,
                                what = "profile",
                                digits = 4L) {

  if(!inherits(model, "llcalist")) {
    stop("model must inherit from class 'llcalist'.")
  }

  result <- lapply(model, FUN = latInspect, what = what, digits = digits)

  class(result) <- "latInspect.llcalist"

  return(result)

}
