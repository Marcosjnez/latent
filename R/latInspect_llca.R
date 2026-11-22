# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 10/09/2026
#'
#' Inspect fitted latent class and hidden Markov models
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
#'     \item{\code{"initial"}, \code{"transition"}}{For a hidden Markov model,
#'       the initial probabilities by trajectory or the common transition matrix.
#'       In the transition matrix, columns are previous states and rows are next states.}
#'     \item{\code{"forward"}, \code{"filtered"}}{Filtered state probabilities
#'       conditional on the observations through the current time.}
#'     \item{\code{"backward"}}{Scaled likelihood of the future observations
#'       for each current state. Its value is one at a trajectory's final time.
#'       Multiplying it by the forward probabilities gives the smoothed probabilities.}
#'     \item{\code{"smoothed"}, \code{"posterior"}}{For a hidden Markov model,
#'       the state probabilities given all observations in the trajectory.}
#'     \item{\code{"posterior_transition"}}{Joint probabilities of adjacent states
#'       given the whole trajectory, returned as an array indexed by next state,
#'       previous state, and destination observation. Trajectory starts have
#'       \code{NA} slices.}
#'     \item{\code{"viterbi"}}{Most probable joint state sequence for each
#'       trajectory. Values are integer state identifiers in the fitted model.}
#'     \item{\code{"loglik_trajectory"}}{Unweighted observed-data log-likelihood
#'       per trajectory. In a temporal model, \code{"loglik_case"} gives the same
#'       result. The sum weighted by trajectory weights gives the model log-likelihood.}
#'   }
#'   Matching is case-insensitive.
#' @param sort Logical. Sort only inspection outputs by decreasing weighted
#'   posterior size. Defaults to TRUE. Class identities and the original
#'   reference class are retained, including in modal assignments.
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
#' For hidden Markov models, observation-level results follow the original
#' input row order. Trajectory-level results follow the order in which IDs first
#' appeared. \code{"fullclasses"} gives the initial state probabilities, and
#' \code{"classes"} averages them using trajectory weights. The class sizes in
#' \code{"profile"} and the default state ordering use smoothed occupancy across
#' all observations, with each trajectory weight applied to its observations.
#' \code{"table"} returns the observations, marginal modal states, and smoothed
#' probabilities. The modal state at each time can differ from the Viterbi path.
#' Use \code{sort = FALSE} to keep the state order from estimation.
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
#' @param x Fitted model or collection of fitted models passed to the method.
#' @param ... Additional arguments to the generic.
#' @method latInspect llca
#' @export
latInspect.llca <- function(x,
                            what = "profile",
                            digits = 4L, sort = TRUE, ..., fit = x) {
  if(!missing(fit) && !missing(x)) return(.match_legacy_model("fit"))

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

  if(!is.logical(sort) || length(sort) != 1L || is.na(sort)) {
    stop("sort must be TRUE or FALSE.")
  }
  what <- tolower(what)
  if(what %in% c("est", "estimates", "parameters", "fixed", "structures", "param",
                 "trans", "labels", "transparameters", "transformed_pars", "se",
                 "standard.errors", "standard_errors", "vcov", "covariance")) {
    return(latInspect.latent(fit, what = what, sort = sort))
  }

  #### Optimization diagnostics ####

  if(what %in% c("elapsed", "timing")) return(fit@Optim$elapsed)
  if(what == "convergence") {
    return(data.frame(Iterations = fit@Optim$iterations,
                      convergence = fit@Optim$convergence,
                      grad.norm = fit@Optim$ng))
  }
  if(what == "gradient") {
    info <- data.frame(name = fit@modelInfo$parameters_labels, gradient = fit@Optim$g,
                       rgradient = fit@Optim$rg, dir = fit@Optim$dir)
    return(list(grad.norm = fit@Optim$ng, gradient.info = info))
  }
  if(what %in% c("loss", "losses", "loglik", "ll", "fit.matrix", "fit_matrix")) {
    fit_matrix <- .llca_fit_matrix(fit)
    if(what %in% c("loss", "losses")) {
      return(fit_matrix[c("loss", "penalized_loss"), "overall", drop = FALSE])
    }
    if(what %in% c("loglik", "ll")) {
      return(fit_matrix[c("loglik", "penalized_loglik"), "overall", drop = FALSE])
    }
    return(fit_matrix)
  }

  #### Observation and trajectory data ####

  d <- fit@dataList
  temporal <- !is.null(d$lengths)
  rows <- if(temporal) d$inverse_ordering else seq_len(nrow(d$data))
  data <- d$data
  row_names <- rownames(data)[rows]
  if(what %in% c("data", "measurement")) {
    if(temporal) data <- data[rows, , drop = FALSE]
    if(what == "measurement") return(data[, d$indicators_names, drop = FALSE])
    return(data)
  }

  if(temporal && what %in% c("pattern", "patterns", "loglik_pattern", "loglik.pattern")) {
    stop("Temporal likelihoods are per trajectory. Use 'loglik_trajectory' or 'table'.")
  }
  if(what %in% c("loglik_case", "loglik.case", "case", "loglik_trajectory",
                 "loglik_pattern", "loglik.pattern")) {
    ll <- as.numeric(fit@Optim$outputs$estimators$vectors[[1L]][[1L]])
    if(temporal) return(setNames(ll, as.character(d$trajectory_ids)))
    if(what == "loglik_trajectory") stop("This request requires a temporal model.")
    if(what %in% c("loglik_pattern", "loglik.pattern")) {
      return(setNames(d$pattern_weights * ll, d$pattern_names))
    }
    return(setNames(ll[d$short2full], row_names))
  }

  if(what %in% c("pattern", "patterns")) {
    patterns <- data[d$full2short, d$variables, drop = FALSE]
    result <- data.frame(patterns, Observed = d$pattern_weights)
    if(ncol(patterns)) result <- result[do.call(order, patterns), , drop = FALSE]
    rownames(result) <- paste0("pattern", seq_len(nrow(result)))
    return(result)
  }

  #### State probabilities ####

  what <- switch(what, filtered = "forward", smoothed = "posterior", what)
  if(!temporal && what %in% c("forward", "backward", "viterbi",
                              "posterior_transition", "transition", "initial")) {
    stop("This request requires a temporal model.")
  }
  if(what == "viterbi") {
    return(setNames(.hmm_quantities(fit, "viterbi")[rows], row_names))
  }

  nclasses <- ncol(fit@transformed_pars$class)
  class_ids <- seq_len(nclasses)
  class_names <- d$class_names
  if(length(class_names) != nclasses) class_names <- paste0("Class", class_ids)
  weights <- .llca_weights(fit)
  need_posterior <- what %in% c("profile", "posterior", "state", "states", "table",
                                "summary", "classification")
  if(sort && what %in% c("backward", "posterior_transition")) {
    recursion <- .hmm_quantities(fit, what, with_smoothed = TRUE)
    posterior <- recursion$smoothed
  } else if(sort || need_posterior) {
    posterior <- .llca_posterior(fit)
  }
  if(sort) class_ids <- .llca_order(posterior, weights)
  class_names <- class_names[class_ids]

  if(what %in% c("forward", "backward")) {
    value <- if(sort && what == "backward") recursion$value else .hmm_quantities(fit, what)
    result <- value[rows, class_ids, drop = FALSE]
    dimnames(result) <- list(row_names, class_names)
    return(result)
  }
  if(what == "posterior_transition") {
    value <- if(sort) recursion$value else .hmm_quantities(fit, what)
    result <- value[class_ids, class_ids, rows, drop = FALSE]
    dimnames(result) <- list(next_state = class_names, previous_state = class_names,
                             observation = row_names)
    return(result)
  }
  if(what == "transition") {
    return(fit@transformed_pars$transition[class_ids, class_ids, drop = FALSE])
  }
  if(need_posterior) {
    expanded <- if(temporal) rows else d$short2full
    # Modal ties retain the original state identity even when sorting the display.
    reported_states <- max.col(posterior, ties.method = "first")
    state_patterns <- match(reported_states, class_ids)
    posterior_patterns <- posterior[, class_ids, drop = FALSE]
    colnames(posterior_patterns) <- if(temporal) class_names else paste0("P(", class_names, "|data)")
    if(what %in% c("state", "states")) return(setNames(reported_states[expanded], row_names))
    if(what == "posterior") {
      result <- posterior_patterns[expanded, , drop = FALSE]
      rownames(result) <- row_names
      return(result)
    }
  }

  if(what %in% c("table", "summary")) {
    if(temporal) {
      variables <- unique(c(d$id, d$time, d$variables))
      return(data.frame(data[rows, variables, drop = FALSE], State = reported_states[rows],
                        Posterior = I(posterior_patterns[rows, , drop = FALSE]),
                        row.names = row_names))
    }
    patterns <- data[d$full2short, d$variables, drop = FALSE]
    ll <- as.numeric(fit@Optim$outputs$estimators$vectors[[1L]][[1L]])
    result <- cbind(patterns, Observed = weights, Estimated = exp(ll)*sum(weights),
                    State = reported_states, Posterior = posterior_patterns,
                    loglik_case = ll, loglik_patterns = weights*ll)
    if(ncol(patterns)) result <- result[do.call(order, patterns), , drop = FALSE]
    rownames(result) <- paste0("pattern", seq_len(nrow(result)))
    return(result)
  }

  if(what %in% c("fullclasses", "initial", "class", "classes", "cluster", "clusters")) {
    priors <- fit@transformed_pars$class[, class_ids, drop = FALSE]
    if(what %in% c("fullclasses", "initial")) return(priors)
    prior_weights <- if(temporal) d$trajectory_weights else weights
    return(colSums(priors * prior_weights) / sum(prior_weights))
  }
  if(what %in% c("beta", "betas", "coef", "coefs", "coefficient", "coefficients")) {
    return(fit@transformed_pars$beta[, class_ids, drop = FALSE])
  }

  if(what == "classification") {

    weighted_posterior <- posterior_patterns * weights

    denominator <- colSums(weighted_posterior)

    class_error_prop <- crossprod(posterior_patterns, weighted_posterior)
    class_error_prop <- sweep(class_error_prop, MARGIN = 1L,
                              STATS = denominator, FUN = "/")

    assigned_modal <- diag(nclasses)[state_patterns, , drop = FALSE]
    weighted_modal <- sweep(assigned_modal, MARGIN = 1L,
                            STATS = weights, FUN = "*")

    class_error_modal <- crossprod(posterior_patterns, weighted_modal)
    class_error_modal <- sweep(class_error_modal, MARGIN = 1L,
                               STATS = denominator, FUN = "/")

    class_error_prop[!is.finite(class_error_prop)] <- NA_real_
    class_error_modal[!is.finite(class_error_modal)] <- NA_real_

    assigned_names <- paste0("Assigned", class_ids)

    rownames(class_error_prop) <- class_names
    colnames(class_error_prop) <- assigned_names
    rownames(class_error_modal) <- class_names
    colnames(class_error_modal) <- assigned_names

    result <- list(class_error_modal = class_error_modal,
                   class_error_prop = class_error_prop)

  } else if(what == "diagnosis") {

    result <- lclass_diag(fit, type = "all", digits = digits)
    if(sort) {
      for(nm in intersect(names(result), c("Mostlikely.Class", "Avg.Mostlikely"))) {
        result[[nm]] <- result[[nm]][class_ids, class_ids, drop = FALSE]
        dimnames(result[[nm]]) <- list(class_names, class_names)
      }
      for(nm in intersect(names(result), c("AvePP", "Sum.Posterior", "Sum.Mostlikely"))) {
        if(nrow(result[[nm]]) == nclasses) {
          result[[nm]] <- result[[nm]][class_ids, , drop = FALSE]
          rownames(result[[nm]]) <- class_names
        }
      }
      for(nm in intersect(names(result), c("OCC", "Misclassification.per.class"))) {
        result[[nm]] <- setNames(result[[nm]][class_ids], class_names)
      }
      attr(result, "class_names") <- class_names
    }
  } else if(what %in% c("profile", "classconditional", "items", "item", "outcome",
                        "outcomes", "respconditional", "probcat")) {

    output <- list(parameters = fit@parameters, transformed_pars = fit@transformed_pars)
    if(sort) output <- .sort_llca_blocks(fit, output, class_ids)
    transformed_pars <- output$transformed_pars
    extract_conditional <- function(variable_names) {

      variable_names <- intersect(variable_names, names(transformed_pars))
      result <- transformed_pars[variable_names]

      gaussian_names_all <- c(d$gaussian$gaussian_names, d$mvgaussian$mvgaussian_names)

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

        } else if(name %in% d$multinomial$multinomial_names) {

          j <- match(name, d$multinomial$multinomial_names)

          if(!is.null(d$multinomial$multinomial_prob_rows)) {
            prob_rows <- d$multinomial$multinomial_prob_rows[[j]]
          } else {
            prob_rows <- seq_len(nrow(x) / 2L)
          }

          result[[name]] <- x[prob_rows, , drop = FALSE]

        }

      }

      return(result)

    }

    ClassConditional <- extract_conditional(d$indicators_names)
    OutcomeConditional <- extract_conditional(d$outcomes_names)
    if(what %in% c("classconditional", "items", "item")) return(ClassConditional)
    if(what %in% c("outcome", "outcomes")) return(OutcomeConditional)
    if(what == "profile") {
      result <- list(class_size = colSums(posterior_patterns * weights) / sum(weights),
                     indicators = ClassConditional)
      if(length(OutcomeConditional)) result$outcomes <- OutcomeConditional
      parameters <- output$parameters
      result$dependencies <- list(gaussian = parameters[grepl("Sigma", names(parameters))],
                                  multinomial = parameters[grepl("log_", names(parameters))])
      return(result)
    }

    RespConditional <- list()
    probCat <- list()

    if(length(d$variable_type) && all(d$variable_type == "multinomial")) {
      classes <- .llca_marginal_classes(fit, class_ids)

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
    result <- if(what == "probcat") probCat else RespConditional

  } else {
    stop("Unknown value of what: '", what,
         "'. See ?latInspect.llca for the supported values.")
  }

  #### Result ####

  return(result)

}

# Average model-implied probabilities over observations, before seeing responses.
.llca_marginal_classes <- function(fit, index) {
  d <- fit@dataList
  priors <- fit@transformed_pars$class
  if(is.null(d$lengths)) {
    return(colSums(priors[, index, drop = FALSE] * d$pattern_weights) / sum(d$pattern_weights))
  }
  total <- numeric(ncol(priors))
  for(q in seq_along(d$lengths)) {
    p <- priors[q, ]
    mass <- p
    if(d$lengths[q] > 1L) {
      for(t in seq_len(d$lengths[q] - 1L)) {
        p <- as.numeric(fit@transformed_pars$transition %*% p)
        mass <- mass + p
      }
    }
    total <- total + d$trajectory_weights[q] * mass
  }
  setNames(total[index] / sum(d$trajectory_weights * d$lengths), colnames(priors)[index])
}

.llca_fit_matrix <- function(fit) {
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

      result <- c(loss             = loss,
                  loss_base        = loss_base,
                  loss_sat         = loss_sat,
                  loglik           = loglik,
                  loglik_base      = loglik_base,
                  loglik_sat       = loglik_sat,
                  penalty          = penalty,
                  penalized_loss   = loss + penalty,
                  penalized_loglik = loglik - penalty)

      return(result)
    },
    FUN.VALUE = numeric(9L)
  )

  estimator_names <- unlist(lapply(fit@modelInfo$control_estimator, FUN = \(x) x$double_names),
                            use.names = FALSE)

  if(length(estimator_names) != ncol(fit_matrix)) {
    estimator_names <- paste0("estimator", seq_len(ncol(fit_matrix)))
  }

  colnames(fit_matrix) <- estimator_names
  fit_matrix <- cbind(fit_matrix, overall = rowSums(fit_matrix))
  fit_matrix
}

#' @rdname latInspect.llca
#' @param model For the \code{"llcalist"} method, a collection containing
#'   \code{"llca"} and/or \code{"llca_sam"} model fits.
#' @method latInspect llcalist
#' @export
latInspect.llcalist <- function(x,
                                what = "profile",
                                digits = 4L, sort = TRUE, ..., model = x) {
  if(!missing(model) && !missing(x)) return(.match_legacy_model("model"))

  if(!inherits(model, "llcalist")) {
    stop("model must inherit from class 'llcalist'.")
  }

  if(!is.logical(sort) || length(sort) != 1L || is.na(sort)) {
    stop("sort must be TRUE or FALSE.")
  }

  result <- lapply(model, FUN = latInspect, what = what, digits = digits, sort = sort)

  class(result) <- "latInspect.llcalist"

  #### Result ####

  return(result)

}
