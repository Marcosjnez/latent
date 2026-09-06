# Author: Marcos Jimenez
# Modification date: 06/09/2026
#
# Sorting is a reporting transformation, not a second optimization. Keep the
# estimation chart intact: fixed markers, equality constraints, penalties,
# multinomial reference categories, and previous-stage identities must not be
# changed by a data-dependent display order.

#' Sort Factors or Latent Classes in a Fitted Model
#'
#' @param fit A fitted lcfa/lrotate/lefa object for sort_factors(), a fitted
#'   llca object for sort_classes(), or a sorted latent object for unsort_latent().
#' @param reorder Logical. For \code{sort_factors()}, reorder factors by their
#'   variance-adjusted sums of squared loadings. If FALSE, preserve the fitted
#'   factor order and only orient signs so the largest absolute loading of each
#'   factor is positive. Defaults to TRUE.
#'
#' @details
#' sort_factors() ranks factors within each group by
#' colSums(lambda^2)*diag(psi), in decreasing order when \code{reorder = TRUE}.
#' With unit factor variances this is the sum of squared loadings. For oblique
#' factors it is a ranking criterion, not an additive allocation of common
#' variance. With \code{reorder = FALSE}, the fitted factor order is retained.
#' In both cases the largest absolute loading in each factor is made positive;
#' ties retain the first item or factor in the original order. Factor names are
#' retained.
#'
#' sort_classes() ranks classes by their frequency-weighted expected posterior
#' sizes. Class names and the original multinomial reference class are retained:
#' the reference class need not occupy the first displayed column. Sorting does
#' not change regression contrasts, penalties, likelihoods, or fitted values.
#'
#' Public parameter matrices, their label templates, and stored sampling
#' covariances use the reported order/sign. Raw optimizer coordinates, Hessians,
#' and free-parameter covariance matrices retain the original fitting chart.
#' Negative factor coordinates are also registered as exact linear transformed
#' parameters, so subsequent delta-method standard errors include the signs.
#' The original fitted object and the reporting map are stored in
#' modelInfo$sorting; unsort_latent() restores the original coordinates without
#' refitting. The sorting functions are idempotent.
#'
#' Standard errors condition on the selected permutation and signs. They do not
#' measure uncertainty about which factor or class ranks first near a tie.
#'
#' @return An object with the same S4 class as fit.
#' @name sort_latent
NULL

#' @rdname sort_latent
#' @export
sort_factors <- function(fit, reorder = TRUE) {

  #### Check inputs ####

  check_sort_flag(reorder)
  check_fitted_sort_latent(fit)

  sorting <- fit@modelInfo$sorting
  current_reorder <- sorting$reorder
  if(is.null(current_reorder) && identical(sorting$type, "factors")) {
    current_reorder <- TRUE
  }

  if(identical(sorting$type, "factors") &&
     identical(current_reorder, reorder)) {
    #### Result ####
    return(fit)
  }

  native <- unsort_latent(fit)
  rotation <- is_rotation_sort_latent(native)

  if(!inherits(native, "lcfa") && !rotation) {
    stop("fit must be a fitted lcfa, lrotate, or lefa object.")
  }

  data_param <- if(rotation) native@modelInfo$data_param else native@dataList$data_param
  groups <- length(data_param$lambda_group)
  if(groups < 1L || length(data_param$psi_group) != groups) {
    stop("The fitted object is missing its factor parameter-block metadata.")
  }
  maps <- list()
  orderings <- signs <- importance <- vector("list", groups)

  #### Factor permutations and signs ####

  for(i in seq_len(groups)) {

    lambda_name <- data_param$lambda_group[i]
    psi_name <- data_param$psi_group[i]
    alpha_name <- data_param$alpha_group[i]
    lambda <- native@transformed_pars[[lambda_name]]
    psi <- native@transformed_pars[[psi_name]]

    if(!is.matrix(lambda) || !is.numeric(lambda) || any(!is.finite(lambda)) ||
       !is.matrix(psi) || !is.numeric(psi) || any(!is.finite(psi)) ||
       ncol(lambda) < 1L || !identical(dim(psi), rep.int(ncol(lambda), 2L))) {
      stop("The fitted factor loading/covariance blocks are invalid for sorting.")
    }

    if(any(diag(psi) < 0)) {
      warning("A negative factor variance was found; sorting does not repair the improper solution.")
    }

    q <- ncol(lambda)
    contribution <- colSums(lambda^2)*diag(psi)
    if(any(!is.finite(contribution))) {
      stop("Factor variance contributions are non-finite.")
    }
    index <- if(reorder) {
      order(-contribution, seq_len(q))
    } else {
      seq_len(q)
    }

    direction <- vapply(seq_len(q), FUN = function(j) {
      value <- lambda[which.max(abs(lambda[, j])), j]
      result <- if(value < 0) -1 else 1
      #### Result ####
      return(result)
    }, FUN.VALUE = numeric(1L))
    direction <- direction[index]

    maps[[lambda_name]] <- sort_block_map(lambda, columns = index,
                                          column_sign = direction)
    maps[[psi_name]] <- sort_block_map(psi, rows = index, columns = index,
                                       row_sign = direction, column_sign = direction)
    if(length(alpha_name) == 1L && alpha_name %in% names(native@transformed_pars)) {
      maps[[alpha_name]] <- sort_block_map(native@transformed_pars[[alpha_name]],
                                           rows = index, row_sign = direction)
    }

    if(rotation) {
      X_name <- data_param$X_group[i]
      Xinv_name <- data_param$Xinv_group[i]
      maps[[X_name]] <- sort_block_map(native@transformed_pars[[X_name]],
                                       columns = index, column_sign = direction)
      if(Xinv_name %in% names(native@transformed_pars)) {
        maps[[Xinv_name]] <- sort_block_map(native@transformed_pars[[Xinv_name]],
                                            rows = index, row_sign = direction)
      }
    } else {
      X_name <- data_param$xpsi_group[i]
      if(length(X_name) == 1L && X_name %in% names(native@transformed_pars)) {
        maps[[X_name]] <- sort_block_map(native@transformed_pars[[X_name]],
                                         columns = index, column_sign = direction)
      }
    }

    orderings[[i]] <- index
    signs[[i]] <- direction
    importance[[i]] <- setNames(contribution[index], colnames(lambda)[index])
  }

  specification <- list(type = "factors", blocks = maps, block_order = NULL,
                         order = orderings, signs = signs, importance = importance,
                         rotation = rotation, reorder = reorder)
  result <- apply_sort_latent(native, specification)

  #### Result ####

  return(result)

}

#' @rdname sort_latent
#' @export
sort_classes <- function(fit) {

  #### Check inputs ####

  check_fitted_sort_latent(fit)
  if(identical(fit@modelInfo$sorting$type, "classes")) {
    #### Result ####
    return(fit)
  }
  native <- unsort_latent(fit)
  if(!inherits(native, "llca")) {
    stop("fit must be a fitted llca object.")
  }

  class_names <- native@dataList$class_names
  nclasses <- length(class_names)
  if(anyNA(class_names) || any(class_names == "") || anyDuplicated(class_names)) {
    stop("Class names must be unique and non-empty.")
  }
  log_posterior <- native@Optim$outputs$estimators$matrices[[1L]][[1L]]
  weights <- as.numeric(native@dataList$pattern_weights)
  if(nclasses < 1L || length(log_posterior) != length(weights)*nclasses ||
     any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
    stop("The fitted posterior probabilities or pattern weights are invalid.")
  }
  posterior <- exp(matrix(log_posterior, nrow = length(weights), ncol = nclasses))
  if(any(!is.finite(posterior)) || any(abs(rowSums(posterior)-1) > 1e-06)) {
    stop("Posterior class probabilities must be finite and sum to one.")
  }
  sizes <- colSums(posterior*weights)
  index <- order(-sizes, seq_len(nclasses))

  #### Class-specific parameter blocks ####

  maps <- list()
  class_blocks <- unique(c("beta", "theta", "class", "loglik",
                            native@dataList$gaussian$gaussian_names,
                            native@dataList$mvgaussian$mvgaussian_names,
                            native@dataList$multinomial$multinomial_names,
                            names(native@dataList$mvmultinomial$indep_pairs_list)))
  for(nm in intersect(class_blocks, names(native@transformed_pars))) {
    x <- native@transformed_pars[[nm]]
    if(!is.matrix(x) || ncol(x) != nclasses) {
      stop("A class-specific parameter block has incompatible dimensions.")
    }
    maps[[nm]] <- sort_block_map(x, columns = index)
  }

  # Residual association matrices without a class axis must not be permuted,
  # even when their number of categories happens to equal the class count.
  block_order <- names(native@modelInfo$trans)
  sigma_names <- native@dataList$mvgaussian$sigma_names
  if(length(sigma_names) == nclasses && all(sigma_names %in% block_order)) {
    block_order[match(sigma_names, block_order)] <- sigma_names[index]
  }

  specification <- list(type = "classes", blocks = maps, block_order = block_order,
                         order = index, sizes = setNames(sizes[index], class_names[index]),
                         reference_class = class_names[1L])
  result <- apply_sort_latent(native, specification)

  #### Result ####

  return(result)

}

#' @rdname sort_latent
#' @export
unsort_latent <- function(fit) {

  if(!inherits(fit, "latent")) {
    stop("fit must inherit from class 'latent'.")
  }
  metadata <- fit@modelInfo$sorting
  result <- fit
  if(!is.null(metadata)) {
    if(!inherits(metadata$native, "latent")) {
      stop("The original fitting coordinates are missing from the sorted object.")
    }
    # Preserve the concrete outer class (for example, lefa wrapping lrotate).
    for(nm in methods::slotNames(metadata$native)) {
      methods::slot(result, nm) <- methods::slot(metadata$native, nm)
    }
    if(!is.null(fit@dataList$lefa_call)) {
      result@dataList$lefa_call <- fit@dataList$lefa_call
    }
  }

  #### Result ####

  return(result)

}

check_sort_flag <- function(sort) {

  if(!is.logical(sort) || length(sort) != 1L || is.na(sort)) {
    stop("sort must be TRUE or FALSE.")
  }

  #### Result ####

  return(invisible(NULL))

}

check_fitted_sort_latent <- function(fit) {

  if(!inherits(fit, "latent") || length(fit@Optim$transparameters) == 0L ||
     length(fit@transformed_pars) == 0L) {
    stop("fit must contain a fitted latent model.")
  }

  #### Result ####

  return(invisible(NULL))

}

is_rotation_sort_latent <- function(fit) {

  x <- fit@modelInfo$data_param
  result <- !is.null(x$X_group) && !is.null(x$lambda_group) &&
    !is.null(x$ulambda_group) && !is.null(fit@dataList$rotation)

  #### Result ####

  return(result)

}

sort_block_map <- function(x, rows = seq_len(nrow(x)), columns = seq_len(ncol(x)),
                            row_sign = rep(1, length(rows)),
                            column_sign = rep(1, length(columns))) {

  if(!is.matrix(x)) stop("Sorting requires matrix parameter blocks.")
  index <- matrix(seq_along(x), nrow = nrow(x), ncol = ncol(x),
                    dimnames = dimnames(x))[rows, columns, drop = FALSE]
  result <- list(index = index, sign = outer(row_sign, column_sign),
                 rows = rows, columns = columns)

  #### Result ####

  return(result)

}

sort_apply_block <- function(x, map, signs = TRUE) {

  result <- x[map$rows, map$columns, drop = FALSE]
  if(signs) result <- result*map$sign

  #### Result ####

  return(result)

}

apply_sort_latent <- function(fit, specification) {

  #### Native estimation coordinates ####

  native <- unsort_latent(fit)
  result <- native
  native_labels <- native@modelInfo$transparameters_labels
  values <- native@Optim$transparameters
  if(length(values) != length(native_labels)) {
    stop("The fitted transformed-parameter vector has incompatible dimensions.")
  }
  names(values) <- native_labels
  templates <- native@modelInfo$trans
  alias_sources <- character(0L)

  for(nm in names(specification$blocks)) {
    if(!(nm %in% names(templates))) stop("A sorted parameter block is missing from the model.")
    map <- specification$blocks[[nm]]
    labels <- sort_apply_block(templates[[nm]], map, signs = FALSE)
    alias_sources <- c(alias_sources, labels[map$sign < 0])
  }
  alias_sources <- unique(alias_sources)
  if(anyNA(match(alias_sources, native_labels))) {
    stop("A factor sign change could not be matched to its native parameter.")
  }
  alias_names <- paste0("sort.neg(", alias_sources, ")")
  while(any(alias_names %in% native_labels)) alias_names <- paste0("sort.", alias_names)
  names(alias_names) <- alias_sources

  #### Reported matrix templates and values ####

  for(nm in names(specification$blocks)) {
    map <- specification$blocks[[nm]]
    labels <- sort_apply_block(templates[[nm]], map, signs = FALSE)
    negative <- map$sign < 0
    labels[negative] <- alias_names[labels[negative]]
    templates[[nm]] <- labels
    result@transformed_pars[[nm]] <- sort_apply_block(native@transformed_pars[[nm]], map)
  }
  if(!is.null(specification$block_order)) {
    block_order <- specification$block_order
    templates <- templates[block_order]
    result@transformed_pars <- result@transformed_pars[block_order]
  }
  result@modelInfo$trans <- templates
  parameter_blocks <- names(native@parameters)
  if(!is.null(specification$block_order)) {
    parameter_blocks <- specification$block_order[specification$block_order %in% parameter_blocks]
  }
  result@parameters <- result@transformed_pars[parameter_blocks]
  result@modelInfo$param <- native@modelInfo$param[parameter_blocks]

  #### Exact negative aliases for subsequent delta-method calculations ####

  if(length(alias_sources) > 0L) {
    input <- match(alias_sources, native_labels)
    output <- length(native_labels)+seq_along(alias_sources)
    names(alias_names) <- NULL
    values <- c(values, setNames(-values[input], alias_names))
    result@modelInfo$transparameters_labels <- names(values)
    result@modelInfo$ntrans <- length(values)
    # The existing column_space transform has an exact constant Jacobian.
    # Bound its workspace: one huge alias block would form a dense m by m
    # Jacobian although the actual transformation is only a sign change.
    chunks <- split(seq_along(input), (seq_along(input)-1L)%/%64L)
    sign_transforms <- lapply(chunks, FUN = function(idx) {
      result <- list(transform = "column_space",
                      indices_in = list(as.integer(input[idx]-1L)),
                      indices_out = list(as.integer(output[idx]-1L)),
                      X = matrix(-1, nrow = 1L, ncol = 1L))
      #### Result ####
      return(result)
    })
    result@modelInfo$control_transform <- c(native@modelInfo$control_transform,
                                            unname(sign_transforms))

    starts <- native@modelInfo$control_optimizer$transparameters
    if(length(starts) > 0L) {
      result@modelInfo$control_optimizer$transparameters <- lapply(starts, FUN = function(x) {
        if(length(x) != length(native_labels)) {
          stop("The original transformed starting values have incompatible dimensions.")
        }
        result <- c(as.numeric(x), -as.numeric(x)[input])
        #### Result ####
        return(result)
      })
    }
    result@modelInfo$control_optimizer$idx_transforms <- NULL
  }
  result@Optim$transparameters <- values

  mapping <- data.frame(label = names(values),
                         source = c(native_labels, alias_sources),
                         sign = c(rep(1, length(native_labels)), rep(-1, length(alias_sources))),
                         stringsAsFactors = FALSE)
  metadata <- c(specification, list(native = native, mapping = mapping,
                                    report_labels = unique(unlist(templates, use.names = FALSE))))
  result@modelInfo$sorting <- metadata

  #### Axis names and cached posterior probabilities ####

  if(specification$type == "factors") {
    for(i in seq_along(specification$order)) {
      result@dataList$factor_label[[i]] <-
        native@dataList$factor_label[[i]][specification$order[[i]]]
    }
  } else {
    index <- specification$order
    result@dataList$class_names <- native@dataList$class_names[index]
    result@dataList$class_reference <- specification$reference_class
    result@dataList$class_order <- index
    sigma_names <- native@dataList$mvgaussian$sigma_names
    if(length(sigma_names) == length(index)) {
      result@dataList$mvgaussian$sigma_names <- sigma_names[index]
    }
    log_posterior <- native@Optim$outputs$estimators$matrices[[1L]][[1L]]
    log_posterior <- matrix(log_posterior, ncol = length(index))[, index, drop = FALSE]
    colnames(log_posterior) <- result@dataList$class_names
    result@Optim$outputs$estimators$matrices[[1L]][[1L]] <- log_posterior
  }

  #### Stored sampling covariance in the reported coordinates ####

  if(length(native@Optim$SE) > 0L) {
    result@Optim$SE <- sort_stored_se_latent(result, native@Optim$SE)
  }

  #### Result ####

  return(result)

}

sort_label_mapping <- function(fit, labels) {

  map <- fit@modelInfo$sorting$mapping
  index <- match(labels, map$label)
  if(anyNA(index)) stop("Unknown parameter requested from a sorted latent object.")
  result <- map[index, , drop = FALSE]

  #### Result ####

  return(result)

}

sort_covariance_latent <- function(V, map) {

  index <- match(map$source, rownames(V))
  if(anyNA(index)) stop("A reported parameter is missing from the source covariance matrix.")
  result <- V[index, index, drop = FALSE]*outer(map$sign, map$sign)
  dimnames(result) <- list(map$label, map$label)

  #### Result ####

  return(result)

}

sort_stored_se_latent <- function(fit, SE) {

  metadata <- fit@modelInfo$sorting
  native <- metadata$native
  V <- SE$VCOV
  if(is.null(V)) V <- SE$vcov
  result <- SE

  if(!is.null(V)) {
    if(!is.matrix(V)) V <- as.matrix(V)
    if(is.null(rownames(V)) && nrow(V) == length(native@modelInfo$parameters_labels)) {
      dimnames(V) <- rep(list(native@modelInfo$parameters_labels), 2L)
    }
    if(is.null(rownames(V)) || !identical(rownames(V), colnames(V))) {
      stop("Stored standard errors cannot be sorted without aligned covariance labels.")
    }
    report <- sort_label_mapping(fit, metadata$report_labels)
    report <- report[report$source %in% rownames(V), , drop = FALSE]
    # Retain native coordinates not represented by any public matrix template.
    extra <- setdiff(rownames(V), report$source)
    if(length(extra) > 0L) report <- rbind(report, sort_label_mapping(fit, extra))
    sorted_V <- sort_covariance_latent(V, report)
    if(!is.null(SE$VCOV)) result$VCOV <- sorted_V
    if(!is.null(SE$vcov)) result$vcov <- sorted_V
    result$se <- sqrt(diag(sorted_V))
    names(result$se) <- report$label
    if(is.null(result$free_VCOV) &&
       setequal(rownames(V), native@modelInfo$parameters_labels)) {
      labs <- native@modelInfo$parameters_labels
      result$free_VCOV <- V[labs, labs, drop = FALSE]
    }
  }

  if(is.list(SE$table_se)) {
    table_se <- SE$table_se
    for(nm in intersect(names(table_se), names(metadata$blocks))) {
      if(is.matrix(table_se[[nm]])) {
        table_se[[nm]] <- sort_apply_block(table_se[[nm]], metadata$blocks[[nm]], signs = FALSE)
      }
    }
    if(!is.null(metadata$block_order)) {
      table_se <- table_se[metadata$block_order[metadata$block_order %in% names(table_se)]]
    }
    result$table_se <- table_se
    blocks <- intersect(names(table_se), names(fit@transformed_pars))
    result$table <- combine_est_se(fit@transformed_pars[blocks], table_se[blocks],
                                    digits = 4L, merge = TRUE)
  } else if(!is.null(result$se)) {
    # A free-coordinate cache does not contain every derived matrix entry.
    # Preserve that requested scope instead of inventing missing matrix SEs.
    result$table_se <- result$se
    estimates <- fit@Optim$transparameters[names(result$se)]
    result$table <- combine_est_se(estimates, result$table_se,
                                    digits = 4L, merge = TRUE)
  }
  if(!is.null(result$jacob)) {
    result$native_jacob <- result$jacob
    result$jacob <- NULL
  }
  result$reporting_coordinates <- "VCOV/se/tables sorted; H/A/B/C/free_VCOV/joint_vcov/steps retain the native fitting chart"

  #### Result ####

  return(result)

}

# Keep native llca charts when reusing a fitted measurement model. Its public
# class ordering is not a request to change fixed labels or regression contrasts.
unsort_model_arguments <- function(model) {

  if(inherits(model, "latent")) {
    result <- unsort_latent(model)
  } else if(is.list(model)) {
    result <- lapply(model, unsort_model_arguments)
  } else {
    result <- model
  }

  #### Result ####

  return(result)

}

# Exact reporting covariance for multistep models. Reconstruction, Hessians,
# KKT constraints and uncertainty propagation are performed in the native chart;
# the selected covariance is then signed and permuted, without approximation.
se_sorted_multistep <- function(fit, parameters = NULL, digits = 4L, ...) {

  metadata <- fit@modelInfo$sorting
  native <- metadata$native
  if(is.null(parameters)) {
    parameters <- fit@modelInfo$trans[names(fit@parameters)]
  }
  labels <- unique(unlist(parameters, use.names = FALSE))
  map <- sort_label_mapping(fit, labels)
  source_labels <- unique(map$source)
  result <- se.multistep(native, parameters = source_labels, digits = digits, ...)
  V <- result$VCOV
  if(is.null(V)) V <- result$vcov
  V <- sort_covariance_latent(as.matrix(V), map)
  result$VCOV <- V
  result$vcov <- V
  result$se <- setNames(sqrt(diag(V)), labels)
  if(is.list(result$table_se) && !is.data.frame(result$table_se)) {
    tables <- result$table_se
    for(nm in intersect(names(tables), names(metadata$blocks))) {
      tables[[nm]] <- sort_apply_block(tables[[nm]], metadata$blocks[[nm]], signs = FALSE)
    }
    result$table_se <- tables
    blocks <- intersect(names(tables), names(fit@transformed_pars))
    result$table <- combine_est_se(fit@transformed_pars[blocks], tables[blocks],
                                    digits = digits, merge = TRUE)
  } else {
    result$table_se <- fill_in(parameters, result$se, miss = NA_real_)
    estimates <- fill_in(parameters, fit@Optim$transparameters, miss = NA_real_)
    result$table <- combine_est_se(estimates, result$table_se, digits = digits, merge = TRUE)
  }
  if(!is.null(result$jacob)) {
    result$native_jacob <- result$jacob
    result$jacob <- NULL
  }
  result$reporting_coordinates <- "VCOV/se/tables sorted; H/A/B/C/free_VCOV/joint_vcov/steps retain the native fitting chart"

  #### Result ####

  return(result)

}

sort_parameter_table_latent <- function(fit, table) {

  if(is.null(fit@modelInfo$sorting)) {
    #### Result ####
    return(table)
  }
  map <- sort_label_mapping(fit, fit@modelInfo$sorting$report_labels)
  map <- map[map$source %in% table$label, , drop = FALSE]
  result <- table[match(map$source, table$label), , drop = FALSE]
  result$label <- map$label
  result$estimate <- result$estimate*map$sign
  if("z" %in% names(result)) result$z <- result$z*map$sign
  rownames(result) <- NULL

  #### Result ####

  return(result)

}

sort_class_diagnostics_latent <- function(fit, diagnostics) {

  index <- fit@modelInfo$sorting$order
  classes <- fit@dataList$class_names
  result <- diagnostics
  for(nm in intersect(names(result), c("Mostlikely.Class", "Avg.Mostlikely"))) {
    result[[nm]] <- result[[nm]][index, index, drop = FALSE]
    dimnames(result[[nm]]) <- list(classes, classes)
  }
  for(nm in intersect(names(result), c("AvePP", "Sum.Posterior", "Sum.Mostlikely"))) {
    if(nrow(result[[nm]]) == length(index)) {
      result[[nm]] <- result[[nm]][index, , drop = FALSE]
      rownames(result[[nm]]) <- classes
    }
  }
  for(nm in intersect(names(result), c("OCC", "Misclassification.per.class"))) {
    result[[nm]] <- setNames(result[[nm]][index], classes)
  }
  attr(result, "class_names") <- classes

  #### Result ####

  return(result)

}

sort_free_labels_latent <- function(fit) {

  free <- fit@modelInfo$parameters_labels
  map <- sort_label_mapping(fit, fit@modelInfo$sorting$report_labels)
  map <- map[map$source %in% free, , drop = FALSE]
  result <- c(map$label, setdiff(free, map$source))

  #### Result ####

  return(result)

}
