# Author: Marcos Jimenez
# Modification date: 13/09/2026

#' Sorted Parameter Structures for Inspection
#'
#' @param fit A fitted llca, lcfa, lefa, or lrotate result.
#'
#' @return A list with four elements: \code{parameters},
#'   \code{transformed_pars}, \code{param}, and \code{trans}. These are copies
#'   of the corresponding estimate and label structures, not a fitted object.
#'   Original labels and factor/class names travel with their entries. The
#'   \code{order} and \code{signs} attributes describe the display permutation
#'   and signs for use by \code{latInspect()}; nothing is stored in fit.
#'
#' @details
#' Classes are ordered by decreasing frequency-weighted posterior size.
#' Factors are ordered by decreasing \code{colSums(lambda^2)*diag(psi)} and
#' oriented so their largest absolute loading is positive. Ties retain the
#' original order. Fixed modeled CFA loadings retain their specified order and
#' signs. Target rotations and components selecting factors retain factor order.
#' The original multinomial reference class and all parameter labels are kept.
#'
#' @export
sort_latent <- function(fit) {

  if(!inherits(fit, "latent")) {
    stop("fit must inherit from class 'latent'.")
  }

  result <- list(parameters = fit@parameters,
                 transformed_pars = fit@transformed_pars,
                 param = fit@modelInfo$param,
                 trans = fit@modelInfo$trans)
  orderings <- NULL
  signs <- list()

  if(length(fit@transformed_pars) == 0L) {
    #### Result ####
    return(result)
  }

  if(inherits(fit, "llca")) {

    q <- ncol(fit@transformed_pars$class)
    weights <- as.numeric(fit@dataList$pattern_weights)
    posterior <- fit@Optim$outputs$estimators$matrices[[1L]][[1L]]
    if(length(posterior) != length(weights)*q || any(!is.finite(weights)) ||
       any(weights < 0) || sum(weights) <= 0) {
      stop("Invalid posterior probabilities or pattern weights for sorting.")
    }
    posterior <- exp(matrix(posterior, nrow = length(weights), ncol = q))
    if(any(!is.finite(posterior)) || any(abs(rowSums(posterior)-1) > 1e-06)) {
      stop("Posterior probabilities must be finite and sum to one.")
    }
    index <- order(-colSums(posterior*weights), seq_len(q))
    orderings <- index
    blocks <- unique(c("beta", "theta", "class", "loglik",
                        fit@dataList$gaussian$gaussian_names,
                        fit@dataList$mvgaussian$mvgaussian_names,
                        fit@dataList$multinomial$multinomial_names,
                        names(fit@dataList$mvmultinomial$indep_pairs_list)))

    for(part in names(result)) {
      for(nm in intersect(blocks, names(result[[part]]))) {
        x <- result[[part]][[nm]]
        if(!is.matrix(x) || ncol(x) != q) {
          stop("A class-specific parameter block has incompatible dimensions.")
        }
        result[[part]][[nm]] <- x[, index, drop = FALSE]
      }
      # Covariance blocks have a class-specific name, not a class axis.
      sigma <- fit@dataList$mvgaussian$sigma_names
      if(length(sigma) == q) {
        block_names <- names(result[[part]])
        present <- intersect(sigma, block_names)
        block_names[match(present, block_names)] <- intersect(sigma[index], present)
        result[[part]] <- result[[part]][block_names]
      }
    }

  } else {

    rotation <- !is.null(fit@modelInfo$data_param$X_group)
    data_param <- if(rotation) fit@modelInfo$data_param else fit@dataList$data_param
    if(is.null(data_param$lambda_group)) {
      #### Result ####
      return(result)
    }
    # Fixed markers are identification choices, not reporting indeterminacies.
    if(!rotation && inherits(fit@dataList$LAV, "lavaan")) {
      tab <- lavaan::parTable(fit@dataList$LAV)
      if(any(tab$op == "=~" & !is.na(tab$free) & tab$free == 0L)) {
        #### Result ####
        return(result)
      }
    }

    reorder <- TRUE
    if(rotation) {
      specification <- fit@dataList$rotation_spec
      if(is.null(specification)) specification <- fit@dataList$rotation
      criteria <- if(is.character(specification)) specification else names(specification)
      reorder <- !any(criteria %in% c("target", "xtarget"))
      components <- if(is.character(specification)) rep(list(list()), length(criteria)) else specification
      for(component in components) {
        factors <- if("factors" %in% names(component)) component$factors else fit@dataList$args$factors
        if(!is.null(factors)) reorder <- FALSE
      }
    }

    orderings <- vector("list", length(data_param$lambda_group))
    for(i in seq_along(data_param$lambda_group)) {
      lambda_name <- data_param$lambda_group[i]
      psi_name <- data_param$psi_group[i]
      alpha_name <- data_param$alpha_group[i]
      mu_name <- data_param$mu_group[i]
      lambda <- fit@transformed_pars[[lambda_name]]
      psi <- fit@transformed_pars[[psi_name]]
      q <- ncol(lambda)
      if(!is.matrix(lambda) || !is.numeric(lambda) || any(!is.finite(lambda)) ||
         !is.matrix(psi) || !is.numeric(psi) || any(!is.finite(psi)) ||
         length(q) != 1L || q < 1L || !identical(dim(psi), c(q, q))) {
        stop("Invalid factor loading/covariance blocks for sorting.")
      }
      if(any(diag(psi) < 0)) {
        warning("A negative factor variance was found; sorting does not repair the solution.")
      }
      contribution <- colSums(lambda^2)*diag(psi)
      if(any(!is.finite(contribution))) stop("Non-finite factor contributions.")
      index <- if(reorder) order(-contribution, seq_len(q)) else seq_len(q)
      largest <- lambda[cbind(max.col(t(abs(lambda)), ties.method = "first"), seq_len(q))]
      direction <- ifelse(largest[index] < 0, -1, 1)
      orderings[[i]] <- index
      X_name <- if(rotation) data_param$X_group[i] else data_param$xpsi_group[i]
      Xinv_name <- if(rotation) data_param$Xinv_group[i] else character(0L)
      blocks <- intersect(c(lambda_name, psi_name, alpha_name, mu_name,
                            X_name, Xinv_name),
                          names(result$transformed_pars))

      for(nm in blocks) {
        x <- result$transformed_pars[[nm]]
        rows <- if(nm %in% c(psi_name, alpha_name, mu_name, Xinv_name)) {
          index
        } else {
          seq_len(nrow(x))
        }
        columns <- if(nm %in% c(lambda_name, psi_name, X_name)) index else seq_len(ncol(x))
        row_sign <- if(nm %in% c(psi_name, alpha_name, mu_name, Xinv_name)) {
          direction
        } else {
          rep(1, length(rows))
        }
        column_sign <- if(nm %in% c(lambda_name, psi_name, X_name)) direction else rep(1, length(columns))
        multiplier <- outer(row_sign, column_sign)
        signs[[nm]] <- multiplier

        for(part in names(result)) {
          if(!(nm %in% names(result[[part]]))) next
          x <- result[[part]][[nm]][rows, columns, drop = FALSE]
          if(part %in% c("parameters", "transformed_pars")) {
            x <- x*multiplier
          } else if(part == "param") {
            # Sign numeric fixed values, never turn labels into new aliases.
            numeric_values <- suppressWarnings(as.numeric(x))
            fixed <- is.finite(numeric_values) & c(multiplier) < 0
            x[fixed] <- -numeric_values[fixed]
          }
          result[[part]][[nm]] <- x
        }
      }
    }

  }

  attr(result, "order") <- orderings
  attr(result, "signs") <- signs

  #### Result ####

  return(result)

}
