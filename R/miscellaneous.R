# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 16/07/2026

# Miscellaneous functions used in latent

allnumeric <- function(lst) {

  # Transform all elements of a list into numeric values
  lst <- rapply(
    lst,
    function(x) {
      if (is.matrix(x) || is.array(x)) { storage.mode(x) <- "double"; x } # keep as matrix or array
      else if (is.factor(x)) as.numeric(as.character(x))             # factors -> numeric values
      else if (is.atomic(x)) as.numeric(x)                           # vectors -> numeric
      else x
    },
    how = "replace",
    classes = c("matrix","array","factor","numeric","integer","logical","character")
  )

  return(lst)

}

insert_object <- function(X, Y) {
  # Ensure matrix type
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  # Get dimnames
  rx <- rownames(X)
  cx <- colnames(X)
  ry <- rownames(Y)
  cy <- colnames(Y)

  # We need names to match positions
  if (is.null(ry) || is.null(cy)) {
    stop("Y must have both rownames and colnames to be inserted by name.")
  }
  if (is.null(rx) || is.null(cx)) {
    stop("X must have both rownames and colnames to insert Y by name.")
  }

  # Match row and column names
  row_idx <- match(ry, rx)
  col_idx <- match(cy, cx)

  # Check that all names from Y exist in X
  if (any(is.na(row_idx))) {
    missing_rows <- ry[is.na(row_idx)]
    stop("These row names from Y are not present in X: ",
         paste(missing_rows, collapse = ", "))
  }
  if (any(is.na(col_idx))) {
    missing_cols <- cy[is.na(col_idx)]
    stop("These column names from Y are not present in X: ",
         paste(missing_cols, collapse = ", "))
  }

  # Insert (overwrite) values
  X[row_idx, col_idx] <- Y

  X
}

fill_in <- function(lst, values, miss = NA) {
  if (!is.atomic(values) || is.null(names(values)) ||
      anyNA(names(values)) || any(names(values) == "")) {
    stop("`values` must be a named atomic vector.")
  }

  make_miss_like <- function(x) {
    out <- values[rep(NA_integer_, length(x))]
    out[] <- miss

    if (!is.null(dim(x))) {
      dim(out) <- dim(x)
      dimnames(out) <- dimnames(x)
    } else {
      names(out) <- names(x)
    }

    out
  }

  fill_leaf <- function(x) {
    if (is.numeric(x)) {
      return(make_miss_like(x))
    }

    if (!is.character(x)) {
      stop("All terminal objects in `lst` must be character or numeric.")
    }

    x_vec <- as.vector(x)
    idx <- match(x_vec, names(values))

    out <- make_miss_like(x)

    hit <- !is.na(idx)
    out[hit] <- unname(values[idx[hit]])

    out
  }

  recurse <- function(x) {
    if (is.list(x) && !is.data.frame(x)) {
      x[] <- lapply(x, recurse)
      return(x)
    }

    if (is.atomic(x)) {
      return(fill_leaf(x))
    }

    stop("Unsupported type found in `lst`.")
  }

  recurse(lst)

}

extract_unique_values <- function(chars, vals, miss = 0) {

  coerce_numeric_leaf <- function(v, miss) {
  }

  align_leaf <- function(ch, vl, miss) {
  }

  walk <- function(ch, vl = NULL) {
  }

  flat <- walk(chars, vals)

  u <- unique(flat$chars)
  out <- rep(miss, length(u))
  names(out) <- u

  idx_split <- split(seq_along(flat$chars), flat$chars)

  for (nm in names(idx_split)) {
    idx <- idx_split[[nm]]
    first_ok <- match(TRUE, flat$ok[idx])
    out[nm] <- if (!is.na(first_ok)) flat$vals[idx[first_ok]] else miss
  }

  out
}

fixed_values_indices <- function(x_list, y_list) {
  # ----- checks -----
  if (!is.list(x_list) || !is.list(y_list)) {
    stop("Both inputs must be lists.")
  }

  nx <- names(x_list)
  ny <- names(y_list)

  if (is.null(nx) || any(nx == "") || anyDuplicated(nx)) {
    stop("The first list must have unique, non-empty names.")
  }
  if (is.null(ny) || any(ny == "") || anyDuplicated(ny)) {
    stop("The second list must have unique, non-empty names.")
  }

  # Only shared objects are used
  shared_names <- intersect(nx, ny)

  if (length(shared_names) == 0) {
    return(list(
      numerical_values = setNames(numeric(0), character(0)),
      indices = integer(0)
    ))
  }

  # Keep the order of the first list
  shared_names <- nx[nx %in% shared_names]

  # Helper: flatten object ignoring names/dimnames
  flat <- function(z) {
    if (is.factor(z)) z <- as.character(z)
    as.vector(unname(z))
  }

  # Helper: detect numeric or numeric-like character
  is_numeric_like <- function(z) {
    if (is.factor(z)) z <- as.character(z)

    if (is.numeric(z)) {
      return(!is.na(z))
    }

    if (is.character(z)) {
      z <- trimws(z)
      return(!is.na(suppressWarnings(as.numeric(z))))
    }

    rep(FALSE, length(z))
  }

  # Helper: convert numeric / numeric-like character to numeric
  to_numeric <- function(z) {
    if (is.factor(z)) z <- as.character(z)

    if (is.numeric(z)) {
      return(as.numeric(z))
    }

    if (is.character(z)) {
      return(suppressWarnings(as.numeric(trimws(z))))
    }

    rep(NA_real_, length(z))
  }

  # Unique characters from the full first list
  x_all_total <- unlist(
    lapply(x_list, function(obj) as.character(flat(obj))),
    use.names = FALSE
  )
  unique_chars_total <- unique(x_all_total)

  # Collect matched labels and numeric values from shared objects only
  x_all_shared <- character(0)
  y_all_shared <- numeric(0)

  for (nm in shared_names) {
    x_obj <- x_list[[nm]]
    y_obj <- y_list[[nm]]

    # Check matching structure only for shared objects
    if (!identical(dim(x_obj), dim(y_obj)) || length(x_obj) != length(y_obj)) {
      stop(sprintf("Shared objects '%s' do not match in structure.", nm))
    }

    x_vec <- as.character(flat(x_obj))
    y_vec <- flat(y_obj)

    num_mask <- is_numeric_like(y_vec)

    x_all_shared <- c(x_all_shared, x_vec[num_mask])
    y_all_shared <- c(y_all_shared, to_numeric(y_vec[num_mask]))
  }

  if (length(x_all_shared) == 0) {
    return(list(
      numerical_values = setNames(numeric(0), character(0)),
      indices = integer(0)
    ))
  }

  # Map label -> numeric value
  split_vals <- split(y_all_shared, x_all_shared)

  # Check that each label corresponds to only one unique numeric value
  inconsistent <- vapply(split_vals, function(v) length(unique(v)) > 1, logical(1))
  if (any(inconsistent)) {
    bad <- names(split_vals)[inconsistent]
    stop(
      sprintf(
        "Some characters are matched to more than one numerical value: %s",
        paste(bad, collapse = ", ")
      )
    )
  }

  # Preserve the order from the unique characters of the first list
  matched_labels <- unique_chars_total[unique_chars_total %in% names(split_vals)]
  matched_values <- vapply(matched_labels, function(ch) split_vals[[ch]][1], numeric(1))
  matched_indices <- match(matched_labels, unique_chars_total)

  list(
    numerical_values = setNames(matched_values, matched_labels),
    indices = matched_indices
  )
}

create_init <- function(trans, param, init_param,
                        control_transform, control) {

  idx_transformed <- unlist(lapply(control_transform,
                                   FUN = \(x) unlist(x$indices_out)+1L))

  # ---------- helpers ----------
  check_named_list <- function(x, arg_name, allow_empty = TRUE) {
    if (!is.list(x)) {
      stop(sprintf("'%s' must be a list.", arg_name))
    }
    if (length(x) == 0) {
      if (allow_empty) return(invisible(NULL))
      stop(sprintf("'%s' cannot be empty.", arg_name))
    }
    nms <- names(x)
    if (is.null(nms) || any(nms == "") || anyDuplicated(nms)) {
      stop(sprintf("'%s' must have unique, non-empty names.", arg_name))
    }
    invisible(NULL)
  }

  flat <- function(z) {
    if (is.factor(z)) z <- as.character(z)
    as.vector(unname(z))
  }

  list_to_char_vector <- function(x) {
    unlist(lapply(x, function(z) as.character(flat(z))), use.names = FALSE)
  }

  is_numeric_like <- function(z) {
    if (is.factor(z)) z <- as.character(z)

    if (is.numeric(z)) {
      return(!is.na(z))
    }

    if (is.character(z)) {
      z <- trimws(z)
      return(!is.na(suppressWarnings(as.numeric(z))))
    }

    rep(FALSE, length(z))
  }

  to_numeric <- function(z) {
    if (is.factor(z)) z <- as.character(z)

    if (is.numeric(z)) {
      return(as.numeric(z))
    }

    if (is.character(z)) {
      return(suppressWarnings(as.numeric(trimws(z))))
    }

    rep(NA_real_, length(z))
  }

  extract_numeric_by_label <- function(reference, values, miss = NA_real_) {
    ref_labels <- unname(unique(list_to_char_vector(reference)))
    out <- setNames(rep(miss, length(ref_labels)), ref_labels)

    if (length(values) == 0) {
      return(out)
    }

    shared_names <- intersect(names(reference), names(values))
    if (length(shared_names) == 0) {
      return(out)
    }
    shared_names <- names(reference)[names(reference) %in% shared_names]

    lbl <- character(0)
    val <- numeric(0)

    for (nm in shared_names) {
      ref_obj <- reference[[nm]]
      val_obj <- values[[nm]]

      if (!identical(dim(ref_obj), dim(val_obj)) || length(ref_obj) != length(val_obj)) {
        stop(sprintf("Objects '%s' do not match in structure between reference and values.", nm))
      }

      ref_vec <- as.character(flat(ref_obj))
      val_vec <- flat(val_obj)

      num_mask <- is_numeric_like(val_vec)
      if (any(num_mask)) {
        lbl <- c(lbl, ref_vec[num_mask])
        val <- c(val, to_numeric(val_vec[num_mask]))
      }
    }

    if (length(lbl) == 0) {
      return(out)
    }

    # Keep only the first numeric match for each label
    first_idx <- match(ref_labels, lbl)
    hit <- !is.na(first_idx)

    out[hit] <- val[first_idx[hit]]
    out
  }

  # ---------- checks ----------
  check_named_list(trans, "trans", allow_empty = FALSE)
  check_named_list(param, "param", allow_empty = TRUE)

  if (!is.list(init_param)) {
    stop("'init_param' must be a list.")
  }

  if (is.null(control$rstarts) || length(control$rstarts) != 1 || !is.numeric(control$rstarts)) {
    stop("'control$rstarts' must be a single numeric value.")
  }
  rstarts <- as.integer(control$rstarts)
  if (rstarts < 1L) {
    stop("'control$rstarts' must be >= 1.")
  }

  if (length(init_param) != rstarts) {
    stop("length(init_param) must be equal to control$rstarts.")
  }

  bad_param_names <- setdiff(names(param), names(trans))
  if (length(bad_param_names) > 0) {
    stop(
      sprintf(
        "'param' contains objects not present in 'trans': %s",
        paste(bad_param_names, collapse = ", ")
      )
    )
  }

  # ---------- labels from trans ----------
  transparameters_labels <- unname(unique(list_to_char_vector(trans)))
  ntrans <- length(transparameters_labels)

  if (length(idx_transformed) > 0) {
    if (!is.numeric(idx_transformed) || any(is.na(idx_transformed))) {
      stop("'idx_transformed' must be a numeric/integer index vector.")
    }
    idx_transformed <- unique(as.integer(idx_transformed))
    if (any(idx_transformed < 1L | idx_transformed > ntrans)) {
      stop("'idx_transformed' contains indices out of range.")
    }
  } else {
    idx_transformed <- integer(0)
  }

  # ---------- fixed values from param ----------
  fixed_full <- extract_numeric_by_label(trans, param, miss = NA_real_)
  fixed_idx <- which(!is.na(fixed_full))

  drop_idx <- sort(unique(c(idx_transformed, fixed_idx)))
  keep_idx <- setdiff(seq_len(ntrans), drop_idx)
  parameters_labels <- transparameters_labels[keep_idx]

  # ---------- build outputs ----------
  parameters <- vector("list", rstarts)
  transparameters <- vector("list", rstarts)

  for (i in seq_len(rstarts)) {
    check_named_list(init_param[[i]], sprintf("init_param[[%d]]", i), allow_empty = TRUE)

    bad_init_names <- setdiff(names(init_param[[i]]), names(trans))
    if (length(bad_init_names) > 0) {
      stop(
        sprintf(
          "'init_param[[%d]]' contains objects not present in 'trans': %s",
          i, paste(bad_init_names, collapse = ", ")
        )
      )
    }

    init_full <- extract_numeric_by_label(trans, init_param[[i]], miss = NA_real_)

    # Start from init values, fill missing with 0, then enforce fixed values
    trans_vec <- init_full
    trans_vec[is.na(trans_vec)] <- 0
    if (length(fixed_idx) > 0) {
      trans_vec[fixed_idx] <- fixed_full[fixed_idx]
    }

    names(trans_vec) <- transparameters_labels
    transparameters[[i]] <- trans_vec

    param_vec <- trans_vec[keep_idx]
    names(param_vec) <- parameters_labels
    parameters[[i]] <- param_vec
  }

  parameters_labels <- names(parameters[[1]])
  nparam <- length(parameters_labels)

  transparameters_labels <- names(transparameters[[1]])
  ntrans <- length(transparameters_labels)

  trans2param <- match(parameters_labels, transparameters_labels)

  list(
    parameters = parameters,
    transparameters = transparameters,
    parameters_labels = parameters_labels,
    transparameters_labels = transparameters_labels,
    nparam = nparam,
    ntrans = ntrans,
    trans2param = trans2param
  )
}

group_lists_by_sublists <- function(...) {
  # Transpose lists of objects with same names
  x <- list(...)
  nm <- names(x[[1]])

  setNames(
    lapply(nm, function(n) lapply(x, `[[`, n)),
    nm
  )
}

combine_est_se <- function(est, se, digits = 2, merge = TRUE) {

  est_suffix <- " (est)"
  se_suffix <- " (se)"

  format_object <- function(x) {

    if(is.null(digits)) {
      result <- x
    } else {
      result <- round(x, digits = digits)
    }

    #### Result ####

    return(result)

  }

  as_parameter_matrix <- function(x) {

    if(is.null(dim(x))) {

      result <- matrix(unclass(x), ncol = 1L,
                       dimnames = list(names(x), NULL))

    } else {

      result <- as.matrix(unclass(x))

    }

    #### Result ####

    return(result)

  }

  align_se <- function(E, S, object_name = NULL) {

    if(!identical(dim(S), dim(E))) {

      if(length(S) == length(E)) {

        S <- matrix(S, nrow = nrow(E), ncol = ncol(E),
                    dimnames = dimnames(E))

      } else if(ncol(S) == 1L && nrow(S) == nrow(E) &&
                ncol(E) > 1L) {

        S <- matrix(rep(S, ncol(E)), nrow = nrow(E), ncol = ncol(E))

      } else if(nrow(S) == 1L && ncol(S) == ncol(E) &&
                nrow(E) > 1L) {

        S <- matrix(rep(S, each = nrow(E)),
                    nrow = nrow(E), ncol = ncol(E))

      } else if(all(dim(S) == rev(dim(E)))) {

        S <- t(S)

      } else {

        object_label <- if(is.null(object_name)) "" else
          paste0(" for object '", object_name, "'")

        stop(sprintf("Dimension mismatch%s: est = %s, se = %s",
                     object_label, paste(dim(E), collapse = "x"),
                     paste(dim(S), collapse = "x")))

      }

    }

    #### Result ####

    return(S)

  }

  combine_single_est_se <- function(E_object, S_object,
                                    object_name = NULL) {

    est_is_vector <- is.null(dim(E_object))

    E <- as_parameter_matrix(E_object)
    S <- as_parameter_matrix(S_object)
    S <- align_se(E, S, object_name = object_name)

    E_fmt <- format_object(E)
    S_fmt <- format_object(S)

    if(merge) {

      M <- matrix(paste0(E_fmt, " (", S_fmt, ")"),
                  nrow = nrow(E), ncol = ncol(E),
                  dimnames = dimnames(E))

      if(est_is_vector) {
        result <- as.vector(M)
        names(result) <- names(E_object)
      } else {
        result <- M
      }

    } else if(est_is_vector) {

      result <- cbind(est = as.vector(E_fmt), se = as.vector(S_fmt))
      rownames(result) <- names(E_object)

    } else {

      nr <- nrow(E)
      nc <- ncol(E)

      result <- matrix(NA, nrow = nr, ncol = 2L*nc,
                       dimnames = list(rownames(E), NULL))

      result[, seq(1L, 2L*nc, by = 2L)] <- E_fmt
      result[, seq(2L, 2L*nc, by = 2L)] <- S_fmt

      cls <- colnames(E)
      if(is.null(cls)) {
        cls <- paste0("col", seq_len(nc))
      }

      colnames(result) <- as.vector(rbind(paste0(cls, est_suffix),
                                          paste0(cls, se_suffix)))

    }

    #### Result ####

    return(result)

  }

  if(is.list(est) || is.list(se)) {

    stopifnot(is.list(est), is.list(se))

    if(is.null(names(est)) || is.null(names(se))) {

      if(length(est) != length(se)) {
        stop("Unnamed lists 'est' and 'se' must have the same length.")
      }

      item_names <- seq_along(est)
      result <- vector("list", length(est))

      for(i in item_names) {
        result[[i]] <- combine_single_est_se(est[[i]], se[[i]])
      }

      names(result) <- names(est)

    } else {

      item_names <- intersect(names(est), names(se))
      result <- vector("list", length(item_names))
      names(result) <- item_names

      for(nm in item_names) {
        result[[nm]] <- combine_single_est_se(est[[nm]], se[[nm]],
                                              object_name = nm)
      }

    }

  } else {

    result <- combine_single_est_se(est, se)

  }

  #### Result ####

  return(result)

}

combine_est_ci <- function(lower, est, upper, digits = 2, merge = TRUE) {

  fmt <- function(x) {

    if(is.null(digits)) {
      result <- as.character(x)
    } else {
      result <- formatC(round(x, digits = digits),
                        format = "f", digits = digits)
    }

    #### Result ####

    return(result)

  }

  as_parameter_matrix <- function(x) {

    if(is.null(dim(x))) {

      result <- matrix(unclass(x), ncol = 1L,
                       dimnames = list(names(x), NULL))

    } else {

      result <- as.matrix(unclass(x))

    }

    #### Result ####

    return(result)

  }

  align_bound <- function(E, B, object_name = NULL) {

    if(!identical(dim(B), dim(E))) {

      if(length(B) == length(E)) {

        B <- matrix(B, nrow = nrow(E), ncol = ncol(E),
                    dimnames = dimnames(E))

      } else if(ncol(B) == 1L && nrow(B) == nrow(E) &&
                ncol(E) > 1L) {

        B <- matrix(rep(B, ncol(E)), nrow = nrow(E), ncol = ncol(E))

      } else if(nrow(B) == 1L && ncol(B) == ncol(E) &&
                nrow(E) > 1L) {

        B <- matrix(rep(B, each = nrow(E)),
                    nrow = nrow(E), ncol = ncol(E))

      } else if(all(dim(B) == rev(dim(E)))) {

        B <- t(B)

      } else {

        object_label <- if(is.null(object_name)) "" else
          paste0(" for object '", object_name, "'")

        stop(sprintf("Dimension mismatch%s: est = %s, bound = %s",
                     object_label, paste(dim(E), collapse = "x"),
                     paste(dim(B), collapse = "x")))

      }

    }

    #### Result ####

    return(B)

  }

  combine_single_ci <- function(L_object, E_object, U_object,
                                object_name = NULL) {

    est_is_vector <- is.null(dim(E_object))

    E <- as_parameter_matrix(E_object)
    L <- as_parameter_matrix(L_object)
    U <- as_parameter_matrix(U_object)

    L <- align_bound(E, L, object_name = object_name)
    U <- align_bound(E, U, object_name = object_name)

    if(merge) {

      result <- matrix(paste0(fmt(E), " (", fmt(L), "|", fmt(U), ")"),
                       nrow = nrow(E), ncol = ncol(E),
                       dimnames = dimnames(E))

    } else {

      result <- matrix(paste0(fmt(E), " [", fmt(L), ", ", fmt(U), "]"),
                       nrow = nrow(E), ncol = ncol(E),
                       dimnames = dimnames(E))

    }

    if(est_is_vector) {
      result <- as.vector(result)
      names(result) <- names(E_object)
    }

    #### Result ####

    return(result)

  }

  if(is.list(est) || is.list(lower) || is.list(upper)) {

    stopifnot(is.list(est), is.list(lower), is.list(upper))

    if(is.null(names(est)) || is.null(names(lower)) ||
       is.null(names(upper))) {

      if(length(est) != length(lower) || length(est) != length(upper)) {
        stop("Unnamed lists 'lower', 'est', and 'upper' must have the same length.")
      }

      result <- vector("list", length(est))

      for(i in seq_along(est)) {
        result[[i]] <- combine_single_ci(lower[[i]], est[[i]], upper[[i]])
      }

      names(result) <- names(est)

    } else {

      item_names <- Reduce(intersect,
                           list(names(est), names(lower), names(upper)))

      result <- vector("list", length(item_names))
      names(result) <- item_names

      for(nm in item_names) {
        result[[nm]] <- combine_single_ci(lower[[nm]], est[[nm]],
                                          upper[[nm]], object_name = nm)
      }

    }

  } else {

    result <- combine_single_ci(lower, est, upper)

  }

  #### Result ####

  return(result)

}

effects_coding <- function(coeffs, vcov) {

  # Transform from zero constraints to effects-coding constraints

  p <- nrow(coeffs)
  q <- ncol(coeffs)

  # Centering matrix:
  C <- diag(q)-matrix(1/q, q, q)

  # Effects-coding parameterization:
  beta_new <- coeffs %*% C

  # Jacobian of coeffs %*% C:
  J <- kronecker(t(C), diag(p))

  # Variance-covariance matrix under effects coding:
  vcov_new <- J %*% vcov %*% t(J)
  se_new <- sqrt(diag(vcov_new))

  table_se <- matrix(se_new, nrow = p, ncol = q,
                     dimnames = dimnames(coeffs))

  #### Result ####

  result <- list(beta = beta_new, se = se_new,
                 vcov = vcov_new, table_se = table_se)

  return(result)

}

move_intercept <- function(coeffs, vcov, column = 1L) {

  # Change the reference column of a zero-constrained parameterization

  p <- nrow(coeffs)
  q <- ncol(coeffs)

  # Current and new free columns:
  current_reference <- which(apply(coeffs, MARGIN = 2,
                                   FUN = \(x) all(x == 0)))[1L]
  old_free <- setdiff(seq_len(q), current_reference)
  new_free <- setdiff(seq_len(q), column)

  # Subtract the coefficient of the new reference class:
  beta_new <- sweep(coeffs, MARGIN = 1,
                    STATS = coeffs[, column], FUN = "-")

  # Transformation from the old free coefficients to the new free coefficients:
  M <- diag(q)
  M[column, ] <- M[column, ]-1

  T <- M[old_free, new_free, drop = FALSE]
  J <- kronecker(t(T), diag(p))

  vcov_new <- J %*% vcov %*% t(J)
  se_new <- sqrt(diag(vcov_new))

  table_se <- matrix(0, nrow = p, ncol = q,
                     dimnames = dimnames(coeffs))
  table_se[, new_free] <- matrix(se_new, nrow = p,
                                 ncol = q-1L)

  #### Result ####

  result <- list(beta = beta_new, se = se_new,
                 vcov = vcov_new, table_se = table_se)

  return(result)

}
