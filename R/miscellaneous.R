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

fill_list_with_vector <- function(lst, values) {

  # Insert a vector of values in an arbitrary list lst

  i <- 1

  assign_recursive <- function(x) {
    if (is.list(x)) {
      lapply(x, assign_recursive)
    } else if (is.matrix(x)) {
      dims <- dim(x)
      n <- prod(dims)
      x[] <- values[i:(i + n - 1)]
      i <<- i + n
      x
    } else if (is.atomic(x)) {
      n <- length(x)
      x[] <- values[i:(i + n - 1)]
      i <<- i + n
      x
    } else {
      stop("Unsupported type")
    }
  }

  assign_recursive(lst)

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

  fill_leaf <- function(x) {
    if (!is.character(x)) {
      stop("All terminal objects in `lst` must be character.")
    }

    x_vec <- as.vector(x)
    idx   <- match(x_vec, names(values))

    # initialize output with the same atomic type as `values`
    out <- values[rep(NA_integer_, length(x_vec))]
    out[] <- miss

    hit <- !is.na(idx)
    out[hit] <- unname(values[idx[hit]])

    # restore original structure
    if (!is.null(dim(x))) {
      dim(out) <- dim(x)
      dimnames(out) <- dimnames(x)
    } else {
      names(out) <- names(x)
    }

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
    # Return a numeric vector. Non-numeric entries become `miss`.
    if (is.null(v) || is.list(v) || !is.atomic(v)) {
      return(NULL)
    }

    if (is.numeric(v)) {
      return(as.vector(v))
    }

    if (is.character(v)) {
      out <- suppressWarnings(as.numeric(v))
      out[is.na(out)] <- miss
      return(out)
    }

    # Anything else does not contribute numeric values
    NULL
  }

  align_leaf <- function(ch, vl, miss) {
    n <- length(ch)
    out <- rep(miss, n)
    ok  <- rep(FALSE, n)

    v <- coerce_numeric_leaf(vl, miss)

    if (!is.null(v)) {
      m <- min(n, length(v))
      if (m > 0L) {
        out[seq_len(m)] <- v[seq_len(m)]
        ok[seq_len(m)]  <- !is.na(v[seq_len(m)])
      }
    }

    list(vals = out, ok = ok)
  }

  walk <- function(ch, vl = NULL) {
    if (is.list(ch) && !is.data.frame(ch)) {
      nms_ch <- names(ch)
      if (is.null(nms_ch)) nms_ch <- as.character(seq_along(ch))

      pieces <- lapply(seq_along(ch), function(i) {
        nm <- nms_ch[i]
        vl_i <- NULL

        if (!is.null(vl) && is.list(vl) && !is.data.frame(vl)) {
          nms_vl <- names(vl)

          if (!is.null(nms_vl)) {
            if (nm %in% nms_vl) vl_i <- vl[[nm]]
          } else if (length(vl) >= i) {
            vl_i <- vl[[i]]
          }
        }

        walk(ch[[i]], vl_i)
      })

      return(list(
        chars = unlist(lapply(pieces, `[[`, "chars"), use.names = FALSE),
        vals  = unlist(lapply(pieces, `[[`, "vals"),  use.names = FALSE),
        ok    = unlist(lapply(pieces, `[[`, "ok"),    use.names = FALSE)
      ))
    }

    if (!is.character(ch)) {
      stop("All leaves in `chars` must be character.")
    }

    tmp <- align_leaf(ch, vl, miss)
    list(
      chars = as.vector(ch),
      vals  = tmp$vals,
      ok    = tmp$ok
    )
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
