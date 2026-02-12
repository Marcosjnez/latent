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
