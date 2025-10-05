combine_est_se <- function(est, se, digits = 2) {

  est_suffix <- " (est)"
  se_suffix <- " (se)"

  format_matrix <- function(mat) {
    if (is.null(digits)) return(mat)
    fmt <- function(x) round(x, digits = digits)
    apply(mat, c(1, 2), fmt)
  }

  # Combine class probabilities
  # if (!is.null(est$classes) && !is.null(se$classes)) {
  #   classes <- rbind(est$classes, se$classes)
  #   rownames(classes) <- c("est", "se")
  #   colnames(classes) <- names(est$classes)
  #   classes <- t(classes)
  #   classes <- format_matrix(classes)
  # } else {
  #   classes <- NULL
  # }

  # Combine item parameters
  combine_item_est_se <- function(est_items, se_items) {
    stopifnot(is.list(est_items), is.list(se_items))
    item_names <- intersect(names(est_items), names(se_items))
    out <- vector("list", length(item_names))
    names(out) <- item_names

    for (nm in item_names) {
      E <- est_items[[nm]]
      S <- se_items[[nm]]

      E <- as.matrix(unclass(E))
      S <- as.matrix(unclass(S))

      nr <- nrow(E)
      nc <- ncol(E)

      if (!identical(dim(S), dim(E))) {
        if (ncol(S) == 1 && ncol(E) > 1) {
          S <- matrix(rep(S, ncol(E)), nrow = nr, ncol = nc)
        } else if (ncol(S) == ncol(E) && nrow(S) == 1) {
          S <- matrix(rep(S, each = nr), nrow = nr, ncol = nc)
        } else if (all(dim(S) == rev(dim(E)))) {
          S <- t(S)
        } else {
          stop(sprintf("Dimension mismatch for item '%s': est = %s, se = %s",
                       nm, paste(dim(E), collapse = "x"), paste(dim(S), collapse = "x")))
        }
      }

      M <- matrix(NA, nrow = nr, ncol = 2 * nc,
                  dimnames = list(rownames(E), NULL))

      E_fmt <- format_matrix(E)
      S_fmt <- format_matrix(S)

      M[, seq(1, 2 * nc, by = 2)] <- E_fmt
      M[, seq(2, 2 * nc, by = 2)] <- S_fmt

      cls <- colnames(E); if (is.null(cls)) cls <- paste0("col", seq_len(nc))
      colnames(M) <- as.vector(rbind(paste0(cls, est_suffix),
                                     paste0(cls, se_suffix)))
      out[[nm]] <- M
    }

    return(out)
  }

  classes <- combine_item_est_se(est[1], se[1])

  items <- combine_item_est_se(est$items, se$items)

  return(list(classes = classes, items = items))
}

combine_est_ci <- function(lower, est, upper, digits = 2) {

  fmt <- function(x) formatC(round(x, digits = digits),
                             format = "f", digits = digits)

  LEU <- function(E, L, U) {
    matrix(paste0(
      fmt(E), " [", fmt(L), ", ", fmt(U), "]"
    ), nrow = nrow(E), ncol = ncol(E),
    dimnames = dimnames(E))
  }

  # Combine class-level estimates
  classes <- LEU(est[[1]], lower[[1]], upper[[1]])

  # Combine item-level estimates
  items <- mapply(LEU, est$items, lower$items, upper$items,
                  SIMPLIFY = FALSE)

  list(classes = classes, items = items)
}
