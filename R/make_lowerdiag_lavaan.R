make_lowerdiag_lavaan <- function(data, nfactors,
                                  factor_prefix = "F",
                                  label_loadings = TRUE) {
  # ---- get variable/item names ----
  if (is.null(dim(data))) stop("`data` must be a matrix/data.frame-like object.")
  is_square <- nrow(data) == ncol(data)

  vars <- colnames(data)
  if (is.null(vars) && is_square) vars <- rownames(data)
  if (is.null(vars)) vars <- paste0("V", seq_len(if (is_square) ncol(data) else ncol(data)))

  p <- length(vars)
  q <- as.integer(nfactors)
  if (q < 1) stop("`nfactors` must be >= 1.")
  if (q > p) stop("`nfactors` cannot exceed the number of items (columns).")

  # ---- build lavaan lines: Fj =~ (items j..p) ----
  lines <- character(q)

  for (j in seq_len(q)) {
    items <- vars[j:p]

    terms <- character(length(items))
    for (k in seq_along(items)) {
      i <- j + k - 1  # item index in 1..p
      lab <- if (label_loadings) paste0("l", i, "_", j, "*") else ""

      if (k == 1) {
        # free the "first" loading for this factor (avoid the default 1-fixing)
        terms[k] <- paste0("NA*", lab, items[k])
      } else {
        terms[k] <- paste0(lab, items[k])
      }
    }

    lines[j] <- paste0(factor_prefix, j, " =~ ", paste(terms, collapse = " + "))
  }

  paste(lines, collapse = "\n")
}
