extract_item_names <- function(x, ngroups = NULL) {
  # x can be:
  # - a lavaan model syntax string
  # - a fitted lavaan object
  #
  # ngroups:
  # - NULL or 1: default behavior
  # - >1: replicate the same item names for each group slot

  if (inherits(x, "lavaan")) {
    pt <- parTable(x)
    items <- unique(pt[pt$op == "=~", c("group", "rhs")])

    if (is.null(ngroups) || ngroups == 1L) {
      item_names <- split(items$rhs, items$group)

      group_labels <- tryCatch(lavInspect(x, "group.label"), error = function(e) NULL)
      if (!is.null(group_labels) && length(group_labels) == length(item_names)) {
        names(item_names) <- group_labels
      } else {
        names(item_names) <- paste0("group", seq_along(item_names))
      }

    } else if (ngroups > 1L) {
      base_items <- unique(items$rhs)
      item_names <- replicate(ngroups, base_items, simplify = FALSE)
      names(item_names) <- paste0("group", seq_len(ngroups))

    } else {
      stop("ngroups must be NULL or a positive integer.")
    }

  } else if (is.character(x) && length(x) == 1L) {
    pt <- lavaan::lavaanify(x)
    base_items <- unique(pt[pt$op == "=~", "rhs"])

    if (is.null(ngroups) || ngroups == 1L) {
      item_names <- list(group1 = base_items)
    } else if (ngroups > 1L) {
      item_names <- replicate(ngroups, base_items, simplify = FALSE)
      names(item_names) <- paste0("group", seq_len(ngroups))
    } else {
      stop("ngroups must be NULL or a positive integer.")
    }

  } else {
    stop("x must be either a lavaan model syntax string or a fitted lavaan object.")
  }

  item_names
}
