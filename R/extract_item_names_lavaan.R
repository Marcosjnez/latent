# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 13/09/2026

extract_item_names_lavaan <- function(x, ngroups = NULL) {

  if(inherits(x, "lavaan")) {

    partable <- lavaan::parTable(x)

    if(is.null(ngroups)) {
      groups <- unique(partable$group[partable$group > 0L])
      ngroups <- length(groups)
      if(ngroups == 0L) ngroups <- 1L
    }

    group_labels <- tryCatch(lavaan::lavInspect(x, "group.label"),
                             error = function(e) NULL)

  } else if(is.character(x) && length(x) == 1L) {

    if(is.null(ngroups)) ngroups <- 1L
    partable <- lavaan::lavaanify(x, ngroups = ngroups)
    group_labels <- NULL

  } else {

    stop("x must be either a lavaan model syntax string or a fitted lavaan object.")

  }

  if(!is.numeric(ngroups) || length(ngroups) != 1L ||
     !is.finite(ngroups) || ngroups < 1L ||
     ngroups != as.integer(ngroups)) {
    stop("ngroups must be NULL or a positive integer.")
  }

  ngroups <- as.integer(ngroups)
  item_names <- vector("list", ngroups)

  for(i in seq_len(ngroups)) {

    group_table <- partable[partable$group == i, , drop = FALSE]
    latent_names <- unique(group_table$lhs[group_table$op == "=~"])
    indicators <- unique(group_table$rhs[group_table$op == "=~"])

    # Latent indicators belong to the beta matrix in lavaan's LISREL
    # representation (for example, higher-order factor loadings), not to the
    # observed-variable covariance matrix used by lcfa.
    item_names[[i]] <- indicators[!(indicators %in% latent_names)]

  }

  if(!is.null(group_labels) && length(group_labels) == ngroups) {
    names(item_names) <- group_labels
  } else {
    names(item_names) <- paste0("group", seq_len(ngroups))
  }

  #### Result ####

  return(item_names)

}
