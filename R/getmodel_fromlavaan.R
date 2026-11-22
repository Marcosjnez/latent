# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 07/04/2026

getmodel_fromlavaan <- function(obj) {

  matrices <- lavaan::lavTech(obj, "free", list.by.group = TRUE,
                              drop.list.single.group = FALSE)
  table <- lavaan::lavMatrixRepresentation(lavaan::partable(obj), representation = "LISREL")

  for(j in seq_along(matrices)) {
    group <- table[table$group == j, , drop = FALSE]
    for(name in names(matrices[[j]])) {
      entries <- group[group$mat == name, , drop = FALSE]
      if(!nrow(entries)) next
      value <- ifelse(entries$label == "", entries$plabel, entries$label)
      fixed <- entries$free == 0 & !is.na(entries$ustart)
      value[fixed] <- entries$ustart[fixed]
      M <- matrices[[j]][[name]]
      M[cbind(entries$row, entries$col)] <- value
      if(name %in% c("psi", "theta")) M[cbind(entries$col, entries$row)] <- value
      matrices[[j]][[name]] <- M
    }
  }

  if(length(matrices) == 1L) matrices[[1L]] else unname(matrices)
}
