# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 22/05/2025

# Recursive function to replace near-zero elements
replace_near_zero <- function(x, tol = .Machine$double.eps) {
  if (is.list(x)) {
    lapply(x, replace_near_zero, tol = tol)
  } else if (is.numeric(x)) {
    x[abs(x) < tol] <- tol
    x
  } else {
    x
  }
}
