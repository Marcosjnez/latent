# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 07/06/2025

cfast_control <- function(control) {

  # Auxiliary function for cfast.R

  # Control input

  if(is.null(control$opt)) {
    control$opt <- "lbfgs"
  }

  if(is.null(control$step_maxit)) {
    control$step_maxit <- 30L
  }

  if(is.null(control$c1)) {
    control$c1 <- 0.5
  }

  if(is.null(control$c2)) {
    control$c2 <- 0.5
  }

  if(is.null(control$step_eps)) {
    control$step_eps <- 1e-09
  }

  if(is.null(control$df_eps)) {
    control$df_eps <- 1e-09
  }

  if(is.null(control$M)) {
    control$M <- 100L
  }

  if(is.null(control$eps)) {
    control$eps <- 1e-06
  }

  if(is.null(control$ss_fac)) {
    control$ss_fac <- 2
  }

  if(is.null(control$maxit)) {
    control$maxit <- 1000L
  }

  if(is.null(control$rstarts)) {
    control$rstarts <- 1L
  }

  if(is.null(control$cores)) {
    control$cores <- 1L
  }

  if(is.null(control$tcg_maxit)) {
    control$tcg_maxit <- 10L
  }

  return(control)

}
