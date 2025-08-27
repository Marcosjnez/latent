# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 24/08/2025

lca_control <- function(control) {

  # Auxiliary function for lca.R

  # Control input

  if(isFALSE(control$penalties)) {

    control$reg <- FALSE

  } else {

    control$reg <- TRUE

    if(!is.list(control$penalties)) {

      control$penalties <- list(
        class = list(alpha = 1),
        prob  = list(alpha = 1),
        sd    = list(alpha = 1)
      )

    }
  }

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
    control$eps <- 1e-05
  }

  if(is.null(control$ss_fac)) {
    control$ss_fac <- 2
  }

  if(is.null(control$maxit)) {
    control$maxit <- 1000L
  }

  if(is.null(control$rstarts)) {
    control$rstarts <- 16L
  }

  if(is.null(control$cores)) {
    control$cores <- 1L
  }

  if(is.null(control$tcg_maxit)) {
    control$tcg_maxit <- 10L
  }

  return(control)

}
