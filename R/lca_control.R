# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 02/09/2025

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

  if(is.null(control$cores)) {
    control$cores <- parallel::detectCores()-1L
  }

  if(is.null(control$tcg_maxit)) {
    control$tcg_maxit <- 10L
  }

  if(is.null(control$opt)) {
    control$opt <- "em-lbfgs"
  }

  if(control$opt == "em-lbfgs") {

    if(is.null(control$em_rstarts)) {
      control$em_rstarts <- 50
    }

    if(is.null(control$maxit_em)) {
      control$maxit_em <- 250L
    }

    if(is.null(control$rstarts)) {
      control$rstarts <- 5L
    }

    # Save number of rstarts to store and compute in lbfgs
    control$pick <- control$rstarts
    # In "em-lbfgs", use at least as many em_rstarts as rstarts in lbgfs
    control$rstarts <- max(control$em_rstarts, control$rstarts)

  }

  if(is.null(control$rstarts)) {
    control$rstarts <- 16L
  }

  return(control)

}
