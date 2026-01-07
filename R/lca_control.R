# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 03/09/2025

lca_control <- function(control) {

  # Auxiliary function for lca.R

  # Control input

  # control$opt <- "lbfgs" # FORCE LBFGS

  if(isFALSE(control$penalties)) {

    control$reg <- FALSE

  } else if(isTRUE(control$penalties)) {

    control$reg <- TRUE

    control$penalties <- list(
      class = list(alpha = 1),
      prob  = list(alpha = 1),
      sd    = list(alpha = 1)
    )

  } else if(is.list(control$penalties)) {

    control$reg <- TRUE

  } else {

    stop("penalties should be TRUE, FALSE, or a list")

  }

  if(is.null(control$step_maxit)) {
    control$step_maxit <- 30L
  } else if(control$step_maxit < 1L) {
    stop("step_maxit must be an integer greater than 0")
  }

  if(is.null(control$c1)) {
    control$c1 <- 0.5
  } else if(control$c1 < 0L) {
    stop("c1 must be a positive number, preferable lower than c2")
  }

  if(is.null(control$c2)) {
    control$c2 <- 0.5
  } else if(control$c2 < 0L) {
    stop("c2 must be a positive number, preferable larger than c1")
  }

  if(is.null(control$step_eps)) {
    control$step_eps <- 1e-09
  } else if(control$step_eps < 0) {
    stop("step_eps must be a positive number, preferable close to 0")
  }

  if(is.null(control$df_eps)) {
    control$df_eps <- 1e-09
  } else if(control$df_eps < 0) {
    stop("df_eps must be a positive number, preferable close to 0")
  }

  if(is.null(control$eps)) {
    control$eps <- 1e-05
  } else if(control$eps < 0) {
    stop("eps must be a positive number, preferable close to 0")
  }

  if(is.null(control$M)) {
    control$M <- 100L
  } else if(control$M < 0L) {
    stop("M must be a positive integer")
  }

  if(is.null(control$ss_fac)) {
    control$ss_fac <- 2
  } else if(control$ss_fac <= 1) {
    stop("ss_fac must be a positive integer larger than 1")
  }

  if(is.null(control$maxit)) {
    control$maxit <- 1000L
  } else if(control$maxit < 0L) {
    stop("maxit must be a positive integer")
  }

  if(is.null(control$cores)) {
    control$cores <- parallel::detectCores()-1L
  } else if(control$cores < 0L) {
    stop("cores must be a positive integer")
  }

  if(is.null(control$tcg_maxit)) {
    control$tcg_maxit <- 10L
  } else if(control$tcg_maxit < 0L) {
    stop("tcg_maxit must be a positive integer")
  }

  if(is.null(control$opt)) {
    control$opt <- "lbfgs"
  }

  if(is.null(control$rstarts)) {
    control$rstarts <- 16L
  } else if(control$rstarts < 0L) {
    stop("rstarts must be a positive integer")
  }

  return(control)

}
