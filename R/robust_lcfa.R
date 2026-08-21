# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 21/08/2026
#'
#' Robust Variance-Covariance Matrix for CFA Models
#'
#' @param fit A fitted object inheriting from class \code{"lcfa"}.
#'
#' @return This method currently stops because casewise/patternwise robust
#'   likelihood scores have not yet been implemented for direct CFA likelihood
#'   models.
#'
#' @method robust lcfa
#' @export
robust.lcfa <- function(fit) {

  stop("Robust standard errors are not yet implemented for ordinary lcfa ",
       "likelihood models. Use information standard errors, or a multistep ",
       "analysis whose parent statistic supplies a robust VCOV matrix.")

  #### Result ####

  return(NULL)

}
