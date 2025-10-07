# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 07/10/2025
#'
#' @title
#' Fitted values for Latent Class Analysis.
#' @description
#'
#' Get fitted class membership probabilities for latent class models.
#'
#' @usage
#'
#' fitted(model)
#'
#' @param model Fitted llca object.
#'
#' @details None.
#'
#' @return Stuff:
#' \item{dof}{Degrees of freedom}
#'
#' @references
#'
#' None yet.
#'
#' @export
fitted.llca <- function(model) {

  X <- model@Optim$data_list$X
  P0 <- predict.llca(model = model, new = X[, -1])

  return(P0)

}
