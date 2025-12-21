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
fitted.llca <- function(model, digits = NULL) {

  X <- model@data_list$X
  P0 <- predict.llca(model = model, new = X[, -1], digits = digits)

  return(P0)

}

#' @method fitted llcalist
#' @export
fitted.llcalist <- function(model, digits = NULL) {

  nmodels <- length(model)
  out <- vector("list", length = nmodels)
  for(i in 1:nmodels) {

    out[[i]] <- fitted.llca(model[[i]], digits = digits)
    names(out)[i] <- paste("nclasses = ", model[[i]]@modelInfo$nclasses,
                           sep = "")

  }

  class(out) <- "fitted.llcalist"

  return(out)

}
