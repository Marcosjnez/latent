# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 27/06/2026
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

  covariates <- model@dataList$design
  P0 <- predict.llca(model = model, new = covariates[, -1])

  return(P0)

}

#' @method fitted llcalist
#' @export
fitted.llcalist <- function(model) {

  nmodels <- length(model)
  out <- vector("list", length = nmodels)
  for(i in 1:nmodels) {

    out[[i]] <- fitted.llca(model[[i]])
    names(out)[i] <- paste("nclasses=",
                           ncol(model[[i]]@modelInfo$param$beta),
                           sep = "")

  }

  class(out) <- "fitted.llcalist"

  return(out)

}
