# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 07/10/2025
#'
#' @title
#' Predicted values for Latent Class Analysis.
#' @description
#'
#' Get predicted class membership probabilities for latent class models.
#'
#' @usage
#'
#' predict(model, new = NULL)
#'
#' @param model Fitted llca object.
#' @param new Matrix or data.frame of new predictors.
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
predict.llca <- function(model, new = NULL, digits = NULL) {

  if(is.null(new) || ncol(new) == 0) {

    X <- model@data_list$X

  } else {

    X <- new
    # Check that X is either a data.frame or a matrix:
    if(!is.data.frame(X) & !is.matrix(X)) {
      stop("data must be a matrix or data.frame")
    }

    # Transform characters into factors:
    X_df <- as.data.frame(X)
    X_df[] <- lapply(X_df, function(x) if (is.character(x)) factor(x) else x)

    # Create the design matrix:
    X <- model.matrix(~ . + 1, X_df)
    # Standardize the variables:
    # X[, -1] <- apply(X[, -1], MARGIN = 2, FUN = \(x) x-mean(x))

    # Put an underscore between the variable names and their level names:
    for (v in names(X_df)[sapply(X_df, is.factor)]) {
      i <- which(startsWith(colnames(X), v))
      colnames(X)[i] <- paste0(v, "_", make.names(levels(X_df[[v]]))[seq_along(i)])
    }

  }

  Z0 <- X %*% model@parameters$beta
  linpreds <- apply(Z0, MARGIN = 1, FUN = latent::soft, a = 1.00)
  if(is.vector(linpreds)) {
    P0 <- matrix(linpreds, ncol = 1)
  } else {
    P0 <- t(linpreds)
  }
  rownames(P0) <- rownames(X)
  colnames(P0) <- colnames(model@transformed_pars$class)

  if(!is.null(digits)) {
    P0 <- round(P0, digits = digits)
  }

  return(data.frame(X[,-1], P0) )

}

#' @method predict llcalist
#' @export
predict.llcalist <- function(model, digits = NULL) {

  nmodels <- length(model)
  out <- vector("list", length = nmodels)
  for(i in 1:nmodels) {

    out[[i]] <- predict.llca(model[[i]], digits = digits)
    names(out)[i] <- paste("nclasses = ", model[[i]]@modelInfo$nclasses,
                           sep = "")

  }

  class(out) <- "predict.llcalist"

  return(out)

}
