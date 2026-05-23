# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 04/05/2026
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
predict.llca <- function(model, new = NULL) {

  if(is.null(new) || ncol(new) == 0) {

    X <- model@dataList$X
    rownames(X) <- rownames(model@dataList$data)

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

  Z0 <- X %*% model@transformed_pars$beta
  linpreds <- apply(Z0, MARGIN = 1, FUN = latent::soft, a = 1.00)
  if(is.vector(linpreds)) {
    P0 <- matrix(linpreds, ncol = 1)
  } else {
    P0 <- t(linpreds)
  }
  colnames(P0) <- colnames(model@transformed_pars$class)

  remove <- -match("(Intercept)", colnames(X))
  return(data.frame(X[, remove], P0))

}

#' @method predict llcalist
#' @export
predict.llcalist <- function(model) {

  nmodels <- length(model)
  out <- vector("list", length = nmodels)
  for(i in 1:nmodels) {

    out[[i]] <- predict.llca(model[[i]])
    names(out)[i] <- paste("nclasses=",
                           ncol(model[[i]]@modelInfo$trans$class),
                           sep = "")

  }

  class(out) <- "predict.llcalist"

  return(out)

}
