# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 18/08/2026
#'
#' Information Variance-Covariance Matrix for Latent Models
#'
#' @param fit A fitted object inheriting from class \code{"latent"}.
#'
#' @return A list containing the Hessian, variance-covariance matrix, and
#'   standard errors of the freely estimated parameters.
#'
#' @method information latent
#' @export
information.latent <- function(fit) {

  if(inherits(fit, "multistep")) {

    result <- information.multistep(fit)

    #### Result ####

    return(result)

  }

  labels <- fit@modelInfo$parameters_labels

  if(length(labels) == 0L) {

    empty <- matrix(numeric(0L), nrow = 0L, ncol = 0L,
                    dimnames = list(labels, labels))

    #### Result ####

    return(list(H = empty, VCOV = empty, se = numeric(0L)))

  }

  stored_VCOV <- tryCatch(fit@Optim$SE$VCOV,
                          error = function(e) NULL)

  # Direct sample-statistic estimators store their already scaled sampling
  # covariance matrix in Optim$SE$VCOV. Reuse it when it has exactly the free
  # parameter dimensions.
  if(!is.null(stored_VCOV) &&
     is.matrix(stored_VCOV) &&
     nrow(stored_VCOV) == length(labels) &&
     ncol(stored_VCOV) == length(labels)) {

    VCOV <- as.matrix(stored_VCOV)

    if(!isSymmetric(VCOV) || any(!is.finite(VCOV))) {
      stop("The stored VCOV matrix is not finite and symmetric.")
    }

    if(!is.null(rownames(VCOV)) && !is.null(colnames(VCOV)) &&
       all(labels %in% rownames(VCOV)) &&
       all(labels %in% colnames(VCOV))) {
      VCOV <- VCOV[labels, labels, drop = FALSE]
    }

    VCOV <- (VCOV+t(VCOV))/2
    rownames(VCOV) <- colnames(VCOV) <- labels

    H <- tryCatch(fit@Optim$SE$H,
                  error = function(e) NULL)

  } else {

    H <- hessian(fit)
    H <- (H+t(H))/2

    eigenvalues <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    tolerance <- sqrt(.Machine$double.eps)*max(1, max(abs(eigenvalues)))

    if(any(eigenvalues <= tolerance)) {
      warning("The Hessian is not positive definite; an approximate inverse was used.")
    }

    VCOV <- approx_Hinv(H)

    if(any(!is.finite(VCOV))) {
      stop("The Hessian could not be inverted.")
    }

    VCOV <- (VCOV+t(VCOV))/2
    rownames(VCOV) <- colnames(VCOV) <- labels

  }

  se <- sqrt(diag(VCOV))
  names(se) <- labels

  #### Result ####

  result <- list(H = H,
                 VCOV = VCOV,
                 se = se)

  return(result)

}
