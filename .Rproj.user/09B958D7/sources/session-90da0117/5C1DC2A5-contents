#' @title
#' Get the default model for Latent Class Analysis.
#' @description
#'
#' Get the default model for Latent Class Analysis.
#'
#' @usage
#'
#' getmodel(data, model = rep("multinomial", ncol(data)), nclasses = 2L)
#'
#' @param data data.frame or matrix of response.
#' @param model Character vector with the model for each item.
#' @param nclasses Number of latent classes.
#'
#' @details \code{getmodel} generates the model for the probability of belonging
#' to the classes and the conditional response probabilities. These models may
#' be modified by the user to set equality constraints or to fix parameters.
#'
#' @return List with the following objects:
#' \item{classes}{The model for the probabilities of the classes.}
#' \item{conditionals}{The model for the conditional response probabilities.}
#'
#' @references
#'
#' None yet.
#'
#' @export
getmodel <- function(data, model = rep("multinomial", ncol(data)), nclasses = 2L,
                     free = FALSE) {

  # This function creates the default LCA model

  # Create the model for the classes:
  classes <- paste("class", 1:nclasses, sep = "")
  if(!free) classes[1] <- "-1" # Set a reference class
  names(classes) <- paste("Class", 1:nclasses, sep = "")

  # Create the model of the conditional parameters for each item:
  nitems <- ncol(data)
  conditionals <- vector("list", length = nitems)
  names(conditionals) <- paste("Item", 1:nitems)

  for(j in 1:nitems) { # Loop across the items

    if(model[j] == "multinomial") { # For categorical items

      # The model of the conditional parameters has as many rows as response
      # categories in the item and as many columns as the number of classes
      ncategories <- max(data[, j])
      catg_v <- rep(1:ncategories, times = nclasses)
      classes_v <- rep(paste("class", 1:nclasses, sep = ""), each = ncategories)
      conditionals[[j]] <- matrix(paste("y", j, catg_v, "|", classes_v, sep = ""),
                                  nrow = ncategories, ncol = nclasses)
      rownames(conditionals[[j]]) <- paste("Category", 1:ncategories, sep = "")
      colnames(conditionals[[j]]) <- paste("Class", 1:nclasses, sep = "")
      if(!free) conditionals[[j]][1, ] <- "-1" # Set a reference category in each item

    } else if(model[j] == "gaussian") { # For continuous items

      # The model of the conditional parameters has 2 rows (one for the mean of
      # the item and another for the standard deviation) and as many columns as
      # the number of classes
      classes_v <- rep(paste("class", 1:nclasses, sep = ""), each = 2)
      conditionals[[j]] <- matrix(paste(c("mean", "std"), j, "|", classes_v, sep = ""),
                                  nrow = 2, ncol = nclasses)
      rownames(conditionals[[j]]) <- c("Means", "Stds")
      colnames(conditionals[[j]]) <- paste("Class", 1:nclasses, sep = "")

    }
  }

  result <- list(classes = classes, conditionals = conditionals)

  return(result)

}
