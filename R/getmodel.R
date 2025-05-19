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
getmodel <- function(data, item_model = rep("multinomial", ncol(data)), nclasses = 2L,
                     model = NULL, constraints = TRUE) {

  # This function creates the default LCA model

  # Create the model of the conditional parameters for each item:
  nitems <- ncol(data)
  conditionals <- vector("list", length = nitems)
  names(conditionals) <- paste("Item", 1:nitems)

  if(is.null(model)) {

    # Create the probability model by default
    # Create the model for the classes:
    classes <- paste("class", 1:nclasses, sep = "")
    names(classes) <- paste("Class", 1:nclasses, sep = "")

    for(j in 1:nitems) { # Loop across the items

      if(item_model[j] == "multinomial") { # For categorical items

        # The model of the conditional parameters has as many rows as response
        # categories in the item and as many columns as the number of classes
        ncategories <- max(data[, j])
        catg_v <- rep(1:ncategories, times = nclasses)
        classes_v <- rep(paste("class", 1:nclasses, sep = ""), each = ncategories)
        conditionals[[j]] <- matrix(paste("y", j, catg_v, "|", classes_v, sep = ""),
                                    nrow = ncategories, ncol = nclasses)
        rownames(conditionals[[j]]) <- paste("Category", 1:ncategories, sep = "")
        colnames(conditionals[[j]]) <- paste("Class", 1:nclasses, sep = "")

      } else if(item_model[j] == "gaussian") { # For continuous items

        # The model of the conditional parameters has 2 rows (one for the mean of
        # the item and another for the standard deviation) and as many columns as
        # the number of classes
        classes_v <- rep(paste("class", 1:nclasses, sep = ""), each = 2)
        conditionals[[j]] <- matrix(paste(c("mean", "log_sd"), j, "|", classes_v, sep = ""),
                                    nrow = 2, ncol = nclasses)
        rownames(conditionals[[j]]) <- c("Means", "log_Sds")
        colnames(conditionals[[j]]) <- paste("Class", 1:nclasses, sep = "")

      }
    }

  } else {

    classes <- model$classes
    conditionals <- model$conditionals

  }

  if(constraints) {

    # Put all the constrainst that are necessary to identify the model
    if(sum(!is.na(suppressWarnings(as.numeric(classes)))) < 1) {
      classes[1] <- "-1" # Set reference class at position 1
    }

    for(j in 1:nitems) { # Loop across the items

      if(item_model[j] == "multinomial") { # For categorical items

        if(sum(!is.na(suppressWarnings(as.numeric(conditionals[[j]][1, ])))) < 1) {
          conditionals[[j]][1, ] <- "-1" # Set a reference category in each item
        }

      }
    }

    model <- list(classes = classes, conditionals = conditionals)

  }

  log_model <- list(classes = model$classes, conditionals = model$conditionals)
  prob_model <- log2prob(log_model, nitems, nclasses)

  result <- list(prob_model = prob_model, log_model = log_model)

  return(result)

}

log2prob <- function(log_model, nitems, nclasses) {

  conditionals <- vector("list", length = nitems)
  names(conditionals) <- paste("Item", 1:nitems)

  # Create the probability model by default
  # Create the model for the classes:
  classes <- paste("P(class", 1:nclasses, ")", sep = "")
  names(classes) <- paste("Class", 1:nclasses, sep = "")

  for(j in 1:nitems) { # Loop across the items

    if(item_model[j] == "multinomial") { # For categorical items

      # The model of the conditional parameters has as many rows as response
      # categories in the item and as many columns as the number of classes
      ncategories <- length(unique(data[, j]))
      classes_v <- rep(paste("class", 1:nclasses, sep = ""), each = ncategories)
      catg_v <- rep(1:ncategories, times = nclasses)
      conditionals[[j]] <- matrix(paste("P(y", j, catg_v, "|", classes_v, ")", sep = ""),
                                  nrow = ncategories, ncol = nclasses)
      rownames(conditionals[[j]]) <- paste("P(Category", 1:ncategories, ")", sep = "")
      colnames(conditionals[[j]]) <- paste("Class", 1:nclasses, sep = "")

    } else if(item_model[j] == "gaussian") { # For continuous items

      # The model of the conditional parameters has 2 rows (one for the mean of
      # the item and another for the standard deviation) and as many columns as
      # the number of classes
      classes_v <- rep(paste("class", 1:nclasses, sep = ""), each = 2)
      conditionals[[j]] <- matrix(paste(c("mean", "sd"), j, "|", classes_v, sep = ""),
                                  nrow = 2, ncol = nclasses)
      rownames(conditionals[[j]]) <- c("Means", "Sds")
      colnames(conditionals[[j]]) <- paste("Class", 1:nclasses, sep = "")

    }
  }

  prob_model <- list(classes = classes, conditionals = conditionals)

  # slots <- slots2 <- vector("list")
  # k <- 1
  # slots[[k]] <- log_model$classes
  # slots2[[k]] <- prob_model$classes
  # for(i in 1:length(log_model$conditionals)) {
  #   for(j in 1:ncol(log_model$conditionals[[i]])) {
  #     k <- k+1
  #     slots[[k]] <- log_model$conditionals[[i]][, j]
  #     slots2[[k]] <- prob_model$conditionals[[i]][, j]
  #   }
  # }
  #
  # nslots <- length(slots)
  # for(i in 1:(nslots-1L)) {
  #
  #   ni <- length(slots[[i]])
  #
  #   for(j in 2:(nslots)) {
  #
  #     nj <- length(slots[[j]])
  #     same_length <- ni == nj
  #
  #     if(same_length) {
  #       if(all(slots[[j]] %in% slots[[i]])) {
  #         mi <- match(slots[[i]], slots[[j]])
  #         slots[[j]] <- slots[[i]][mi]
  #         slots2[[j]] <- slots2[[i]][mi]
  #       }
  #     }
  #
  #   }
  # }
  #
  # prob_model <- fill_list_with_vector(prob_model, unlist(slots2))

  return(prob_model)

}

fill_list_with_vector <- function(lst, values) {
  i <- 1

  assign_recursive <- function(x) {
    if (is.list(x)) {
      lapply(x, assign_recursive)
    } else if (is.matrix(x)) {
      dims <- dim(x)
      n <- prod(dims)
      x[] <- values[i:(i + n - 1)]
      i <<- i + n
      x
    } else if (is.atomic(x)) {
      n <- length(x)
      x[] <- values[i:(i + n - 1)]
      i <<- i + n
      x
    } else {
      stop("Unsupported type")
    }
  }

  assign_recursive(lst)
}
