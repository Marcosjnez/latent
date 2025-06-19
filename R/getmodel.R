# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 07/06/2025
#'
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
#' @param item Character vector with the model for each item.
#' @param nclasses Number of latent classes.
#' @param model List of parameter labels. See 'details' for more information.
#' @param constraints Should the model be checked for identification? Defaults to TRUE.
#'
#' @details \code{getmodel} generates the model for the probability of belonging
#' to the classes and the conditional response probabilities. These models may
#' be modified by the user to set equality constraints or to fix parameters.
#'
#' @return List with the following objects:
#' \item{none}{.}
#' \item{none}{.}
#'
#' @references
#'
#' None yet.
#'
#' @export
getmodel <- function(data, item = rep("multinomial", ncol(data)), nclasses = 2L,
                     model = NULL, constraints = TRUE) {

  # This function creates and checks the LCA model in logarithm and probability
  # scales

  # Initialize the logarithm model for the items:
  nitems <- ncol(data)
  conditionals <- vector("list", length = nitems)
  names(conditionals) <- paste("Item", 1:nitems)

  # If the user did not supply a custom model, create a default model:
  if(is.null(model)) {

    # Create the model for the classes:
    classes <- paste("class", 1:nclasses, sep = "")
    names(classes) <- paste("Class", 1:nclasses, sep = "")

    # Create the model for the items:
    for(j in 1:nitems) { # Loop across all the items

      if(item[j] == "multinomial") { # For categorical items

        # The model of the items has as many rows as there are response
        # categories in an item and as many columns as the number of classes
        ncategories <- max(data[, j], na.rm = TRUE)
        catg_v <- rep(1:ncategories, times = nclasses)
        classes_v <- rep(classes, each = ncategories)
        conditionals[[j]] <- matrix(paste("y", j, catg_v, "|", classes_v, sep = ""),
                                    nrow = ncategories, ncol = nclasses)
        rownames(conditionals[[j]]) <- paste("Category", 1:ncategories, sep = "")
        colnames(conditionals[[j]]) <- paste("Class", 1:nclasses, sep = "")

      } else if(item[j] == "gaussian") { # For continuous items

        # The model of the items has as 2 rows (one for the mean of the item and
        # another for the standard deviation) and as many columns as the number
        # of classes
        classes_v <- rep(classes, each = 2)
        conditionals[[j]] <- matrix(paste(c("mean", "log_sd"), j, "|", classes_v, sep = ""),
                                    nrow = 2, ncol = nclasses)
        rownames(conditionals[[j]]) <- c("Means", "log_Sds")
        colnames(conditionals[[j]]) <- paste("Class", 1:nclasses, sep = "")

      }
    }

  } else {

    # If the user supplied a custom model, check that there is a model for the
    # classes and another for the items

    if(is.null(model$classes)) stop("No model was found for the classes")
    if(is.null(model$conditionals)) stop("No model was found for the items")
    classes <- model$classes
    conditionals <- model$conditionals

  }

  # Set constraints to get an identifiable model:
  if(constraints) {

    # Set a reference category for the classes:
    if(sum(!is.na(suppressWarnings(as.numeric(classes)))) < 1) {
      classes[1] <- "0" # Set reference class at position 1
    }

    for(j in 1:nitems) { # Loop across the items

      # For categorical items, set a reference category:
      if(item[j] == "multinomial") {

        if(sum(!is.na(suppressWarnings(as.numeric(conditionals[[j]][1, ])))) < 1) {
          conditionals[[j]][1, ] <- "0" # Set a reference category in each item
        }

      }
    }

    model <- list(classes = classes, conditionals = conditionals)

  }

  # Build the logarithm model:
  log_model <- list(classes = model$classes, conditionals = model$conditionals)
  # Generate the probability model from the logarithm model:
  prob_model <- log2prob(log_model, data, item, nitems, nclasses)

  # Return the model in the logarithm and probability scales:
  result <- list(prob_model = prob_model, log_model = log_model)

  return(result)

}

log2prob <- function (log_model, data, item, nitems, nclasses) {
  conditionals <- vector("list", length = nitems)
  names(conditionals) <- paste("Item", 1:nitems)
  classes <- paste("P(class", 1:nclasses, ")", sep = "")
  names(classes) <- paste("Class", 1:nclasses, sep = "")
  for (j in 1:nitems) {
    if (item[j] == "multinomial") {
      uniques <- unique(data[, j])
      ncategories <- sum(!is.na(uniques))
      classes_v <- rep(paste("class", 1:nclasses, sep = ""),
                       each = ncategories)
      catg_v <- rep(1:ncategories, times = nclasses)
      conditionals[[j]] <- matrix(paste("P(y", j, catg_v,
                                        "|", classes_v, ")", sep = ""), nrow = ncategories,
                                  ncol = nclasses)
      rownames(conditionals[[j]]) <- paste("P(Category",
                                           1:ncategories, ")", sep = "")
      colnames(conditionals[[j]]) <- paste("Class", 1:nclasses,
                                           sep = "")
    }
    else if (item[j] == "gaussian") {
      classes_v <- rep(paste("class", 1:nclasses, sep = ""),
                       each = 2)
      conditionals[[j]] <- matrix(paste(c("mean", "sd"),
                                        j, "|", classes_v, sep = ""), nrow = 2, ncol = nclasses)
      rownames(conditionals[[j]]) <- c("Means", "Sds")
      colnames(conditionals[[j]]) <- paste("Class", 1:nclasses,
                                           sep = "")
    }
  }
  prob_model <- list(classes = classes, conditionals = conditionals)
  slots <- slots2 <- vector("list")
  k <- 1
  slots[[k]] <- log_model$classes
  slots2[[k]] <- prob_model$classes
  for (i in 1:length(log_model$conditionals)) {
    for (j in 1:ncol(log_model$conditionals[[i]])) {
      k <- k + 1
      slots[[k]] <- log_model$conditionals[[i]][, j]
      slots2[[k]] <- prob_model$conditionals[[i]][, j]
    }
  }
  nslots <- length(slots)
  for (i in 1:(nslots - 1L)) {
    ni <- length(slots[[i]])
    for (j in 2:(nslots)) {
      nj <- length(slots[[j]])
      same_length <- ni == nj
      if (same_length) {
        if (all(slots[[j]] %in% slots[[i]])) {
          mi <- match(slots[[i]], slots[[j]])
          slots[[j]] <- slots[[i]][mi]
          slots2[[j]] <- slots2[[i]][mi]
        }
      }
    }
  }
  prob_model <- fill_list_with_vector(prob_model, unlist(slots2))
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
