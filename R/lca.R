# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 10/08/2025
#'
#' @title
#' Latent Class Analysis.
#' @description
#'
#' Estimate latent class models with gaussian and multinomial item models.
#'
#' @usage
#'
#' lca(data, item = rep("gaussian", ncol(data)), nclasses = 2L, model = NULL,
#'     do.fit = TRUE, control = list(opt = "lbfgs", rstarts = 30L, cores = 1L))
#'
#' @param data data frame or matrix.
#' @param nclasses Number of latent classes.
#' @param item Character vector with the model for each item (i.e., "gaussian" or "multinomial"). Defaults to "gaussian" for all the items.
#' @param model List of parameter labels. See 'details' for more information.
#' @param do.fit TRUE to fit the model and FALSE to return only the model setup. Defaults to TRUE.
#' @param control List of control parameters for the optimization algorithm. See 'details' for more information.
#'
#' @details \code{lca} estimates models with categorical and continuous data.
#'
#' @return List with the following objects:
#' \item{parameters}{The model for the logarithm probabilities of the classes.}
#' \item{f}{Logarithm likelihood at the maximum.}
#'
#' @references
#'
#' None yet.
#'
#' @export
lca <- function(data, nclasses = 2L, item = rep("gaussian", ncol(data)),
                model = NULL, control = NULL, do.fit = TRUE) {

  #### Initial input checks ####

  # Check that data is either a data.frame or a matrix:
  if(!is.data.frame(data) & !is.matrix(data)) {
    stop("data must be a matrix or data.frame")
  }

  # Number of items:
  nitems <- ncol(data)

  # Check that item is a character vector with a string for each column of data:
  if(!is.character(item) || length(item) != nitems) {

    stop("item must be a character vector with as many elements as columns in data")

  }

  #### Process the data ####

  # Number of subjects:
  nobs <- nrow(data)
  # Convert data to a data.table object:
  dt <- data.table::as.data.table(data)
  ## Collect some information from the data ##
  counts_dt <- dt[, .(index = .I[1], count = .N), by = names(dt)]
  # Data matrix with the unique response patterns:
  patterns <- as.matrix(counts_dt[, names(dt), with = FALSE])
  # Number of unique response patterns:
  npatterns <- nrow(counts_dt)
  # Counts of each response pattern:
  weights <- counts_dt$count
  # Indices to map the original data to the matrix of unique patterns:
  full2short <- counts_dt$index
  # Indices to map the matrix of unique patterns to the original data:
  short2full <- match(do.call(paste, dt),
                      do.call(paste, counts_dt[, -c("index", "count"), with = FALSE]))

  # Put in a list the objects generated form the data:
  data_list <- vector("list")
  data_list$dt <- dt
  data_list$patterns <- patterns
  data_list$npatterns <- npatterns
  data_list$nitems <- nitems
  data_list$weights <- weights
  data_list$full2short <- full2short
  data_list$short2full <- short2full

  #### Initialize objects to store all the models ####

  model0 <- model
  NCLASSES <- nclasses
  nmodels <- length(NCLASSES)

  llca_list <- list()

  model0 <- model
  NCLASSES <- nclasses
  nmodels <- length(NCLASSES)

  llca_list <- list()
  # Loop to create and fit the models for different nclasses:
  for(NK in 1:nmodels) {

    print(paste0("Model nclasses = ", NCLASSES[NK]) )

    # Ensure that nclasses is an integer:
    nclasses <- as.integer(NCLASSES[NK])
    # Check that nclasses is positive:
    if(nclasses < 1L) {
      stop("nclasses must be a positive integer")
    }

    ## store original call
    mc  <- match.call()

    #### Create the model ####

    # Get the model specification:
    full_model <- get_full_lca_model(data_list = data_list, item = item,
                                     nclasses = nclasses, model = model,
                                     control = control)
    list2env(full_model, envir = environment())

    # Get the short model specification (in logarithm and probability scale) with
    # labels for each parameter:
    short_model <- get_short_lca_model(data = data, item = item, nclasses = nclasses,
                                       lca_trans = lca_trans, model = model)
    list2env(short_model, envir = environment())

    #### Create the structures ####

    # Generate the structures for optimization:
    structures <- get_lca_structures(data_list = data_list,
                                     full_model = full_model)
    list2env(structures, envir = environment())

    #### Fit the model ####

    result <- vector("list") # Initialize the object to be returned

    # Model information:
    modelInfo <- list(item = item,
                      nobs = nobs,
                      npatterns = npatterns,
                      nparam = nparam,
                      df = npatterns - nparam,
                      ntrans = ntrans,
                      model = log_model,
                      prob_model = prob_model)

    # Data for the optimization algorithms:
    Optim <- list(data = data,
                  data_list = data_list,
                  control_manifold = control_manifold,
                  control_transform = control_transform,
                  control_estimator = control_estimator,
                  control = control)

    # Fit the model or just get the model specification:
    if(!do.fit) {

      llca_list[[NK]] <- new("llca",
                             version            = as.character( packageVersion('latent') ),
                             call               = mc, # matched call
                             timing             = numeric(), # timing information
                             modelInfo          = modelInfo, # modelInfo
                             Optim              = Optim, # Optim
                             parameters         = list(),
                             transformed_pars   = list(),
                             posterior          = matrix(),
                             state              = vector(),
                             loglik             = numeric(), # loglik values and info
                             loglik_case        = numeric(),
                             summary_table      = list(),
                             ClassConditional   = list(),
                             RespConditional    = list(),
                             probCat            = list()
      )

      next

    }

    # Perform the optimization (fit the model):
    x <- optimizer(control_manifold = control_manifold,
                   control_transform = control_transform,
                   control_estimator = control_estimator,
                   control_optimizer = control)

    # Collect all the information about the optimization:

    Optim <- list(data = data,
                  data_list = data_list,
                  control_manifold = control_manifold,
                  control_transform = control_transform,
                  control_estimator = control_estimator,
                  control = control,
                  opt = x)

    #### Process the outputs ####

    ## Collect information from the fitted model ##
    # Logarithm likelihood:
    loglik <- -x$f
    # Logarithm likelihood of each response pattern:
    loglik_case <- x$outputs$estimators$vectors[[1]][[3]]
    # Sum of logarithm likelihoods by response pattern:
    loglik_pattern <- weights * loglik_case

    ## Summary table with information for each response pattern ##
    # Estimated counts for each response pattern:
    estimated <- round(exp(loglik_case) * nobs, 0)
    # Posterior:
    posterior <- exp(matrix(x$outputs$estimators$matrices[[1]][[2]],
                            nrow = npatterns, ncol = nclasses))
    colnames(posterior) <- paste("P(", "Class", 1:nclasses, "|Y)", sep = "")
    # Posterior classification:
    state <- apply(posterior, MARGIN = 1, FUN = which.max)
    # Data table of response patterns:
    summary_table <- cbind(Pattern = patterns + 1,
                           Observed = weights,
                           Estimated = estimated,
                           Posterior = posterior,
                           State = state,
                           loglik_case = loglik_case,
                           loglik_pattern = loglik_pattern)
    summary_table <- as.data.frame(summary_table)
    # Sort the patterns by increasing order:
    summary_table <- summary_table[do.call(order, as.data.frame(summary_table)), ]
    rownames(summary_table) <- paste("pattern", 1:nrow(summary_table), sep = "")

    # Check the existence of gaussian items:
    gauss <- "gaussian" %in% item
    # Check the existence of multinomial items:
    multin <- "multinomial" %in% item

    # Create the parameter table:
    selection <- match(unlist(prob_model), transparameters_labels)
    out <- x$transparameters[selection]
    user_model <- fill_list_with_vector(prob_model, out)
    user_model <- allnumeric(user_model)

    selection <- match(unlist(log_model), transparameters_labels)
    out <- x$transparameters[selection]
    raw_model <- fill_list_with_vector(log_model, out)
    raw_model <- allnumeric(raw_model)

    ClassConditional <- RespConditional <- probCat <- list()

    # Additional outputs only for multinomial models:
    # if(multin & !gauss) {
    #
    #   classes <- x$outputs$estimators$vectors[[2]][[1]]
    #
    #   # For multinomial items:
    #   indices <- which(item == "multinomial")
    #   # Initialize the outputs:
    #   nitems <- length(indices) # Number of items in the data
    #   ClassConditional <- RespConditional <- vector("list", length = nitems)
    #   names(ClassConditional) <- names(RespConditional) <- paste("Item", 1:nitems)
    #   multinomial_conditionals <- lapply(x$outputs$estimators$list_matrices[[1]][[1]], FUN = \(x) {
    #     # Extract P(y|X):
    #     ncategories <- length(x)/nclasses
    #     result <- matrix(x, nrow = ncategories, ncol = nclasses)
    #     colnames(result) <- paste("Class", 1:nclasses, sep = "")
    #     rownames(result) <- paste("Category", 1:ncategories, sep = "")
    #     return(result)
    #   })
    #   names(multinomial_conditionals) <- paste("Item", 1:nitems)
    #
    #   ClassConditional[indices] <- multinomial_conditionals
    #
    #   probCat <- lapply(multinomial_conditionals, FUN = \(mat) {
    #     # Calculate P(y|X)*P(X), the joint probability:
    #     jointp <- t(mat) * classes
    #     # Calculate P(y), the denominator of the posterior:
    #     probCat <- colSums(jointp)
    #     return(probCat)
    #   })
    #
    #   RespConditional[indices] <- lapply(multinomial_conditionals, FUN = \(mat) {
    #     # Calculate P(y|X)*P(X), the joint probability:
    #     jointp <- t(mat) * classes
    #     # Calculate P(y), the denominator of the posterior:
    #     probCat <- colSums(jointp)
    #     # Calculate P(X|y) = P(y|X)*P(X)/P(y), the posterior:
    #     posterior <- t(jointp) / probCat
    #     return(posterior)
    #   })
    #
    #   result$ClassConditional <- multinomial_conditionals
    #   result$RespConditional <- RespConditional[indices]
    #   names(result$RespConditional) <- paste("Item", 1:nitems)
    #   result$probCat <- probCat
    #
    # }

    #### Return ####

    loglik_case <- loglik_case[short2full]
    posterior <- posterior[short2full, ]
    state <- state[short2full]
    rownames(posterior) <- names(state) <-
      names(loglik_case) <- rownames(data)

    llca_list[[NK]] <- new("llca",
                           version            = as.character( packageVersion('latent') ),
                           call               = mc, # matched call
                           timing             = x$elapsed, # timing information
                           modelInfo          = modelInfo, # modelInfo
                           Optim              = Optim, # Optim
                           parameters         = raw_model,
                           transformed_pars   = user_model,
                           posterior          = posterior,
                           state              = state,
                           loglik             = loglik, # loglik values and info
                           loglik_case        = loglik_case,
                           summary_table      = summary_table,
                           ClassConditional   = ClassConditional,
                           RespConditional    = RespConditional,
                           probCat            = probCat
    )


  }

  class(llca_list) <- "llca.list"

  if(nmodels == 1){
    return(llca_list[[1]])
  }else{
    return(llca_list)
  }

  return(result)

}

# Andres Chull, fit indices
# input output lavaan
# missing data CFA EFA LCA
# 1:nclasses
# Estimated values

# parameters package R
