# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 07/03/2026
#'
#' @title
#' Fit a Confirmatory Factor Analysis (CFA) model with lavaan syntax.
#'
#' @usage
#'
#' lcfa(data, model = NULL, estimator = "ml",
#' ordered = FALSE, group = NULL,
#' sample.cov = NULL, nobs = NULL,
#' positive = FALSE, penalties = TRUE,
#' missing = "pairwise.complete.obs",
#' std.lv = FALSE, do.fit = TRUE,
#' message = FALSE, mimic = 'latent',
#' control = NULL, ...)
#'
#' @param data data frame or matrix.
#' @param model lavaan's model syntax.
#' @param estimator Available estimators: "ml", "uls", and "dwls". Defaults to "ml".
#' @param ordered Logical. Defaults to TRUE.
#' @param group String. Name of the variable that splits the data in different groups.
#' @param sample.cov Covariance matrix between the items. Defaults to NULL.
#' @param nobs Number of observations. Defaults to NULL.
#' @param positive Force a positive-definite solution. Defaults to FALSE.
#' @param penalties list of penalty terms for the parameters.
#' @param missing Method to handle missing data.
#' @param std.lv Provide the parameters of the standardized model.
#' @param acov String. "standard" or "robust". Default is "standard".
#' @param do.fit TRUE to fit the model and FALSE to return only the model setup. Defaults to TRUE.
#' @param message Logical. Defaults to FALSE.
#' @param mimic String. Choose the output you want to obtain. Defaults to 'latent'.
#' @param control List of control parameters for the optimization algorithm. See 'details' for more information.
#' @param ... Additional lavaan arguments. See ?lavaan for more information.
#'
#' @details \code{lcfa} estimates confirmatory factor models.
#'
#' @return List with the following objects:
#' \item{version}{Version number of 'latent' when the model was estimated.}
#' \item{call}{Code used to estimate the model.}
#' \item{ModelInfo}{Model information.}
#' \item{Optim}{Output of the optimizer.}
#' \item{parameters}{Structure with all model parameters.}
#' \item{transparameters}{Structure with all transformed model parameters.}
#' \item{loglik}{Logarithm likelihood of the model.}
#' \item{penalized_loglik}{Logarithm likelihood + logarithm priors of the model.}
#'
#' @examples
#'
#' \dontrun{
#' # The famous Holzinger and Swineford (1939) example
#' HS.model <- ' visual  =~ x1 + x2 + x3
#'               textual =~ x4 + x5 + x6
#'               speed   =~ x7 + x8 + x9 '
#'
#' fit <- lcfa(model = HS.model, data = HolzingerSwineford1939)
#' summary(fit, digits = 3L)
#'}
#'
#' @export
lcfa <- function(data, model = NULL, estimator = "ml",
                 ordered = FALSE, group = NULL,
                 sample.cov = NULL, nobs = NULL,
                 positive = FALSE, penalties = TRUE,
                 missing = "pairwise.complete.obs",
                 std.lv = TRUE, acov = "standard",
                 do.fit = TRUE, message = FALSE,
                 mimic = "latent", control = NULL,
                 ...) {

  if(ordered) {
    cor <- "poly"
    control$deltaparam <- TRUE
  } else {
    cor <- "pearson"
  }

  # Check the arguments to control_optimizer and create defaults:
  control$positive <- positive
  control$penalties <- penalties
  control$estimator <- tolower(estimator)
  control <- lcfa_control(control)

  # stop("ahh")
  if(is.null(group)) {
    ngroups <- 1L
  } else {
    group_label <- unique(data[[group]])
    ngroups <- length(group_label)
  }

  item_names <- extract_item_names(model, ngroups = ngroups)

  # Estimate the correlation matrix for each group::
  X <- correl <- nobs <- vector("list", length = ngroups)
  if(ngroups > 1) {
    names(X) <- names(correl) <- group_label
  }

  for(i in 1:ngroups) {

    if(ngroups > 1) {
      X[[i]] <- data[data[[group]] == group_label[i], item_names[[i]]]
    } else {
      X[[i]] <- data[, item_names[[i]]]
    }
    nobs[[i]] <- nrow(X[[i]])

    if(is.null(sample.cov)) {

      correl[[i]] <- correlation(data = X[[i]],
                                 item_names = item_names[[i]],
                                 cor = cor,
                                 estimator = estimator,
                                 acov = acov,
                                 nobs = nobs[[i]],
                                 missing = missing)

    } else {

      correl[[i]]$R <- sample.cov
      p <- nrow(sample.cov)
      # ACOV?
      correl[[i]]$W <- matrix(1, nrow = p, ncol = p)

    }

  }

  # Extract the lavaan model:
  model_syntax <- model
  # stop("ahhh")
  LAV <- lavaan::cfa(model = model_syntax,
                     data = data,
                     sample.cov = sample.cov,
                     sample.nobs = nobs,
                     group = group,
                     # estimator = estimator,
                     # ordered = ordered,
                     # parameterization = "theta",
                     std.lv = std.lv,
                     do.fit = FALSE,
                     ...)

  LAV@Options$positive <- positive

  # extract slots from dummy lavaan object
  lavpartable    <- LAV@ParTable
  lavmodel       <- LAV@Model
  lavdata        <- LAV@Data
  lavoptions     <- LAV@Options
  lavsamplestats <- LAV@SampleStats
  lavcache       <- LAV@Cache
  timing         <- LAV@timing
  ngroups        <- LAV@Model@ngroups
  item_names     <- LAV@Data@ov.names
  nobs           <- LAV@Data@nobs
  group_label    <- LAV@Data@group.label
  item_label     <- LAV@Data@ov.names
  factor_label   <- replicate(ngroups, list(LAV@Model@dimNames[[1]][[2]]))
  X              <- LAV@Data@X

  # Rename columns:
  for(i in 1:ngroups) {
    colnames(X[[i]]) <- item_names[[i]]
  }

  # Model for the parameters:
  model <- getmodel_fromlavaan(LAV)

  if(ngroups == 1) model <- list(model)

  ## store original call
  mc  <- match.call()
  args <- as.list(match.call(expand.dots = TRUE))[-1]

  # Data and structure information:
  nitems <- as.list(lavmodel@nvar)
  npatterns <- lapply(nitems, FUN = \(p) 0.5*p*(p+1))
  nfactors <- lapply(model, FUN = \(x) ncol(x$lambda))

  data_list <- vector("list")
  data_list$ngroups <- ngroups
  data_list$data <- data
  data_list$data_per_group <- X
  data_list$nobs <- nobs
  data_list$nitems <- nitems
  data_list$npatterns <- npatterns
  data_list$nfactors <- nfactors
  data_list$correl <- correl
  data_list$positive <- positive
  data_list$estimator <- estimator
  data_list$cor <- cor
  data_list$group_label <- group_label
  data_list$item_label <- item_label
  data_list$factor_label <- factor_label
  data_list$args <- args

  #### Create the model ####

  # Get the model specification:
  full_model <- get_full_cfa_model(data_list = data_list,
                                   model = model,
                                   control = control)
  list2env(full_model, envir = environment())

  #### Create the structures ####

  # Generate the structures for optimization:
  structures <- get_cfa_structures(data_list = data_list,
                                   full_model = full_model,
                                   control = control)
  list2env(structures, envir = environment())

  #### Collect all the model information ####

  # Model information:
  modelInfo <- list(nobs = nobs,
                    nparam = nparam - rest,
                    npatterns = npatterns,
                    dof = sum(unlist(npatterns)) - nparam + rest,
                    ntrans = ntrans,
                    parameters_labels = parameters_labels,
                    transparameters_labels = transparameters_labels,
                    param = param,
                    trans = trans,
                    control_manifold = control_manifold,
                    control_transform = control_transform,
                    control_estimator = control_estimator,
                    control = control)

  #### Fit the model ####

  if(!do.fit) {

    lcfa_list <- new("lcfa",
                     version            = as.character( packageVersion('latent') ),
                     call               = mc, # matched call
                     timing             = numeric(), # timing information
                     data_list          = data_list,
                     modelInfo          = modelInfo,
                     Optim              = list(),
                     parameters         = list(),
                     transformed_pars   = list(),
                     loglik             = numeric(), # loglik values
                     penalized_loglik   = numeric(),
                     loss               = numeric(),
                     penalized_loss     = numeric()
    )

    return(lcfa_list)

  }

  if(message) {
    msg <- "Fitting the model"
    w <- nchar(msg) + 4
    cat("\n", "+", strrep("-", w), "+\n",
        "|  ", msg, "  |\n",
        "+", strrep("-", w), "+\n\n", sep = "")
  }

  control$cores <- min(control$rstarts, control$cores)
  # Fit the model:
  Optim <- optimizer(control_manifold = control_manifold,
                     control_transform = control_transform,
                     control_estimator = control_estimator,
                     control_optimizer = control)
  names(Optim$parameters) <- modelInfo$parameters_labels
  names(Optim$transparameters) <- modelInfo$transparameters_labels

  # Collect all the information about the optimization:

  elapsed <- Optim$elapsed

  #### Estimated model structures ####

  # Create the structures of transformed parameters:
  indices_trans <- match(unlist(modelInfo$trans),
                         modelInfo$transparameters_labels)
  values <- Optim$transparameters[indices_trans]
  transformed_pars <- fill_list_with_vector(modelInfo$trans, values)
  transformed_pars <- allnumeric(transformed_pars)

  # Create the structures of untransformed parameters:
  parameters <- transformed_pars[names(modelInfo$param)]

  #### Process the fit information ####

  # Get the indices of the estimator structures "cfa_dwls" and "cfa_ml":
  all_estimators <- unlist(lapply(modelInfo$control_estimator, FUN = \(x) x$estimator))
  indices_cfa <- which(all_estimators == "cfa_dwls" |
                         all_estimators == "cfa_ml" |
                         all_estimators == "cfa_fml")

  # Get the indices of the estimator structures "logdetmat" (penalties):
  indices_logdetmat <- which(all_estimators == "logdetmat")

  # Initialize the objects to be returned:
  loss <- penalized_loss <- loglik <- penalized_loglik <- penalty <-
    vector("list", length = ngroups)

  # For each group, extract the loss, penalized loss, loglik and penalized loglik
  for(i in 1:ngroups) {

    k <- indices_cfa[i]

    loss[[i]] <- c(Optim$outputs$estimators$doubles[[k]][[1]])
    loglik[[i]] <- c(Optim$outputs$estimators$doubles[[k]][[2]])

    # If there are penalties, add the penalties to the loss or loglik:
    if(length(indices_logdetmat) > 0) {

      l <- indices_logdetmat[i]
      penalty[[i]] <- c(Optim$outputs$estimators$doubles[[l]][[1]])
      penalized_loss[[i]] <- loss[[i]] + penalty[[i]]
      penalized_loglik[[i]] <- loglik[[i]] + penalty[[i]]

    } else {

      penalized_loss[[i]] <- loss[[i]]
      penalized_loglik[[i]] <- loglik[[i]]

    }

  }

  loss <- sum(unlist(loss))
  penalized_loss <- sum(unlist(penalized_loss))
  loglik <- sum(unlist(loglik))
  penalized_loglik <- sum(unlist(penalized_loglik))

  #### latent object ####

  result <- new("lcfa",
                version            = as.character( packageVersion('latent') ),
                call               = mc, # matched call
                timing             = elapsed, # timing information
                data_list          = data_list,
                modelInfo          = modelInfo,
                Optim              = Optim,
                parameters         = parameters,
                transformed_pars   = transformed_pars,
                loglik             = loglik, # loglik values
                penalized_loglik   = penalized_loglik,
                loss               = loss,
                penalized_loss     = penalized_loss
  )

  if(message) {
    msg <- "Computing standard errors"
    w <- nchar(msg) + 4
    cat("\n", "+", strrep("-", w), "+\n",
        "|  ", msg, "  |\n",
        "+", strrep("-", w), "+\n\n", sep = "")
  }

  # Standard errors:
  Optim$SE <- se(result, type = "standard", digits = 9)

  result@Optim <- Optim

  #### Return ####

  if(mimic == "lavaan") {

    ## for lavaan like object
    ## need to fill in these objects into the lav dummy objects
    ## will take me longer to find how to match these slots


    ###
    ### cecks for latent things that cannot be mimimc == lavaan
    ###

    ## optiosn
    lavoptions$se <- "standard" ## se type argument?
    lavoptions$test <- "standard" ## test type argument?
    lavoptions$do.fit <- do.fit

    ## parameters
    parsdf <- data.frame(plabel = names(Optim$parameters),
                         est = Optim$parameters)

    sedf <- data.frame(plabel = modelInfo$parameters_labels,
                       se = Optim$SE$se)

    parsdf <- merge(parsdf, sedf)
    #parsdf

    ## partable
    lavpartable$est <- NULL
    lavpartable <- merge(lavpartable, parsdf, all =T)
    lavpartable$est <- ifelse(lavpartable$free == 0,
                              lavpartable$start,
                              lavpartable$est)
    lavpartable$se <- ifelse(lavpartable$free == 0, 0,
                             lavpartable$se)
    lavpartable <- lavpartable[order(lavpartable$id),]
    lavpartable <- as.list(lavpartable)
    #lp

    ## pta
    pta <- LAV@pta
    pta$names <- names(lavpartable)

    ## lavmodel
    ## only for 1 group for now
    ## how to adjust multiple group labels?
    psi <- transformed_pars$psi.group1
    psi[upper.tri(psi)]<-t(psi)[upper.tri(psi)]

    lavmodel@GLIST$lambda <- transformed_pars$lambda.group1
    lavmodel@GLIST$psi <- psi
    lavmodel@GLIST$theta <- transformed_pars$theta.group1

    ## implied
    implied <- lav_model_implied(lavmodel)


    ## lavvcov
    ### order vcov from 1  up .. in names
    lavvcov <- LAV@vcov
    lavvcov$se <- lavoptions$se
    lavvcov$vcov <- Optim$SE$vcov

    ## test
    chi2 <- -2*(loglik - Optim$outputs$estimators$doubles[[1]][5])
    df <- modelInfo$dof
    pv <- 1 - pchisq(chi2, df)

    TEST <- list(standard =
                   list(test = lavoptions$test,
                        stat = chi2,
                        stat.group = chi2,
                        df = df,
                        refdistr = "chisq",
                        pvalue = pv))
    ## I am copying here as text, we will need to retrieve this information
    attributes(TEST) <- list(names = lavoptions$test,
                             info = list(ngroups = ngroups,
                                         group.label = group_label,
                                         information = c("expected", "expected"),
                                         h1.information = c("structured", "structured"),
                                         observed.information = c("hessian", "hessian") ))

    ### FIT
    lavfit <- lat_model_fit(lavpartable = lavpartable,
                            lavmodel = lavmodel,
                            x = lavpartable$est[lavpartable$free != 0],
                            VCOV = Optim$SE$vcov,
                            TEST = TEST,
                            loss = c(loss = loss, loglik = loglik),
                            Optim = Optim,
                            modelInfo = modelInfo)

    ##
    ## lavoptim
    optnames <- c('x','npar','iterations','converged','fx','fx.group','logl.group',
                  'logl','control')
    lavoptim <- lapply(optnames, function(x) slot(lavfit, x))
    names(lavoptim) <- optnames

    ## h1
    h1 <- list(implied = implied,
               logl = list(loglik = loglik,
                           loglik.group = loglik))

    ## add baseline

    ####
    ####
    result <- new("lavaan",
                  version      = as.character( packageVersion('latent') ),
                  call         = mc,                  # match.call
                  timing       = timing,              # list
                  Options      = lavoptions,          # list *
                  ParTable     = lavpartable,         # list *
                  pta          = pta,             # list
                  Data         = lavdata,             # S4 class
                  SampleStats  = lavsamplestats,      # S4 class
                  Model        = lavmodel,            # S4 class *
                  Cache        = lavcache,            # list
                  Fit          = lavfit,              # S4 class * blav_model_fit
                  boot         = list(),
                  optim        = lavoptim,
                  implied      = implied,          # list *
                  vcov         = lavvcov,             #*
                  test         = TEST,                #* blav_model_test
                  h1           = h1,
                  baseline     = list(),
                  internal     = list(),
                  external     = list()              # can add extra info from latent
    )

  }

  return(result)

}


### helper function for lavaan like objects
### edited from blavaan internal functions

lat_model_fit <- function (lavpartable = NULL,
                           lavmodel = NULL,
                           x = NULL,
                           VCOV = NULL,
                           TEST = NULL,
                           loss = NULL,
                           Optim = NULL,
                           modelInfo = NULL){
  stopifnot(is.list(lavpartable), inherits(lavmodel, c("Model",
                                                       "lavModel")))

  iterations <- Optim$iterations
  converged <- Optim$convergence
  fx <- loss["loss"]
  fx.group <- loss["loss"]
  logl.group <- loss["loglik"]
  logl <- loss["loglik"]
  control <- modelInfo$control
  attributes(fx) <- NULL
  x.copy <- x
  attributes(x.copy) <- NULL
  est <- lavpartable$est
  se <- lavpartable$se

  if (is.null(TEST)) {
    test <- list()
  }
  else {
    test <- TEST
  }
  implied <- lav_model_implied(lavmodel)

  if (lavmodel@conditional.x) {
    names(implied) <- c("cov", "mean", "slopes", "th", "group.w")
  }
  if (!is.null(attr(x, "partrace"))) {
    PARTRACE <- attr(x, "partrace")
  }
  else {
    PARTRACE <- matrix(0, 0L, 0L)
  }
  new("Fit",
      npar = as.integer(max(lavpartable$free)),
      x = x.copy,
      partrace = PARTRACE,
      start = lavpartable$start,
      est = est,
      se = se,
      fx = fx,
      fx.group = fx.group,
      logl = logl,
      logl.group = logl.group,
      iterations = as.integer(iterations),
      converged = converged,
      control = control,
      Sigma.hat = implied$cov,
      Mu.hat = implied$mean,
      TH = implied$th,
      test = test)
}

