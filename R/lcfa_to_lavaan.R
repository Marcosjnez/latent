# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 21/08/2026
#'
#' Convert a Fitted lcfa Object to lavaan
#'
#' @param object A fitted object inheriting from class \code{"lcfa"}.
#' @param ... Additional arguments used only when rebuilding the lavaan
#'   scaffold.
#'
#' @return A \code{"lavaan"} object populated with estimates and covariance
#'   information from the fitted \code{"lcfa"} object.
#'
#' @details
#' This converter uses a lavaan scaffold to preserve lavaan's parameter-table
#' and model-matrix conventions. It depends on non-exported lavaan constructors,
#' so compatibility is checked at run time and a descriptive error is returned
#' when the installed lavaan version is incompatible.
#'
#' @export
lcfa_to_lavaan <- function(object, ...) {

  #### Check inputs ####

  if(!inherits(object, "lcfa")) {
    stop("object must inherit from class 'lcfa'.")
  }

  if(!requireNamespace("lavaan", quietly = TRUE)) {
    stop("The lavaan package is required.")
  }

  if(length(object@Optim$transparameters) == 0L) {
    stop("The lcfa object has not been fitted.")
  }

  if(!is.null(object@dataList$likelihood) &&
     !identical(tolower(object@dataList$likelihood), "normal")) {
    stop("lcfa_to_lavaan() currently supports likelihood = 'normal' only.")
  }

  #### Lavaan scaffold ####

  if(isTRUE(object@dataList$direct_fiml)) {

    args <- object@dataList$args

    scaffold_args <- list(
      model = args$model,
      data = object@dataList$data,
      estimator = "ML",
      missing = "ml",
      std.lv = isTRUE(args$std.lv),
      # The working data have already been standardized by lcfa when std.ov
      # was requested.
      std.ov = FALSE,
      meanstructure = TRUE,
      do.fit = FALSE,
      warn = FALSE
    )

    if(object@dataList$ngroups > 1L) {
      scaffold_args$group <- object@dataList$group
    }

    if(!is.null(object@dataList$likelihood)) {
      scaffold_args$likelihood <- object@dataList$likelihood
    }

    scaffold_args <- utils::modifyList(
      scaffold_args,
      list(...)
    )

    scaffold <- do.call(
      lavaan::cfa,
      scaffold_args
    )

  } else {

    scaffold <- object@dataList$LAV

    if(!inherits(scaffold, "lavaan")) {
      stop("The lavaan scaffold is missing from the lcfa object.")
    }

  }

  lavpartable <- scaffold@ParTable
  lavmodel <- scaffold@Model
  lavsamplestats <- scaffold@SampleStats
  lavdata <- scaffold@Data
  lavoptions <- scaffold@Options
  lavcache <- scaffold@Cache
  ngroups <- object@dataList$ngroups

  #### Parameter estimates ####

  transparameters <- object@Optim$transparameters

  if(is.null(names(transparameters))) {
    names(transparameters) <-
      object@modelInfo$transparameters_labels
  }

  est <- lavpartable$start
  fixed <- lavpartable$free == 0L

  for(i in seq_along(est)) {

    plabel <- lavpartable$plabel[i]

    if(!is.na(plabel) &&
       nzchar(plabel) &&
       plabel %in% names(transparameters)) {
      est[i] <- transparameters[[plabel]]
    }

  }

  est[fixed] <- lavpartable$start[fixed]

  #### Standard errors and free-parameter covariance ####

  se_vector <- rep(0, length(est))
  all_free_rows <- which(lavpartable$free > 0L)
  free_ids <- sort(unique(
    lavpartable$free[all_free_rows]
  ))
  free_rows <- match(
    free_ids,
    lavpartable$free
  )
  free_plabels <- lavpartable$plabel[free_rows]

  stored_se <- tryCatch(
    object@Optim$SE$se,
    error = function(e) NULL
  )
  stored_vcov <- tryCatch(
    object@Optim$SE$VCOV,
    error = function(e) NULL
  )

  if(!is.null(stored_se)) {

    if(is.null(names(stored_se)) &&
       length(stored_se) == length(free_plabels)) {
      names(stored_se) <- object@modelInfo$parameters_labels
    }

    all_free_plabels <- lavpartable$plabel[all_free_rows]
    matched <- match(
      all_free_plabels,
      names(stored_se)
    )
    available <- !is.na(matched)

    se_vector[all_free_rows[available]] <-
      stored_se[matched[available]]

  }

  lcfa_vcov <- NULL

  if(!is.null(stored_vcov)) {

    stored_vcov <- validate_covariance_matrix(
      stored_vcov,
      object_name = "stored lcfa variance-covariance matrix"
    )

    if(!is.null(rownames(stored_vcov)) &&
       all(free_plabels %in% rownames(stored_vcov)) &&
       all(free_plabels %in% colnames(stored_vcov))) {

      lcfa_vcov <- stored_vcov[
        free_plabels,
        free_plabels,
        drop = FALSE
      ]

    } else if(nrow(stored_vcov) == length(free_plabels)) {

      lcfa_vcov <- stored_vcov
      rownames(lcfa_vcov) <-
        colnames(lcfa_vcov) <- free_plabels

    } else {

      warning("The stored lcfa VCOV could not be aligned with lavaan's ",
              "free-parameter order and was omitted.")

    }

  }

  lavpartable$est <- est
  lavpartable$se <- se_vector
  lavpartable$start <- est

  #### Rebuild the lavaan model ####

  required_internal <- c(
    "lav_model",
    "lav_model_x2glist",
    "lav_h1_implied_logl",
    "lav_model_loglik",
    "lav_model_test",
    "lav_model_fit",
    "lav_partable_attributes",
    "lav_object_independence"
  )

  namespace <- asNamespace("lavaan")

  unavailable <- required_internal[
    !vapply(
      required_internal,
      exists,
      FUN.VALUE = logical(1L),
      envir = namespace,
      inherits = FALSE
    )
  ]

  if(length(unavailable) > 0L) {
    stop("The installed lavaan version does not provide required internal ",
         "constructor(s): ",
         paste(unavailable, collapse = ", "))
  }

  lav_model <- getFromNamespace(
    "lav_model",
    "lavaan"
  )
  lav_model_x2glist <- getFromNamespace(
    "lav_model_x2glist",
    "lavaan"
  )

  lavmodel <- lav_model(
    lavpartable = lavpartable,
    lavoptions = lavoptions,
    th.idx = lavmodel@th.idx
  )

  x <- est[free_rows]

  lavmodel@GLIST <- lav_model_x2glist(
    lavmodel = lavmodel,
    x = x,
    type = "free"
  )

  lavimplied <- lavaan::lav_model_implied(lavmodel)

  lavoptions$se <- if(is.null(lcfa_vcov)) {
    "none"
  } else {
    "standard"
  }
  lavoptions$test <- "standard"
  lavoptions$do.fit <- TRUE

  #### Log-likelihood and test statistics ####

  lav_h1_implied_logl <- getFromNamespace(
    "lav_h1_implied_logl",
    "lavaan"
  )
  lav_model_loglik <- getFromNamespace(
    "lav_model_loglik",
    "lavaan"
  )
  lav_model_test <- getFromNamespace(
    "lav_model_test",
    "lavaan"
  )
  lav_model_fit <- getFromNamespace(
    "lav_model_fit",
    "lavaan"
  )

  lavh1 <- tryCatch(
    lav_h1_implied_logl(
      lavdata = lavdata,
      lavsamplestats = lavsamplestats,
      lavoptions = lavoptions
    ),
    error = function(e) list()
  )

  lavloglik <- tryCatch(
    lav_model_loglik(
      lavdata = lavdata,
      lavsamplestats = lavsamplestats,
      lavimplied = lavimplied,
      lavmodel = lavmodel,
      lavoptions = lavoptions,
      lavh1 = lavh1
    ),
    error = function(e) list()
  )

  fit_matrix <- lcfa_fit_matrix(object)

  fit_indices <- getfit.lcfa(
    object,
    digits = 12L,
    fit_matrix = fit_matrix
  )

  loss_group <- fit_matrix[
    "penalized_loss",
    seq_len(ngroups),
    drop = TRUE
  ]
  loglik_group <- fit_matrix[
    "penalized_loglik",
    seq_len(ngroups),
    drop = TRUE
  ]

  fx <- if(is.finite(fit_indices[["chisq"]])) {

    fit_indices[["chisq"]]/
      sum(unlist(object@dataList$nobs))

  } else {

    fit_matrix["penalized_loss", "overall"]

  }

  attr(x, "fx") <- fx
  attr(x, "iterations") <-
    as.integer(object@Optim$iterations)
  attr(x, "converged") <-
    isTRUE(object@Optim$convergence)
  attr(x, "control") <- list()

  lavtest <- tryCatch(
    lav_model_test(
      lavobject = NULL,
      lavmodel = lavmodel,
      lavpartable = lavpartable,
      lavsamplestats = lavsamplestats,
      lavimplied = lavimplied,
      lavh1 = lavh1,
      lavoptions = lavoptions,
      x = x,
      VCOV = lcfa_vcov,
      lavcache = lavcache,
      lavdata = lavdata,
      lavloglik = if(is.list(lavloglik)) lavloglik else NULL
    ),
    error = function(e) list()
  )

  lavfit <- lav_model_fit(
    lavpartable = lavpartable,
    lavmodel = lavmodel,
    lavimplied = lavimplied,
    x = x,
    VCOV = lcfa_vcov,
    TEST = lavtest
  )

  lavoptim <- list(
    x = x,
    npar = length(x),
    iterations = as.integer(object@Optim$iterations),
    converged = isTRUE(object@Optim$convergence),
    fx = as.numeric(fx),
    fx.group = as.numeric(loss_group),
    logl.group = as.numeric(loglik_group),
    logl = as.numeric(
      fit_matrix["penalized_loglik", "overall"]
    ),
    control = list()
  )

  lavvcov <- list(
    se = if(is.null(lcfa_vcov)) "none" else "standard"
  )

  if(!is.null(lcfa_vcov)) {
    lavvcov$vcov <- lcfa_vcov
  }

  lav_partable_attributes <- getFromNamespace(
    "lav_partable_attributes",
    "lavaan"
  )
  lavpta <- lav_partable_attributes(
    lavpartable
  )

  timing <- list(
    total = if(length(object@timing) >= 1L) {
      object@timing[[1L]]
    } else {
      0
    },
    optim = if(length(object@timing) >= 2L) {
      object@timing[[2L]]
    } else {
      0
    }
  )

  result <- new(
    "lavaan",
    version = as.character(
      utils::packageVersion("lavaan")
    ),
    call = object@call,
    timing = timing,
    Options = lavoptions,
    ParTable = lavpartable,
    pta = lavpta,
    Data = lavdata,
    SampleStats = lavsamplestats,
    Model = lavmodel,
    Cache = lavcache,
    Fit = lavfit,
    boot = list(),
    optim = lavoptim,
    implied = lavimplied,
    loglik = if(is.list(lavloglik)) lavloglik else list(),
    vcov = lavvcov,
    test = if(is.list(lavtest)) lavtest else list(),
    h1 = lavh1,
    baseline = list(),
    internal = list(latent.version = object@version),
    external = list()
  )

  lav_object_independence <- getFromNamespace(
    "lav_object_independence",
    "lavaan"
  )

  result@baseline <- tryCatch({

    baseline <- lav_object_independence(
      object = result
    )

    list(
      partable = baseline@ParTable,
      test = baseline@test
    )

  }, error = function(e) list())

  #### Result ####

  return(result)

}
