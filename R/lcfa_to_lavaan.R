#' Convert a fitted lcfa object (latent) to a lavaan object
#'
#' @param object A fitted \code{lcfa} object from the \code{latent} package.
#' @param ... Additional arguments passed to \code{lavaan::cfa()}.
#'
#' @return A \code{lavaan} object populated with estimates from the \code{lcfa} fit.
#'
#' @details
#' This function reconstructs a \code{lavaan} S4 object from a fitted \code{lcfa}
#' object, allowing the use of lavaan post-processing functions such as
#' \code{fitMeasures()}, \code{lavResiduals()}, and \code{modindices()}.
#'
#' The conversion works by:
#' \enumerate{
#'   \item Building a lavaan scaffold (un-fitted) with the same model syntax, data,
#'         estimator, and options.
#'   \item Injecting lcfa parameter estimates into the lavaan parameter table via
#'         plabel matching.
#'   \item Rebuilding the lavaan Model object and implied moments from the updated
#'         parameter table.
#'   \item Computing test statistics, log-likelihood, baseline model, and VCOV
#'         using standard lavaan internals.
#' }
#'
#' @examples
#' \dontrun{
#' library(latent)
#' library(lavaan)
#'
#' HS.model <- ' visual  =~ x1 + x2 + x3
#'               textual =~ x4 + x5 + x6
#'               speed   =~ x7 + x8 + x9 '
#'
#' fit_latent <- lcfa(model = HS.model, data = HolzingerSwineford1939,
#'                    estimator = "ml", std.lv = TRUE, meanstructure = TRUE)
#' fit_lavaan <- lcfa_to_lavaan(fit_latent)
#' fitMeasures(fit_lavaan)
#' lavResiduals(fit_lavaan)
#' modindices(fit_lavaan, sort = TRUE)
#' }
#'
#' @export
lcfa_to_lavaan <- function(object, ...) {

  # ------------------------------------------------------------------
  # 0. Validate input
  # ------------------------------------------------------------------
  if (!inherits(object, "lcfa")) {
    stop("'object' must be a fitted lcfa object from the latent package.")
  }
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("The 'lavaan' package is required but not installed.")
  }

  # ------------------------------------------------------------------
  # 1. Recover original call arguments from the lcfa object
  # ------------------------------------------------------------------
  mc        <- object@call
  dataList <- object@dataList

  model_syntax <- eval(mc$model)
  ngroups      <- dataList$ngroups
  group_labels <- dataList$group_label   # e.g. "g1" or c("Pasteur","Grant")

  # Raw data / sample covariance
  raw_data    <- tryCatch(eval(mc$data,        envir = parent.frame()), error = function(e) NULL)
  sample_cov  <- tryCatch(eval(mc$sample.cov,  envir = parent.frame()), error = function(e) NULL)
  nobs_val    <- tryCatch(eval(mc$nobs,         envir = parent.frame()), error = function(e) NULL)

  # Group variable
  group_var <- tryCatch(eval(mc$group, envir = parent.frame()), error = function(e) NULL)

  # Estimator
  latent_est <- tolower(as.character(
    if (!is.null(mc$estimator)) eval(mc$estimator) else "ml"
  ))
  lavaan_estimator <- switch(latent_est,
                             "ml"         = "ML",
                             "fml"        = "ML",
                             "means_fml"  = "ML",
                             "uls"        = "ULS",
                             "means_uls"  = "ULS",
                             "dwls"       = "DWLS",
                             "means_dwls" = "DWLS",
                             "ML"
  )

  # ordered?
  ordered_val <- tryCatch(eval(mc$ordered, envir = parent.frame()), error = function(e) FALSE)
  if (is.null(ordered_val)) ordered_val <- FALSE

  # std.lv
  std_lv <- tryCatch(eval(mc$std.lv, envir = parent.frame()), error = function(e) TRUE)
  if (is.null(std_lv)) std_lv <- TRUE

  # std.ov
  std_ov <- tryCatch(eval(mc$std.ov, envir = parent.frame()), error = function(e) FALSE)
  if (is.null(std_ov)) std_ov <- FALSE

  # meanstructure
  meanstructure <- tryCatch(eval(mc$meanstructure, envir = parent.frame()), error = function(e) NULL)
  if (is.null(meanstructure)) {
    if (isTRUE(ordered_val)) {
      meanstructure <- FALSE
    } else {
      missing_arg <- tryCatch(eval(mc$missing, envir = parent.frame()), error = function(e) NULL)
      if (!is.null(missing_arg) && tolower(missing_arg) == "fiml") {
        meanstructure <- TRUE
      } else {
        meanstructure <- FALSE
      }
    }
  }

  # likelihood
  likelihood_arg <- tryCatch(eval(mc$likelihood, envir = parent.frame()), error = function(e) NULL)

  # ------------------------------------------------------------------
  # 2. Build the lavaan scaffold (do.fit = FALSE)
  # ------------------------------------------------------------------
  cfa_args <- list(
    model         = model_syntax,
    estimator     = lavaan_estimator,
    std.lv        = std_lv,
    std.ov        = std_ov,
    meanstructure = TRUE,   # always TRUE to ensure correct Delta / A1 dimensions
    do.fit        = FALSE
  )

  if (!is.null(raw_data)) {
    cfa_args$data <- raw_data
  } else if (!is.null(sample_cov) && !is.null(nobs_val)) {
    cfa_args$sample.cov  <- sample_cov
    cfa_args$sample.nobs <- nobs_val
  } else {
    stop("Cannot recover data or sample covariance from the lcfa call.")
  }

  if (!is.null(group_var) && nzchar(as.character(group_var))) {
    cfa_args$group <- as.character(group_var)
  }
  if (isTRUE(ordered_val)) {
    cfa_args$ordered <- TRUE
  }
  if (!is.null(likelihood_arg)) {
    cfa_args$likelihood <- as.character(likelihood_arg)
  }

  dots <- list(...)
  cfa_args <- c(cfa_args, dots)

  lav_scaffold <- do.call(lavaan::cfa, cfa_args)

  # ------------------------------------------------------------------
  # 3. Extract key scaffold components
  # ------------------------------------------------------------------
  lavpartable    <- lav_scaffold@ParTable
  lavmodel       <- lav_scaffold@Model
  lavsamplestats <- lav_scaffold@SampleStats
  lavdata        <- lav_scaffold@Data
  lavoptions     <- lav_scaffold@Options
  lavcache       <- lav_scaffold@Cache

  npar_lav <- length(lavpartable$lhs)

  # ------------------------------------------------------------------
  # 4.  Build est_vec by matching plabels to lcfa named parameter vector
  #
  #     lcfa stores all *transformed* parameters in:
  #       object@Optim$transparameters   (named by modelInfo$transparameters_labels)
  #     and the free (estimated) parameters in:
  #       object@Optim$parameters        (named by modelInfo$parameters_labels)
  #
  #     The names are lavaan plabels (".p1.", ".p2.", ...) so we can match
  #     directly on lavpartable$plabel.
  # ------------------------------------------------------------------
  tp    <- object@Optim$transparameters   # full named vector incl. fixed
  pars  <- object@Optim$parameters        # free parameters only

  est_vec <- lavpartable$start            # initialise with scaffold start values
  se_vec  <- rep(0, npar_lav)

  # Match by plabel (most reliable: plabels are unique per free parameter)
  if (!is.null(lavpartable$plabel) && !is.null(names(tp))) {
    for (r in seq_len(npar_lav)) {
      pl <- lavpartable$plabel[r]
      if (!is.null(pl) && nzchar(pl) && pl %in% names(tp)) {
        est_vec[r] <- tp[[pl]]
      }
    }
  }

  # For fixed parameters, always keep the fixed (ustart) value
  fixed_mask <- lavpartable$free == 0L
  est_vec[fixed_mask] <- lavpartable$start[fixed_mask]

  # ------------------------------------------------------------------
  # 5.  Fill SE from lcfa by plabel matching
  # ------------------------------------------------------------------
  lcfa_se   <- tryCatch(object@Optim$SE$se,   error = function(e) NULL)
  lcfa_vcov <- tryCatch(object@Optim$SE$vcov, error = function(e) NULL)

  if (!is.null(lcfa_se) && !is.null(names(lcfa_se))) {
    for (r in seq_len(npar_lav)) {
      pl <- lavpartable$plabel[r]
      if (!is.null(pl) && nzchar(pl) && pl %in% names(lcfa_se)) {
        se_vec[r] <- lcfa_se[[pl]]
      }
    }
  } else if (!is.null(lcfa_se)) {
    # Fallback: positional matching on free parameters
    free_idx <- which(!fixed_mask)
    if (length(lcfa_se) == length(free_idx)) {
      se_vec[free_idx] <- as.numeric(lcfa_se)
    }
  }

  # ------------------------------------------------------------------
  # 6. Update the parameter table
  # ------------------------------------------------------------------
  lavpartable$est   <- est_vec
  lavpartable$se    <- se_vec
  lavpartable$start <- est_vec

  # ------------------------------------------------------------------
  # 7. Rebuild lavModel from the updated parameter table so that all
  #    internal index arrays (m.free.idx, x.free.idx, Delta rows, ...)
  #    stay consistent with the actual matrix contents.
  # ------------------------------------------------------------------
  lavmodel <- lavaan:::lav_model(
    lavpartable = lavpartable,
    lavoptions  = lavoptions,
    th.idx      = lavmodel@th.idx
  )

  # Extract the free-parameter vector in lavaan's canonical order
  free_rows  <- which(lavpartable$free > 0L)
  free_order <- order(lavpartable$free[free_rows])   # sort by free index (1,2,3,...)
  x_ordered  <- est_vec[free_rows][free_order]

  lavmodel@GLIST <- lavaan:::lav_model_x2glist(
    lavmodel = lavmodel,
    x        = x_ordered,
    type     = "free"
  )

  # ------------------------------------------------------------------
  # 8. Model-implied moments
  # ------------------------------------------------------------------
  lavimplied <- lavaan::lav_model_implied(lavmodel)

  # ------------------------------------------------------------------
  # 9. Options tweaks required by downstream lavaan functions
  # ------------------------------------------------------------------
  lavoptions$se     <- "standard"
  lavoptions$test   <- "standard"
  lavoptions$do.fit <- TRUE

  # Reset NACOV so lavaan recomputes it internally (avoids dimension
  # mismatches that arise when latent used a different ACOV structure)
  lavsamplestats@NACOV <- vector("list", ngroups)

  # ------------------------------------------------------------------
  # 10. Saturated (H1) model
  # ------------------------------------------------------------------
  lavh1 <- tryCatch(
    lavaan:::lav_h1_implied_logl(
      lavdata        = lavdata,
      lavsamplestats = lavsamplestats,
      lavoptions     = lavoptions
    ),
    error = function(e) list()
  )

  # ------------------------------------------------------------------
  # 11. Log-likelihood (lavaan style)
  # ------------------------------------------------------------------
  lavloglik <- tryCatch(
    lavaan:::lav_model_loglik(
      lavdata        = lavdata,
      lavsamplestats = lavsamplestats,
      lavimplied     = lavimplied,
      lavmodel       = lavmodel,
      lavoptions     = lavoptions,
      lavh1          = lavh1
    ),
    error = function(e) list()
  )

  # ------------------------------------------------------------------
  # 12. Optimizer result bookkeeping
  # ------------------------------------------------------------------
  x_est <- lavaan::lav_model_get_parameters(lavmodel, type = "free")
  npar  <- length(x_est)

  iterations <- tryCatch(as.integer(object@Optim$iterations), error = function(e) 0L)
  if (length(iterations) == 0L || is.na(iterations)) iterations <- 0L

  converged <- tryCatch(isTRUE(object@Optim$convergence), error = function(e) TRUE)

  # ------------------------------------------------------------------
  # fx / loss scaling convention:
  #
  # lavaan internally uses fx = F (the full loss), where for ML:
  #   F = sum_g { w_g * F_g }
  # For a single group model lcfa stores F = 2 * (log-likelihood difference),
  # i.e. it matches  nobs * F_ml  before dividing, so the stored @loss
  # already equals what lavaan would call 2*fx for a single group.
  # Dividing by 2 therefore recovers lavaan's fx for single-group models.
  #
  # For multiple groups lcfa sums the per-group losses with group weights
  # w_g = n_g / N, which already matches lavaan's multigroup fx directly
  # (no additional factor of 2 needed).
  # ------------------------------------------------------------------
  loss_val <- object@loss
  if (length(loss_val) == 0L || is.na(loss_val)) loss_val <- 0

  fx_val <- loss_val / 2

  # Per-group fx and logl
  fit_by_group <- tryCatch(
    latInspect(object, what = "fit"),
    error = function(e) NULL
  )

  if (!is.null(fit_by_group)) {
    if (ngroups == 1L) {
      fx_group <- vapply(fit_by_group,
                         function(x) sum(unlist(x[["loss"]])) / 2,
                         numeric(1))
    } else {
      fx_group <- vapply(fit_by_group,
                         function(x) sum(unlist(x[["loss"]])),
                         numeric(1))
    }
    logl_group <- vapply(fit_by_group,
                         function(x) sum(unlist(x[["loglik"]])),
                         numeric(1))
  } else {
    fx_group   <- rep(fx_val  / ngroups, ngroups)
    logl_group <- rep(as.numeric(NA),    ngroups)
  }

  logl <- object@loglik
  if (length(logl) == 0L || is.na(logl)) logl <- as.numeric(NA)

  attr(fx_val, "fx.group") <- fx_group
  attr(x_est,  "fx")         <- fx_val
  attr(x_est,  "iterations") <- iterations
  attr(x_est,  "converged")  <- converged
  attr(x_est,  "control")    <- list()

  lavoptim <- list(
    x          = x_est,
    npar       = npar,
    iterations = iterations,
    converged  = converged,
    fx         = fx_val,
    fx.group   = fx_group,
    logl.group = logl_group,
    logl       = logl,
    control    = list()
  )

  # ------------------------------------------------------------------
  # 13. Test statistics
  # ------------------------------------------------------------------
  lavtest <- tryCatch(
    lavaan:::lav_model_test(
      lavobject      = NULL,
      lavmodel       = lavmodel,
      lavpartable    = lavpartable,
      lavsamplestats = lavsamplestats,
      lavimplied     = lavimplied,
      lavh1          = lavh1,
      lavoptions     = lavoptions,
      x              = x_est,
      VCOV           = lcfa_vcov,
      lavcache       = lavcache,
      lavdata        = lavdata,
      lavloglik      = if (is.list(lavloglik)) lavloglik else NULL
    ),
    error = function(e) list()
  )

  # ------------------------------------------------------------------
  # 14. Fit S4 object
  # ------------------------------------------------------------------
  lavfit <- tryCatch(
    lavaan:::lav_model_fit(
      lavpartable = lavpartable,
      lavmodel    = lavmodel,
      lavimplied  = lavimplied,
      x           = x_est,
      VCOV        = lcfa_vcov,
      TEST        = lavtest
    ),
    error = function(e) {
      new("Fit",
          npar       = as.integer(npar),
          x          = as.numeric(x_est),
          partrace   = matrix(0, 0L, 0L),
          start      = lavpartable$start,
          est        = lavpartable$est,
          se         = se_vec,
          fx         = as.numeric(fx_val),
          fx.group   = as.numeric(fx_group),
          logl       = if (!is.na(logl)) as.numeric(logl) else as.numeric(NA),
          logl.group = if (!anyNA(logl_group)) as.numeric(logl_group)
          else rep(as.numeric(NA), ngroups),
          iterations = iterations,
          converged  = converged,
          control    = list(),
          Sigma.hat  = if (!is.null(lavimplied$cov))  lavimplied$cov  else list(),
          Mu.hat     = if (!is.null(lavimplied$mean)) lavimplied$mean else list(),
          TH         = if (!is.null(lavimplied$th))   lavimplied$th   else list(),
          test       = if (is.list(lavtest)) lavtest else list()
      )
    }
  )

  # ------------------------------------------------------------------
  # 15. VCOV structure
  # ------------------------------------------------------------------
  lavvcov <- list(se = "standard")
  if (!is.null(lcfa_vcov)) {
    lavvcov$vcov <- lcfa_vcov
  }

  # ------------------------------------------------------------------
  # 16. Partial table attributes
  # ------------------------------------------------------------------
  lavpta <- lavaan::lav_partable_attributes(lavpartable)

  # ------------------------------------------------------------------
  # 17. Timing
  # ------------------------------------------------------------------
  timing <- list(
    total = if (length(object@timing) >= 1L) object@timing[[1L]] else 0,
    optim = if (length(object@timing) >= 2L) object@timing[[2L]] else 0
  )

  # ------------------------------------------------------------------
  # 18. Assemble the lavaan object (without baseline first)
  # ------------------------------------------------------------------
  result <- new("lavaan",
                version     = as.character(utils::packageVersion("lavaan")),
                call        = mc,
                timing      = timing,
                Options     = lavoptions,
                ParTable    = lavpartable,
                pta         = lavpta,
                Data        = lavdata,
                SampleStats = lavsamplestats,
                Model       = lavmodel,
                Cache       = lavcache,
                Fit         = lavfit,
                boot        = list(),
                optim       = lavoptim,
                implied     = lavimplied,
                loglik      = if (is.list(lavloglik)) lavloglik else list(),
                vcov        = lavvcov,
                test        = if (is.list(lavtest)) lavtest else list(),
                h1          = lavh1,
                baseline    = list(),
                internal    = list(latent.version = object@version),
                external    = list()
  )

  # ------------------------------------------------------------------
  # 19. Baseline (independence) model
  # ------------------------------------------------------------------
  lavbaseline <- tryCatch({
    tt <- lavaan:::lav_object_independence(object = result)
    list(partable = tt@ParTable, test = tt@test)
  }, error = function(e) list())

  result@baseline <- lavbaseline

  return(result)
}
