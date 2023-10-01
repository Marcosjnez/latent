#' @title
#' Compute fit measures for exploratory factor models.
#' @description
#'
#' Compute fit measures for exploratory factor models.
#'
#' @usage
#'
#' fit(model, nobs = NULL)
#'
#' @param model Object of class efa or cfa.
#' @param nobs Sample size. Defaults to NULL.
#'
#' @details \code{fit}... to be explained
#'
#' @return Vector of fit measures.
#'
#' @author
#'
#' Vithor R. Franco & Marcos Jiménez
#'
#' @export
fit <- function(efa, nobs = NULL, digits = 3) {
  # Check if nobs was provided
  if(is.null(nobs)) {
    if(is.null(efa$modelInfo$nobs) | isTRUE(efa$modelInfo$nobs == 0L)) {
      warning("Sample size was not provided. Some Chi-squared-based statistics will not be computed.")
      nobs <- NA
    } else {
      nobs <- efa$modelInfo$nobs
    }
  }

  # Basic measures
  ObjFn <- efa$efa$f
  ObjFn_null <- efa$modelInfo$f_null
  residuals <- efa$efa$residuals
  p <- efa$modelInfo$nvars
  q <- efa$modelInfo$nfactors
  correction <- if(is.null(nobs)) NULL else nobs-1-1/6*(2*p+5)-2/3*q
  df_null <- efa$modelInfo$df_null
  df <- efa$modelInfo$df
  t <- df_null - df
  chisq_null <- if(is.null(nobs)) NULL else {nobs - 1} * ObjFn_null
  chisq_null.corrected <- if(is.null(nobs)) NULL else correction * ObjFn_null
  chisq <- if(is.null(nobs)) NULL else {nobs - 1} * ObjFn
  pvalue <- if(is.null(nobs)) NULL else 1 - stats::pchisq(chisq, df)
  chisq.corrected <- if(is.null(nobs)) NULL else correction * ObjFn
  pvalue.corrected <- if(is.null(nobs)) NULL else 1 - stats::pchisq(chisq.corrected, df)

  # Incremental fit indices
  cfi <- if(is.null(nobs)) NULL else {
    # 1 - max(chisq - df, 0) / max(chisq_null - df_null, 0)
    {max(chisq_null-df_null,0) - max(chisq - df, 0)} / max(chisq_null-df_null, 0)
  }
  cfi.corrected <- if(is.null(nobs)) NULL else {
    # 1 - max(chisq.corrected - df, 0) / max(chisq_null.corrected - df_null, 0)
    {max(chisq_null.corrected - df_null, 0) - max(chisq.corrected - df, 0)} /
      max(chisq_null.corrected - df_null, 0)
  }
  tli <- if(is.null(nobs)) NULL else {
    {chisq_null/df_null - chisq/df}/{chisq_null/{df_null-1}}
  }
  tli.corrected <- if(is.null(nobs)) NULL else {
    {chisq_null.corrected / df_null - chisq.corrected / df} /
      {chisq_null.corrected / {df_null-1}}
  }
  nfi <- {ObjFn_null - ObjFn} / ObjFn_null

  # Absolute fit indices
  rmsea <- if(is.null(nobs)) NULL else {
    sqrt(max(chisq - df, 0) / {df * {nobs-1}})
  }
  rmsea.corrected <- if(is.null(nobs)) NULL else {
    sqrt(max(chisq.corrected - df, 0) / {df * {nobs-1}})
  }
  srmr <- sqrt(sum(residuals[lower.tri(residuals,diag=T)]^2)/
                 {{efa$modelInfo$nvars*{efa$modelInfo$nvars+1}}/2})
  lavsrc <- max(abs(residuals))

  # Comparative fit indices
  aic           <- if(is.null(nobs)) NULL else chisq + {2 * t}
  aic.corrected  <- if(is.null(nobs)) NULL else chisq.corrected + {2 * t}
  bic           <- if(is.null(nobs)) NULL else chisq + {log(nobs) * t}
  bic.corrected  <- if(is.null(nobs)) NULL else chisq.corrected + {log(nobs) * t}
  hq            <- if(is.null(nobs)) NULL else chisq + {2 * log(log(nobs)) * t}
  hq.corrected   <- if(is.null(nobs)) NULL else chisq.corrected + {2 * log(log(nobs)) * t}
  ecvi          <- if(is.null(nobs)) NULL else {chisq/{nobs-1}} + {2*{t/{nobs-1}}}
  ecvi.corrected <- if(is.null(nobs)) NULL else {chisq.corrected/{nobs-1}} + {2*{t/{nobs-1}}}

  indices <- cbind(c(chisq, df, pvalue, chisq_null, df_null, rmsea, srmr, cfi, tli, nfi, lavsrc),
                   c(chisq.corrected, df, pvalue.corrected, chisq_null.corrected, df_null,
                     rmsea.corrected, srmr, cfi.corrected, tli.corrected, nfi, lavsrc))
  indices <- round(indices, digits)
  rownames(indices) <- c("Chi-square", "DF", "p-value",
                         "Chi-square (null)", "DF (null)",
                         "RMSEA", "SRMR", "CFI", "TLI", "NFI", "Max Res")
  colnames(indices) <- c("Unscaled", "Barlett's correction")
  indices <- as.data.frame(indices)
  name.width <- max(sapply(names(indices), nchar))
  names(indices) <- format(names(indices), width = name.width, justify = "centre")
  cat("\n", "Fit Indices \n", sep=""); print(format(indices,
                                                    justify = "centre", width = name.width))

  information <- cbind(c(aic, bic, hq, ecvi),
                       c(aic.corrected, bic.corrected, hq.corrected, ecvi.corrected))
  information <- round(information, digits)
  rownames(information) <- c("AIC", "BIC", "HQ", "ECVI")
  colnames(information) <- c("Unscaled", "Barlett's correction")
  information <- as.data.frame(information)
  name.width <- max(sapply(names(information), nchar))
  names(information) <- format(names(information), width = name.width, justify = "centre")
  cat("\n", "Information Criteria \n", sep=""); print(format(information,
                                                             justify = "centre", width = name.width))

  # Results
  Results <- list(indices = indices, information = information)
  return(invisible(Results))
}
