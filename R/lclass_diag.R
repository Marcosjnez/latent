# Author: Mauricio Garnier-Villarreal
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 10/09/2026

#' Classification Diagnostics for Latent Class Analysis
#'
#' Computes various classification diagnostics from a fitted latent class model.
#'
#' @param x An object of class `"llca"` from the \pkg{latent} package.
#' @param digits Number of decimal places to use for rounding (default = 4).
#' @param type Character vector specifying which diagnostics to compute.
#'   Available options:
#'   \itemize{
#'     \item \code{"Entropy"} – Entropy‑based R² (classification uncertainty).
#'     \item \code{"AvePP"} – Average posterior probability for each class.
#'     \item \code{"OCC"} – Odds of Correct Classification.
#'     \item \code{"Overall.Misclassification"} – Overall misclassification rate.
#'     \item \code{"Misclassification.per.class"} – Misclassification rate per true class.
#'     \item \code{"Sum.Posterior"} – Class proportions based on posterior sums.
#'     \item \code{"Sum.Mostlikely"} – Class proportions based on modal assignment.
#'     \item \code{"Mostlikely.Class"} – Classification probabilities (true × assigned).
#'     \item \code{"Avg.Mostlikely"} – Average posterior probabilities by assigned class.
#'     \item \code{"all"} – A convenience shortcut to request all of the above.
#'   }
#' @param ... Additional arguments passed to \code{latInspect} or \code{getfit}.
#'
#' @return An object of class \code{"lclassd"}, a list containing the requested diagnostics.
#'
#' @examples
#' \dontrun{
#'   # Fit a model
#'   library(latent)
#'   # Fit a 3‑class model with multinomial indicators
#'   fit <- lca(data = gss82, nclasses = 3,
#'              multinomial = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"))
#'
#'   # Get default diagnostics (Entropy, AvePP, Mostlikely.Class, Sum.Mostlikely)
#'   diag <- lclass_diag(fit)
#'   print(diag)
#'
#'   # Request all diagnostics
#'   diag_all <- lclass_diag(fit, type = "all")
#'   print(diag_all)
#' }
#' @export
lclass_diag <- function(x, digits = 4,
                        type = c("Entropy", "AvePP", "Mostlikely.Class", "Sum.Mostlikely"),
                        ...) {

  if(!inherits(x, "llca") || !length(x@Optim)) {
    stop("x must be a fitted llca object.")
  }
  if(!is.null(x@dataList$lengths)) {
    stop("Use latInspect(x, what = 'classification') for hidden Markov models.")
  }

  # ---- Argument validation ----
  valid_types <- c("Entropy", "AvePP", "OCC", "Overall.Misclassification",
                   "Misclassification.per.class", "Sum.Posterior",
                   "Sum.Mostlikely", "Mostlikely.Class", "Avg.Mostlikely")

  if ("all" %in% type) {
    type <- valid_types
  }

  type <- match.arg(type, choices = valid_types, several.ok = TRUE)

  # Work with response-pattern frequencies, including any observation weights.
  post_probs <- .llca_posterior(x)
  weights <- x@dataList$pattern_weights
  classes <- seq_len(ncol(post_probs))
  class_names <- paste0("P(", x@dataList$class_names, "|data)")
  state <- max.col(post_probs, ties.method = "first")
  assigned <- diag(length(classes))[state, , drop = FALSE]
  counts <- crossprod(assigned, post_probs * weights)
  modal_count <- rowSums(counts)
  posterior_count <- colSums(counts)
  total <- sum(weights)

  # Rows of the first matrix are assigned classes; rows of the second are true classes.
  avg_most <- sweep(counts, 1L, modal_count, "/")
  avg_most[!is.finite(avg_most)] <- NA_real_
  dimnames(avg_most) <- list(paste0("assigned.", classes), paste0("meanprob.", class_names))
  class_most <- sweep(t(counts), 1L, posterior_count, "/")
  class_most[!is.finite(class_most)] <- NA_real_
  dimnames(class_most) <- list(paste0("avgprob.", classes), paste0("assigned.", classes))
  misclass <- setNames(1 - diag(class_most), class_names)

  if("AvePP" %in% type) {
    AV <- matrix(NA_real_, length(classes), 5L,
                 dimnames = list(NULL, c("mean", "sd", "median", "min", "max")))
    for(i in classes) {
      keep <- state == i
      if(!any(keep)) next
      p <- post_probs[keep, i]
      w <- weights[keep]
      mean_p <- avg_most[i, i]
      sd_p <- if(sum(w) > 1) sqrt(sum(w * (p-mean_p)^2)/(sum(w)-1)) else NA_real_
      order_p <- order(p)
      cumulative <- cumsum(w[order_p])
      middle <- which(cumulative >= sum(w)/2)[1L]
      median_p <- p[order_p[middle]]
      if(cumulative[middle] == sum(w)/2 && middle < length(p)) {
        median_p <- (median_p + p[order_p[middle + 1L]])/2
      }
      AV[i, ] <- c(mean_p, sd_p, median_p, min(p), max(p))
    }
    AV <- as.data.frame(round(AV, digits))
  }

  if("OCC" %in% type) {
    prior <- latInspect(x, what = "classes", sort = FALSE, ...)
    mean_p <- diag(avg_most)
    OCC <- (mean_p / (1 - mean_p)) / (prior / (1 - prior))
    OCC[is.nan(OCC)] <- NA_real_
  }

  size_table <- function(count) {
    data.frame(class = class_names, count = round(count, digits),
               proportion = round(count/total, digits), row.names = NULL)
  }

  out <- lapply(type, function(thetype) {
    switch(thetype,
           "Mostlikely.Class"           = round(class_most, digits),
           "Avg.Mostlikely"             = round(avg_most, digits),
           "Sum.Posterior"              = size_table(posterior_count),
           "Sum.Mostlikely"             = size_table(modal_count),
           "AvePP"                      = AV,
           "Entropy"                    = getfit(x, digits = digits)[["R2_entropy"]],
           "OCC"                        = round(OCC, digits),
           "Misclassification.per.class" = round(misclass, digits),
           "Overall.Misclassification"  = round(1 - sum(diag(counts))/total, digits))
  })

  names(out) <- type
  class(out) <- "lclassd"
  return(out)
}

#' Print Method for Latent Class Diagnostics
#'
#' @param x An object of class \code{"lclassd"} returned by \code{lclass_diag()}.
#' @param ... Additional arguments passed to \code{print.default} or \code{print.data.frame}.
#'
#' @method print lclassd
#' @export
print.lclassd <- function(x, ...) {
  displayed_classes <- attr(x, "class_names", exact = TRUE)
  cat("\n\u2014\u2014\u2014 Latent Class Classification Diagnostics \u2014\u2014\u2014\n")
  cat(rep("\u2500", 50), "\n\n", sep = "")

  # If "Mostlikely.Class" is present, display number of classes
  if (!is.null(x[["Mostlikely.Class"]])) {
    nclass <- ncol(x[["Mostlikely.Class"]])
    cat("Number of latent classes:", nclass, "\n\n")
  }

  for (i in seq_along(x)) {
    nm <- names(x)[i]
    obj <- x[[i]]

    # ---- Custom formatting for each known component ----
    if (nm == "Entropy") {
      cat("\u25c6 Entropy R\u00b2:\n")
      cat("  ", obj, "\n\n")

    } else if (nm == "AvePP") {
      cat("\u25c6 Average posterior probabilities (AvePP) by true class:\n")
      df <- obj
      colnames(df) <- c("Mean", "SD", "Median", "Min", "Max")
      rownames(df) <- if(is.null(displayed_classes)) paste("Class", seq_len(nrow(df))) else displayed_classes
      print(df, ...)
      cat("\n")

    } else if (nm == "OCC") {
      cat("\u25c6 Odds of Correct Classification (OCC):\n")
      occ_vec <- obj
      names(occ_vec) <- if(is.null(displayed_classes)) paste("Class", seq_along(occ_vec)) else displayed_classes
      print(occ_vec, ...)
      cat("\n")

    } else if (nm == "Overall.Misclassification") {
      cat("\u25c6 Overall misclassification rate:\n")
      cat("  ", obj, "\n\n")

    } else if (nm == "Misclassification.per.class") {
      cat("\u25c6 Misclassification rate by true class:\n")
      mc_vec <- obj
      names(mc_vec) <- if(is.null(displayed_classes)) paste("Class", seq_along(mc_vec)) else displayed_classes
      print(mc_vec, ...)
      cat("\n")

    } else if (nm == "Mostlikely.Class") {
      cat("\u25c6 Classification probabilities (rows = true class, columns = assigned class):\n")
      mat <- obj
      rownames(mat) <- if(is.null(displayed_classes)) paste("True class", seq_len(nrow(mat))) else paste("True", displayed_classes)
      colnames(mat) <- if(is.null(displayed_classes)) paste("Assigned", seq_len(ncol(mat))) else paste("Assigned", displayed_classes)
      print(mat, ...)
      cat("\n")

    } else if (nm == "Avg.Mostlikely") {
      cat("\u25c6 Average posterior probabilities by assigned class:\n")
      mat <- obj
      rownames(mat) <- if(is.null(displayed_classes)) paste("Assigned class", seq_len(nrow(mat))) else paste("Assigned", displayed_classes)
      colnames(mat) <- if(is.null(displayed_classes)) paste("Prob for class", seq_len(ncol(mat))) else paste("Prob for", displayed_classes)
      print(mat, ...)
      cat("\n")

    } else if (nm == "Sum.Posterior") {
      cat("\u25c6 Class proportions based on posterior sums:\n")
      df <- obj
      print(df, row.names = FALSE, ...)
      cat("\n")

    } else if (nm == "Sum.Mostlikely") {
      cat("\u25c6 Class proportions based on modal assignment:\n")
      df <- obj
      print(df, row.names = FALSE, ...)
      cat("\n")

    } else {
      # Fallback for any unknown components
      cat("[[", i, "]] ", nm, ":\n", sep = "")
      cat(rep("\u2500", nchar(nm) + 6), "\n", sep = "")
      if (is.matrix(obj) || is.data.frame(obj)) {
        print(obj, ...)
      } else {
        print(obj, ...)
      }
      cat("\n")
    }
  }

  invisible(x)
}

