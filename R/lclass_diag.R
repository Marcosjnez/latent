#' Classification Diagnostics for Latent Class Analysis
#'
#' Computes various classification diagnostics from a fitted latent class model.
#'
#' @param x An object of class `"lca"` from the \pkg{latent} package.
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
#'   fit <- lca(data = gss82, nclasses = 3,item = rep("multinomial", ncol(gss82)))
#'
#'   # Get default diagnostics (Entropy, AvePP, Mostlikely.Class, Sum.Mostlikely)
#'   diag <- lclass_diag(fit)
#'   print(diag)
#'
#'   # Request all diagnostics
#'   diag_all <- lclass_diag(fit, type = "all")
#'   print(diag_all)
#' }
#' @importFrom psych describeBy
#' @export
lclass_diag <- function(x, digits = 4,
                        type = c("Entropy", "AvePP", "Mostlikely.Class", "Sum.Mostlikely"),
                        ...) {

  # ---- Argument validation ----
  valid_types <- c("Entropy", "AvePP", "OCC", "Overall.Misclassification",
                   "Misclassification.per.class", "Sum.Posterior",
                   "Sum.Mostlikely", "Mostlikely.Class", "Avg.Mostlikely")

  if ("all" %in% type) {
    type <- valid_types
  }

  type <- match.arg(type, choices = valid_types, several.ok = TRUE)

  # ---- Extract basic quantities from model ----
  post_probs <- latInspect(x, what = "posterior", digits = digits, ...)
  state      <- latInspect(x, what = "state", ...)
  class_prop <- latInspect(x, what = "profile", digits = digits, ...)$class

  # ---- Helper functions (kept inside for portability) ----
  avgprobs_mostlikely <- function(post_prob, class = NULL) {
    if (is.null(dim(post_prob))) return(1)
    if (is.null(class)) class <- apply(post_prob, 1, which.max)
    tab <- t(sapply(1:ncol(post_prob), function(i) {
      colMeans(post_prob[class == i, , drop = FALSE])
    }))
    rownames(tab) <- paste0("assigned.", 1:nrow(tab))
    colnames(tab) <- paste0("meanprob.", colnames(tab))
    return(tab)
  }

  classification_probs_mostlikely <- function(post_prob, class = NULL) {
    if (is.null(dim(post_prob))) return(1)
    if (is.null(class)) class <- apply(post_prob, 1, which.max)
    avg_probs <- avgprobs_mostlikely(post_prob, class)
    avg_probs[is.na(avg_probs)] <- 0
    class_counts <- as.integer(table(ordered(class, levels = 1:ncol(post_prob))))
    tab <- avg_probs * class_counts
    tab <- tab / matrix(colSums(avg_probs * class_counts),
                        ncol = ncol(tab), nrow = nrow(tab), byrow = TRUE)
    rownames(tab) <- paste0("assigned.", 1:nrow(tab))
    colnames(tab) <- paste0("avgprob.", 1:nrow(tab))
    return(t(tab))   # rows = true class, cols = assigned class
  }

  avepp <- function(post_prob, state, digits = 4) {
    desc_k <- psych::describeBy(post_prob, group = state, digits = digits, mat = TRUE)
    unq <- sort(unique(state))
    out <- NULL
    for (j in seq_along(unq)) {
      temp0 <- subset(desc_k, group1 == unq[j])
      temp <- temp0[j, c("mean", "sd", "median", "min", "max")]
      out <- rbind(out, temp)
    }
    return(out)
  }

  sum_mostlikely <- function(post_prob, state, digits = 4) {
    classif <- table(state)
    mix_names <- colnames(post_prob)
    data.frame(class = mix_names,
               count = as.vector(classif),
               proportion = round(as.vector(prop.table(classif)), digits))
  }

  sum_postprob <- function(post_prob, class_pp) {
    mix_names <- colnames(post_prob)
    numObs <- nrow(post_prob)
    data.frame(class = mix_names,
               count = as.vector(class_pp * numObs),
               proportion = as.vector(class_pp))
  }

  # ---- Pre‑compute commonly used matrices (avoids redundant calculations) ----
  avg_most   <- avgprobs_mostlikely(post_probs)          # for "Avg.Mostlikely"
  class_most <- classification_probs_mostlikely(post_probs)  # for "Mostlikely.Class"

  # ---- AvePP ----
  AV <- avepp(post_probs, state, digits = digits)

  # ---- Odds of Correct Classification (OCC) ----
  OCC_k <- (AV$mean / (1 - AV$mean)) / (class_prop / (1 - class_prop))

  # ---- Misclassification per true class ----
  misclass_by_true <- 1 - diag(class_most)
  names(misclass_by_true) <- colnames(post_probs)

  # ---- Overall misclassification (correctly weighted by true class proportions) ----
  overall_misclass <- sum(class_prop * misclass_by_true)

  # ---- Sum.Mostlikely (modal assignment proportions) ----
  sum_ml <- sum_mostlikely(post_probs, state, digits)

  # ---- Assemble requested diagnostics ----
  out <- lapply(type, function(thetype) {
    switch(thetype,
           "Mostlikely.Class"           = round(class_most, digits),
           "Avg.Mostlikely"              = round(avg_most, digits),
           "Sum.Posterior"                = sum_postprob(post_probs, class_prop),
           "Sum.Mostlikely"                = sum_ml,
           "AvePP"                          = AV,
           "Entropy"                        = getfit(x, digits = digits)[["R2_entropy"]],
           "OCC"                            = round(OCC_k, digits),
           "Misclassification.per.class"    = round(misclass_by_true, digits),
           "Overall.Misclassification"      = round(overall_misclass, digits)
    )
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
#' @export
print.lclassd <- function(x, ...) {
  cat("\n——— Latent Class Classification Diagnostics ———\n")
  cat(rep("─", 50), "\n\n", sep = "")

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
      cat("◆ Entropy R²:\n")
      cat("  ", obj, "\n\n")

    } else if (nm == "AvePP") {
      cat("◆ Average posterior probabilities (AvePP) by true class:\n")
      df <- obj
      colnames(df) <- c("Mean", "SD", "Median", "Min", "Max")
      rownames(df) <- paste("Class", seq_len(nrow(df)))
      print(df, ...)
      cat("\n")

    } else if (nm == "OCC") {
      cat("◆ Odds of Correct Classification (OCC):\n")
      occ_vec <- obj
      names(occ_vec) <- paste("Class", seq_along(occ_vec))
      print(occ_vec, ...)
      cat("\n")

    } else if (nm == "Overall.Misclassification") {
      cat("◆ Overall misclassification rate:\n")
      cat("  ", obj, "\n\n")

    } else if (nm == "Misclassification.per.class") {
      cat("◆ Misclassification rate by true class:\n")
      mc_vec <- obj
      names(mc_vec) <- paste("Class", seq_along(mc_vec))
      print(mc_vec, ...)
      cat("\n")

    } else if (nm == "Mostlikely.Class") {
      cat("◆ Classification probabilities (rows = true class, columns = assigned class):\n")
      mat <- obj
      rownames(mat) <- paste("True class", seq_len(nrow(mat)))
      colnames(mat) <- paste("Assigned", seq_len(ncol(mat)))
      print(mat, ...)
      cat("\n")

    } else if (nm == "Avg.Mostlikely") {
      cat("◆ Average posterior probabilities by assigned class:\n")
      mat <- obj
      rownames(mat) <- paste("Assigned class", seq_len(nrow(mat)))
      colnames(mat) <- paste("Prob for class", seq_len(ncol(mat)))
      print(mat, ...)
      cat("\n")

    } else if (nm == "Sum.Posterior") {
      cat("◆ Class proportions based on posterior sums:\n")
      df <- obj
      print(df, row.names = FALSE, ...)
      cat("\n")

    } else if (nm == "Sum.Mostlikely") {
      cat("◆ Class proportions based on modal assignment:\n")
      df <- obj
      print(df, row.names = FALSE, ...)
      cat("\n")

    } else {
      # Fallback for any unknown components
      cat("[[", i, "]] ", nm, ":\n", sep = "")
      cat(rep("─", nchar(nm) + 6), "\n", sep = "")
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

