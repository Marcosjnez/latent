#' Plot method for LCA model fit comparisons (scree plot)
#'
#' Plots the requested information criteria across latent class solutions.
#'
#' @param x An object of class \code{"getfit.llcalist"}, typically from \code{getfit()}.
#' @param indices Character vector of information criteria to plot. If \code{NULL},
#'   defaults to \code{c("AICp","BICp")} when the object is penalised, otherwise
#'   \code{c("AIC","BIC")}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{ggplot} object.
#'
#' @method plot getfit.llcalist
#' @export
#'
#' @examples
#' \dontrun{
#' fit_ind_cont <- getfit(fit)
#' plot(fit_ind_cont)
#' plot(fit_ind_cont, indices = c("AIC", "BIC", "KIC"))
#' }
plot.getfit.llcalist <- function(x, indices = NULL, ...) {

  penalized <- attr(x, "penalized")

  if (is.null(indices)) {
    if (penalized) {
      indices <- c("AICp", "BICp")
    } else {
      indices <- c("AIC", "BIC")
    }
  }

  # Extract required columns (nclasses + chosen indices)
  wide_df <- data.frame(x[, c("nclasses", indices)])

  # Stack into long format using base R
  long_df <- do.call(rbind, lapply(indices, function(idx) {
    data.frame(
      nclasses = wide_df$nclasses,
      IC       = idx,
      Value    = wide_df[[idx]],
      stringsAsFactors = FALSE
    )
  }))

  ggplot(long_df, aes(x = nclasses, y = Value, group = IC)) +
    geom_point(aes(color = IC)) +
    geom_line(aes(color = IC)) +
    labs(
      title = "Model fit by number of classes",
      x     = "Class Solution",
      y     = "Value"
    ) +
    theme_minimal() +
    scale_color_discrete(name = "Information Criterion")
}

