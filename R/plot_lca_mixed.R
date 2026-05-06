#' Plot mixed LCA results (both Gaussian and multinomial indicators)
#'
#' @param mixed_data Output from \code{reshape_lca_mixed}.
#' @param ... Extra arguments passed to \code{plot_lca_gaussian} and \code{plot_lca_multinomial}.
#' @return A list of ggplot objects, named \code{gaussian} and \code{multinomial}.
#' @export
plot_lca_mixed <- function(mixed_data, ...) {

  plots <- list()

  # Continuous
  if (!is.null(mixed_data$gaussian) && nrow(mixed_data$gaussian) > 0) {
    plots$gaussian <- plot_lca_gaussian(mixed_data$gaussian, ...)
  } else {
    warning("No gaussian indicators to plot.")
    plots$gaussian <- NULL
  }

  # Categorical
  if (!is.null(mixed_data$multinomial) && nrow(mixed_data$multinomial) > 0) {
    plots$multinomial <- plot_lca_multinomial(mixed_data$multinomial, ...)
  } else {
    warning("No multinomial indicators to plot.")
    plots$multinomial <- NULL
  }

  plots
}
