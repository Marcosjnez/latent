#' Plot mean and standard deviation profiles for Gaussian LCA indicators
#'
#' @param data Data frame from \code{reshape_lca_continuous}.
#' @param title Plot title (default \code{NULL}).
#' @param ylab y-axis label.
#' @param caption Plot caption.
#' @param nrow,ncol Number of rows/columns in facet wrap.
#' @param point_size Size of the mean points.
#' @param err_width Width of error bars.
#' @param base_size Base font size.
#' @param angle Angle for x-axis labels.
#' @param ... Other arguments (currently ignored).
#' @return A ggplot object.
#' @export
plot_lca_gaussian <- function(
    data,
    title       = NULL,
    ylab        = "Mean score (\u00b11 SD)",
    caption     = "Error bars indicate \u00b11 standard deviation",
    nrow        = NULL,
    ncol        = NULL,
    point_size  = 2.5,
    err_width   = 0.2,
    base_size   = 12,
    angle       = 45,
    ...
) {
  # Ensure variable order is as provided
  if (!is.factor(data$variable)) {
    data$variable <- factor(data$variable, levels = unique(data$variable))
  }

  p <- ggplot(data, aes(x = variable, y = mean, group = 1)) +
    geom_point(size = point_size, color = "steelblue") +
    geom_line(color = "steelblue") +
    geom_errorbar(
      aes(ymin = mean - sd, ymax = mean + sd),
      width = err_width,
      color = "grey40"
    ) +
    facet_wrap(~ class, nrow = nrow, ncol = ncol) +
    scale_y_continuous(limits = c(NA, NA)) +
    labs(
      x       = NULL,
      y       = ylab,
      title   = title,
      caption = caption
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.x   = element_text(angle = angle, hjust = 1),
      panel.spacing = unit(1, "lines")
    )

  p
}

