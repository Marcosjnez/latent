#' Plot response probabilities for multinomial LCA indicators
#'
#' @param data Data frame from \code{reshape_lca_multinomial}.
#' @param variables Optional character vector of variables to include.
#' @param bars Variable to map on x-axis (default \code{"variable"}).
#' @param facet Faceting variable (default \code{"class"}).
#' @param bw Logical; if \code{TRUE}, use greyscale instead of colour.
#' @param ... Other arguments (currently ignored).
#' @return A ggplot object.
#' @export
plot_lca_multinomial <- function(
    data,
    variables = NULL,
    bars = "variable",
    facet = "class",
    bw = FALSE,
    angle = 45,
    ...
) {
  df <- data

  if (!is.null(variables)) {
    df <- subset(df, variable %in% variables)
  }

  resp_vals <- unique(df$response)
  if (all(!is.na(suppressWarnings(as.integer(as.character(resp_vals)))))) {
    df$response <- ordered(df$response)
  }

  if (nrow(df) < 2) stop("Insufficient data to plot.")

  p <- ggplot(df) +
    aes(
      x     = .data[[bars]],
      y     = .data[["probability"]],
      fill  = .data[["response"]]
    ) +
    geom_bar(stat = "identity", position = position_fill(reverse = TRUE)) +
    facet_wrap(facet, scales = "free_x") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), name = NULL) +
    theme_bw() +
    scale_fill_viridis_d(guide = guide_legend(reverse = TRUE)) +
    theme(
      axis.text.x   = element_text(angle = angle, hjust = 1),
      panel.spacing = unit(1, "lines")
    )

  if (bw) {
    p <- p + scale_fill_grey(guide = guide_legend(reverse = TRUE))
  }
  p
}
