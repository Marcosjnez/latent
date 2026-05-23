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

  p <- ggplot2::ggplot(df) +
    ggplot2::aes(
      x     = .data[[bars]],
      y     = .data[["probability"]],
      fill  = .data[["response"]]
    ) +
    ggplot2::geom_bar(stat = "identity",
                      position = ggplot2::position_fill(reverse = TRUE)) +
    ggplot2::facet_wrap(facet, scales = "free_x") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), name = NULL) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_viridis_d(guide = ggplot2::guide_legend(reverse = TRUE)) +
    ggplot2::theme(
      axis.text.x   = ggplot2::element_text(angle = angle, hjust = 1),
      panel.spacing = ggplot2::unit(1, "lines")
    )

  if (bw) {
    p <- p + ggplot2::scale_fill_grey(guide = ggplot2::guide_legend(reverse = TRUE))
  }
  p
}
