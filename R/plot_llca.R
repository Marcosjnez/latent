# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 12/07/2026

#' Determine a panel layout
#'
#' @param npanels Number of panels.
#' @param nrow,ncol Optional requested numbers of panel rows and columns.
#'
#' @return An integer vector containing the number of rows and columns.
#'
#' @noRd
plot_panel_layout <- function(npanels, nrow = NULL, ncol = NULL) {

  if(length(npanels) != 1L || !is.numeric(npanels) || is.na(npanels) ||
     npanels < 1L) {
    stop("npanels must be a positive integer")
  }
  npanels <- as.integer(npanels)

  if(!is.null(nrow)) {
    if(length(nrow) != 1L || !is.numeric(nrow) || is.na(nrow) || nrow < 1L) {
      stop("nrow must be a positive integer or NULL")
    }
    nrow <- as.integer(nrow)
  }
  if(!is.null(ncol)) {
    if(length(ncol) != 1L || !is.numeric(ncol) || is.na(ncol) || ncol < 1L) {
      stop("ncol must be a positive integer or NULL")
    }
    ncol <- as.integer(ncol)
  }

  if(is.null(nrow) && is.null(ncol)) {
    ncol <- ceiling(sqrt(npanels))
    nrow <- ceiling(npanels/ncol)
  } else if(is.null(nrow)) {
    nrow <- ceiling(npanels/ncol)
  } else if(is.null(ncol)) {
    ncol <- ceiling(npanels/nrow)
  }

  if(nrow*ncol < npanels) {
    stop("nrow multiplied by ncol must provide at least npanels panels")
  }

  result <- c(nrow = nrow, ncol = ncol)

  #### Result ####

  return(result)

}

#' Draw rotated horizontal-axis labels
#'
#' @param at Numeric positions of the labels.
#' @param labels Character vector of labels.
#' @param angle Rotation angle in degrees.
#' @param cex Character expansion factor.
#'
#' @return The labels, invisibly.
#'
#' @noRd
plot_axis_labels <- function(at, labels, angle = 45, cex = 0.85) {

  usr <- graphics::par("usr")
  offset <- 0.045*diff(usr[3:4])
  y <- usr[3] - offset

  if(angle == 0) {
    graphics::text(at, y, labels = labels, xpd = NA, adj = c(0.5, 1),
                   cex = cex)
  } else {
    graphics::text(at, y, labels = labels, srt = angle, xpd = NA,
                   adj = c(1, 1), cex = cex)
  }

  #### Result ####

  return(invisible(labels))

}

#' Sort response levels for plotting
#'
#' @param values Response values.
#'
#' @return A character vector containing the response levels.
#'
#' @noRd
plot_response_levels <- function(values) {

  result <- unique(as.character(values))
  numeric_values <- suppressWarnings(as.numeric(result))
  if(length(result) > 1L && all(!is.na(numeric_values))) {
    result <- result[order(numeric_values)]
  }

  #### Result ####

  return(result)

}

#' Plot fit indices across latent class solutions
#'
#' Creates a scree-style base R plot comparing information criteria across the
#' models contained in a fitted latent class enumeration.
#'
#' @param x An object of class \code{"getfit.llcalist"}, normally returned by
#'   \code{getfit()} for an \code{llcalist} object.
#' @param indices Character vector containing the names of the fit indices to
#'   plot. When \code{NULL}, the function uses \code{AICp} and \code{BICp}
#'   for penalized models and \code{AIC} and \code{BIC} otherwise.
#' @param title Character string with the plot title.
#' @param xlab Character string with the horizontal-axis label.
#' @param ylab Character string with the vertical-axis label.
#' @param base_size Numeric value controlling the overall text size.
#' @param colors Optional character vector of line colors. Colors are recycled
#'   when fewer colors than indices are supplied.
#' @param pch Plotting symbol used for the fitted values.
#' @param lwd Width of the lines connecting fitted values.
#' @param ... Additional arguments. Currently unused.
#'
#' @details
#' Each selected information criterion is represented by a separate line. The
#' horizontal axis gives the number of latent classes and the vertical axis
#' gives the value of the selected criterion.
#'
#' @return The plotted data, invisibly.
#'
#' @examples
#' \dontrun{
#' fits <- lca(data = gss82, nclasses = 1:5,
#'             multinomial = c("PURPOSE", "ACCURACY"))
#' fit_indices <- getfit(fits)
#' plot(fit_indices)
#' plot(fit_indices, indices = c("AIC", "BIC", "SABIC"))
#' }
#'
#' @method plot getfit.llcalist
#' @export
plot.getfit.llcalist <- function(x, indices = NULL,
                                  title = "Model fit by number of classes",
                                  xlab = "Class solution", ylab = "Value",
                                  base_size = 12, colors = NULL, pch = 19,
                                  lwd = 2, ...) {

  if(!is.matrix(x) && !is.data.frame(x)) {
    stop("x must be a matrix or data.frame containing model-fit indices")
  }
  if(is.null(colnames(x)) || !"nclasses" %in% colnames(x)) {
    stop("x must contain a column named 'nclasses'")
  }
  if(length(base_size) != 1L || !is.numeric(base_size) || is.na(base_size) ||
     base_size <= 0) {
    stop("base_size must be a single positive number")
  }

  penalized <- isTRUE(attr(x, "penalized"))
  if(is.null(indices)) {
    if(penalized && all(c("AICp", "BICp") %in% colnames(x))) {
      indices <- c("AICp", "BICp")
    } else {
      indices <- c("AIC", "BIC")
    }
  }

  if(!is.character(indices) || length(indices) == 0L || anyNA(indices)) {
    stop("indices must be a non-empty character vector")
  }
  indices <- unique(indices)
  missing_indices <- setdiff(indices, colnames(x))
  if(length(missing_indices) > 0L) {
    stop("Unknown fit index or indices: ", paste(missing_indices, collapse = ", "))
  }

  nclasses <- as.numeric(x[, "nclasses"])
  values <- as.matrix(x[, indices, drop = FALSE])
  storage.mode(values) <- "double"
  finite_values <- values[is.finite(values)]
  if(length(finite_values) == 0L) stop("The selected fit indices contain no finite values")

  if(is.null(colors)) {
    colors <- grDevices::hcl.colors(length(indices), palette = "Dark 3")
  } else {
    colors <- rep(colors, length.out = length(indices))
  }

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  graphics::par(cex = base_size/12)

  y_range <- grDevices::extendrange(finite_values, f = 0.04)
  graphics::plot(nclasses, values[, 1L], type = "n", xaxt = "n",
                 ylim = y_range, main = title, xlab = xlab, ylab = ylab)
  graphics::axis(1, at = sort(unique(nclasses)))
  graphics::grid(nx = NA, ny = NULL, col = "grey90", lty = 1)
  for(i in seq_along(indices)) {
    complete <- is.finite(nclasses) & is.finite(values[, i])
    graphics::lines(nclasses[complete], values[complete, i], col = colors[i],
                    lwd = lwd)
    graphics::points(nclasses[complete], values[complete, i], col = colors[i],
                     pch = pch)
  }
  graphics::legend("topright", legend = indices, col = colors, lty = 1,
                   lwd = lwd, pch = pch, bty = "n",
                   title = "Information criterion")

  result <- data.frame(nclasses = nclasses, values, check.names = FALSE)

  #### Result ####

  return(invisible(result))

}

#' Reshape class-conditional indicator profiles
#'
#' Divides class-conditional indicator profiles into multinomial and Gaussian
#' components and converts them into long data frames suitable for plotting.
#'
#' @param item_output Named list of class-conditional indicator matrices,
#'   normally returned by \code{latInspect(model, what = "item")}.
#' @param item_types Named list containing the character vectors
#'   \code{multinomial} and \code{gaussian}. An optional \code{mvgaussian}
#'   component is combined with \code{gaussian}.
#'
#' @return A list with the following components:
#' \item{multinomial}{A data frame with columns \code{variable},
#'   \code{response}, \code{class}, and \code{probability}.}
#' \item{gaussian}{A data frame with columns \code{variable}, \code{class},
#'   \code{mean}, and \code{sd}.}
#'
#' @export
reshape_lca_mixed <- function(item_output, item_types) {

  if(!is.list(item_output) || is.null(names(item_output))) {
    stop("item_output must be a named list of class-conditional profiles")
  }
  if(!is.list(item_types)) stop("item_types must be a list")

  multinomial_names <- item_types$multinomial
  gaussian_names <- c(item_types$gaussian, item_types$mvgaussian)
  if(is.null(multinomial_names)) multinomial_names <- character(0L)
  if(is.null(gaussian_names)) gaussian_names <- character(0L)

  multinomial_names <- intersect(unique(multinomial_names), names(item_output))
  gaussian_names <- intersect(unique(gaussian_names), names(item_output))

  result <- list(multinomial = reshape_lca_multinomial(item_output[multinomial_names]),
                 gaussian = reshape_lca_continuous(item_output[gaussian_names]))

  #### Result ####

  return(result)

}

#' Reshape multinomial class-conditional probabilities
#'
#' @param object Named list of multinomial class-conditional probability
#'   matrices.
#'
#' @return A long data frame with one row per indicator, response category, and
#'   latent class.
#'
#' @noRd
reshape_lca_multinomial <- function(object) {

  if(length(object) == 0L) {
    result <- data.frame(variable = character(0L), response = character(0L),
                         class = character(0L), probability = numeric(0L),
                         stringsAsFactors = FALSE)
  } else {
    if(is.null(names(object)) || any(names(object) == "")) {
      stop("Multinomial profiles must be supplied in a named list")
    }

    output <- vector("list", length(object))
    for(i in seq_along(object)) {
      variable <- names(object)[i]
      profile <- as.matrix(object[[i]])
      if(nrow(profile) == 0L || ncol(profile) == 0L) {
        stop("The profile for '", variable, "' is empty")
      }

      response_names <- rownames(profile)
      class_names <- colnames(profile)
      if(is.null(response_names)) response_names <- as.character(seq_len(nrow(profile)))
      if(is.null(class_names)) class_names <- paste("Class", seq_len(ncol(profile)))

      output[[i]] <- expand.grid(response = response_names, class = class_names,
                                 stringsAsFactors = FALSE)
      output[[i]]$probability <- as.vector(profile)
      output[[i]]$variable <- variable
      output[[i]] <- output[[i]][, c("variable", "response", "class",
                                    "probability")]
    }

    result <- do.call(rbind, output)
    rownames(result) <- NULL
  }

  #### Result ####

  return(result)

}

#' Reshape Gaussian class-conditional profiles
#'
#' @param object Named list of Gaussian class-conditional profile matrices.
#'
#' @return A long data frame with one row per indicator and latent class,
#'   containing the class-specific mean and standard deviation.
#'
#' @noRd
reshape_lca_continuous <- function(object) {

  if(length(object) == 0L) {
    result <- data.frame(variable = character(0L), class = character(0L),
                         mean = numeric(0L), sd = numeric(0L),
                         stringsAsFactors = FALSE)
  } else {
    if(is.null(names(object)) || any(names(object) == "")) {
      stop("Gaussian profiles must be supplied in a named list")
    }

    output <- vector("list", length(object))
    for(i in seq_along(object)) {
      variable <- names(object)[i]
      profile <- as.matrix(object[[i]])
      if(is.null(rownames(profile)) ||
         !all(c("mean", "stdv") %in% rownames(profile))) {
        stop("The Gaussian profile for '", variable,
             "' must contain rows named 'mean' and 'stdv'")
      }

      class_names <- colnames(profile)
      if(is.null(class_names)) class_names <- paste("Class", seq_len(ncol(profile)))

      output[[i]] <- data.frame(variable = variable, class = class_names,
                                mean = as.numeric(profile["mean", ]),
                                sd = as.numeric(profile["stdv", ]),
                                stringsAsFactors = FALSE)
    }

    result <- do.call(rbind, output)
    rownames(result) <- NULL
  }

  #### Result ####

  return(result)

}

#' Plot Gaussian indicator profiles
#'
#' Displays the class-specific means and standard deviations of Gaussian or
#' multivariate-Gaussian indicators using base R graphics.
#'
#' @param data Data frame produced by \code{reshape_lca_continuous()}.
#' @param variables Optional character vector selecting the indicators to plot.
#' @param title Optional character string with the plot title.
#' @param ylab Character string with the vertical-axis label.
#' @param caption Optional character string with the plot caption.
#' @param nrow,ncol Optional numbers of rows and columns in the panel layout.
#' @param point_size Numeric expansion factor for the points representing means.
#' @param err_width Numeric half-width of the horizontal caps on the standard-
#'   deviation bars.
#' @param base_size Numeric value controlling the overall text size.
#' @param angle Numeric angle of the horizontal-axis labels.
#' @param color Color used for class-profile lines and points.
#' @param ... Additional arguments. Currently unused.
#'
#' @details
#' A separate panel is drawn for every latent class. Error bars represent one
#' class-specific standard deviation below and above each mean.
#'
#' @return The plotted profile data, invisibly.
#'
#' @export
plot_lca_gaussian <- function(data, variables = NULL, title = NULL,
                              ylab = "Mean score (±1 SD)",
                              caption = "Error bars indicate ±1 standard deviation",
                              nrow = NULL, ncol = NULL, point_size = 1.2,
                              err_width = 0.12, base_size = 12, angle = 45,
                              color = "steelblue", ...) {

  required_names <- c("variable", "class", "mean", "sd")
  if(!is.data.frame(data) || !all(required_names %in% names(data))) {
    stop("data must contain the columns: ", paste(required_names, collapse = ", "))
  }
  if(!is.null(variables)) data <- data[data$variable %in% variables, , drop = FALSE]
  if(nrow(data) == 0L) stop("There are no Gaussian indicator profiles to plot")
  if(length(base_size) != 1L || !is.numeric(base_size) || is.na(base_size) ||
     base_size <= 0) {
    stop("base_size must be a single positive number")
  }
  if(length(point_size) != 1L || !is.numeric(point_size) ||
     is.na(point_size) || point_size <= 0) {
    stop("point_size must be a single positive number")
  }
  if(length(err_width) != 1L || !is.numeric(err_width) ||
     is.na(err_width) || err_width < 0) {
    stop("err_width must be a single non-negative number")
  }
  if(length(angle) != 1L || !is.numeric(angle) || is.na(angle)) {
    stop("angle must be a single number")
  }

  variable_levels <- unique(as.character(data$variable))
  class_levels <- unique(as.character(data$class))
  layout_dimensions <- plot_panel_layout(length(class_levels), nrow = nrow,
                                          ncol = ncol)

  lower <- data$mean - data$sd
  upper <- data$mean + data$sd
  finite_values <- c(lower[is.finite(lower)], upper[is.finite(upper)])
  if(length(finite_values) == 0L) stop("The Gaussian profiles contain no finite values")
  y_range <- range(finite_values)
  if(diff(y_range) == 0) y_range <- y_range + c(-0.5, 0.5)

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::layout(matrix(1))
    graphics::par(oldpar)
  }, add = TRUE)

  bottom_margin <- if(angle == 0) 4.5 else 7
  outer_bottom <- if(is.null(caption) || !nzchar(caption)) 0.5 else 2
  outer_top <- if(is.null(title) || !nzchar(title)) 0.5 else 2
  graphics::layout(matrix(seq_len(prod(layout_dimensions)),
                          nrow = layout_dimensions["nrow"], byrow = TRUE))
  graphics::par(mar = c(bottom_margin, 4.2, 2.8, 1.2),
                oma = c(outer_bottom, 0.5, outer_top, 0.5),
                cex = base_size/12)

  for(panel in seq_len(prod(layout_dimensions))) {
    if(panel > length(class_levels)) {
      graphics::plot.new()
      next
    }

    class_name <- class_levels[panel]
    class_data <- data[as.character(data$class) == class_name, , drop = FALSE]
    class_data <- class_data[match(variable_levels, as.character(class_data$variable)),
                             , drop = FALSE]
    x <- seq_along(variable_levels)

    graphics::plot(x, class_data$mean, type = "n", xaxt = "n",
                   xlim = c(0.5, length(x) + 0.5), ylim = y_range,
                   xlab = "", ylab = ylab, main = class_name)
    graphics::abline(h = graphics::axTicks(2), col = "grey92", lty = 1)
    graphics::axis(2)
    graphics::box(col = "grey70")

    complete <- is.finite(class_data$mean) & is.finite(class_data$sd)
    graphics::segments(x[complete], class_data$mean[complete] - class_data$sd[complete],
                       x[complete], class_data$mean[complete] + class_data$sd[complete],
                       col = "grey40", lwd = 1.5)
    graphics::segments(x[complete] - err_width,
                       class_data$mean[complete] - class_data$sd[complete],
                       x[complete] + err_width,
                       class_data$mean[complete] - class_data$sd[complete],
                       col = "grey40", lwd = 1.5)
    graphics::segments(x[complete] - err_width,
                       class_data$mean[complete] + class_data$sd[complete],
                       x[complete] + err_width,
                       class_data$mean[complete] + class_data$sd[complete],
                       col = "grey40", lwd = 1.5)
    if(sum(is.finite(class_data$mean)) > 1L) {
      graphics::lines(x, class_data$mean, col = color, lwd = 2)
    }
    graphics::points(x, class_data$mean, pch = 21, bg = color, col = "white",
                     cex = point_size, lwd = 0.8)
    plot_axis_labels(x, variable_levels, angle = angle,
                     cex = min(0.9, base_size/13))
  }

  if(!is.null(title) && nzchar(title)) {
    graphics::mtext(title, side = 3, outer = TRUE, line = 0.5,
                    font = 2, cex = 1.1)
  }
  if(!is.null(caption) && nzchar(caption)) {
    graphics::mtext(caption, side = 1, outer = TRUE, line = 0.3,
                    cex = 0.8, col = "grey35")
  }

  result <- data

  #### Result ####

  return(invisible(result))

}

#' Draw indicator-specific multinomial legends
#'
#' @param response_levels Named list of response levels by indicator.
#' @param colors Shared response-position colors.
#' @param legend_ncol Number of legend columns.
#' @param base_size Numeric value controlling text size.
#'
#' @return The response levels, invisibly.
#'
#' @noRd
plot_multinomial_legends <- function(response_levels, colors,
                                     legend_ncol = 1L, base_size = 12) {

  variables <- names(response_levels)
  nvariables <- length(variables)
  legend_ncol <- max(1L, min(as.integer(legend_ncol), nvariables))
  legend_nrow <- ceiling(nvariables/legend_ncol)

  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

  cell_width <- 1/legend_ncol
  cell_height <- 1/legend_nrow
  for(i in seq_along(variables)) {
    column <- ceiling(i/legend_nrow)
    row <- i - (column - 1L)*legend_nrow
    x_left <- (column - 1L)*cell_width + 0.04*cell_width
    x_text <- x_left + 0.13*cell_width
    y_top <- 1 - (row - 1L)*cell_height - 0.08*cell_height
    levels_variable <- response_levels[[variables[i]]]
    line_height <- 0.78*cell_height/(length(levels_variable) + 1L)
    cex_legend <- min(0.9, base_size/13,
                      max(0.5, 6/(legend_nrow*(length(levels_variable) + 1L))))

    graphics::text(x_left, y_top, labels = variables[i], adj = c(0, 1),
                   font = 2, cex = cex_legend, xpd = NA)
    for(j in seq_along(levels_variable)) {
      y <- y_top - j*line_height
      box_width <- 0.075*cell_width
      box_height <- min(0.022, 0.28*line_height)
      graphics::rect(x_left, y - box_height, x_left + box_width,
                     y + box_height, col = colors[j], border = "grey50",
                     xpd = NA)
      graphics::text(x_text, y, labels = levels_variable[j], adj = c(0, 0.5),
                     cex = cex_legend, xpd = NA)
    }
  }

  #### Result ####

  return(invisible(response_levels))

}

#' Plot multinomial response probabilities
#'
#' Displays class-specific response probabilities for multinomial indicators as
#' stacked bar charts using base R graphics.
#'
#' @param data Data frame produced by \code{reshape_lca_multinomial()}.
#' @param variables Optional character vector selecting the indicators to plot.
#' @param bars Character string naming the variable mapped to the horizontal
#'   axis. The default is \code{"variable"}.
#' @param facet Character string naming the variable defining the plot panels.
#'   The default is \code{"class"}.
#' @param bw Logical value indicating whether a greyscale palette should be used.
#' @param title Optional character string with the plot title.
#' @param ylab Character string with the vertical-axis label.
#' @param base_size Numeric value controlling the overall text size.
#' @param angle Numeric angle of the horizontal-axis labels.
#' @param nrow,ncol Optional numbers of rows and columns in the panel layout.
#' @param legend_ncol Optional number of columns used in the legend area. When
#'   \code{NULL}, it is selected from the number of indicators.
#' @param legend_width Relative width assigned to the legend area.
#' @param border Color of the borders separating stacked response categories.
#' @param ... Additional arguments. Currently unused.
#'
#' @details
#' A separate legend is drawn for every indicator. Each legend displays the
#' actual response levels of that indicator. Colors are assigned according to
#' response-category position, so the first category uses the same color for
#' every indicator, the second category uses the same color for every indicator,
#' and so forth. No package other than base R is required.
#'
#' @return The plotted probability data, invisibly.
#'
#' @export
plot_lca_multinomial <- function(data, variables = NULL, bars = "variable",
                                 facet = "class", bw = FALSE, title = NULL,
                                 ylab = "Probability", base_size = 12,
                                 angle = 45, nrow = NULL, ncol = NULL,
                                 legend_ncol = NULL, legend_width = 0.8,
                                 border = "white", ...) {

  required_names <- c("variable", "response", "class", "probability")
  if(!is.data.frame(data) || !all(required_names %in% names(data))) {
    stop("data must contain the columns: ", paste(required_names, collapse = ", "))
  }
  if(!is.character(bars) || length(bars) != 1L || !bars %in% names(data)) {
    stop("bars must name one column in data")
  }
  if(!is.character(facet) || length(facet) != 1L || !facet %in% names(data)) {
    stop("facet must name one column in data")
  }
  if(!is.logical(bw) || length(bw) != 1L || is.na(bw)) {
    stop("bw must be TRUE or FALSE")
  }
  if(!is.numeric(data$probability)) stop("probability must be numeric")
  if(!is.null(variables)) data <- data[data$variable %in% variables, , drop = FALSE]
  if(nrow(data) == 0L) stop("There are no multinomial indicator profiles to plot")
  if(any(!is.finite(data$probability))) {
    stop("The multinomial probabilities must be finite")
  }
  if(length(base_size) != 1L || !is.numeric(base_size) || is.na(base_size) ||
     base_size <= 0) {
    stop("base_size must be a single positive number")
  }
  if(length(angle) != 1L || !is.numeric(angle) || is.na(angle)) {
    stop("angle must be a single number")
  }
  if(!is.null(legend_ncol)) {
    if(length(legend_ncol) != 1L || !is.numeric(legend_ncol) ||
       is.na(legend_ncol) || legend_ncol < 1L) {
      stop("legend_ncol must be a positive integer or NULL")
    }
    legend_ncol <- as.integer(legend_ncol)
  }

  plot_variables <- unique(as.character(data$variable))
  response_levels <- lapply(plot_variables, FUN = function(variable) {

    result <- plot_response_levels(
      data$response[as.character(data$variable) == variable])

    #### Result ####

    return(result)

  })
  names(response_levels) <- plot_variables

  max_responses <- max(lengths(response_levels))
  if(bw) {
    colors <- grDevices::gray.colors(max_responses, start = 0.25, end = 0.85)
  } else {
    colors <- grDevices::hcl.colors(max_responses, palette = "viridis")
  }

  facet_levels <- unique(as.character(data[[facet]]))
  bar_levels <- unique(as.character(data[[bars]]))
  bar_variables <- character(length(bar_levels))
  for(i in seq_along(bar_levels)) {
    rows <- as.character(data[[bars]]) == bar_levels[i]
    variables_bar <- unique(as.character(data$variable[rows]))
    if(length(variables_bar) != 1L) {
      stop("Each value of bars must correspond to exactly one indicator")
    }
    bar_variables[i] <- variables_bar
  }

  layout_dimensions <- plot_panel_layout(length(facet_levels), nrow = nrow,
                                          ncol = ncol)
  if(is.null(legend_ncol)) {
    legend_ncol <- max(1L, min(3L, ceiling(length(plot_variables)/6L)))
  }
  if(length(legend_width) != 1L || !is.numeric(legend_width) ||
     is.na(legend_width) || legend_width <= 0) {
    stop("legend_width must be a single positive number")
  }

  panel_count <- prod(layout_dimensions)
  legend_panel <- panel_count + 1L
  panel_matrix <- matrix(seq_len(panel_count), nrow = layout_dimensions["nrow"],
                         byrow = TRUE)
  layout_matrix <- cbind(panel_matrix,
                         rep(legend_panel, layout_dimensions["nrow"]))

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::layout(matrix(1))
    graphics::par(oldpar)
  }, add = TRUE)

  bottom_margin <- if(angle == 0) 4.5 else 7
  outer_top <- if(is.null(title) || !nzchar(title)) 0.5 else 2
  graphics::layout(layout_matrix,
                   widths = c(rep(1, layout_dimensions["ncol"]), legend_width))
  graphics::par(mar = c(bottom_margin, 4.2, 2.8, 1.2),
                oma = c(0.5, 0.5, outer_top, 0.5), cex = base_size/12)

  for(panel in seq_len(panel_count)) {
    if(panel > length(facet_levels)) {
      graphics::plot.new()
      next
    }

    facet_value <- facet_levels[panel]
    x <- seq_along(bar_levels)
    graphics::plot(x, rep(0, length(x)), type = "n", xaxt = "n",
                   xlim = c(0.4, length(x) + 0.6), ylim = c(0, 1),
                   xlab = "", ylab = ylab, main = facet_value, yaxs = "i")
    graphics::abline(h = seq(0, 1, by = 0.25), col = "grey92", lty = 1)
    graphics::axis(2, at = seq(0, 1, by = 0.25), las = 1)
    graphics::box(col = "grey70")

    for(i in seq_along(bar_levels)) {
      variable <- bar_variables[i]
      levels_variable <- response_levels[[variable]]
      bottom <- 0
      for(j in seq_along(levels_variable)) {
        rows <- as.character(data[[facet]]) == facet_value &
          as.character(data[[bars]]) == bar_levels[i] &
          as.character(data$response) == levels_variable[j]
        probability <- sum(data$probability[rows], na.rm = TRUE)
        top <- bottom + probability
        graphics::rect(x[i] - 0.36, bottom, x[i] + 0.36, top,
                       col = colors[j], border = border)
        bottom <- top
      }
    }

    plot_axis_labels(x, bar_levels, angle = angle,
                     cex = min(0.9, base_size/13))
  }

  graphics::par(mar = c(0.5, 0.5, 0.5, 0.5))
  plot_multinomial_legends(response_levels, colors,
                           legend_ncol = legend_ncol, base_size = base_size)

  if(!is.null(title) && nzchar(title)) {
    graphics::mtext(title, side = 3, outer = TRUE, line = 0.5,
                    font = 2, cex = 1.1)
  }

  result <- data

  #### Result ####

  return(invisible(result))

}

#' Plot mixed indicator profiles
#'
#' Draws the appropriate base R profile plot for each available indicator
#' family.
#'
#' @param mixed_data List produced by \code{reshape_lca_mixed()}.
#' @param variables Optional character vector selecting indicators to plot.
#' @param ... Additional arguments passed to \code{plot_lca_gaussian()} and
#'   \code{plot_lca_multinomial()}.
#'
#' @details
#' If only one indicator family is available, one plot is drawn. For mixed
#' models, the Gaussian and multinomial plots are drawn sequentially as two
#' pages on the active graphics device.
#'
#' @return An object of class \code{"llca_plots"}, invisibly. The object stores
#'   the data and plotting arguments and can be printed to redraw the plots.
#'
#' @export
plot_lca_mixed <- function(mixed_data, variables = NULL, ...) {

  if(!is.list(mixed_data) ||
     !all(c("gaussian", "multinomial") %in% names(mixed_data))) {
    stop("mixed_data must be the result of reshape_lca_mixed()")
  }
  if(!is.data.frame(mixed_data$gaussian) ||
     !is.data.frame(mixed_data$multinomial)) {
    stop("The gaussian and multinomial components must be data frames")
  }

  dots <- list(...)
  result <- list()
  if(nrow(mixed_data$gaussian) > 0L) {
    result$gaussian <- list(data = mixed_data$gaussian,
                            variables = variables, arguments = dots)
  }
  if(nrow(mixed_data$multinomial) > 0L) {
    result$multinomial <- list(data = mixed_data$multinomial,
                               variables = variables, arguments = dots)
  }
  if(length(result) == 0L) stop("There are no indicator profiles to plot")

  class(result) <- c("llca_plots", "list")
  print.llca_plots(result)

  #### Result ####

  return(invisible(result))

}

#' Print mixed latent class profile plots
#'
#' Redraws the base R plots stored in a mixed latent class profile object.
#'
#' @param x An object of class \code{"llca_plots"}.
#' @param ... Additional arguments. Currently unused.
#'
#' @return The input object, invisibly.
#'
#' @method print llca_plots
#' @export
print.llca_plots <- function(x, ...) {

  if(!is.null(x$gaussian)) {
    arguments <- c(list(data = x$gaussian$data,
                        variables = x$gaussian$variables),
                   x$gaussian$arguments)
    do.call(plot_lca_gaussian, arguments)
  }
  if(!is.null(x$multinomial)) {
    arguments <- c(list(data = x$multinomial$data,
                        variables = x$multinomial$variables),
                   x$multinomial$arguments)
    do.call(plot_lca_multinomial, arguments)
  }

  #### Result ####

  return(invisible(x))

}

#' Plot regression coefficients from an LCA model
#'
#' Computes confidence intervals for the latent-class regression coefficients
#' and displays them using \code{forestplot()}.
#'
#' @param fit A fitted object of class \code{"llca"}.
#' @param se_type Character string specifying the type of standard errors passed
#'   to \code{se()}.
#' @param what Character string specifying whether to plot odds ratios
#'   (\code{"OR"}) or log odds ratios (\code{"log"}).
#' @param effects Character string specifying whether coefficients are displayed
#'   using effects coding (\code{"coding"}) or dummy coding (\code{"dummy"}).
#' @param confidence Numeric confidence level between zero and one.
#' @param predictors Optional character vector selecting predictors. A supplied
#'   value selects every coefficient name containing that value. When
#'   \code{NULL}, all available predictors are included.
#' @param intercept Logical value indicating whether intercept coefficients are
#'   included.
#' @param show_est_ci Logical value indicating whether coefficient estimates and
#'   confidence intervals are appended to the vertical-axis labels.
#' @param est_ci_header_cex Numeric expansion factor for the estimate-and-
#'   confidence-interval header.
#' @param cex_y Numeric expansion factor for coefficient labels.
#' @param mfrow Numeric vector of length two passed to \code{par(mfrow = ...)}.
#' @param ... Additional arguments passed to \code{forestplot()}.
#'
#' @details
#' The confidence intervals use the same asymptotic-normal calculation as the
#' former \code{plot_coeffs()} function. For \code{what = "OR"}, estimates and
#' interval limits are exponentiated before plotting.
#'
#' @return A list containing the plotted standard errors, variance-covariance
#'   matrix, confidence limits, and coefficient estimates, invisibly.
#'
#' @noRd
plot_lca_coefficients <- function(fit, se_type = "standard",
                                  what = c("OR", "log"),
                                  effects = c("coding", "dummy"),
                                  confidence = 0.95, predictors = NULL,
                                  intercept = TRUE, show_est_ci = TRUE,
                                  est_ci_header_cex = 0.5, cex_y = 0.5,
                                  mfrow = c(1, 1), ...) {

  if(!inherits(fit, "llca")) stop("fit must inherit from class 'llca'")
  if(length(fit@transformed_pars) == 0L) {
    stop("fit must be a fitted llca object")
  }
  if(is.null(fit@transformed_pars$beta)) {
    stop("The fitted model does not contain latent-class regression coefficients")
  }
  if(!is.character(se_type) || length(se_type) != 1L || is.na(se_type)) {
    stop("se_type must be a single character string")
  }
  if(!is.numeric(confidence) || length(confidence) != 1L ||
     is.na(confidence) || confidence <= 0 || confidence >= 1) {
    stop("confidence must be a single number between zero and one")
  }
  if(!is.null(predictors) &&
     (!is.character(predictors) || anyNA(predictors))) {
    stop("predictors must be a character vector or NULL")
  }
  if(!is.logical(intercept) || length(intercept) != 1L || is.na(intercept)) {
    stop("intercept must be TRUE or FALSE")
  }
  if(!is.logical(show_est_ci) || length(show_est_ci) != 1L ||
     is.na(show_est_ci)) {
    stop("show_est_ci must be TRUE or FALSE")
  }
  if(!is.numeric(mfrow) || length(mfrow) != 2L || anyNA(mfrow) ||
     any(mfrow < 1L)) {
    stop("mfrow must contain two positive integers")
  }

  what <- match.arg(tolower(what), choices = c("or", "log"))
  effects <- match.arg(effects)

  SE <- se(fit = fit, type = se_type, digits = 9)

  beta_labels <- as.vector(fit@modelInfo$trans$beta)

  #### log(odds ratio) ####

  if(effects == "coding") {
    EF <- effects_coding(fit@transformed_pars$beta,
                         as.matrix(SE$vcov[beta_labels, beta_labels, drop = FALSE]))
    betas <- as.matrix(EF$beta)
    vcov <- as.matrix(EF$vcov)
    colnames(betas) <- colnames(fit@transformed_pars$beta)
  } else {
    if(ncol(fit@transformed_pars$beta) < 2L) {
      stop("Dummy-coded coefficients require at least two latent classes")
    }
    betas <- fit@transformed_pars$beta[, -1L, drop = FALSE]
    vcov <- SE$vcov[select_betas, select_betas, drop = FALSE]
  }

  if(is.null(rownames(betas))) {
    rownames(betas) <- paste0("Predictor", seq_len(nrow(betas)))
  }
  if(is.null(colnames(betas))) {
    colnames(betas) <- paste0("Class", seq_len(ncol(betas)))
  }
  if(nrow(vcov) != length(betas) || ncol(vcov) != length(betas)) {
    stop("The coefficient matrix and its variance-covariance matrix have incompatible dimensions")
  }

  allnames <- paste(rep(rownames(betas), times = ncol(betas)),
                    rep(colnames(betas), each = nrow(betas)), sep = "|")
  parameters <- c(betas)
  names(parameters) <- allnames
  rownames(vcov) <- colnames(vcov) <- allnames

  group <- factor(rep(colnames(betas), each = nrow(betas)),
                  levels = colnames(betas))
  coef_names <- rep(rownames(betas), times = ncol(betas))

  standard_errors <- sqrt(diag(vcov))
  names(standard_errors) <- allnames
  critical_value <- sqrt(stats::qchisq(confidence, df = 1))
  lower <- parameters - critical_value*standard_errors
  upper <- parameters + critical_value*standard_errors

  confidence_label <- paste0(format(100*confidence, trim = TRUE,
                                    scientific = FALSE), "% CI")
  if(what == "log") {
    xlab <- paste0("log(odds ratio) (", confidence_label, ")")
    est_ci_header <- xlab
    refline <- 0
  } else {
    xlab <- paste0("odds ratio (", confidence_label, ")")
    est_ci_header <- xlab
    refline <- 1

    parameters <- exp(parameters)
    J <- diag(parameters)
    rownames(J) <- colnames(J) <- allnames
    vcov <- J %*% vcov %*% t(J)
    standard_errors <- sqrt(diag(vcov))
    names(standard_errors) <- allnames
    lower <- exp(lower)
    upper <- exp(upper)
  }

  available_predictors <- unique(rownames(betas))
  if(is.null(predictors)) predictors <- available_predictors
  if(intercept) {
    predictors <- unique(c("(Intercept)", predictors))
  } else {
    predictors <- setdiff(predictors, "(Intercept)")
  }

  if(length(predictors) == 0L) {
    stop("No predictors were selected")
  }
  selected <- Reduce(`|`, lapply(predictors, FUN = function(predictor) {

    result <- grepl(predictor, coef_names, fixed = TRUE)

    #### Result ####

    return(result)

  }))
  selection <- which(selected)
  if(length(selection) == 0L) {
    stop("None of the requested predictors were found. Available predictors are: ",
         paste(available_predictors, collapse = ", "))
  }

  result <- list(se = standard_errors[selection],
                 vcov = vcov[selection, selection, drop = FALSE],
                 upper = upper[selection], lower = lower[selection],
                 parameters = parameters[selection])

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  graphics::par(mfrow = as.integer(mfrow))
  forestplot(parm = parameters[selection], lower = lower[selection],
             upper = upper[selection], group = group[selection],
             labels = coef_names[selection], refline = refline, xlab = xlab,
             show_est_ci = show_est_ci, est_ci_header = est_ci_header,
             est_ci_header_cex = est_ci_header_cex, cex_y = cex_y, ...)

  #### Result ####

  return(invisible(result))

}

#' Plot results from an LCA model
#'
#' Displays class-conditional measurement profiles or latent-class regression
#' coefficients from a fitted latent class model using base R graphics.
#'
#' @param x A fitted object of class \code{"llca"}.
#' @param type Character string specifying the result to plot. Available values
#'   are \code{"all"}, \code{"gaussian"}, \code{"multinomial"}, and
#'   \code{"coefficients"}. Partial matching is supported.
#' @param variables Optional character vector selecting measurement indicators
#'   for profile plots. This argument is ignored when
#'   \code{type = "coefficients"}; use \code{predictors} through \code{...}
#'   instead.
#' @param ... Additional arguments passed to \code{plot_lca_mixed()} for
#'   profile plots or to the coefficient-plotting routine for
#'   \code{type = "coefficients"}. Coefficient options include
#'   \code{se_type}, \code{what}, \code{effects}, \code{confidence},
#'   \code{predictors}, \code{intercept}, and the graphical arguments of
#'   \code{forestplot()}.
#'
#' @details
#' Gaussian and multivariate-Gaussian indicators are displayed using their
#' class-specific means and standard deviations. Multinomial indicators are
#' displayed using their class-specific response probabilities. Response colors
#' are reused across indicators according to category position, while each
#' indicator has its own legend containing its actual response levels. Outcomes
#' and covariates are not included in the measurement-profile plots.
#'
#' For mixed models with \code{type = "all"}, the Gaussian and multinomial
#' plots are drawn sequentially as separate pages on the active graphics device.
#' With \code{type = "coefficients"}, the latent-class regression coefficients
#' are displayed in forest-plot form. Use \code{what = "OR"} for odds ratios or
#' \code{what = "log"} for log odds ratios.
#'
#' @return For profile plots, an object of class \code{"llca_plots"},
#'   invisibly. For coefficient plots, a list containing the displayed
#'   coefficient estimates, standard errors, variance-covariance matrix, and
#'   confidence limits, invisibly.
#'
#' @examples
#' \dontrun{
#' fit <- lca(data = gss82, nclasses = 3,
#'            multinomial = c("PURPOSE", "ACCURACY"))
#' plot(fit)
#' plot(fit, type = "multinomial", bw = TRUE)
#'
#' structural_fit <- lca(data = gss82, nclasses = 3,
#'                       multinomial = c("PURPOSE", "ACCURACY"),
#'                       covariates = c("RACE", "SEX"))[["structural"]]
#' plot(structural_fit, type = "coefficients", what = "OR",
#'      effects = "coding", predictors = "RACE")
#' }
#'
#' @method plot llca
#' @export
plot.llca <- function(x,
                      type = c("all", "gaussian", "multinomial",
                               "coefficients"),
                      variables = NULL, ...) {

  if(!inherits(x, "llca")) stop("x must inherit from class 'llca'")
  if(length(x@transformed_pars) == 0L) stop("x must be a fitted llca object")

  type <- match.arg(type)
  if(type == "coefficients") {
    if(!is.null(variables)) {
      warning("variables is ignored for coefficient plots; use predictors instead")
    }
    result <- plot_lca_coefficients(fit = x, ...)
  } else {
    item_output <- latInspect(x, what = "item")
    item_types <- list(multinomial = x@dataList$multinomial$multinomial_names,
                       gaussian = c(x@dataList$gaussian$gaussian_names,
                                    x@dataList$mvgaussian$mvgaussian_names))

    item_profiles <- reshape_lca_mixed(item_output, item_types)
    if(type == "gaussian") {
      item_profiles$multinomial <- item_profiles$multinomial[0, ]
    }
    if(type == "multinomial") {
      item_profiles$gaussian <- item_profiles$gaussian[0, ]
    }

    result <- plot_lca_mixed(item_profiles, variables = variables, ...)
  }

  #### Result ####

  return(invisible(result))

}

#' Plot results from a collection of LCA models
#'
#' Selects the appropriate fitted model from an \code{llcalist} object and
#' displays measurement profiles or latent-class regression coefficients using
#' the corresponding \code{plot.llca()} method.
#'
#' @param x An object of class \code{"llcalist"}.
#' @param type Character string specifying the result to plot. Available values
#'   are \code{"all"}, \code{"gaussian"}, \code{"multinomial"}, and
#'   \code{"coefficients"}.
#' @param variables Optional character vector selecting measurement indicators
#'   for profile plots. This argument is passed to \code{plot.llca()}.
#' @param ... Additional arguments passed to \code{plot.llca()}.
#'
#' @details
#' For a structural-after-measurement result containing components named
#' \code{measurement} and \code{structural}, profile plots use the measurement
#' model and coefficient plots use the structural model. This allows the full
#' object returned by \code{lca()} to be plotted directly.
#'
#' For other \code{llcalist} objects, a single model is plotted directly. When
#' several models are present, profile plots are produced sequentially for all
#' models. Coefficient plots are produced only for models that include
#' covariates.
#'
#' @return The result returned by \code{plot.llca()}, invisibly. When several
#'   models are plotted, a named list of results is returned invisibly.
#'
#' @examples
#' \dontrun{
#' fit <- lca(data = empathy, nclasses = 4,
#'            gaussian = c("ec1", "ec2", "ec3", "ec4", "ec5", "ec6"),
#'            covariates = c("pt1", "pt2", "pt3", "pt4"))
#' plot(fit)
#' plot(fit, type = "coefficients", what = "OR")
#' }
#'
#' @method plot llcalist
#' @export
plot.llcalist <- function(x,
                          type = c("all", "gaussian", "multinomial",
                                   "coefficients"),
                          variables = NULL, ...) {

  if(!inherits(x, "llcalist")) stop("x must inherit from class 'llcalist'")
  if(length(x) == 0L) stop("x must contain at least one fitted llca object")

  valid_models <- vapply(x, FUN = inherits, FUN.VALUE = logical(1L),
                         what = "llca")
  if(!all(valid_models)) {
    stop("Every element of x must inherit from class 'llca'")
  }

  type <- match.arg(type)
  model_names <- names(x)
  if(is.null(model_names)) model_names <- rep("", length(x))

  if("measurement" %in% model_names || "structural" %in% model_names) {

    if(type == "coefficients") {
      if(!"structural" %in% model_names) {
        stop("The llcalist object does not contain a structural model")
      }
      result <- plot.llca(x[["structural"]], type = type,
                          variables = variables, ...)
    } else {
      if("measurement" %in% model_names) {
        selected_model <- x[["measurement"]]
      } else {
        selected_model <- x[["structural"]]
      }
      result <- plot.llca(selected_model, type = type,
                          variables = variables, ...)
    }

  } else if(length(x) == 1L) {

    result <- plot.llca(x[[1L]], type = type, variables = variables, ...)

  } else {

    if(type == "coefficients") {
      has_covariates <- vapply(x, FUN = function(model) {

        covariates_names <- model@dataList$covariates_names
        result <- !is.null(covariates_names) && length(covariates_names) > 0L

        #### Result ####

        return(result)

      }, FUN.VALUE = logical(1L))
      selected <- which(has_covariates)
      if(length(selected) == 0L) {
        stop("None of the fitted models contains covariates")
      }
    } else {
      selected <- seq_along(x)
    }

    result <- vector("list", length(selected))
    result_names <- model_names[selected]
    empty_names <- is.na(result_names) | result_names == ""
    if(any(empty_names)) {
      result_names[empty_names] <- vapply(x[selected[empty_names]],
                                         FUN = function(model) {

        nclasses <- ncol(model@transformed_pars$class)
        result <- paste0("nclasses=", nclasses)

        #### Result ####

        return(result)

      }, FUN.VALUE = character(1L))
    }
    names(result) <- make.unique(result_names)

    for(i in seq_along(selected)) {
      result[[i]] <- plot.llca(x[[selected[i]]], type = type,
                               variables = variables, ...)
    }

  }

  #### Result ####

  return(invisible(result))

}
