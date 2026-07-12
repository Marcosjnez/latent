# Author: Marcos Jimenez (using ChatGPT 5.2)
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 12/07/2026
#'
#' Forest plot for grouped parameter estimates
#'
#' Draws one base R forest plot for each group of parameter estimates. Each
#' panel displays point estimates, confidence intervals, and a reference line.
#'
#' @param parm Numeric vector containing parameter estimates.
#' @param lower,upper Numeric vectors containing the lower and upper confidence
#'   limits. They must have the same length as \code{parm}.
#' @param group Optional vector defining groups of parameters. One plot is drawn
#'   for each group. When \code{NULL}, all parameters are placed in one group.
#' @param labels Optional character vector containing coefficient labels. When
#'   \code{NULL}, names from \code{parm} are used when available.
#' @param refline Numeric value indicating the vertical reference line.
#' @param xlab Character string with the horizontal-axis label.
#' @param xlim Optional numeric vector of length two defining the horizontal-axis
#'   limits. When \code{NULL}, limits are calculated from the confidence
#'   intervals and reference line.
#' @param order_within_group Character string indicating whether coefficients
#'   retain their supplied order (\code{"as_is"}) or are ordered by their
#'   estimate within each group (\code{"by_parm"}).
#' @param point_pch Plotting symbol used for parameter estimates.
#' @param point_cex Numeric expansion factor for parameter-estimate symbols.
#' @param ci_lwd Width of confidence-interval segments.
#' @param ref_lty Line type used for the reference line.
#' @param mar Numeric vector passed to \code{par(mar = ...)}.
#' @param show_est_ci Logical value indicating whether estimates and confidence
#'   intervals are appended to coefficient labels.
#' @param digits Number of decimal places used when formatting estimates and
#'   confidence intervals.
#' @param cex_y,cex_x Numeric expansion factors for vertical- and horizontal-axis
#'   labels.
#' @param cex_lab Numeric expansion factor for the horizontal-axis title.
#' @param cex_main Numeric expansion factor for group titles.
#' @param est_ci_header Optional character string displayed above the coefficient
#'   labels when \code{show_est_ci = TRUE}. By default, \code{xlab} is used.
#' @param est_ci_header_cex Optional expansion factor for the estimate-and-
#'   confidence-interval header.
#' @param ... Additional graphical arguments passed to \code{plot.default()}.
#'
#' @details
#' Invalid rows containing non-finite estimates or confidence limits, or missing
#' groups, are removed with a warning. All groups use the same horizontal-axis
#' limits so their estimates can be compared directly.
#'
#' @return A list containing the horizontal-axis limits and plotted group names,
#'   invisibly.
#'
#' @examples
#' \dontrun{
#' forestplot(parm = c(0.2, -0.1, 0.5),
#'            lower = c(0.05, -0.3, 0.2),
#'            upper = c(0.35, 0.1, 0.8),
#'            group = c("Class 1", "Class 1", "Class 2"),
#'            labels = c("Age", "Sex", "Age"))
#' }
#'
#' @export
forestplot <- function(parm, lower, upper, group = NULL, labels = NULL,
                       refline = 0, xlab = "Estimate", xlim = NULL,
                       order_within_group = c("as_is", "by_parm"),
                       point_pch = 16, point_cex = 1, ci_lwd = 2,
                       ref_lty = 2, mar = c(5, 9, 4, 2) + 0.1,
                       show_est_ci = FALSE, digits = 2, cex_y = 1,
                       cex_x = 1, cex_lab = 1, cex_main = 1,
                       est_ci_header = NULL, est_ci_header_cex = NULL, ...) {

  if(!is.numeric(parm) || !is.numeric(lower) || !is.numeric(upper)) {
    stop("parm, lower, and upper must be numeric vectors")
  }
  if(is.null(group)) group <- rep("All", length(parm))
  if(length(parm) != length(lower) || length(parm) != length(upper) ||
     length(parm) != length(group)) {
    stop("parm, lower, upper, and group must have the same length")
  }
  if(!is.null(labels) && length(labels) != length(parm)) {
    stop("labels must have the same length as parm or be NULL")
  }
  if(!is.numeric(refline) || length(refline) != 1L || is.na(refline)) {
    stop("refline must be a single numeric value")
  }
  if(!is.null(xlim) &&
     (!is.numeric(xlim) || length(xlim) != 2L || anyNA(xlim) ||
      !all(is.finite(xlim)) || xlim[1L] >= xlim[2L])) {
    stop("xlim must be an increasing numeric vector of length two or NULL")
  }
  if(!is.logical(show_est_ci) || length(show_est_ci) != 1L ||
     is.na(show_est_ci)) {
    stop("show_est_ci must be TRUE or FALSE")
  }
  if(!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
     digits < 0) {
    stop("digits must be a single non-negative number")
  }

  ok <- is.finite(parm) & is.finite(lower) & is.finite(upper) & !is.na(group)
  if(!all(ok)) {
    warning(sprintf("Dropping %d row(s) with non-finite estimates or missing groups",
                    sum(!ok)))
    parm <- parm[ok]
    lower <- lower[ok]
    upper <- upper[ok]
    group <- group[ok]
    if(!is.null(labels)) labels <- labels[ok]
  }
  if(length(parm) == 0L) {
    stop("No valid rows remain after removing incomplete estimates")
  }

  if(is.factor(group)) {
    group <- droplevels(group)
  } else {
    group <- factor(group, levels = unique(group))
  }
  groups <- levels(group)
  if(length(groups) == 0L) stop("No groups remain to plot")

  order_within_group <- match.arg(order_within_group)
  if(is.null(xlim)) {
    plot_range <- range(c(lower, upper, refline), finite = TRUE)
    padding <- 0.04*diff(plot_range)
    if(!is.finite(padding) || padding == 0) padding <- 1
    xlim <- c(plot_range[1L] - padding, plot_range[2L] + padding)
  }

  format_values <- function(x) {

    result <- formatC(x, digits = as.integer(digits), format = "f")

    #### Result ####

    return(result)

  }

  get_coefficient_labels <- function(indices) {

    if(!is.null(labels)) {
      result <- as.character(labels[indices])
    } else if(!is.null(names(parm))) {
      result <- as.character(names(parm)[indices])
    } else {
      result <- paste0("p", indices)
    }

    #### Result ####

    return(result)

  }

  if(is.null(est_ci_header_cex)) est_ci_header_cex <- cex_y

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)

  for(group_name in groups) {
    indices <- which(group == group_name)
    coefficient_labels <- get_coefficient_labels(indices)

    if(order_within_group == "by_parm") {
      ordering <- order(parm[indices], decreasing = FALSE)
      indices <- indices[ordering]
      coefficient_labels <- coefficient_labels[ordering]
    }

    if(show_est_ci) {
      coefficient_labels <- sprintf("%s %s (%s, %s)", coefficient_labels,
                                    format_values(parm[indices]),
                                    format_values(lower[indices]),
                                    format_values(upper[indices]))
    }

    ncoefficients <- length(indices)
    if(show_est_ci) {
      y_header <- ncoefficients + 1
      y_coefficients <- ncoefficients:1
      ylim <- c(0.5, ncoefficients + 1.5)
    } else {
      y_header <- NA_real_
      y_coefficients <- ncoefficients:1
      ylim <- c(0.5, ncoefficients + 0.5)
    }

    graphics::par(mar = mar, cex.axis = cex_x, cex.lab = cex_lab)
    graphics::plot(NA, xlim = xlim, ylim = ylim, yaxt = "n", ylab = "",
                   xlab = xlab, main = "", ...)
    graphics::title(main = as.character(group_name), cex.main = cex_main)
    graphics::abline(v = refline, lty = ref_lty, col = "gray50")
    graphics::segments(lower[indices], y_coefficients, upper[indices],
                       y_coefficients, lwd = ci_lwd)
    graphics::points(parm[indices], y_coefficients, pch = point_pch,
                     cex = point_cex)
    graphics::axis(2, at = y_coefficients, labels = coefficient_labels,
                   las = 1, tick = FALSE, cex.axis = cex_y)

    if(show_est_ci) {
      header <- if(is.null(est_ci_header)) xlab else est_ci_header
      xpd <- graphics::par(xpd = NA)
      graphics::mtext(header, side = 2, at = y_header,
                      line = graphics::par("mgp")[2L], las = 1, adj = 1,
                      cex = est_ci_header_cex, font = 2)
      graphics::par(xpd)
    }

    graphics::box()
  }

  result <- list(xlim = xlim, groups = groups)

  #### Result ####

  return(invisible(result))

}
