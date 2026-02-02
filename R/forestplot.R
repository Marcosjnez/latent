# Author: Marcos Jimenez (using ChatGPT 5.2)
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 29/01/2026
#'
#' @export
forestplot <- function(parm, lower, upper,
                       group = NULL,
                       labels = NULL,
                       refline = 0,
                       xlab = "Estimate",
                       xlim = NULL,
                       order_within_group = c("as_is", "by_parm"),
                       point_pch = 16,
                       point_cex = 1,
                       ci_lwd = 2,
                       ref_lty = 2,
                       mar = c(5, 9, 4, 2) + 0.1,
                       show_est_ci = FALSE,
                       digits = 2,
                       cex_y = 1,
                       cex_x = 1,
                       cex_lab = 1,
                       cex_main = 1,
                       est_ci_header = NULL,
                       est_ci_header_cex = NULL,
                       ...) {

  # ---- checks ----
  if (missing(group) || is.null(group)) group <- rep("All", length(parm))
  if (length(parm) != length(lower) || length(parm) != length(upper) || length(parm) != length(group)) {
    stop("parm, lower, upper, and group must have the same length.")
  }

  ok <- is.finite(parm) & is.finite(lower) & is.finite(upper) & !is.na(group)
  if (!all(ok)) {
    warning(sprintf("Dropping %d row(s) with non-finite parm/lower/upper or NA group.", sum(!ok)))
    parm <- parm[ok]; lower <- lower[ok]; upper <- upper[ok]; group <- group[ok]
    if (!is.null(labels)) labels <- labels[ok]
  }
  if (length(parm) == 0L) stop("No valid rows remain after filtering finite parm/lower/upper and non-NA group.")

  group <- as.factor(group)
  levs <- levels(group)
  if (length(levs) == 0L) stop("No groups to plot (group has zero levels after filtering).")

  order_within_group <- match.arg(order_within_group)

  # shared xlim
  if (is.null(xlim)) {
    rng <- range(c(lower, upper, refline), finite = TRUE)
    pad <- 0.04 * diff(rng)
    if (!is.finite(pad) || pad == 0) pad <- 1
    xlim <- c(rng[1] - pad, rng[2] + pad)
  }

  fmt <- function(x) formatC(x, digits = digits, format = "f")

  get_coef_labels <- function(idx) {
    if (!is.null(labels)) return(as.character(labels[idx]))
    if (!is.null(names(parm))) return(as.character(names(parm)[idx]))
    paste0("p", idx)
  }

  if (is.null(est_ci_header_cex)) est_ci_header_cex <- cex_y

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  # ---- one forest plot per group, sequentially ----
  for (g in levs) {
    idx <- which(group == g)

    lab_g <- get_coef_labels(idx)

    if (order_within_group == "by_parm") {
      o <- order(parm[idx], decreasing = FALSE)
      idx <- idx[o]
      lab_g <- lab_g[o]
    }

    # Build y-label strings (optionally include estimates/CI)
    if (isTRUE(show_est_ci)) {
      lab_g <- sprintf("%s %s (%s, %s)",
                       lab_g, fmt(parm[idx]), fmt(lower[idx]), fmt(upper[idx]))
    }

    n <- length(idx)

    # y positions: reserve a first row for the header when show_est_ci=TRUE
    if (isTRUE(show_est_ci)) {
      y_header <- n + 1
      y_coef   <- n:1
      ylim_use <- c(0.5, n + 1.5)
    } else {
      y_header <- NA
      y_coef   <- n:1
      ylim_use <- c(0.5, n + 0.5)
    }

    par(mar = mar, cex.axis = cex_x, cex.lab = cex_lab)

    plot(NA, xlim = xlim, ylim = ylim_use,
         yaxt = "n", ylab = "", xlab = xlab, main = "", ...)

    title(main = as.character(g), cex.main = cex_main)

    abline(v = refline, lty = ref_lty, col = "gray50")

    # draw CI + points (ONLY for coefficient rows)
    segments(lower[idx], y_coef, upper[idx], y_coef, lwd = ci_lwd)
    points(parm[idx],  y_coef, pch = point_pch, cex = point_cex)

    # y-axis labels (ONLY for coefficient rows)
    axis(2, at = y_coef, labels = lab_g, las = 1, tick = FALSE, cex.axis = cex_y)

    # header as an extra "first row" (ONLY when show_est_ci=TRUE)
    if (isTRUE(show_est_ci)) {
      hdr <- if (is.null(est_ci_header)) xlab else est_ci_header
      op <- par(xpd = NA)
      mtext(hdr, side = 2, at = y_header, line = par("mgp")[2],
            las = 1, adj = 1, cex = est_ci_header_cex, font = 2)
      par(op)
    }

    box()
  }

  invisible(list(xlim = xlim, groups = levs))
}
