# Forest plot for grouped parameter estimates

Draws one base R forest plot for each group of parameter estimates. Each
panel displays point estimates, confidence intervals, and a reference
line.

## Usage

``` r
forestplot(
  parm,
  lower,
  upper,
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
  ...
)
```

## Arguments

- parm:

  Numeric vector containing parameter estimates.

- lower, upper:

  Numeric vectors containing the lower and upper confidence limits. They
  must have the same length as `parm`.

- group:

  Optional vector defining groups of parameters. One plot is drawn for
  each group. When `NULL`, all parameters are placed in one group.

- labels:

  Optional character vector containing coefficient labels. When `NULL`,
  names from `parm` are used when available.

- refline:

  Numeric value indicating the vertical reference line.

- xlab:

  Character string with the horizontal-axis label.

- xlim:

  Optional numeric vector of length two defining the horizontal-axis
  limits. When `NULL`, limits are calculated from the confidence
  intervals and reference line.

- order_within_group:

  Character string indicating whether coefficients retain their supplied
  order (`"as_is"`) or are ordered by their estimate within each group
  (`"by_parm"`).

- point_pch:

  Plotting symbol used for parameter estimates.

- point_cex:

  Numeric expansion factor for parameter-estimate symbols.

- ci_lwd:

  Width of confidence-interval segments.

- ref_lty:

  Line type used for the reference line.

- mar:

  Numeric vector passed to `par(mar = ...)`.

- show_est_ci:

  Logical value indicating whether estimates and confidence intervals
  are appended to coefficient labels.

- digits:

  Number of decimal places used when formatting estimates and confidence
  intervals.

- cex_y, cex_x:

  Numeric expansion factors for vertical- and horizontal-axis labels.

- cex_lab:

  Numeric expansion factor for the horizontal-axis title.

- cex_main:

  Numeric expansion factor for group titles.

- est_ci_header:

  Optional character string displayed above the coefficient labels when
  `show_est_ci = TRUE`. By default, `xlab` is used.

- est_ci_header_cex:

  Optional expansion factor for the estimate-and- confidence-interval
  header.

- ...:

  Additional graphical arguments passed to
  [`plot.default()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

A list containing the horizontal-axis limits and plotted group names,
invisibly.

## Details

Invalid rows containing non-finite estimates or confidence limits, or
missing groups, are removed with a warning. All groups use the same
horizontal-axis limits so their estimates can be compared directly.

## Examples

``` r
if (FALSE) { # \dontrun{
forestplot(parm = c(0.2, -0.1, 0.5),
           lower = c(0.05, -0.3, 0.2),
           upper = c(0.35, 0.1, 0.8),
           group = c("Class 1", "Class 1", "Class 2"),
           labels = c("Age", "Sex", "Age"))
} # }
```
