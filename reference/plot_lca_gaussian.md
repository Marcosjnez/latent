# Plot Gaussian indicator profiles

Displays the class-specific means and standard deviations of Gaussian or
multivariate-Gaussian indicators using base R graphics.

## Usage

``` r
plot_lca_gaussian(
  data,
  variables = NULL,
  title = NULL,
  ylab = "Mean score (<U+00B1>1 SD)",
  caption = "Error bars indicate <U+00B1>1 standard deviation",
  nrow = NULL,
  ncol = NULL,
  point_size = 1.2,
  err_width = 0.12,
  base_size = 12,
  angle = 45,
  color = "steelblue",
  ...
)
```

## Arguments

- data:

  Data frame produced by `reshape_lca_continuous()`.

- variables:

  Optional character vector selecting the indicators to plot.

- title:

  Optional character string with the plot title.

- ylab:

  Character string with the vertical-axis label.

- caption:

  Optional character string with the plot caption.

- nrow, ncol:

  Optional numbers of rows and columns in the panel layout.

- point_size:

  Numeric expansion factor for the points representing means.

- err_width:

  Numeric half-width of the horizontal caps on the standard- deviation
  bars.

- base_size:

  Numeric value controlling the overall text size.

- angle:

  Numeric angle of the horizontal-axis labels.

- color:

  Color used for class-profile lines and points.

- ...:

  Additional arguments. Currently unused.

## Value

The plotted profile data, invisibly.

## Details

A separate panel is drawn for every latent class. Error bars represent
one class-specific standard deviation below and above each mean.
