# Plot mean and standard deviation profiles for Gaussian LCA indicators

Plot mean and standard deviation profiles for Gaussian LCA indicators

## Usage

``` r
plot_lca_gaussian(
  data,
  title = NULL,
  ylab = "Mean score (±1 SD)",
  caption = "Error bars indicate ±1 standard deviation",
  nrow = NULL,
  ncol = NULL,
  point_size = 2.5,
  err_width = 0.2,
  base_size = 12,
  angle = 45,
  ...
)
```

## Arguments

- data:

  Data frame from `reshape_lca_continuous`.

- title:

  Plot title (default `NULL`).

- ylab:

  y-axis label.

- caption:

  Plot caption.

- nrow, ncol:

  Number of rows/columns in facet wrap.

- point_size:

  Size of the mean points.

- err_width:

  Width of error bars.

- base_size:

  Base font size.

- angle:

  Angle for x-axis labels.

- ...:

  Other arguments (currently ignored).

## Value

A ggplot object.
