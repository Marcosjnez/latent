# Plot response probabilities for multinomial LCA indicators

Plot response probabilities for multinomial LCA indicators

## Usage

``` r
plot_lca_multinomial(
  data,
  variables = NULL,
  bars = "variable",
  facet = "class",
  bw = FALSE,
  angle = 45,
  ...
)
```

## Arguments

- data:

  Data frame from `reshape_lca_multinomial`.

- variables:

  Optional character vector of variables to include.

- bars:

  Variable to map on x-axis (default `"variable"`).

- facet:

  Faceting variable (default `"class"`).

- bw:

  Logical; if `TRUE`, use greyscale instead of colour.

- ...:

  Other arguments (currently ignored).

## Value

A ggplot object.
