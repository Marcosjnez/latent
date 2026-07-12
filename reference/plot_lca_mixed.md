# Plot mixed indicator profiles

Draws the appropriate base R profile plot for each available indicator
family.

## Usage

``` r
plot_lca_mixed(mixed_data, variables = NULL, ...)
```

## Arguments

- mixed_data:

  List produced by
  [`reshape_lca_mixed()`](https://marcosjnez.github.io/latent/reference/reshape_lca_mixed.md).

- variables:

  Optional character vector selecting indicators to plot.

- ...:

  Additional arguments passed to
  [`plot_lca_gaussian()`](https://marcosjnez.github.io/latent/reference/plot_lca_gaussian.md)
  and
  [`plot_lca_multinomial()`](https://marcosjnez.github.io/latent/reference/plot_lca_multinomial.md).

## Value

An object of class `"llca_plots"`, invisibly. The object stores the data
and plotting arguments and can be printed to redraw the plots.

## Details

If only one indicator family is available, one plot is drawn. For mixed
models, the Gaussian and multinomial plots are drawn sequentially as two
pages on the active graphics device.
