# Plot mixed LCA results (both Gaussian and multinomial indicators)

Plot mixed LCA results (both Gaussian and multinomial indicators)

## Usage

``` r
plot_lca_mixed(mixed_data, ...)
```

## Arguments

- mixed_data:

  Output from `reshape_lca_mixed`.

- ...:

  Extra arguments passed to `plot_lca_gaussian` and
  `plot_lca_multinomial`.

## Value

A list of ggplot objects, named `gaussian` and `multinomial`.
