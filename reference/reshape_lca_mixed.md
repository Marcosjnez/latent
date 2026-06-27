# Split and reshape mixed (multinomial + Gaussian) LCA item output

Split and reshape mixed (multinomial + Gaussian) LCA item output

## Usage

``` r
reshape_lca_mixed(item_output, item_types)
```

## Arguments

- item_output:

  Result of `latInspect(fit, "item")`.

- item_types:

  Character vector of indicator types (e.g., `fit@dataList$item`).

## Value

A list with components `multinomial` and `gaussian`, each a data frame.
