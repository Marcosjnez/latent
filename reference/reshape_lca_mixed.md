# Reshape class-conditional indicator profiles

Divides class-conditional indicator profiles into multinomial and
Gaussian components and converts them into long data frames suitable for
plotting.

## Usage

``` r
reshape_lca_mixed(item_output, item_types)
```

## Arguments

- item_output:

  Named list of class-conditional indicator matrices, normally returned
  by `latInspect(model, what = "item")`.

- item_types:

  Named list containing the character vectors `multinomial` and
  `gaussian`. An optional `mvgaussian` component is combined with
  `gaussian`.

## Value

A list with the following components:

- multinomial:

  A data frame with columns `variable`, `response`, `class`, and
  `probability`.

- gaussian:

  A data frame with columns `variable`, `class`, `mean`, and `sd`.
