# Sorted Parameter Structures for Inspection

Sorted Parameter Structures for Inspection

## Usage

``` r
sort_latent(fit)
```

## Arguments

- fit:

  A fitted llca, lcfa, lefa, or lrotate result.

## Value

A list with four elements: `parameters`, `transformed_pars`, `param`,
and `trans`. These are copies of the corresponding estimate and label
structures, not a fitted object. Original labels and factor/class names
travel with their entries. The `order` and `signs` attributes describe
the display permutation and signs for use by
[`latInspect()`](https://marcosjnez.github.io/latent/reference/latInspect.md);
nothing is stored in fit.

## Details

Classes are ordered by decreasing frequency-weighted posterior size.
Factors are ordered by decreasing `colSums(lambda^2)*diag(psi)` and
oriented so their largest absolute loading is positive. Ties retain the
original order. Fixed modeled CFA loadings retain their specified order
and signs. Target rotations and components selecting factors retain
factor order. The original multinomial reference class and all parameter
labels are kept.
