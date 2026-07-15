# Standard errors for a collection of latent class models

Applies `se()` to every fitted model in an `"llcalist"` object.

## Usage

``` r
# S3 method for class 'llcalist'
se(model, type = "standard", parameters = NULL, digits = 4L, ...)
```

## Arguments

- model:

  An object of class `"llcalist"` containing fitted `"llca"` objects.

- type:

  Character string indicating the standard-error estimator. See
  [`se.llca()`](https://marcosjnez.github.io/latent/reference/se.llca.md).

- digits:

  Non-negative integer indicating the number of decimal places used in
  the formatted tables, or `NULL` to avoid rounding.

- ...:

  Additional arguments passed to
  [`se.llca()`](https://marcosjnez.github.io/latent/reference/se.llca.md).

## Value

A named list with one standard-error result per fitted model and class
`"se.llcalist"`.

## Details

Existing names are preserved. Consequently, class-enumeration results
retain names such as `"nclasses=2"`, whereas multiple-step models retain
names such as `"measurement"` and `"structural"`. Unnamed elements are
labelled according to their number of latent classes.

## Examples

``` r
if (FALSE) { # \dontrun{
fits <- lca(data = empathy, nclasses = 2:4,
            gaussian = c("ec1", "ec2", "ec3"))
se(fits)
} # }
```
