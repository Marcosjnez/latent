# Confidence intervals for a collection of latent class models

Applies `ci()` to every fitted model in an `"llcalist"` object.

## Usage

``` r
# S3 method for class 'llcalist'
ci(model, type = "standard", confidence = 0.95, digits = 3L, ...)
```

## Arguments

- model:

  An object of class `"llcalist"` containing fitted `"llca"` objects.

- type:

  Character string indicating the standard-error estimator. See
  [`ci.llca()`](https://marcosjnez.github.io/latent/reference/ci.llca.md).

- confidence:

  Numeric scalar strictly between zero and one specifying the confidence
  level.

- digits:

  Non-negative integer indicating the number of decimal places used in
  the formatted tables.

- ...:

  Additional arguments passed to
  [`ci.llca()`](https://marcosjnez.github.io/latent/reference/ci.llca.md).

## Value

A named list with one confidence-interval result per fitted model and
class `"ci.llcalist"`.

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
ci(fits, confidence = 0.90)
} # }
```
