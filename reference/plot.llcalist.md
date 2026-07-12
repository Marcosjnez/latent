# Plot results from a collection of LCA models

Selects the appropriate fitted model from an `llcalist` object and
displays measurement profiles or latent-class regression coefficients
using the corresponding
[`plot.llca()`](https://marcosjnez.github.io/latent/reference/plot.llca.md)
method.

## Usage

``` r
# S3 method for class 'llcalist'
plot(
  x,
  type = c("all", "gaussian", "multinomial", "coefficients"),
  variables = NULL,
  ...
)
```

## Arguments

- x:

  An object of class `"llcalist"`.

- type:

  Character string specifying the result to plot. Available values are
  `"all"`, `"gaussian"`, `"multinomial"`, and `"coefficients"`.

- variables:

  Optional character vector selecting measurement indicators for profile
  plots. This argument is passed to
  [`plot.llca()`](https://marcosjnez.github.io/latent/reference/plot.llca.md).

- ...:

  Additional arguments passed to
  [`plot.llca()`](https://marcosjnez.github.io/latent/reference/plot.llca.md).

## Value

The result returned by
[`plot.llca()`](https://marcosjnez.github.io/latent/reference/plot.llca.md),
invisibly. When several models are plotted, a named list of results is
returned invisibly.

## Details

For a structural-after-measurement result containing components named
`measurement` and `structural`, profile plots use the measurement model
and coefficient plots use the structural model. This allows the full
object returned by
[`lca()`](https://marcosjnez.github.io/latent/reference/lca.md) to be
plotted directly.

For other `llcalist` objects, a single model is plotted directly. When
several models are present, profile plots are produced sequentially for
all models. Coefficient plots are produced only for models that include
covariates.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lca(data = empathy, nclasses = 4,
           gaussian = c("ec1", "ec2", "ec3", "ec4", "ec5", "ec6"),
           covariates = c("pt1", "pt2", "pt3", "pt4"))
plot(fit)
plot(fit, type = "coefficients", what = "OR")
} # }
```
