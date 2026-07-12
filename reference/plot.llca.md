# Plot results from an LCA model

Displays class-conditional measurement profiles or latent-class
regression coefficients from a fitted latent class model using base R
graphics.

## Usage

``` r
# S3 method for class 'llca'
plot(
  x,
  type = c("all", "gaussian", "multinomial", "coefficients"),
  variables = NULL,
  ...
)
```

## Arguments

- x:

  A fitted object of class `"llca"`.

- type:

  Character string specifying the result to plot. Available values are
  `"all"`, `"gaussian"`, `"multinomial"`, and `"coefficients"`. Partial
  matching is supported.

- variables:

  Optional character vector selecting measurement indicators for profile
  plots. This argument is ignored when `type = "coefficients"`; use
  `predictors` through `...` instead.

- ...:

  Additional arguments passed to
  [`plot_lca_mixed()`](https://marcosjnez.github.io/latent/reference/plot_lca_mixed.md)
  for profile plots or to the coefficient-plotting routine for
  `type = "coefficients"`. Coefficient options include `se_type`,
  `what`, `effects`, `confidence`, `predictors`, `intercept`, and the
  graphical arguments of
  [`forestplot()`](https://marcosjnez.github.io/latent/reference/forestplot.md).

## Value

For profile plots, an object of class `"llca_plots"`, invisibly. For
coefficient plots, a list containing the displayed coefficient
estimates, standard errors, variance-covariance matrix, and confidence
limits, invisibly.

## Details

Gaussian and multivariate-Gaussian indicators are displayed using their
class-specific means and standard deviations. Multinomial indicators are
displayed using their class-specific response probabilities. Response
colors are reused across indicators according to category position,
while each indicator has its own legend containing its actual response
levels. Outcomes and covariates are not included in the
measurement-profile plots.

For mixed models with `type = "all"`, the Gaussian and multinomial plots
are drawn sequentially as separate pages on the active graphics device.
With `type = "coefficients"`, the latent-class regression coefficients
are displayed in forest-plot form. Use `what = "OR"` for odds ratios or
`what = "log"` for log odds ratios.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lca(data = gss82, nclasses = 3,
           multinomial = c("PURPOSE", "ACCURACY"))
plot(fit)
plot(fit, type = "multinomial", bw = TRUE)

structural_fit <- lca(data = gss82, nclasses = 3,
                      multinomial = c("PURPOSE", "ACCURACY"),
                      covariates = c("RACE", "SEX"))[["structural"]]
plot(structural_fit, type = "coefficients", what = "OR",
     effects = "coding", predictors = "RACE")
} # }
```
