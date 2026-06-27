# Convert a fitted lcfa object (latent) to a lavaan object

Convert a fitted lcfa object (latent) to a lavaan object

## Usage

``` r
lcfa_to_lavaan(object, ...)
```

## Arguments

- object:

  A fitted `lcfa` object from the `latent` package.

- ...:

  Additional arguments passed to
  [`lavaan::cfa()`](https://rdrr.io/pkg/lavaan/man/cfa.html).

## Value

A `lavaan` object populated with estimates from the `lcfa` fit.

## Details

This function reconstructs a `lavaan` S4 object from a fitted `lcfa`
object, allowing the use of lavaan post-processing functions such as
[`fitMeasures()`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html),
[`lavResiduals()`](https://rdrr.io/pkg/lavaan/man/lavResiduals.html),
and
[`modindices()`](https://rdrr.io/pkg/lavaan/man/modificationIndices.html).

The conversion works by:

1.  Building a lavaan scaffold (un-fitted) with the same model syntax,
    data, estimator, and options.

2.  Injecting lcfa parameter estimates into the lavaan parameter table
    via plabel matching.

3.  Rebuilding the lavaan Model object and implied moments from the
    updated parameter table.

4.  Computing test statistics, log-likelihood, baseline model, and VCOV
    using standard lavaan internals.

## Examples

``` r
if (FALSE) { # \dontrun{
library(latent)
library(lavaan)

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit_latent <- lcfa(model = HS.model, data = HolzingerSwineford1939,
                   estimator = "ml", std.lv = TRUE, meanstructure = TRUE)
fit_lavaan <- lcfa_to_lavaan(fit_latent)
fitMeasures(fit_lavaan)
lavResiduals(fit_lavaan)
modindices(fit_lavaan, sort = TRUE)
} # }
```
