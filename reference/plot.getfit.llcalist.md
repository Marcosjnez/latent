# Plot method for LCA model fit comparisons (scree plot)

Plots the requested information criteria across latent class solutions.

## Usage

``` r
# S3 method for class 'getfit.llcalist'
plot(x, indices = NULL, ...)
```

## Arguments

- x:

  An object of class `"getfit.llcalist"`, typically from `getfit()`.

- indices:

  Character vector of information criteria to plot. If `NULL`, defaults
  to `c("AICp","BICp")` when the object is penalised, otherwise
  `c("AIC","BIC")`.

- ...:

  Additional arguments (currently unused).

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) { # \dontrun{
fit_ind_cont <- getfit(fit)
plot(fit_ind_cont)
plot(fit_ind_cont, indices = c("AIC", "BIC", "KIC"))
} # }
```
