# Plot fit indices across latent class solutions

Creates a scree-style base R plot comparing information criteria across
the models contained in a fitted latent class enumeration.

## Usage

``` r
# S3 method for class 'getfit.llcalist'
plot(
  x,
  indices = NULL,
  title = "Model fit by number of classes",
  xlab = "Class solution",
  ylab = "Value",
  base_size = 12,
  colors = NULL,
  pch = 19,
  lwd = 2,
  ...
)
```

## Arguments

- x:

  An object of class `"getfit.llcalist"`, normally returned by
  `getfit()` for an `llcalist` object.

- indices:

  Character vector containing the names of the fit indices to plot. When
  `NULL`, the function uses `AICp` and `BICp` for penalized models and
  `AIC` and `BIC` otherwise.

- title:

  Character string with the plot title.

- xlab:

  Character string with the horizontal-axis label.

- ylab:

  Character string with the vertical-axis label.

- base_size:

  Numeric value controlling the overall text size.

- colors:

  Optional character vector of line colors. Colors are recycled when
  fewer colors than indices are supplied.

- pch:

  Plotting symbol used for the fitted values.

- lwd:

  Width of the lines connecting fitted values.

- ...:

  Additional arguments. Currently unused.

## Value

The plotted data, invisibly.

## Details

Each selected information criterion is represented by a separate line.
The horizontal axis gives the number of latent classes and the vertical
axis gives the value of the selected criterion.

## Examples

``` r
if (FALSE) { # \dontrun{
fits <- lca(data = gss82, nclasses = 1:5,
            multinomial = c("PURPOSE", "ACCURACY"))
fit_indices <- getfit(fits)
plot(fit_indices)
plot(fit_indices, indices = c("AIC", "BIC", "SABIC"))
} # }
```
