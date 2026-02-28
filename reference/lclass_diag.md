# Classification Diagnostics for Latent Class Analysis

Computes various classification diagnostics from a fitted latent class
model.

## Usage

``` r
lclass_diag(
  x,
  digits = 4,
  type = c("Entropy", "AvePP", "Mostlikely.Class", "Sum.Mostlikely"),
  ...
)
```

## Arguments

- x:

  An object of class \`"lca"\` from the latent package.

- digits:

  Number of decimal places to use for rounding (default = 4).

- type:

  Character vector specifying which diagnostics to compute. Available
  options:

  - `"Entropy"` – Entropy‑based R² (classification uncertainty).

  - `"AvePP"` – Average posterior probability for each class.

  - `"OCC"` – Odds of Correct Classification.

  - `"Overall.Misclassification"` – Overall misclassification rate.

  - `"Misclassification.per.class"` – Misclassification rate per true
    class.

  - `"Sum.Posterior"` – Class proportions based on posterior sums.

  - `"Sum.Mostlikely"` – Class proportions based on modal assignment.

  - `"Mostlikely.Class"` – Classification probabilities (true ×
    assigned).

  - `"Avg.Mostlikely"` – Average posterior probabilities by assigned
    class.

  - `"all"` – A convenience shortcut to request all of the above.

- ...:

  Additional arguments passed to `latInspect` or `getfit`.

## Value

An object of class `"lclassd"`, a list containing the requested
diagnostics.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Fit a model
  library(latent)
  # Fit a 3‑class model with multinomial indicators
  fit <- lca(data = gss82, nclasses = 3,item = rep("multinomial", ncol(gss82)))

  # Get default diagnostics (Entropy, AvePP, Mostlikely.Class, Sum.Mostlikely)
  diag <- lclass_diag(fit)
  print(diag)

  # Request all diagnostics
  diag_all <- lclass_diag(fit, type = "all")
  print(diag_all)
} # }
```
