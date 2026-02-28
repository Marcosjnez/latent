# Local Bivariate Residuals for Latent Class Analysis

Computes bivariate residuals (BVR) and related diagnostics for a fitted
latent class model with mixed indicators (continuous and categorical).
The function returns a matrix of residuals, p-values, a correlation-like
effect size, and a summary table. The root mean square residual (RMSR)
provides a global fit measure.

## Usage

``` r
lbvr(model, digits = 4)
```

## Arguments

- model:

  An object of class \`"lca"\` from the latent package.

- digits:

  Integer; number of decimal places to use in rounding the output
  matrices and the summary table. Default is 4.

## Value

A list of class \`"lbvr"\` with components:

- residual_matrix:

  A symmetric matrix of residual measures:

  - For continuous–continuous pairs: Pearson correlation of residuals.

  - For categorical–categorical pairs: maximum absolute standardized
    residual.

  - For continuous–categorical pairs: maximum absolute Z‑score of mean
    residuals per category.

- pvalue_matrix:

  Symmetric matrix of p-values corresponding to the residual measures
  (using correlation test, chi‑square test, or ANOVA F‑test).

- r_mat:

  Symmetric matrix of effect sizes comparable to a correlation:

  - Continuous–continuous: Pearson \\r\\.

  - Categorical–categorical: difference of Cramér's V (\\V\_{obs} -
    V\_{exp}\\).

  - Continuous–categorical: \\\eta\\ (square root of eta‑squared from
    ANOVA).

- rmsr:

  Root mean square of the upper triangle of \`r_mat\`; a global fit
  index.

- res_tab:

  A data frame with one row per variable pair, containing: `pair`,
  `statistic` (the residual measure), `p_value`, `r` (the
  correlation‑like effect size), `var_types`, and `test` (a descriptive
  label). The table is sorted by absolute `r`.

## Details

The function extracts data, variable types, posterior probabilities,
class‑specific means, and response probabilities from a fitted latent
model. It then computes residuals for each pair of indicators using only
pairwise complete observations. Missing data are handled by scaling
expected frequencies to the number of complete cases for the pair.

For categorical variables, the original data (which may be zero‑based)
are incremented by 1 to ensure proper table construction.

## Examples

``` r
if (FALSE) { # \dontrun{
library(latent)
# Fit a 3‑class model with multinomial indicators
fit <- lca(data = gss82, nclasses = 3,item = rep("multinomial", ncol(gss82)))
# Compute bivariate residuals
lbvr(fit)
} # }
```
