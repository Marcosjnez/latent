# Local bivariate residuals for latent class analysis

Computes pairwise residual diagnostics for the indicators of a fitted
latent class model. Different diagnostics are used for pairs of Gaussian
indicators, pairs of multinomial indicators, and mixed
Gaussian-multinomial pairs.

## Usage

``` r
# S3 method for class 'llca'
lbvr(x, digits = 4L, ...)
```

## Arguments

- x:

  A fitted object of class `"llca"`.

- digits:

  A non-negative integer indicating the number of decimal places used to
  round the returned matrices and summary table. The default is `4L`.

- ...:

  Additional arguments reserved for other methods.

## Value

An object of class `"lbvr"`, consisting of a list with the following
components:

- residual_matrix:

  A symmetric matrix containing the standardized pairwise diagnostic.

- raw_residual_matrix:

  A symmetric matrix containing the raw residual score. Mixed
  Gaussian-multinomial pairs contain `NA_real_`.

- pvalue_matrix:

  A symmetric matrix containing the corresponding p-values.

- r_mat:

  A symmetric matrix containing the correlation-like effect sizes.

- rmsr:

  The root mean square of the upper-triangular values of `r_mat`.

- res_tab:

  A data frame with one row per indicator pair and columns `pair`,
  `residual`, `statistic`, `p_value`, `r`, `var_types`, and `test`. Rows
  are ordered by the absolute value of `r`.

## Details

The function evaluates every pair of measurement indicators. Covariates
and distal outcomes are not included. Pairwise complete observations are
used for each diagnostic.

For two Gaussian indicators, the standardized diagnostic is the Pearson
correlation between their posterior-weighted residuals. The raw residual
is obtained with the LatentGold-style score calculation implemented in
[`latentgold_cont_cont_bvr()`](https://marcosjnez.github.io/latent/reference/latentgold_cont_cont_bvr.md).

For two multinomial indicators, the function compares the observed
cross-classification table with the model-implied table. It reports the
maximum absolute standardized residual, its chi-squared p-value, the
difference between the observed and expected Cramer's V values, and the
raw residual score used by the package.

For a Gaussian and a multinomial indicator, the function reports the
maximum absolute category-specific Z-score of the Gaussian residuals,
the p-value from the ANOVA calculation, and eta as a correlation-like
effect size. A raw residual is not defined for these mixed pairs and is
returned as `NA_real_`.

Indicator pairs already included as residual dependencies in the fitted
model are assigned a residual of zero, a p-value of one, and an effect
size of zero.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lca(data = gss82, nclasses = 3L,
           multinomial = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"))

residuals <- lbvr(fit)
residuals$res_tab
} # }
```
