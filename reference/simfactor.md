# Simulate factor structures with misspecification errors.

Simulate factor and bifactor structures with crossloadings, correlated
factors, and more.

## Usage

``` r
simfactor(nfactors = 5, nitems = 6, loadings = "medium",
crossloadings = 0, correlations = 0,
estimator = "minres", fit = "rmsr", misfit = 0,
error_method = "cudeck", efa = FALSE,
ngenerals = 0, loadings_g = "medium", correlations_g = 0,
pure = FALSE,
lambda = NULL, Phi = NULL, Psi = NULL)
```

## Arguments

- nfactors:

  Number of factors.

- nitems:

  Number of items per factor.

- loadings:

  Loadings' magnitude on the factors: "low", "medium" or "high".
  Defaults to "medium".

- crossloadings:

  Magnitude of the cross-loadings among the group factors. Defaults to
  0.

- correlations_g:

  Correlation among the general factors. Defaults to 0.

- correlations:

  Correlation among the factors. Defaults to 0.

- estimator:

  estimator used to generate population error: "minres" or "ml".

- fit:

  Fit index to control the population error.

- misfit:

  Misfit value to generate population error.

- error_method:

  Method used to control population error: c("yuan", "cudeck"). Defaults
  to "cudeck".

- efa:

  Reproduce the error with EFA or CFA. Defaults to FALSE (CFA).

- ngenerals:

  Number of general factors.

- loadings_g:

  Loadings' magnitude on the general factors: "low", "medium" or "high".
  Defaults to "medium".

- pure:

  Fix a pure item on each general factor. Defaults to FALSE.

- lambda:

  Custom loading matrix. If Phi is NULL, then all the factors will be
  correlated at the value given in correlations.

- Phi:

  Custom Phi matrix. If lambda is NULL, then Phi should be conformable
  to the loading matrix specified with the above arguments.

- Psi:

  Custom Psi matrix.

## Value

List with the following objects:

- lambda:

  Population loading matrix.

- Phi:

  Population factor correlation matrix.

- Psi:

  Population covariance matrix between the errors.

- R:

  Model correlation matrix.

- R_error:

  Model correlation matrix with misspecification errors.

- uniquenesses:

  Population uniquenesses.

- delta:

  Minimum of the loss function that correspond to the misfit value.

## Details

`simfactor` generates bi-factor and generalized bifactor patterns with
cross-loadings, pure items and correlations among the general and group
factors. When `crossloading` is different than 0, one cross-loading is
introduced for an item pertaining to each group factor. When `pure` is
TRUE, one item loading of each group factor is removed so that the item
loads entirely on the general factor. To maintain the item communalities
constant upon these modifications, the item loading on the other factors
may shrunk (if adding cross-loadings) or increase (if setting pure
items).

Loading magnitudes may range between 0.3-0.5 ("low"), 0.4-0.6 ("medium")
and 0.5-0.7 ("high"). Custom ranges can be supplied as vectors (i.e.,
c(0.2, 0.5))

## References

Cudeck, R., & Browne, M. W. (1992). Constructing a covariance matrix
that yields a specified minimizer and a specified minimum discrepancy
function value. *Psychometrika, 57*(3), 357–369.
[doi:10.1007/BF02295424](https://doi.org/10.1007/BF02295424)

Jiménez, M., Abad, F. J., Garcia-Garzon, E., & Garrido, L. E. (2023).
Exploratory Bi-factor Analysis with Multiple General Factors.
*Multivariate Behavioral Research, 58*(6), 1072–1089.
[doi:10.1080/00273171.2023.2189571](https://doi.org/10.1080/00273171.2023.2189571)

Jiménez, M., Abad, F. J., Garcia-Garzon, E., Golino, H., Christensen, A.
P., & Garrido, L. E. (2023). Dimensionality assessment in bifactor
structures with multiple general factors: A network psychometrics
approach. *Psychological Methods*. Advance online publication.
[doi:10.1037/met0000590](https://doi.org/10.1037/met0000590)

Yuan, K.-H., & Hayashi, K. (2003). Bootstrap approach to inference and
power analysis based on three test statistics for covariance structure
models. *British Journal of Mathematical and Statistical Psychology,
56*(1), 93–110.
[doi:10.1348/000711003321645368](https://doi.org/10.1348/000711003321645368)

## Examples

``` r
# Simulate data:
sim <- simfactor(nfactors = 3, nitems = 4, correlations = 0.40,
                 crossloadings = 0.30)
sim$lambda
#>               S1        S2        S3
#> item1  0.4807076 0.0000000 0.3000000
#> item2  0.4127323 0.0000000 0.0000000
#> item3  0.4777403 0.0000000 0.0000000
#> item4  0.5951096 0.0000000 0.0000000
#> item5  0.3000000 0.4579785 0.0000000
#> item6  0.0000000 0.5356761 0.0000000
#> item7  0.0000000 0.5470639 0.0000000
#> item8  0.0000000 0.4391913 0.0000000
#> item9  0.0000000 0.3000000 0.5961079
#> item10 0.0000000 0.0000000 0.5483043
#> item11 0.0000000 0.0000000 0.4102893
#> item12 0.0000000 0.0000000 0.5060425
sim$Phi
#>      [,1] [,2] [,3]
#> [1,]  1.0  0.4  0.4
#> [2,]  0.4  1.0  0.4
#> [3,]  0.4  0.4  1.0
scores <- MASS::mvrnorm(1e3, rep(0, nrow(sim$R_error)), Sigma = sim$R_error)
s <- cor(scores)
```
