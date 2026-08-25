# Standard Errors for Multistep Latent Models

Compute standard errors for deterministic multistep estimators by
propagating the joint variance-covariance matrix of parameters estimated
in previous steps.

## Usage

``` r
# S3 method for class 'multistep'
se(fit, type = NULL, parameters = NULL, digits = 4L, ...)

# S3 method for class 'multistep_lcfa'
se(fit, ...)
```

## Arguments

- fit:

  A fitted object inheriting from class `"multistep"`.

- type:

  Optional compatibility argument. If supplied, it is validated but does
  not replace the sample-statistic covariance method selected when the
  multistep object was fitted.

- parameters:

  Optional parameter specification identifying the parameters or
  transformed parameters for which standard errors should be returned.
  If `NULL`, the parameter structures of the top-level model are used.

- digits:

  Non-negative integer specifying the number of decimal places used in
  formatted parameter tables. If `NULL`, values are not rounded.

- ...:

  Additional arguments reserved for future multistep methods.

## Value

A list containing parameter tables, standard errors, the selected
variance-covariance matrix, the final joint variance-covariance matrix,
and the Hessian and cross-derivative matrices used in the final step.

## Details

A multistep estimator is deterministic conditional on the parameter
estimates passed from preceding steps. Consequently, no inverse-Hessian
variance is added for the new step.

Let \\A\\ be the joint variance-covariance matrix of parameters
estimated in all preceding steps, \\P_2\\ the inverse-Hessian operator
for parameters estimated in the current step, and \\C\\ the
cross-Hessian between preceding and current parameters in the
corresponding unrestricted model. The covariance matrix is enlarged at
every step as \$\$ V = \begin{pmatrix} A & -A C P_2 \\ -P_2^\top C^\top
A & P_2 C^\top A C P_2^\top \end{pmatrix}. \$\$

By default, \\P_2=H_2^{-1}\\. If the current-step control contains
`se_method = "KKT"` and active equality constraints are available,
\\P_2\\ is the parameter block of the inverse KKT matrix.

The final joint covariance matrix is passed to
[`vcov.latent()`](https://marcosjnez.github.io/latent/reference/vcov.latent.md),
which propagates it through all requested parameter transformations.
