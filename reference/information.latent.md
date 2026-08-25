# Information Variance-Covariance Matrix for Latent Models

Compute the variance-covariance matrix of the freely estimated
parameters of a fitted latent-variable model.

## Usage

``` r
# S3 method for class 'latent'
information(fit)
```

## Arguments

- fit:

  A fitted object inheriting from class `"latent"`.

## Value

A list containing the Hessian, variance-covariance matrix, standard
errors, and covariance-method metadata. Multistep objects additionally
return the propagation matrices and joint covariance matrix.

## Details

Three covariance calculations are handled directly. Deterministic
multistep models use the covariance propagated from preceding estimation
steps. Source estimators such as
[`lmean()`](https://marcosjnez.github.io/latent/reference/lmean.md),
[`lpearson()`](https://marcosjnez.github.io/latent/reference/lpearson.md),
[`lpoly()`](https://marcosjnez.github.io/latent/reference/lpoly.md),
[`lyule()`](https://marcosjnez.github.io/latent/reference/lyule.md), and
[`lmvnorm()`](https://marcosjnez.github.io/latent/reference/lmvnorm.md)
use their analytic information covariance. Remaining latent models use
the inverse Hessian.

If `fit@modelInfo$control_optimizer$se_method = "KKT"`, the covariance
of the freely estimated parameters is obtained from the parameter block
of the inverse bordered KKT matrix. The upper-left block is the Hessian
of the Lagrangian, including the second derivatives of active
constraints. When no active constraint is found, the ordinary Hessian
inverse is used.

For ordinary CFA maximum-likelihood models, the Hessian is evaluated
with the corresponding `cfa_ml` estimator because `cfa_fml` represents
the fitted discrepancy function rather than the information contribution
used for the covariance matrix.
