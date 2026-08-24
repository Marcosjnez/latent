# Hessian Matrix for Latent Models

The Hessian is computed with respect to the freely estimated parameters
on their optimization scale. Transformed parameters are therefore not
included as separate rows or columns.

## Usage

``` r
# S3 method for class 'latent'
hessian(fit)
```

## Arguments

- fit:

  A fitted object inheriting from class `"latent"`.

- riemannian:

  Logical. If `FALSE`, return the Euclidean Hessian. If `TRUE`, return
  the ambient constrained inverse operator \\P=T H_R^{-1}T^\top\\, where
  \\H_R\\ is the tangent-coordinate Riemannian Hessian and \\T\\ is the
  tangent-space basis.

## Value

A symmetric numeric matrix with one row and one column for each freely
estimated parameter. Row and column names correspond to
`fit@modelInfo$parameters_labels`. When `riemannian = FALSE`, the matrix
is the Euclidean Hessian. When `riemannian = TRUE`, it is the ambient
constrained inverse operator \\P\\.

## Details

Compute the Hessian matrix of the objective function evaluated at the
parameter estimates of a fitted `"latent"` object.

The Hessian is evaluated at the parameter estimates stored in
`fit@Optim`. For models estimated by an alternative optimization
criterion, such as expectation-maximization, the estimator stored after
fitting determines the objective whose derivatives are evaluated.

The Hessian is used internally by
[`vcov.latent`](https://marcosjnez.github.io/latent/reference/vcov.latent.md)
to obtain covariance matrices and standard errors.

## See also

[`hessian`](https://marcosjnez.github.io/latent/reference/hessian.md),
[`vcov.latent`](https://marcosjnez.github.io/latent/reference/vcov.latent.md),
[`jacobian.latent`](https://marcosjnez.github.io/latent/reference/jacobian.latent.md)
