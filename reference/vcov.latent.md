# Variance-Covariance Matrix for Latent Models

`vcov.latent()` computes the variance-covariance matrix of freely
estimated or transformed parameters of a fitted object inheriting from
class `"latent"`.

For ordinary one-step models, covariance matrices are obtained from the
inverse Hessian and transformed to the requested parameterization using
the delta method.

When the fitted model contains a previous `"latent"` object among the
model specifications stored in `fit@modelInfo$control_optimizer$model`,
uncertainty from the previous estimation step is propagated to the
current parameter estimates. Nested fitted models are handled
recursively.

## Usage

``` r
# S3 method for class 'latent'
vcov(fit, parameters = NULL, H = NULL)
```

## Arguments

- fit:

  A fitted object inheriting from class `"latent"`.

- parameters:

  Optional parameter specification identifying the parameters for which
  the covariance matrix should be returned. Labels must occur in
  `fit@modelInfo$transparameters_labels`. If `NULL`, the parameter
  blocks corresponding to the fitted model parameters are used.

- H:

  Optional Hessian or equivalent information matrix for the freely
  estimated parameters. If `NULL`, the Hessian is computed from `fit`.
  Supplying `H` allows alternative covariance estimators, such as
  sandwich estimators, to use the same transformation and sequential
  uncertainty machinery.

## Value

A list containing:

- `vcov`:

  Variance-covariance matrix for the selected parameters.

- `se`:

  Standard errors for the selected parameters.

- `jacob`:

  Jacobian matrix used to propagate covariance through parameter
  transformations.

- `H`:

  Hessian or information matrix used in the covariance calculation.

- `B`:

  For sequential models, the additional structural uncertainty component
  \\C^\top A C\\. For ordinary one-step models, an empty matrix.

- `newH`:

  For sequential models, the joint corrected precision matrix used to
  obtain the final covariance matrix.

## Details

Compute covariance matrices and standard errors for fitted latent
variable models, including uncertainty propagated across sequential
estimation steps.

For a one-step model, let \\H\\ denote the Hessian of the objective
function. The covariance matrix of the freely estimated parameters is
based on \\H^{-1}\\. Covariances of transformed parameters are
subsequently obtained with the delta method.

For a sequential model, suppose that measurement or nuisance parameters
\\\theta_M\\ are estimated in an earlier step and subsequently treated
as fixed while structural parameters \\\theta_S\\ are estimated. Let
\\A\\ denote the covariance matrix of the earlier parameter estimates,
\\H_2\\ the Hessian for the structural parameters, and \\C\\ the
measurement-structural cross-derivative matrix.

The propagated structural uncertainty contains the additional term
\$\$H_2^{-1} C^\top A C H_2^{-1}.\$\$

The method constructs the corresponding joint precision matrix \$\$
\begin{pmatrix} A^{-1} + C H_2^{-1} C^\top & C \\ C^\top & H_2
\end{pmatrix} \$\$ so that covariance between parameters estimated in
different stages is retained.

If the earlier fitted model itself contains a previous `"latent"`
object, [`vcov()`](https://rdrr.io/r/stats/vcov.html) is called
recursively before uncertainty is propagated to the current stage.
Consequently, chains of sequential plug-in estimation steps can be
handled recursively.

At most one nested `"latent"` object is currently supported at each
estimation level.

## References

Bakk, Z., and Kuha, J. Two-step estimation of models between latent
classes and external variables.

## See also

[`hessian.latent`](https://marcosjnez.github.io/latent/reference/hessian.latent.md),
[`jacobian.latent`](https://marcosjnez.github.io/latent/reference/jacobian.latent.md),
[`constraints_derivs.latent`](https://marcosjnez.github.io/latent/reference/constraints_derivs.latent.md),
[`se.llca`](https://marcosjnez.github.io/latent/reference/se.llca.md)
