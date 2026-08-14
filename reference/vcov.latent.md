# Variance-Covariance Matrix for Latent Models

`vcov.latent()` computes the variance-covariance matrix of freely
estimated or transformed parameters of a fitted object inheriting from
class `"latent"`.

For ordinary one-step models, covariance matrices are obtained from the
inverse Hessian, or from an alternative information matrix supplied
through `H`, and transformed to the requested parameterization using the
delta method.

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

  Optional information specification for the freely estimated
  parameters. If `NULL`, the ordinary Hessian is used. A numeric matrix
  can be supplied as the information matrix for the current model. A
  function can also be supplied; it must take a fitted model as its
  first argument and return an information matrix. Function-valued
  specifications are propagated recursively so that the information
  matrix is recomputed separately for each nested model.

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

  Information matrix used in the final covariance calculation.

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

For a one-step model, let \\P\\ denote the information matrix used for
inference. With ordinary Hessian-based inference, \\P = H\\, where \\H\\
is the Hessian of the objective function. The covariance matrix of the
freely estimated parameters is based on \\P^{-1}\\. Covariances of
transformed parameters are subsequently obtained with the delta method.

For a sequential model, suppose that measurement or nuisance parameters
\\\theta_M\\ are estimated in an earlier step and subsequently treated
as fixed while structural parameters \\\theta_S\\ are estimated. Let
\\A\\ denote the covariance matrix of the earlier parameter estimates,
\\H_2\\ the ordinary Hessian for the structural parameters, \\P_2\\ the
information matrix whose inverse gives the structural covariance when
the earlier parameters are treated as known, and \\C\\ the
measurement-structural cross-derivative matrix.

The propagated structural uncertainty contains the additional term
\$\$H_2^{-1} C^\top A C H_2^{-1}.\$\$

The corresponding joint covariance matrix is \$\$ \begin{pmatrix} A & -A
C H_2^{-1} \\ -H_2^{-1} C^\top A & P_2^{-1} + H_2^{-1} C^\top A C
H_2^{-1} \end{pmatrix}. \$\$

Rather than constructing and inverting this covariance matrix directly,
the method constructs its inverse and passes it through the ordinary
covariance machinery. When \\P_2 = H_2\\, this reduces to the standard
sequential pseudo-likelihood correction.

If `H` is supplied as a function and the earlier fitted model itself
contains a previous `"latent"` object, the same function is passed
recursively and evaluated independently for each nested fit. This
allows, for example, robust information matrices to be used at every
estimation stage.

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
