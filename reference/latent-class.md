# Base Class for Fitted Latent Models

The `"latent"` class provides the common representation used by fitted
latent-variable models that share the generic derivative and covariance
infrastructure of latent.

## Slots

- `version`:

  Character string containing the package version used to fit the model.

- `call`:

  Original model-fitting call.

- `timing`:

  Numeric vector containing optimization timing information.

- `dataList`:

  List containing processed data and model-specific data objects.

- `modelInfo`:

  List containing parameter structures, transformations, estimators,
  manifolds, and optimizer controls.

- `Optim`:

  List containing raw optimization results.

- `parameters`:

  List containing fitted model parameters.

- `transformed_pars`:

  List containing transformed parameter estimates.

- `extra`:

  List reserved for additional model-specific information.

## See also

[`hessian`](https://marcosjnez.github.io/latent/reference/hessian.md),
[`jacobian`](https://marcosjnez.github.io/latent/reference/jacobian.md),
[`constraints_derivs`](https://marcosjnez.github.io/latent/reference/constraints_derivs.md),
[`vcov`](https://rdrr.io/r/stats/vcov.html)
