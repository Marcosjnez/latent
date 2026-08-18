# Standard Errors for Multistage Latent Models

Compute standard errors for a multistage latent model by propagating the
joint variance-covariance matrix of parameters estimated in previous
stages.

## Usage

``` r
# S3 method for class 'multistage'
se(fit, type = "information", parameters = NULL, digits = 4L, ...)
```

## Arguments

- fit:

  A fitted object inheriting from class `"multistage"`.

- type:

  Character string specifying the covariance estimator used at each
  estimation stage. Available options are `"information"` and
  `"robust"`.

- parameters:

  Optional parameter specification identifying the parameters or
  transformed parameters for which standard errors should be returned.
  If `NULL`, the freely estimated parameters of the top-level model are
  used.

- digits:

  Non-negative integer specifying the number of decimal places used in
  the formatted parameter tables. If `NULL`, values are not rounded.

- ...:

  Additional arguments reserved for future extensions.

## Value

A list containing the formatted estimates and standard errors,
variance-covariance matrix, Hessian of the top-level estimated
parameters, final propagation matrix `B`, transformation Jacobian, and
the final joint variance-covariance matrix before propagation to
transformed parameters.

## Details

Previous fitted models are obtained from the `extra` slot. The sequence
of models is ordered from the ultimate ancestor to the top-level model
and processed iteratively.

Let \\A\\ be the joint variance-covariance matrix of the parameters
fixed from previous stages, \\H_2\\ the Hessian of the parameters
estimated in the current stage, \\C\\ the cross-Hessian between
previous-stage and current-stage parameters obtained from the
corresponding unrestricted model, and \\V_2\\ the conditional
variance-covariance matrix of the current-stage estimator. At each stage
the joint covariance matrix is enlarged as \$\$ V = \begin{pmatrix} A &
-A C H_2^{-1} \\ -H_2^{-1} C^\top A & V_2 + H_2^{-1} C^\top A C H_2^{-1}
\end{pmatrix}. \$\$

If `type = "robust"`, robust covariance matrices are used for \\V_2\\
whenever a class-specific robust method is available. The ordinary
Hessian \\H_2\\ is still used for propagation of uncertainty between
estimation stages.

The top-level unrestricted model must contain all parameter labels
carried forward from previous stages. This is required to obtain the
complete cross-Hessian between the previously estimated parameters and
the parameters estimated in the new stage.
