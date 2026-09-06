# Rotate factor loading and covariance matrices

`lrotate` rotates the factor loading and factor covariance matrices
supplied directly or extracted from a fitted `lcfa` object using an
orthogonal or oblique projection and a selected rotation criterion.

## Usage

``` r
lrotate(fit = NULL, lambda = NULL, psi = NULL,
        projection = "oblq", rotation = "oblimin",
        se = TRUE, do.fit = TRUE, control = NULL, sort = TRUE, ...)
```

## Arguments

- fit:

  Optional fitted object inheriting from class `"lcfa"`.

- lambda:

  Optional loading matrix or list of loading matrices. This is an
  alternative to supplying `fit`.

- psi:

  Optional factor covariance matrix or list of factor covariance
  matrices corresponding to `lambda`. If omitted, identity matrices are
  used. This argument cannot be used together with `fit`.

- projection:

  Character string. Available projections are `"orth"`, `"oblq"`, and
  `"poblq"`.

- rotation:

  Character string identifying the rotation criterion.

- se:

  Logical. If `TRUE` and `fit` is supplied, propagate standard errors
  from the fitted `lcfa` model to the rotated parameters. Standard
  errors are not available when matrices are supplied directly.

- do.fit:

  Logical. If `TRUE`, fit the rotation. If `FALSE`, return the model
  specification. With a fitted `lcfa` input, the unrestricted
  specification used for derivative calculations is returned.

- control:

  List of optimization-control arguments.

- sort:

  Logical. By default, orient every factor so its largest absolute
  loading is positive and, except for target rotations, sort factors by
  decreasing variance-adjusted sums of squared loadings. With
  `rotation = "target"` or `"xtarget"`, factor order is always retained
  because the target defines the intended factor positions; signs are
  still oriented when `sort = TRUE`. See
  [`sort_factors()`](https://marcosjnez.github.io/latent/reference/sort_latent.md).
  FALSE retains both the fitted order and signs. The fitted criterion
  and its constraints retain their native coordinate system.

- ...:

  Additional arguments required by the selected projection or rotation
  criterion. If omitted (or NULL), `weight` defaults to `1-target` and
  `psiweight` to `1-psitarget`. Explicit weights, including zero
  matrices, are preserved. Group-specific lists are supported.

## Value

An object inheriting from class `"latent"`.

## Details

Exactly one of `fit` and `lambda` must be supplied. Let \\X\\ be the
rotation matrix and let \\\Lambda_0\\, \\\Psi_0\\, and \\\alpha_0\\
denote the unrotated factor loadings, factor covariance matrix, and
factor means. The rotated quantities are
\$\$\Lambda_r=\Lambda_0X^{-T},\$\$ \$\$\Psi_r=X^T\Psi_0X,\$\$ and
\$\$\alpha_r=X^T\alpha_0.\$\$ For an orthogonal projection,
\\X^{-T}=X\\. If \\\Psi_0\\ is a fixed identity matrix, \\\Psi_r\\ is
computed as \\X^TX\\.

With `projection = "poblq"`, either `constraints` or `oblique` must be
supplied through `...`. The former uses arbitrary structural
constraints, whereas the latter gives the sizes of consecutive oblique
blocks; any remaining factors form one orthogonal block. They cannot be
used together.

When `fit` is supplied, the returned object inherits from `"multistep"`
and the fitted `lcfa` object is stored in `extra`. When matrices are
supplied directly, the returned object inherits only from `"latent"`;
sampling uncertainty is not propagated because no fitted source model is
available.

## Examples

``` r
if (FALSE) { # \dontrun{
fit_cfa <- lcfa(data = HolzingerSwineford1939,
                model = model,
                std.lv = TRUE)
fit_rotation <- lrotate(fit = fit_cfa,
                        projection = "oblq",
                        rotation = "oblimin")

direct_rotation <- lrotate(lambda = lambda,
                           psi = psi,
                           projection = "oblq",
                           rotation = "oblimin")
} # }
```
