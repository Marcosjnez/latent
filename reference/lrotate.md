# Rotate factor loading and covariance matrices

`lrotate` rotates the factor loading and factor covariance matrices
supplied directly or extracted from a fitted `lcfa` object using an
orthogonal, oblique, or orthoblique projection and one or more rotation
criteria.

## Usage

``` r
lrotate(fit = NULL, lambda = NULL, psi = NULL,
        projection = "oblq", rotation = "oblimin",
        se = TRUE, do.fit = TRUE, control = NULL, ...)
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

  A criterion name, a character vector of criterion names, or a named
  list of parameter lists. Vector entries apply to the full rotated
  loading matrix by default and their losses are summed. For a named
  list, names identify the criteria (duplicate names are allowed), and
  each list contains that component's arguments, optionally including
  `items` and `factors` to select rows and columns. Missing selectors
  use all rows or columns. See Details for defaults and target
  subsetting.

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

- ...:

  Additional projection or rotation arguments shared by the components
  to which they apply. Component-specific arguments take precedence.
  Missing (or NULL) loading weights default to `1-target`, and
  covariance weights to `1-psitarget`. Explicit weights, including zero
  matrices, are preserved. Group-specific lists are supported.

## Value

An object inheriting from class `"latent"`.

## Details

Parameter labels are generated from their block names, just as in the
CFA parameter constructor. Single-group labels have no group suffix, for
example `X[1,1]` and `lambda_rotated[1,1]`. Multiple groups use the
corresponding group names to keep parameter labels distinct.

Fitted objects retain the estimated factor order and signs. Sorting is
available only through `latInspect(fit, sort = TRUE)` and does not
modify the fitted object or its standard-error calculations.

All components are optimized simultaneously over the same rotation
matrix. With row sets \\I_s\\ and factor sets \\J_s\\, the total
criterion is \$\$Q(\Lambda,\Psi)=\sum_s
Q_s(\Lambda\_{I_s,J_s},\Psi\_{J_s,J_s}).\$\$ Only `xtarget` uses the
selected factor covariance matrix; other criteria use the loading
submatrix only. Overlapping selections are allowed and their objective,
gradient, and Hessian contributions are added. There is no automatic
rescaling or averaging of the component criteria.

`items` and `factors` accept positive integer positions, matrix
row/column names, or logical vectors of the full corresponding length.
Duplicated or empty selections are rejected. Selector order is retained.
A target or weight matrix may have full loading-matrix dimensions, in
which case it is subset automatically, or the selected submatrix's
dimensions, in which case it is used in the supplied selector order.
When both sizes coincide, it is treated as a full matrix. For `xtarget`,
the same rule applies to the principal factor submatrices of `psitarget`
and `psiweight`; `items` affects only the loading part.

Defaults are resolved independently for each component and group.
`geomin` defaults to `epsilon = 0.01`; `oblimin` defaults to
`gamma = 0`. `alpha` is an alias for oblimin's `gamma`, not a component
weight. Conflicting values supplied at the same level are rejected. The
existing required arguments of other criteria remain required, including
`k` for `cf`, `epsilon` for `lclf`, and `w` for `xtarget`. A local NULL
requests the criterion default rather than inheriting a shared value.
Group-specific selectors or parameters may be supplied as lists with one
entry per group, in the input group order.

For target criteria the loss uses squared weighted residuals. To
multiply a target loss by \\c\\, multiply its weight matrix by
\\\sqrt{c}\\. For `xtarget`, scale both weight matrices to scale the
entire term; its `w` remains the relative covariance-target weight.
Sparse selections must jointly provide a sufficiently identified
rotation for standard errors; selecting a submatrix alone does not
guarantee this.

`dataList$rotation` is a printable criterion label, while
`dataList$rotation_spec` retains the vector or named-list specification.
`modelInfo$rotation_components` records the resolved row/factor
positions and parameters for each component in each group.

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

mixed_rotation <- lrotate(lambda = lambda,
                           rotation = c("oblimin", "target", "geomin"),
                           target = target, weight = 1-target)

subset_rotation <- lrotate(lambda = lambda,
                            rotation = list(
                              oblimin = list(alpha = 0, items = 1:5, factors = 1:3),
                              oblimin = list(alpha = 0.5, items = 6:15, factors = 4:5),
                              target = list(target = target, weight = 1-target),
                              geomin = list(items = 16:20)))
} # }
```
