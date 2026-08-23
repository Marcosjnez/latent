# Rotate a fitted CFA model

`lrotate` rotates the factor loading matrices of a fitted `lcfa` object
using an orthogonal or oblique projection and a selected rotation
criterion.

## Usage

``` r
lrotate(fit, projection = "oblq", rotation = "oblimin",
        do.fit = TRUE, control = NULL, ...)
```

## Arguments

- fit:

  A fitted object inheriting from class `"lcfa"`.

- projection:

  Character string. Available projections are `"orth"`, `"oblq"`, and
  `"poblq"`.

- rotation:

  Character string identifying the rotation criterion.

- do.fit:

  Logical. If `TRUE`, fit the rotation. If `FALSE`, return the
  unrestricted rotation specification used for derivative calculations.

- control:

  List of optimization-control arguments.

- ...:

  Additional arguments required by the selected projection or rotation
  criterion.

## Value

An object of class `"multistep"`. The fitted `lcfa` object is stored in
`extra` as the preceding estimation step.

## Details

Let \\X\\ be the rotation matrix and let \\\Lambda_0\\, \\\Psi_0\\, and
\\\alpha_0\\ denote the unrotated factor loadings, factor covariance
matrix, and factor means. The rotated quantities are
\$\$\Lambda_r=\Lambda_0X^{-T},\$\$ \$\$\Psi_r=X^T\Psi_0X,\$\$ and
\$\$\alpha_r=X^T\alpha_0.\$\$ For an orthogonal projection,
\\X^{-T}=X\\. If \\\Psi_0\\ is a fixed identity matrix, \\\Psi_r\\ is
computed as \\X^TX\\.

## Examples

``` r
if (FALSE) { # \dontrun{
fit_cfa <- lcfa(data = HolzingerSwineford1939,
                model = model,
                std.lv = TRUE)
fit_rotation <- lrotate(fit = fit_cfa,
                        projection = "oblq",
                        rotation = "oblimin")
} # }
```
