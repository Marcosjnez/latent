# Inspect Latent Model Results

Inspect Latent Model Results

## Usage

``` r
latInspect(x, ..., sort = TRUE)

# S3 method for class 'latent'
latInspect(fit, what = "est", sort = TRUE)
```

## Arguments

- x:

  A fitted latent model or collection of latent class models.

- ...:

  Arguments passed to the inspection method, including `what`.

- sort:

  Logical. Apply factor ordering/sign orientation or class ordering to
  the requested output. Defaults to TRUE. FALSE returns native outputs.

- fit:

  A fitted latent model. This method also handles lrotate results, which
  inherit from `"latent"` or `"multistep"`.

- what:

  Character string identifying the requested component.

## Value

The requested model component, without altering `x`.

## Details

Sorting never modifies the supplied object, its optimizer coordinates,
or its model information. All uncertainty calculations use the native
fit. `what = "structures"` returns the estimate and template structures;
`"param"` and `"trans"` return the corresponding label templates.
`"labels"` has the layout of `"est"` with original labels. The actual
transformed-estimate slot is named `transformed_pars`. See
[`sort_latent()`](https://marcosjnez.github.io/latent/reference/sort_latent.md)
for the ordering and sign conventions. Unrotated requests keep native
source coordinates. Standard errors are never negated; inspected
covariance matrices receive both axis signs. A shared label with
opposite display signs can occur twice in an inspected covariance
matrix, without inventing new labels. Previously saved objects using
constructor-side sorting should be refitted.
