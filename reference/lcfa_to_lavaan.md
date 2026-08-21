# Convert a Fitted lcfa Object to lavaan

Convert a Fitted lcfa Object to lavaan

## Usage

``` r
lcfa_to_lavaan(object, ...)
```

## Arguments

- object:

  A fitted object inheriting from class `"lcfa"`.

- ...:

  Additional arguments used only when rebuilding the lavaan scaffold.

## Value

A `"lavaan"` object populated with estimates and covariance information
from the fitted `"lcfa"` object.

## Details

This converter uses a lavaan scaffold to preserve lavaan's
parameter-table and model-matrix conventions. It depends on non-exported
lavaan constructors, so compatibility is checked at run time and a
descriptive error is returned when the installed lavaan version is
incompatible.
