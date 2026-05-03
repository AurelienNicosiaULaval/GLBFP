# Predict from GLBFP fit objects

Prediction helper for fitted GLBFP objects.

## Usage

``` r
# S3 method for class 'glbfp_fit'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  A fitted object of class `"glbfp_fit"` or `"glbfp_grid"`.

- newdata:

  Optional matrix/data frame with points where prediction is requested.
  For `"glbfp_fit"`, `newdata` is not supported.

- ...:

  Additional arguments (unused).

## Value

Numeric vector of predicted densities.
