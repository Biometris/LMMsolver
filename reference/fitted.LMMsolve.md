# Fitted values of an LMMsolve object.

Obtain the fitted values from a mixed model fitted using LMMSolve.

## Usage

``` r
# S3 method for class 'LMMsolve'
fitted(object, ...)
```

## Arguments

- object:

  an object of class LMMsolve

- ...:

  some methods for this generic require additional arguments. None are
  used in this method.

## Value

A vector of fitted values.

## Examples

``` r
## Fit model on oats data
data(oats.data)

## Fit simple model with only fixed effects.
LMM1 <- LMMsolve(fixed = yield ~ rep + gen,
                data = oats.data)

## Obtain fitted values.
fitted1 <- fitted(LMM1)
```
