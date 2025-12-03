# Coefficients from the mixed model equations of an LMMsolve object.

Obtain the coefficients from the mixed model equations of an LMMsolve
object.

## Usage

``` r
# S3 method for class 'LMMsolve'
coef(object, se = FALSE, ...)
```

## Arguments

- object:

  an object of class LMMsolve

- se:

  calculate standard errors, default FALSE.

- ...:

  some methods for this generic require additional arguments. None are
  used in this method.

## Value

A list of vectors, containing the estimated effects for each fixed
effect and the predictions for each random effect in the defined linear
mixed model.

## Examples

``` r
## Fit model on oats data
data(oats.data)

## Fit simple model with only fixed effects.
LMM1 <- LMMsolve(fixed = yield ~ rep + gen,
                data = oats.data)

## Obtain coefficients.
coefs1 <- coef(LMM1)

## Obtain coefficients with standard errors.
coefs2 <- coef(LMM1, se = TRUE)
```
