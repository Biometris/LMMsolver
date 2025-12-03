# Family Object for Multinomial Model

The Multinomial model is not part of the standard family. The
implementation is based on Chapter 6 in Fahrmeir et al. (2013).

## Usage

``` r
multinomial()
```

## Value

An object of class `familyLMMsolver` with the following components:

- family:

  character string with the family name.

- linkfun:

  the link function.

- linkinv:

  the inverse of the link function.

- dev.resids:

  function giving the deviance for each observation as a function of (y,
  mu, wt)

## References

Fahrmeir, Ludwig, Thomas Kneib, Stefan Lang, Brian Marx, Regression
models. Springer Berlin Heidelberg, 2013.
