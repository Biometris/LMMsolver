# Construct design matrix for B-Splines

Construct design matrix for B-Splines.

## Usage

``` r
Bsplines(knots, x, deriv = 0)
```

## Arguments

- knots:

  A numerical vector of knot positions.

- x:

  a numeric vector of values at which to evaluate the B-spline functions
  or derivatives.

- deriv:

  A numerical value. The derivative of the given order is evaluated at
  the x positions.
