# Create a data frame for use in making predictions.

Constructs a grid of values spanning the ranges of the spline covariates
stored in an `LMMsolver` object.

## Usage

``` r
makeGrid(object, grid)
```

## Arguments

- object:

  An `LMMsolver` object.

- grid:

  A numeric vector specifying the number of grid points for each spline
  dimension. Its length must equal the number of spline variables.

## Value

A data frame containing all combinations of grid values for the spline
covariates.

## Details

For each spline variable, equally spaced values are generated between
the minimum and maximum values of the B-splines. The Cartesian product
of these sequences is returned using
[`expand.grid`](https://rdrr.io/r/base/expand.grid.html).

## Examples

``` r
if (FALSE) { # \dontrun{
## Create a 200 x 300 grid for a two-dimensional spline term
grd <- makeGrid(fit, grid = c(200, 300))

head(grd)
} # }
```
