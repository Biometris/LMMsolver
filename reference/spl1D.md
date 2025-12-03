# Fit P-splines

Fit multi dimensional P-splines using sparse implementation.

## Usage

``` r
spl1D(
  x,
  nseg,
  pord = 2,
  degree = 3,
  cyclic = FALSE,
  scaleX = TRUE,
  xlim = range(x),
  cond = NULL,
  level = NULL
)

spl2D(
  x1,
  x2,
  nseg,
  pord = 2,
  degree = 3,
  cyclic = c(FALSE, FALSE),
  scaleX = TRUE,
  x1lim = range(x1),
  x2lim = range(x2),
  cond = NULL,
  level = NULL
)

spl3D(
  x1,
  x2,
  x3,
  nseg,
  pord = 2,
  degree = 3,
  scaleX = TRUE,
  x1lim = range(x1),
  x2lim = range(x2),
  x3lim = range(x3)
)
```

## Arguments

- x, x1, x2, x3:

  The variables in the data containing the values of the `x` covariates.

- nseg:

  The number of segments

- pord:

  The order of penalty, default `pord = 2`

- degree:

  The degree of B-spline basis, default `degree = 3`

- cyclic:

  Cyclic or linear B-splines; default `cyclic=FALSE`

- scaleX:

  Should the fixed effects be scaled.

- xlim, x1lim, x2lim, x3lim:

  A numerical vector of length 2 containing the domain of the
  corresponding x covariate where the knots should be placed. Default
  set to `NULL`, when the covariate range will be used.

- cond:

  Conditional factor: splines are defined conditional on the level.
  Default `NULL`.

- level:

  The level of the conditional factor. Default `NULL`.

## Value

A list with the following elements:

- `X` - design matrix for fixed effect. The intercept is not included.

- `Z` - design matrix for random effect.

- `lGinv` - a list of precision matrices

- `knots` - a list of vectors with knot positions

- `dim.f` - the dimensions of the fixed effect.

- `dim.r` - the dimensions of the random effect.

- `term.labels.f` - the labels for the fixed effect terms.

- `term.labels.r` - the labels for the random effect terms.

- `x` - a list of vectors for the spline variables.

- `pord` - the order of the penalty.

- `degree` - the degree of the B-spline basis.

- `scaleX` - logical indicating if the fixed effects are scaled.

- `EDnom` - the nominal effective dimensions.

## Functions

- `spl2D()`: 2-dimensional splines

- `spl3D()`: 3-dimensional splines

## See also

[`LMMsolve`](https://biometris.github.io/LMMsolver/index.html/reference/LMMsolve.md)

## Examples

``` r
## Fit model on oats data
data(oats.data)

## Fit a model with a 1-dimensional spline at the plot level.
LMM1_spline <- LMMsolve(fixed = yield ~ rep + gen,
                       spline = ~spl1D(x = plot, nseg = 20),
                       data = oats.data)

summary(LMM1_spline)
#> Table with effective dimensions and penalties: 
#> 
#>         Term Effective Model Nominal Ratio Penalty
#>  (Intercept)      1.00     1       1  1.00    0.00
#>          rep      2.00     2       2  1.00    0.00
#>          gen     23.00    23      23  1.00    0.00
#>    lin(plot)      1.00     1       1  1.00    0.00
#>      s(plot)      3.99    23      21  0.19 3310.21
#>     residual     41.01    72      45  0.91   13.21
#> 
#>  Total Effective Dimension: 72 

## Fit model on US precipitation data from spam package.
data(USprecip, package = "spam")

## Only use observed data
USprecip <- as.data.frame(USprecip)
USprecip <- USprecip[USprecip$infill == 1, ]

## Fit a model with a 2-dimensional P-spline.
LMM2_spline <- LMMsolve(fixed = anomaly ~ 1,
                       spline = ~spl2D(x1 = lon, x2 = lat, nseg = c(41, 41)),
                       data = USprecip)

summary(LMM2_spline)
#> Table with effective dimensions and penalties: 
#> 
#>           Term Effective Model Nominal Ratio Penalty
#>    (Intercept)      1.00     1       1  1.00    0.00
#>  lin(lon, lat)      3.00     3       3  1.00    0.00
#>         s(lon)    302.60  1936    1932  0.16    0.26
#>         s(lat)    409.09  1936    1932  0.21    0.08
#>       residual   5190.31  5906    5902  0.88   13.53
#> 
#>  Total Effective Dimension: 5906 
```
