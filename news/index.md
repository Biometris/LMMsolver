# Changelog

## LMMsolver 1.0.12

- First derivatives for `predict` using `deriv` argument now also
  implemented for `spl2D` and `spl3D`.
- function
  [`effDim()`](https://biometris.github.io/LMMsolver/index.html/reference/effDim.md)
  added to get data.frame with effective dimensions.
- In the vignette, an example added how the generalized heritability can
  be calculated.
- Improved code coverage \> 95%.
- Data sets `barley.uniformity.trial` and `oats.data` added.
- All data included in the package that are needed for tests.

## LMMsolver 1.0.11

CRAN release: 2025-08-20

- New function
  [`mLogLik()`](https://biometris.github.io/LMMsolver/index.html/reference/mLogLik.md)
  for the calculations of the log-likelihood and first derivatives as
  function of precision parameters `theta`.
- A new argument `deriv` added to
  [`predict()`](https://rdrr.io/r/stats/predict.html) to calculate the
  first derivatives for
  [`spl1D()`](https://biometris.github.io/LMMsolver/index.html/reference/spl1D.md)
  functions.
- Two examples in vignette updated with predictions of derivatives and
  corresponding standard errors.
- bug fixed for `theta` argument of
  [`LMMsolve()`](https://biometris.github.io/LMMsolver/index.html/reference/LMMsolve.md).

## LMMsolver 1.0.10

CRAN release: 2025-05-14

- Cyclic B-splines models added for
  [`spl1D()`](https://biometris.github.io/LMMsolver/index.html/reference/spl1D.md)
  and
  [`spl2D()`](https://biometris.github.io/LMMsolver/index.html/reference/spl1D.md)
  functions.
- Third order differences (`pord=3`) added for `splxD()` functions.
- New argument `type = c("response", "link")` for
  [`predict()`](https://rdrr.io/r/stats/predict.html) function.
- bug fixed for GLMM models if weights are close to zero.

## LMMsolver 1.0.9

CRAN release: 2025-01-14

- Binomial response can now also be modelled as
  `fixed = cbind(failure, succes)`
- Categorical response using `family = multinomial()`
- Vignette updated, with separate section for GLMM.
- doi-link added for `LMMsolver`.
- argument `offset` can be defined as numeric or (new) as column name in
  data frame.
- example added to [`predict()`](https://rdrr.io/r/stats/predict.html)
  function.
- problem with calculation of standard errors fixed, because of minor
  change in `spam`.
- bug fixed related to convergence for GLMM.

## LMMsolver 1.0.8

CRAN release: 2024-08-26

- Vignette has been rewritten, with a new introduction section.
- The function `predict.LMMsolve` added.
- Extension of gam models, combining different `splxD()` is possible
  now.  
- Correction of upper bound nominal effective dimension for large data
  sets.
- new 2D example Sea Surface Temperature added.
- Issue with product of two large matrices fixed.
- Improved efficiency initialization for large datasets.
- Bug in `grpTheta` argument of
  [`LMMsolve()`](https://biometris.github.io/LMMsolver/index.html/reference/LMMsolve.md)
  fixed.
- Deviance function changes, with extra argument `relative`, giving the
  relative conditional deviance as defined in McCullagh and Nelder. The
  default is `relative=TRUE`, for `relative=FALSE` it returns
  `-2*logLik(obj)`

## LMMsolver 1.0.7

CRAN release: 2024-04-16

- Improved efficiency for models where the `residual` argument of
  [`LMMsolve()`](https://biometris.github.io/LMMsolver/index.html/reference/LMMsolve.md)
  is used.
- A data.frame `trace` with convergence sequence for log-likelihood and
  effective dimensions, added as extra output returned by
  [`LMMsolve()`](https://biometris.github.io/LMMsolver/index.html/reference/LMMsolve.md).
- Bug in v1.0.6 for GLMM models fixed.
- Coefficients for three way interactions with one factor and two
  non-factors are now labelled correctly.
- Standard errors in function
  [`obtainSmoothTrend()`](https://biometris.github.io/LMMsolver/index.html/reference/obtainSmoothTrend.md)
  for GLMM models are now calculated.

## LMMsolver 1.0.6

CRAN release: 2023-11-27

- A new argument `grpTheta` for
  [`LMMsolve()`](https://biometris.github.io/LMMsolver/index.html/reference/LMMsolve.md)
  to give components in the model the same penalty.
- The dependency package `sp` is replaced by `sf`.
- A small bug for models with more than 10.000 observations and only a
  numeric variable in the random part of the model is fixed.
- Weights are now checked for missing values after removing observations
  with missing values in response. This prevents spurious errors when
  both response and weight are missing.

## LMMsolver 1.0.5

CRAN release: 2023-04-14

- Small bugs in assignment of names to fixed model coefficients when
  columns were dropped from the model are fixed.  
- Calculation of standard errors for coefficients, with
  `coef(obj, se = TRUE)`.
- Implementation of Generalized Linear Mixed Models (GLMM) with
  additional argument `family` in `LMMsolve` function.
- Variance components and splines can be conditional on a factor. For
  variance components, this is implemented in the `cf(var, cond, level)`
  function. For 1D and 2D splines, additional arguments `cond` and
  `level` are added.
- Several small bugs fixed.

## LMMsolver 1.0.4

CRAN release: 2022-12-15

- Improved computation time for calculation of standard errors.
  Implementation in C++ and using the ‘sparse inverse’.
- Row-wise Kronecker product for `spam` matrices implemented in C++.
  Important for tensor product P-splines with improved computation time
  and memory allocation.

## LMMsolver 1.0.3

CRAN release: 2022-08-19

- Improved computation time and memory allocation, especially important
  for big data with many observations (the number of rows in the data
  frame).
- Replaced the default `model.matrix` function by
  [`Matrix::sparse.model.matrix`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html)
  to generate sparse design matrices.
- In function `obtainSmoothTrend` the standard errors are only
  calculated if `includeIntercept = TRUE`.
- Several small bugs fixed.

## LMMsolver 1.0.2

CRAN release: 2022-04-21

- First and second order derivatives are now calculated correctly.
- Several small bugs fixed.
- Updated tests to pass checks on macM1.

## LMMsolver 1.0.1

CRAN release: 2022-03-28

- `weights` argument in LMMsolve function added
- Function `obtainSmoothTrend` returns in addition to the predictions
  the standard errors.
- Generalized Additive Model (GAM) added for one-dimensional splines,
  i.e. more
  [`spl1D()`](https://biometris.github.io/LMMsolver/index.html/reference/spl1D.md)
  components can be added to the `spline` argument of LMMsolve function
- Improved efficiency of calculating the sparse inverse using
  super-nodes.
- Replaced the original P-splines penalty `D'D` with a scaled version
  which is far more stable if there are many knots.  
- Several bugs fixed.

## LMMsolver 1.0.0

CRAN release: 2021-11-02

- Initial CRAN version
