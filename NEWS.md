# LMMsolver 1.0.10

-   Cyclic B-splines models added for `spl1D()` and `spl2D()` functions. 
-   Third order differences (`pord=3`) added for `splxD()` functions.
-   New argument `type = c("response", "link")` for `predict()` function.
-   bug fixed for GLMM models if weights are close to zero. 

# LMMsolver 1.0.9

-   Binomial response can now also be modelled as `fixed = cbind(failure, succes)`
-   Categorial response using `family = multinomial()`
-   Vignette updated, with separate section for GLMM.
-   doi-link added for `LMMsolver`. 
-   argument `offset` can be defined as numeric or (new) as column name in data frame. 
-   example added to `predict()` function.
-   problem with calculation of standard errors fixed, because of minor change in `spam`.
-   bug fixed related to convergence for GLMM.

# LMMsolver 1.0.8

-   Vignette has been rewritten, with a new introduction section. 
-   The function `predict.LMMsolve` added.
-   Extension of gam models, combining different `splxD()` is possible now.   
-   Correction of upper bound nominal effective dimension for large data sets.
-   new 2D example Sea Surface Temperature added.
-   Issue with product of two large matrices fixed. 
-   Improved efficiency initialization for large datasets.
-   Bug in `grpTheta` argument of `LMMsolve()` fixed. 
-   Deviance function changes, with extra argument `relative`, giving the relative conditional deviance as defined in McCullagh and Nelder. The default is `relative=TRUE`, for `relative=FALSE` it returns `-2*logLik(obj)`

# LMMsolver 1.0.7

-   Improved efficiency for models where the `residual` argument of `LMMsolve()` is used.
-   A data.frame `trace` with convergence sequence for log-likelihood and effective dimensions, added as extra output returned by `LMMsolve()`.
-   Bug in v1.0.6 for GLMM models fixed.
-   Coefficients for three way interactions with one factor and two non-factors are now labelled correctly.
-   Standard errors in function `obtainSmoothTrend()` for GLMM models are now calculated.

# LMMsolver 1.0.6

-   A new argument `grpTheta` for `LMMsolve()` to give components in the model the same penalty. 
-   The dependency package `sp` is replaced by `sf`.
-   A small bug for models with more than 10.000 observations and only a numeric variable in the random part of the model is fixed.
-   Weights are now checked for missing values after removing observations with missing values in response. This prevents spurious errors when both response and weight are missing.

# LMMsolver 1.0.5

-   Small bugs in assignment of names to fixed model coefficients when columns were dropped from the model are fixed.  
-   Calculation of standard errors for coefficients, with `coef(obj, se = TRUE)`.
-   Implementation of Generalized Linear Mixed Models (GLMM) with additional argument `family` in `LMMsolve` function.
-   Variance components and splines can be conditional on a factor. For variance components, this is implemented in the `cf(var, cond, level)` function. For 1D and 2D splines, additional arguments `cond` and `level` are added. 
-   Several small bugs fixed. 

# LMMsolver 1.0.4

-   Improved computation time for calculation of standard errors. Implementation in C++ and using the 'sparse inverse'. 
-   Row-wise Kronecker product for `spam` matrices implemented in C++. Important for tensor product P-splines with improved computation time and memory allocation.  

# LMMsolver 1.0.3

-   Improved computation time and memory allocation, especially important for big data with many observations (the number of rows in the data frame).
-   Replaced the default `model.matrix` function by `Matrix::sparse.model.matrix` to generate sparse design matrices.
-   In function `obtainSmoothTrend` the standard errors are only calculated if `includeIntercept = TRUE`. 
-   Several small bugs fixed.

# LMMsolver 1.0.2

-   First and second order derivatives are now calculated correctly.
-   Several small bugs fixed.
-   Updated tests to pass checks on macM1.

# LMMsolver 1.0.1

-   `weights` argument in LMMsolve function added
-   Function `obtainSmoothTrend` returns in addition to the predictions the standard errors.
-   Generalized Additive Model (GAM) added for one-dimensional splines, i.e. more `spl1D()` components can be added to the `spline` argument of LMMsolve function
-   Improved efficiency of calculating the sparse inverse using super-nodes.
-   Replaced the original P-splines penalty `D'D` with a scaled version which is far more stable if there are many knots.    
-   Several bugs fixed.

# LMMsolver 1.0.0

-   Initial CRAN version
