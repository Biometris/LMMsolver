# LMMsolver 1.0.4

-   Improved computation time for calculation of standard errors. Implementation in C++ and using the 'sparse inverse'. 
-   Row-wise Kronecker product for `spam` matrices implemented in C++. Important for tensor product P-splines with improved computation time and memory allocation.  

# LMMsolver 1.0.3

-   Improved computation time and memory allocation, especially important for big data with many observations (the number of rows in the data frame).
-   Replaced the default `model.matrix` function by `Matrix::sparse.model.matrix` to generate sparse design matrices.
-   In function `obtainSmoothTrend` the standard errors are only calculated if `includeIntercept=TRUE`. 
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
