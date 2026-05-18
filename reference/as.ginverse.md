# Construct a ginverse Object from Precision Matrices

Creates a `ginverse` object from a named list of precision (inverse
covariance) matrices. These matrices are typically used to specify the
inverse of covariance structures for random effects in `LMMsolve`.

## Usage

``` r
as.ginverse(precisionMatrices, tol = 1e-10)
```

## Arguments

- precisionMatrices:

  A named list of square matrices (base `matrix` or objects inheriting
  from `Matrix`). Each element represents a precision matrix
  corresponding to a random effect. The names of the list must match the
  variable names used in the `random` argument of `LMMsolve`.

- tol:

  A numeric tolerance used for numerical stability (e.g. during
  inversion or eigenvalue truncation). Stored as an attribute of the
  resulting object.

## Value

An object of class `"ginverse"` (a named list) containing the supplied
precision matrices, with attribute `"tol"`.

## Details

Each matrix must have identical row and column names corresponding to
the levels of the associated random effect. Alignment with the data is
checked internally within `LMMsolve`.

The function performs basic validation:

- `precisionMatrices` must be a named list.

- Each matrix must be square with identical row and column names.

- Row and column names are used later to align matrices with factor
  levels in the data.

No reordering or alignment with the data is performed at this stage.
This is handled internally by `LMMsolve`.

## See also

[`LMMsolve`](https://biometris.github.io/LMMsolver/index.html/reference/LMMsolve.md)

## Examples

``` r
library(Matrix)

# Create a simple precision matrix
K <- Diagonal(5)
rownames(K) <- colnames(K) <- as.character(1:5)

# Construct ginverse object
g <- as.ginverse(list(id = K))

g
#> $id
#> 5 x 5 diagonal matrix of class "ddiMatrix"
#>   1 2 3 4 5
#> 1 1 . . . .
#> 2 . 1 . . .
#> 3 . . 1 . .
#> 4 . . . 1 .
#> 5 . . . . 1
#> 
#> attr(,"class")
#> [1] "ginverse" "list"    
#> attr(,"tol")
#> [1] 1e-10
```
