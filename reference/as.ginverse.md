# Construct a ginverse Object from Precision Matrices

Creates a `ginverse` object from a named list of precision (inverse
covariance) matrices. These matrices are typically used to specify the
inverse of covariance structures for random effects in `LMMsolve`.

## Usage

``` r
as.ginverse(precisionMatrices, levels = NULL, tol = 1e-10)
```

## Arguments

- precisionMatrices:

  A named list of square matrices. Each element must be a base `matrix`,
  an object inheriting from `Matrix`, or a `spam` object. Each element
  represents a precision matrix corresponding to a random effect.

- levels:

  An optional named list giving the levels corresponding to the rows and
  columns of the precision matrices. This is required for `spam`
  objects, which do not have row and column names. For `matrix` and
  `Matrix` objects, levels are obtained from the row names; if supplied,
  they are checked for consistency with the row and column names.

- tol:

  A numeric tolerance used for numerical stability (e.g. during
  inversion or eigenvalue truncation). Stored as an attribute of the
  resulting object.

## Value

An object of class `"ginverse"` (a named list) containing the supplied
precision matrices, with attributes `"levels"` and `"tol"`.

## Details

Each precision matrix must be square. For `matrix` and `Matrix` objects,
the row and column names define the corresponding levels. For `spam`
objects, which do not use row and column names, the corresponding levels
must be supplied through `levels`.

The function performs basic validation:

- `precisionMatrices` must be a named list.

- If supplied, `levels` must be a named list with matching names.

- Each precision matrix must be square.

- For `matrix` and `Matrix` objects, row and column names must be
  present and identical.

- For `spam` objects, `levels` must be supplied.

- Levels must be a character vector with length equal to the
  corresponding matrix dimension and contain no duplicates.

- If `levels` is supplied for a `matrix` or `Matrix` object, it must
  agree with its row and column names.

No reordering or alignment with the data is performed at this stage.
This is handled internally by `LMMsolve`.

## See also

[`LMMsolve`](https://biometris.github.io/LMMsolver/index.html/reference/LMMsolve.md)

## Examples

``` r
K <- diag(1, 5)
dimnames(K) <- list(as.character(1:5), as.character(1:5))

# Construct ginverse object from a matrix with names
g <- as.ginverse(list(id = K))
g
#> $id
#>   1 2 3 4 5
#> 1 1 0 0 0 0
#> 2 0 1 0 0 0
#> 3 0 0 1 0 0
#> 4 0 0 0 1 0
#> 5 0 0 0 0 1
#> 
#> attr(,"class")
#> [1] "ginverse" "list"    
#> attr(,"levels")
#> attr(,"levels")$id
#> [1] "1" "2" "3" "4" "5"
#> 
#> attr(,"tol")
#> [1] 1e-10

# A spam matrix requires levels to be supplied
# Kspam <- spam::as.spam(K)
# g <- as.ginverse(list(id = Kspam),
#                  levels = list(id = as.character(1:5)))
```
