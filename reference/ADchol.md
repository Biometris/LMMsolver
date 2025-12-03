# construct object for Automated Differentiation Cholesky decomposition

Construct object for reverse Automated Differentiation of Cholesky
decomposition, with as input a list of semi-positive symmetric sparse
matrices \\P_i\\, each of dimension \\q \times q\\. The function
`ADchol` calculates the matrix \\C\\, the sum the precision matrices
\\P_i\\: \\C = \sum\_{i} P_i\\. Next, it calculates the Cholesky
Decomposition using the multiple minimum degree (MMD) algorithm of the
`spam` package.

## Usage

``` r
ADchol(lP)
```

## Arguments

- lP:

  a list of symmetric matrices of class spam, each of dimension \\q
  \times q\\, and with sum of the matrices assumed to be positive
  definite.

## Value

An object of class `ADchol`. This object is used to calculate the
partial partial derivatives of \\log\|C\|\\ in an efficient way.

## References

Furrer, R., & Sain, S. R. (2010). spam: A sparse matrix R package with
emphasis on MCMC methods for Gaussian Markov random fields. Journal of
Statistical Software, 36, 1-25.
