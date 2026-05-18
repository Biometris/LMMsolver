# Generalized heritability of a random term

Computes the generalized heritability of a random-effect term from a
fitted linear mixed model. By default, a single scalar heritability
value is returned. Optionally, a spectral decomposition is provided that
reveals how genetic signal is distributed across estimable genetic
directions.

## Usage

``` r
getHeritability(obj, geno.term, type = c("scalar", "spectral"), tol = 1e-08)
```

## Arguments

- obj:

  An object of class `"LMMsolve"`.

- geno.term:

  A character string giving the name of the genetic random-effect term.

- type:

  Character string specifying the output: `"scalar"` (default) returns a
  single numeric heritability value; `"spectral"` returns a data frame
  with the spectral decomposition.

- tol:

  Numerical tolerance used to determine the estimable genetic space.

## Value

If `type = "scalar"`, a numeric value giving the generalized
heritability. If `type = "spectral"`, a data frame with columns:

- component:

  Index of the spectral component

- lambda:

  Canonical heritability for the component

- w:

  Weight of the component (genetic capacity)

- h2_comp:

  Contribution of the component to total heritability

In the spectral case, the scalar generalized heritability is also
available as the attribute `"h2_G"`.

## Details

Generalized heritability is defined as the proportion of estimable
genetic signal retained by the design relative to the available genetic
capacity, accounting for the genetic covariance structure.

For independent genotypes, this definition reduces to classical
generalized heritability measures based on effective dimension (Cullis,
Oakey, Rodríguez-Álvarez). When genotypes are correlated, the spectral
decomposition reveals anisotropy in information retention across genetic
directions.
