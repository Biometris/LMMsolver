# Solve Linear Mixed Models

Solve Linear Mixed Models using REML.

## Usage

``` r
LMMsolve(
  fixed,
  random = NULL,
  spline = NULL,
  group = NULL,
  ginverse = NULL,
  weights = NULL,
  data,
  residual = NULL,
  family = gaussian(),
  offset = 0,
  tolerance = 1e-06,
  trace = FALSE,
  maxit = 250,
  theta = NULL,
  grpTheta = NULL
)
```

## Arguments

- fixed:

  A formula for the fixed part of the model. Should be of the form
  "response ~ pred"

- random:

  A formula for the random part of the model. Should be of the form "~
  pred".

- spline:

  A formula for the spline part of the model. Should be of the form "~
  spl1D()", ~ spl2D()" or "~spl3D()". Generalized Additive Models (GAMs)
  can also be used, for example "~ spl1D() + spl2D()"

- group:

  A named list where each component is a numeric vector specifying
  contiguous fields in data that are to be considered as a single term.

- ginverse:

  A named list with each component a symmetric matrix, the precision
  matrix of a corresponding random term in the model. The row and column
  order of the precision matrices should match the order of the levels
  of the corresponding factor in the data.

- weights:

  A character string identifying the column of data to use as relative
  weights in the fit. Default value NULL, weights are all equal to one.

- data:

  A data.frame containing the modeling data.

- residual:

  A formula for the residual part of the model. Should be of the form "~
  pred".

- family:

  An object of class `family` or `familyLMMsolver` specifying the
  distribution and link function. See class
  [`family`](https://rdrr.io/r/stats/family.html) and and
  [`multinomial`](https://biometris.github.io/LMMsolver/index.html/reference/multinomial.md)
  for details.

- offset:

  An a priori known component to be included in the linear predictor
  during fitting. `Offset` be a numeric vector, or a character string
  identifying the column of data. Default `offset = 0`.

- tolerance:

  A numerical value. The convergence tolerance for the modified
  Henderson algorithm to estimate the variance components.

- trace:

  Should the progress of the algorithm be printed? Default
  `trace = FALSE`.

- maxit:

  A numerical value. The maximum number of iterations for the algorithm.
  Default `maxit = 250`.

- theta:

  initial values for penalty or precision parameters. Default `NULL`,
  all precision parameters set equal to 1.

- grpTheta:

  a vector to give components the same penalty. Default `NULL`, all
  components have a separate penalty.

## Value

An object of class `LMMsolve` representing the fitted model. See
[`LMMsolveObject`](https://biometris.github.io/LMMsolver/index.html/reference/LMMsolveObject.md)
for a full description of the components in this object.

## Details

A Linear Mixed Model (LMM) has the form \$\$y = X \beta + Z u + e, u
\sim N(0,G), e \sim N(0,R)\$\$ where \\y\\ is a vector of observations,
\\\beta\\ is a vector with the fixed effects, \\u\\ is a vector with the
random effects, and \\e\\ a vector of random residuals. \\X\\ and \\Z\\
are design matrices.

LMMsolve can fit models where the matrices \\G^{-1}\\ and \\R^{-1}\\ are
a linear combination of precision matrices \\Q\_{G,i}\\ and
\\Q\_{R,i}\\: \$\$G^{-1} = \sum\_{i} \psi_i Q\_{G,i} \\, R^{-1} =
\sum\_{i} \phi_i Q\_{R,i}\$\$ where the precision parameters \\\psi_i\\
and \\\phi_i\\ are estimated using REML. For most standard mixed models
\\1/{\psi_i}\\ are the variance components and \\1/{\phi_i}\\ the
residual variances. We use a formulation in terms of precision
parameters to allow for non-standard mixed models using tensor product
splines.

## See also

[`LMMsolveObject`](https://biometris.github.io/LMMsolver/index.html/reference/LMMsolveObject.md),
[`spl1D`](https://biometris.github.io/LMMsolver/index.html/reference/spl1D.md),
[`spl2D`](https://biometris.github.io/LMMsolver/index.html/reference/spl1D.md),
[`spl3D`](https://biometris.github.io/LMMsolver/index.html/reference/spl1D.md)

## Examples

``` r
## Fit models on oats.data
data(oats.data)

## Fit simple model with only fixed effects.
LMM1 <- LMMsolve(fixed = yield ~ rep + gen,
                data = oats.data)

## Fit the same model with genotype as random effect.
LMM1_rand <- LMMsolve(fixed = yield ~ rep,
                     random = ~gen,
                     data = oats.data)

## Fit the model with a 1-dimensional spline at the plot level.
LMM1_spline <- LMMsolve(fixed = yield ~ rep + gen,
                       spline = ~spl1D(x = plot, nseg = 20),
                       data = oats.data)

## Fit models on multipop data included in the package.
data(multipop)

## The residual variances for the two populations can be different.
## Allow for heterogeneous residual variances using the residual argument.
LMM2 <- LMMsolve(fixed = pheno ~ cross,
                residual = ~cross,
                data = multipop)

## QTL-probabilities are defined by the columns pA, pB, pC.
## They can be included in the random part of the model by specifying the
## group argument and using grp() in the random part.

# Define groups by specifying columns in data corresponding to groups in a list.
# Name used in grp() should match names specified in list.
lGrp <- list(QTL = 3:5)
LMM2_group <- LMMsolve(fixed = pheno ~ cross,
                      group = lGrp,
                      random = ~grp(QTL),
                      residual = ~cross,
                      data = multipop)
```
