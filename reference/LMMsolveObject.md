# Fitted LMMsolve Object

An object of class `LMMsolve` returned by the LMMsolve function,
representing a fitted linear mixed model. Objects of this class have
methods for the generic functions coef, fitted, residuals, loglik and
deviance.

## Value

An object of class `LMMsolve` contains the following components:

- logL:

  The restricted log-likelihood at convergence

- sigma2e:

  The residual error

- tau2e:

  The estimated variance components

- EDdf:

  The effective dimensions

- varPar:

  The number of variance parameters for each variance component

- VarDf:

  The table with variance components

- theta:

  The precision parameters

- coefMME:

  A vector with all the estimated effects from mixed model equations

- ndxCoefficients:

  The indices of the coefficients with the names

- yhat:

  The fitted values

- residuals:

  The residuals

- nIter:

  The number of iterations for the mixed model to converge

- y:

  Response variable

- X:

  The design matrix for the fixed part of the mixed model

- Z:

  The design matrix for the random part of the mixed model

- lGinv:

  List with precision matrices for the random terms

- lRinv:

  List with precision matrices for the residual

- C:

  The mixed model coefficient matrix after last iteration

- cholC:

  The cholesky decomposition of coefficient matrix C

- constantREML:

  The REML constant

- dim:

  The dimensions for each of the fixed and random terms in the mixed
  model

- term.labels.f:

  The names of the fixed terms in the mixed model

- term.labels.r:

  The names of the random terms in the mixed model

- respVar:

  The name(s) of the response variable(s).

- splRes:

  An object with definition of spline argument

- deviance:

  The relative deviance

- family:

  An object of class family specifying the distribution and link
  function

- trace:

  A data.frame with the convergence sequence for the log likelihood and
  effective dimensions

.
