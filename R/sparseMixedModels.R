
# mixed model coefficient matrix...

# calcalate the weighted sum of list of matrices.
linearSum  <- function(theta,matrixList)
{
  C = Reduce('+', Map(function(x,y) x*y, theta, matrixList))
}

# calculate log determinant of matrix M,
logDet <- function(cholM) { as.double(determinant(cholM)$modulus) }

# calculate quadratic form xAx (efficient):
quadForm <- function(x, A, y=x) { sum(x*(A %*% y))}

REMLlogL <- function(ADcholRinv, ADcholGinv, ADcholC,phi,psi,theta, yPy)
{
  # calculute logLikelihood...
  logdetR = -logdet(ADcholRinv, phi)
  logdetG =  ifelse(!is.null(ADcholGinv), -logdet(ADcholGinv, psi), 0)
  logdetC =  logdet(ADcholC, theta)

  # See e.g. Smith 1995..
  logL = -0.5*(logdetR + logdetG + logdetC + yPy)
}

calcSumSquares <- function(lRinv,Q, r, a, Nvarcomp)
{
  SSr = sapply(lRinv, FUN= function(X) { quadForm(r,X) })
  if (Nvarcomp > 0)
  {
    SSa = sapply(Q, FUN = function(X) { quadForm(a,X) })
    SS_all = c(SSr, SSa)
  } else {
    SS_all = SSr
  }
}

solveMME <- function(cholC, listC, lWtRinvY, phi, theta)
{
  C = linearSum(theta, listC)
  cholC <-update(cholC, C)
  WtRinvy = as.vector(linearSum(phi, lWtRinvY))
  a = backsolve.spam(cholC, forwardsolve.spam(cholC, WtRinvy))
  a
}

sparseMixedModels <- function(y, X, Z, lGinv, lRinv, maxiter=100, eps=1.0e-6, display=FALSE ,
                              monitor = FALSE) {
  Ntot = length(y)
  p = ncol(X)
  q = ncol(Z)
  Nres = length(lRinv)
  Nvarcomp = length(lGinv)
  NvarcompTot = Nres + Nvarcomp
  dimMME = p + q

  W = as.spam(cbind(X,Z))
  Wt = t(W)

  lWtRinvW = lapply(lRinv, FUN = function(x) { Wt %*% x %*% W})
  lWtRinvY = lapply(lRinv, FUN = function(x) { Wt %*% (x %*% y)})
  # extend lGinv with extra zero's for fixed effect...
  lQ = lapply(lGinv, FUN = function(x)
                          { zero = spam(0, ncol=p, nrow=p)
                            bdiag.spam(zero,x)
                          })
  listC = c(lWtRinvW, lQ)

  # remove some extra zero's....
  lWtRinvW <- lapply(lWtRinvW, FUN = function(X) {cleanup(X)} )
  lQ       <- lapply(lQ,       FUN = function(X) {cleanup(X)} )
  lGinv    <- lapply(lGinv,  FUN = function(X) {cleanup(X)} )
  listC    <- lapply(listC,  FUN=function(X) {cleanup(X)})

  phi = rep(1.0, Nres)
  psi = rep(1.0, Nvarcomp)
  theta = c(psi, phi)
  #cholRinv = chol(Rinv)
  #if (Nvarcomp > 0.0) {cholGinv = chol(linearSum(psi, lGinv))}
  cholC = chol(linearSum(theta, listC))

  if (display) {
    C = linearSum(theta, listC)
    display(C)
    L = t(as.spam(cholC))
    display(L)
  }

  # make ADchol for Rinv, Ginv and C:
  ADcholRinv = ADchol(lRinv)
  if (Nvarcomp >0) {
    ADcholGinv = ADchol(lGinv)
  } else {
    ADcholGinv = NULL
  }
  ADcholC = ADchol(listC)

  logLprev <- 1.0e20
  if (monitor) { cat("iter logLik\n") }

  for (it in 1:maxiter)
  {
    if (Nvarcomp > 0)
    {
      phi = theta[c(1:length(phi))]
      psi = theta[-c(1:length(phi))]
    } else {
      phi = theta
      psi = NULL
    }

    # calculate effective dimensions....
    EDmax_phi = phi*dlogdet(ADcholRinv, phi)
    if (!is.null(ADcholGinv))
    {
      EDmax_psi = psi*dlogdet(ADcholGinv, psi)
    } else {
      EDmax_psi = NULL
    }
    EDmax = c(EDmax_phi, EDmax_psi)
    ED = EDmax - theta*dlogdet(ADcholC, theta)

    # solve mixed model equations and calculate residuals...
    a <- solveMME(cholC, listC, lWtRinvY, phi, theta)
    r = y - W %*% a

    SS_all <- calcSumSquares(lRinv, lQ, r, a, Nvarcomp)
    Rinv = linearSum(phi,lRinv)
    # here we use Johnson and Thompson 1995:
    yPy = quadForm(y,Rinv,r)
    # yPy2 = sum(theta*SS_all) (should give same results?)

    logL = REMLlogL(ADcholRinv, ADcholGinv, ADcholC,phi,psi,theta, yPy)

    if (monitor) { cat(sprintf("%4d %8.4f\n", it, logL)) }
    if (abs(logLprev-logL) < eps) { break }

    theta = ED / (SS_all + 1.0e-20)
    logLprev = logL

  }
  names(phi) = names(lRinv)
  names(psi) = names(lGinv)
  EDnames = c(names(lRinv),names(lGinv))
  L = list(logL=logL, sigma2e=1.0/phi, tau2e = 1.0/psi, ED = ED, EDmax = EDmax,
           EDnames = EDnames, a=a,
           yhat = W%*%a)
  L
}


# generate list of Ginv matrices, based on dimensions.
generateGinv <- function(dim,namesVarComp)
{
  N <- length(dim)
  M <- diag(N)
  K <- list()
  for (i in 1:N)
  {
    L <- mapply(diag.spam,M[i,],dim)
    K[[i]]<-do.call("bdiag.spam",L)
  }
  names(K) <- namesVarComp
  K
}

LMMsolve <- function(fixed, random = NULL, randomMatrices = NULL, lGinverse=NULL, data, residualterm=NULL, eps=1.0e-6,monitor=FALSE,display=FALSE,
                     maxiter=250) {
  ## make random part:
  if (!is.null(random)) {
    mf <- model.frame(random, data, drop.unused.levels = TRUE, na.action = NULL)
    mt <- terms(mf)
    f.terms <- all.vars(mt)[attr(mt,"dataClasses") == "factor"]
    Z1 <- model.matrix(mt, data = mf,
                    contrasts.arg = lapply(mf[,f.terms, drop = FALSE], contrasts, contrasts = FALSE))
    dim1.r <- table(attr(Z1,"assign"))[-1]
    term1.labels.r <- attr(mt,"term.labels")
    Z1 <- Z1[,-1]
  } else {
    dim1.r = NULL
    Z1 = NULL
    term1.labels.r = NULL
  }

  if (!is.null(randomMatrices))
  {
    ndx <- unlist(randomMatrices)
    Z2 = data[, ndx]
    dim2.r <- sapply(randomMatrices, length)
    term2.labels.r <- names(randomMatrices)
  } else {
    dim2.r = NULL
    Z2 = NULL
    term2.labels.r = NULL
  }
  if (!(is.null(random) & is.null(randomMatrices)))
  {
    if (is.null(random)) {
      Z = Z2
    } else if (is.null(randomMatrices)) {
      Z = Z1
    } else {
      Z = cbind(Z1,Z2)
    }
    Z = as.matrix(Z)
    dim.r = c(dim1.r,dim2.r)
    term.labels.r = c(term1.labels.r,term2.labels.r)

    e <- cumsum(dim.r)
    s <- e - dim.r + 1

    lGinv <- list()
    for(i in 1:length(dim.r)) {
      if (term.labels.r[i] %in% names(lGinverse))
      {
        lGinv[[i]] <- lGinverse[[term.labels.r]]
      } else {
        tmp <- rep(0, sum(dim.r))
        tmp[s[i]:e[i]] <- 1
        lGinv[[i]] <- diag.spam(tmp)
      }
    }
    names(lGinv) = term.labels.r
  } else {
    Z = NULL
    lGinv=NULL
    dim.r=NULL
    term.labels.r=NULL
  }

  ## make fixed part:
  mf <- model.frame(fixed, data)
  mt <- terms(mf)
  f.terms <- all.vars(mt)[attr(mt,"dataClasses") == "factor"]
  X = model.matrix(mt, data=mf, contrasts.arg = lapply(mf[,f.terms, drop = FALSE],
                                                       contrasts, contrasts = TRUE))
  dim.f <- table(attr(X,"assign"))
  term.labels.f <- attr(mt,"term.labels")

  # add intercept....
  if (attr(mt,"intercept")==1)
    term.labels.f <- c("(Intercept)", term.labels.f)


  if (!is.null(residualterm))
  {
    lRinv <- makeRlist(df=data, column=residualterm)
  } else {
    lRinv <- list()
    n = nrow(data)
    lRinv[[1]] <- diag.spam(1,n)
    names(lRinv) = "residual"
  }
  y <- mf[,1]
  obj <- sparseMixedModels(y, X, Z, lGinv, lRinv, eps, monitor, display=display,maxiter=maxiter)

  dim <- as.numeric(c(dim.f,dim.r))
  dim
  term.labels <- c(term.labels.f,term.labels.r)
  term.labels
  obj$dim <- dim
  obj$term.labels <- term.labels
  class(obj) = "LMMsolve"
  obj
}

coef.LMMsolve <- function(obj)
{
  result <- list()
  dim <- obj$dim
  e <- cumsum(dim)
  s <- e - dim + 1

  for(i in 1:length(dim)) {
    result[[i]] <- obj$a[s[i]:e[i]]
  }
  names(result) <- obj$term.labels
  result
}



