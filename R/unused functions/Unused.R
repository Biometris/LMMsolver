# calculate log determinant of matrix M,
logDet <- function(cholM) { as.double(determinant(cholM)$modulus) }

# generate list of Ginv matrices, based on dimensions.
generateGinv <- function(dim,namesVarComp)
{
  N <- length(dim)
  M <- diag(N)
  K <- list()
  for (i in 1:N)
  {
    L <- mapply(spam::diag.spam,M[i,],dim)
    K[[i]]<-do.call("spam::bdiag.spam",L)
  }
  names(K) <- namesVarComp
  K
}

#' @importFrom stats model.frame terms model.matrix
makeZmatrix <- function(ran,dat)
{
  # random part of model, see implementation in SpATS....
  mf <- model.frame(ran, dat, drop.unused.levels = TRUE, na.action = NULL)
  mt <- terms(mf)
  f.terms <- all.vars(mt)[attr(mt,"dataClasses") == "factor"]
  Z <- model.matrix(mt, data = mf, contrasts.arg = lapply(mf[,f.terms, drop = FALSE], contrasts, contrasts = FALSE))
  Z = Z[,-1, drop = FALSE]
  Z
}

#' @importFrom stats model.frame terms
makeGlist <- function(ran,dat)
{
  mf <- model.frame(ran, dat, drop.unused.levels = TRUE, na.action = NULL)
  mt <- terms(mf)
  f.terms <- all.vars(mt)[attr(mt,"dataClasses") == "factor"]
  nlevelsRandom = sapply(mf,nlevels)

  # components of Ginv
  e <- cumsum(nlevelsRandom)
  s <- e - nlevelsRandom + 1
  lGinv = list()
  for (i in 1:length(s))
  {
    tmp <- rep(0,ncol(Z))
    tmp[(s[i]:e[i])] = 1.0
    lGinv[[i]] = spam::cleanup(spam::diag.spam(tmp))
  }
  names(lGinv) = names(mf)
  lGinv
}

# Automatic Differentiation of the Cholesky Algorithm
partial.deriv.logdet = function(G,dG, bandwidth,border)
{
  # add require spam?
  n = dim(G)[1]
  dlogdet_border_banded(spam::triplet(G[n:1,n:1],tri=TRUE),triplet(dG[n:1,n:1],tri=TRUE),bandwidth,border,n)
}

# example with two parameters
gradient.deriv.logdet = function(G,Px,Py)
{
  # add require spam?
  n = dim(G)[1]
  b = spam::bandwidth(G)[1]
  #dlogC2(triplet(G,tri=TRUE),triplet(Px,tri=TRUE),triplet(Py,tri=TRUE),b,n)
}