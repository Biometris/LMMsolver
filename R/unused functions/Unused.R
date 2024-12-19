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

#' Construct index matrix
#'
#' Construct index matrix.
#'
#' @param df A data.frame.
#' @param lZ A list of matrices.
#' @param names A character vector of names.
#'
#' @returns The index matrix.
#'
#' @keywords internal
ndxMatrix <- function(df,
                      lZ,
                      names) {
  n <- length(lZ)
  dim <- sapply(X = lZ, FUN = ncol)
  e <- cumsum(dim) + ncol(df)
  s <- e - dim + 1
  lM <- list()
  for (i in 1:n) {
    lM[[i]] <- c(s[i]:e[i])
  }
  names(lM) <- names
  return(lM)
}

# //' Spectral decomposition of D'D
# //'
# //' Spectral decomposition of D'D
# //'
# //' @param q A numeric value.
# //' @param ord A numeric value.
# //'
# //' @returns A q x (q-ord) matrix, the spectral decomposition D'D.
# //'
# //' @export
# // [[Rcpp::export]]
# arma::mat calcUsc(const double& q,
#                   const double& ord) {
#   arma::mat D = arma::diff(arma::eye(q, q), ord);
#   arma::mat DtD = D.t() * D;
#
#   arma::vec eigVals( D.n_cols );
#   arma::mat eigVecs( size(D) );
#   arma::eig_sym(eigVals, eigVecs, DtD);
#   // fliplr and reverse are needed because eigenVals and eigenVecs are returned
#   // sorted in ascending order.
#   eigVecs = fliplr( eigVecs.tail_cols(q - ord) );
#   arma::mat S = arma::diagmat( 1 / sqrt(reverse(eigVals.tail(q - ord))) );
#   return eigVecs * S;
# }

# //' Row-wise kronecker product
# //'
# //' Row-wise kronecker product
# //'
# //' @param X1 A matrix.
# //' @param X2 A matrix.
# //'
# //' @returns The row-wise kronecker product of X1 and X2.
# //'
# //' @export
# // [[Rcpp::export]]
# arma::mat RowKronecker(const arma::mat& X1,
#                        const arma::mat& X2) {
#   arma::mat one1 = arma::ones(1, X1.n_cols);
#   arma::mat one2 = arma::ones(1, X2.n_cols);
#   arma::mat rowKron = arma::kron(X1, one2) % arma::kron(one1, X2);
#   return rowKron;
# }
