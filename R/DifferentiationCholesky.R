# Automatic Differentiation of the Cholesky Algorithm
partial.deriv.logdet = function(G,dG, bandwidth,border)
{
  # add require spam?
  n = dim(G)[1]
  dlogdet_border_banded(triplet(G[n:1,n:1],tri=TRUE),triplet(dG[n:1,n:1],tri=TRUE),bandwidth,border,n)
}

# example with two parameters
gradient.deriv.logdet = function(G,Px,Py)
{
  # add require spam?
  n = dim(G)[1]
  b = bandwidth(G)[1]
  #dlogC2(triplet(G,tri=TRUE),triplet(Px,tri=TRUE),triplet(Py,tri=TRUE),b,n)
}

#
#' @export
setClass("ADchol", slots=c(colpointers = 'numeric', 
                           rowindices  = 'numeric',
                           P           = 'matrix'))

#
# Automatic differentiation Cholesky, ZtZ and P spam matrices,
# with C = lambda[1]*P1 + lambda[2]*P2 + .....
ADchol = function(P_list)
{
  nelem = length(P_list)
  nCol = ncol(P_list[[1]])
  
  lambda = rep(1.0,nelem)

  C = spam(0,nrow=nCol, ncol=nCol)
  for (i in 1:nelem) {
    C = C + lambda[i]*P_list[[i]]
  }
  cholC = chol(C)
  L = construct_ADchol_Rcpp(cholC,P_list)
  new("ADchol", colpointers = L$colpointers,
                rowindices  = L$rowindices,
                P           = L$P)
}

