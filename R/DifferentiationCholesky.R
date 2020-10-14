#' @export
setClass("ADchol",
         slots = c(colpointers = "numeric",
                   rowindices = "numeric",
                   P = "matrix"))


#' Automatic differentiation Cholesky, ZtZ and P spam matrices,
#' with C = lambda[1]*P1 + lambda[2]*P2 + .....
#'
#' @importFrom methods new
ADchol <- function(P_list) {
  nelem <- length(P_list)
  nCol <- ncol(P_list[[1]])

  lambda <- rep(1.0, nelem)

  C <- spam::spam(0, nrow = nCol, ncol = nCol)
  for (i in 1:nelem) {
    C <- C + lambda[i] * P_list[[i]]
  }
  cholC <- chol(C)
  L <- construct_ADchol_Rcpp(cholC, P_list)
  new("ADchol", colpointers = L$colpointers, rowindices = L$rowindices,
      P = L$P)
}

