#' @keywords internal
setClass("ADchol",
         slots = c(supernodes = "numeric",
                   rowpointers = "numeric",
                   colpointers = "numeric",
                   rowindices = "numeric",
                   entries = "numeric",
                   ADentries  = "numeric",
                   P = "matrix"))


#' Automatic differentiation Cholesky, ZtZ and P spam matrices,
#' with C = lambda[1]*P1 + lambda[2]*P2 + .....
#'
#' @importFrom methods new
#' @keywords internal
ADchol <- function(P_list) {
  C <- Reduce(`+`, P_list)
  opt <- summary(C)
  cholC <- chol(C, memory = list(nnzR = 8 * opt$nnz,
                                 nnzcolindices = 4 * opt$nnz))
  L <- construct_ADchol_Rcpp(cholC, P_list)
  new("ADchol",
      supernodes = L$supernodes,
      rowpointers = L$rowpointers,
      colpointers = L$colpointers,
      rowindices = L$rowindices,
      entries = L$entries,
      ADentries = L$ADentries,
      P = L$P)
}

#' This function saves result of partial derivatives of Cholesky to a
#' a spam matrix, and is used to calculate standard errors and for predictions.
#' @keywords internal
DerivCholesky <- function(cholC) {
  cholC@entries <- partialDerivCholesky(cholC)
  A <- spam::as.spam(cholC)
  ## reordering, can this be done in more efficient way?
  A <- A[cholC@invpivot, cholC@invpivot]
  return(A)
}

