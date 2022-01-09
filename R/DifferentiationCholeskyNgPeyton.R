#' @keywords internal
setClass("ADcholnew",
         slots = c(supernodes = "numeric",
                   rowpointers = "numeric",
                   colpointers = "numeric",
                   rowindices = "numeric",
                   P = "matrix"))


#' Automatic differentiation Cholesky, ZtZ and P spam matrices,
#' with C = lambda[1]*P1 + lambda[2]*P2 + .....
#'
#' @importFrom methods new
#' @keywords internal
ADcholnew <- function(P_list) {
  C <- Reduce(`+`, P_list)
  opt <- summary(C)
  cholC <- chol(C, memory = list(nnzR = 8 * opt$nnz,
                                 nnzcolindices = 4 * opt$nnz))
  L <- construct_ADchol_Rcpp_new(cholC, P_list)
  new("ADcholnew",
      supernodes = L$supernodes,
      rowpointers = L$rowpointers,
      colpointers = L$colpointers,
      rowindices = L$rowindices,
      P = L$P)
}


