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

#
# This function saves result of partial derivatives of Cholesky to a
# a spam matrix. It is assumed that cholC and objAD have same values, i.e.:
# - cholC has been updated with:
#     C <- LMMsolver:::linearSum(theta, listP)
#     cholC <- update(cholC, C)
# - objAD has been calculated with:
#     objAD <- LMMsolver:::ADchol(listP)
#     x <- LMMsolver:::dlogdet(objAD, theta)
#
# this function can be useful for example
# for calculating standard errors.
#' @keywords internal
DerivCholesky <- function(cholC, objAD)
{
  # save the results of AD in the numeric vector of spam object:
  cholC@entries <- objAD@ADentries
  A <- as.spam(cholC)
  # put matrix A in the original order
  A <- A[cholC@invpivot, cholC@invpivot]
  A
}


