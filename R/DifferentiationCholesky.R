#' @keywords internal
setClass("ADchol",
         slots = c(supernodes = "numeric",
                   rowpointers = "numeric",
                   colpointers = "numeric",
                   rowindices = "numeric",
                   entries = "numeric",
                   ADentries  = "numeric",
                   P = "matrix"))


#' construct object for Automated Differentiation Cholesky decomposition
#'
#' Construct object for reverse Automated Differentiation of Cholesky decomposition,
#' with as input a list of semi-positive symmetric sparse matrices \eqn{P_i}, each of
#' dimension \eqn{q \times q}. The function `ADchol` calculates the matrix \eqn{C}, the sum
#' the precision matrices \eqn{P_i}: \eqn{C = \sum_{i}  P_i}. Next, it calculates the Cholesky
#' Decomposition using the multiple minimum degree (MMD) algorithm
#' of the `spam` package.
#'
#' @param P_list a list of symmetric matrices of class spam, each of dimension \eqn{q \times q},
#' and with sum of the matrices assumed to be positive definite.
#
#' @return An object of class `ADchol`.
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


