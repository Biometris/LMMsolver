#' @keywords internal
setClass("ADchol",
         slots = c(supernodes = "numeric",
                   rowpointers = "numeric",
                   colpointers = "numeric",
                   rowindices = "numeric",
                   pivot = "numeric",
                   invpivot = "numeric",
                   entries = "numeric",
                   ADentries  = "numeric",
                   P = "matrix"))


#' construct object for Automated Differentiation Cholesky decomposition
#'
#' Construct object for reverse Automated Differentiation of Cholesky decomposition,
#' with as input a list of semi-positive symmetric sparse matrices \eqn{P_i}, each of
#' dimension \eqn{q \times q}. The function \code{ADchol} calculates the matrix \eqn{C}, the sum
#' the precision matrices \eqn{P_i}: \eqn{C = \sum_{i}  P_i}. Next, it calculates the Cholesky
#' Decomposition using the multiple minimum degree (MMD) algorithm
#' of the \code{spam} package.
#'
#' @param lP a list of symmetric matrices of class spam, each of dimension \eqn{q \times q},
#' and with sum of the matrices assumed to be positive definite.
#
#' @returns An object of class \code{ADchol}. This object is used to calculate the partial
#' partial derivatives of \eqn{log|C|} in an efficient way.
#'
#' @references
#' Furrer, R., & Sain, S. R. (2010). spam: A sparse matrix R package with emphasis
#' on MCMC methods for Gaussian Markov random fields.
#' Journal of Statistical Software, 36, 1-25.
#'
#' @importFrom methods new
#' @keywords internal
ADchol <- function(lP) {
  C <- Reduce(`+`, lP)
  opt <- summary(C)
  cholC <- chol(C, memory = list(nnzR = 8 * opt$nnz,
                                 nnzcolindices = 4 * opt$nnz))
  # reorder the matrices in list lP by double transpose, row-permutations are much faster
  # than column permutations (see help permutation() function in spam library)
  lQ <- lapply(lP, function(x) {
    z <- x[cholC@pivot,]
    tz <- spam::t(z)
    tz <- tz[cholC@pivot,]
    return(spam::t(tz)) })
  L <- construct_ADchol_Rcpp(cholC, lQ)
  new("ADchol",
      supernodes = L$supernodes,
      rowpointers = L$rowpointers,
      colpointers = L$colpointers,
      rowindices = L$rowindices,
      pivot = L$pivot,
      invpivot = L$invpivot,
      entries = L$entries,
      ADentries = L$ADentries,
      P = L$P)
}

