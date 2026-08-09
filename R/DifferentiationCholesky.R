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
                   P = "ANY",
                   user_def = "ANY",
                   mode = "character"))


# Reorder a spam matrix according to the Cholesky pivot.
#
# The permutation is applied to both rows and columns. The implementation
# uses two row permutations (via transpose) because these are much faster
# than a direct column permutation for spam matrices.
#
# @keywords internal
reorderSpam <- function(x, permutation)
{
  z <- x[permutation, ]
  z <- spam::t(z)
  z <- z[permutation, ]
  spam::t(z)
}

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
  lQ <- lapply(lP, reorderSpam, permutation = cholC@pivot)
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
      P = L$P,
      mode = "linear",
      user_def = NULL)
}


convertSparseMatrices <- function(lX, obj)
{
  do.call(cbind,
          lapply(lX,
                 convertSparseMatrix_Rcpp,
                 ADobj = obj))
}

dlogdet <- function(obj, theta, b = NULL)
{
  if (obj@mode == "linear") {
    return(dlogdet_cpp_linear(obj, theta, b))
  }

  if (obj@mode == "nonlinear") {

    model_eval <- obj@user_def(theta)

    # use R-indexed pivot:
    pivot_R <- obj@pivot + 1
    ## reorder C and derivatives
    C <- reorderSpam(model_eval$C, pivot_R)
    dC <- lapply(model_eval$dC,
                 reorderSpam,
                 permutation = pivot_R)

    entries <- convertSparseMatrix_Rcpp(C, obj)
    derivatives <- convertSparseMatrices(dC, obj)

    #stop("Need entries(C)")

    return(dlogdet_cpp_general(obj,
                               entries,
                               derivatives,
                               b))
  }

  stop("Unknown ADchol mode.")
}

#
# new LMMsolver.chol class
#

#' @keywords internal
setClass("LMMsolver.chol",
         slots = c(supernodes = "numeric",
                   rowpointers = "numeric",
                   colpointers = "numeric",
                   rowindices = "numeric",
                   pivot = "numeric",
                   invpivot = "numeric",
                   entries = "numeric",
                   ADentries  = "numeric"))

SparseCholesky <- function(C) {
  opt <- summary(C)
  cholC <- chol(C, memory = list(nnzR = 8 * opt$nnz,
                                 nnzcolindices = 4 * opt$nnz))
  N_entries <- length(cholC@entries)

  # Exchange row and columns compared to spam object, as in Ng and Peyton 1993
  # LMMsolver.chol uses C-index (0) instead of R-index (1)
  obj <- new("LMMsolver.chol",
            supernodes = cholC@supernodes - 1,
            colpointers = cholC@rowpointers - 1,
            rowpointers = cholC@colpointers - 1,
            rowindices = cholC@colindices - 1,
            pivot = cholC@pivot - 1,
            invpivot = cholC@invpivot - 1,
            entries = rep(0, N_entries),
            ADentries = rep(0, N_entries))

  obj@entries <- vec(obj, C)
  L <- update_Rcpp_fun(obj)
  obj@entries <- L$entries
  obj@ADentries <- L$ADentries
  return(obj)
}

#' @export
setMethod("update", "LMMsolver.chol",
          function(object, C, ...) {
            object@entries <- vec(object, C)
            L <- update_Rcpp_fun(object)
            object@entries <- L$entries
            object@ADentries <- L$ADentries
            object
          })

updateLinear <- function(object, V, theta) {
   object@entries <- as.vector(V %*% theta)
   L <- update_Rcpp_fun(object)
   object@entries <- L$entries
   object@ADentries <- L$ADentries
   object
}

logdet <- function(object) { 2*sum(object@entries)}

dlogdetLinear <- function(obj, V, theta) {
  g <- as.vector(crossprod(obj@ADentries, V))

  n <- length(obj@pivot)

  # Correct numerical error using the homogeneity identity
  # sum(theta * g) = n.
  n * g / sum(theta * g)
}

dlogdet <- function(obj, dC) {
  dF <- obj@ADentries
  v <- vec(obj, dC)
  sum(dF * v)
}

vecList <- function(obj, x) {
  do.call(cbind, lapply(x, function(dC) vec(obj, dC)))
}

#' @export
setMethod("solve", "LMMsolver.chol",
          function(a, b, ...) {
            solve_Rcpp_fun(a, b)
          })

logdetGradient <- function(obj, C) {
  # use R-indexed pivot:
  pivot_R <- obj@pivot + 1
  ## reorder C and derivatives
  C <- reorderSpam(C, pivot_R)
  entries <- convertSparseMatrix_Rcpp(C, obj)
  return(dlogdetVector_Rcpp(obj, entries))
}


