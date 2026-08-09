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

setMethod("update", "LMMsolver.chol",
          function(object, C, ...) {
            object@entries <- vec(object, C)
            L <- update_Rcpp_fun(object)
            object@entries <- L$entries
            object@ADentries <- L$ADentries
            object
          })

setMethod("solve", "LMMsolver.chol",
          function(a, b, ...) {
            solve_Rcpp_fun(a, b)
          })

updateLinear <- function(object, V, theta) {
   object@entries <- as.vector(V %*% theta)
   L <- update_Rcpp_fun(object)
   object@entries <- L$entries
   object@ADentries <- L$ADentries
   object
}

logdet <- function(object) {
  logdet_Rcpp_fun(object)
}

dlogdet <- function(obj, dC) {
  dF <- obj@ADentries
  v <- vec(obj, dC)
  sum(dF * v)
}

dlogdetLinear <- function(obj, V, theta) {
  g <- as.vector(crossprod(obj@ADentries, V))

  n <- length(obj@pivot)

  # Correct numerical error using the homogeneity identity
  # sum(theta * g) = n.
  n * g / sum(theta * g)
}

vecList <- function(obj, x) {
  do.call(cbind, lapply(x, function(dC) vec(obj, dC)))
}


