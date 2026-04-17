#' Construct a ginverse Object from Precision Matrices
#'
#' Creates a \code{ginverse} object from a named list of precision (inverse covariance)
#' matrices. These matrices are typically used to specify the inverse of covariance
#' structures for random effects in \code{LMMsolve}.
#'
#' Each matrix must have identical row and column names corresponding to the levels
#' of the associated random effect. Alignment with the data is checked internally
#' within \code{LMMsolve}.
#'
#' @param precisionMatrices A named list of square matrices (base \code{matrix} or
#'   objects inheriting from \code{Matrix}). Each element represents a precision
#'   matrix corresponding to a random effect. The names of the list must match the
#'   variable names used in the \code{random} argument of \code{LMMsolve}.
#' @param tol A numeric tolerance used for numerical stability (e.g. during inversion
#'   or eigenvalue truncation). Stored as an attribute of the resulting object.
#'
#' @details
#' The function performs basic validation:
#' \itemize{
#'   \item \code{precisionMatrices} must be a named list.
#'   \item Each matrix must be square with identical row and column names.
#'   \item Row and column names are used later to align matrices with factor levels
#'         in the data.
#' }
#'
#' No reordering or alignment with the data is performed at this stage. This is
#' handled internally by \code{LMMsolve}.
#'
#' @return
#' An object of class \code{"ginverse"} (a named list) containing the supplied
#' precision matrices, with attribute \code{"tol"}.
#'
#' @seealso \code{\link{LMMsolve}}
#'
#' @examples
#' library(Matrix)
#'
#' # Create a simple precision matrix
#' K <- Diagonal(5)
#' rownames(K) <- colnames(K) <- as.character(1:5)
#'
#' # Construct ginverse object
#' g <- as.ginverse(list(id = K))
#'
#' g
#'
#' @export
as.ginverse <- function(precisionMatrices, tol = 1e-10) {

  if (!is.list(precisionMatrices)) {
    stop("precisionMatrices must be a list")
  }

  if (is.null(names(precisionMatrices))) {
    stop("precisionMatrices must be a named list")
  }

  for (nm in names(precisionMatrices)) {
    K <- precisionMatrices[[nm]]

    if (!inherits(K, c("matrix", "Matrix"))) {
      stop(sprintf("'%s' must be a matrix or Matrix", nm))
    }

    if (is.null(rownames(K)) || is.null(colnames(K))) {
      stop(sprintf("'%s' must have row and column names", nm))
    }

    if (!identical(rownames(K), colnames(K))) {
      stop(sprintf("'%s' row and column names must be identical", nm))
    }
  }

  structure(
    precisionMatrices,
    class = c("ginverse", "list"),
    tol = tol
  )
}

checkGinverseAgainstData <- function(ginverse, data, random_terms) {

  if (!inherits(ginverse, "ginverse")) {
    stop("ginverse must be created with as.ginverse()")
  }

  out <- list()

  for (nm in names(ginverse)) {

    if (!nm %in% random_terms) {
      stop(sprintf("ginverse '%s' not present in random effects", nm))
    }

    if (!nm %in% names(data)) {
      stop(sprintf("Column '%s' not found in data", nm))
    }

    # Ensure factor
    data[[nm]] <- as.factor(data[[nm]])
    data[[nm]] <- droplevels(data[[nm]])

    ids <- levels(data[[nm]])
    K <- ginverse[[nm]]

    # Check coverage
    if (!all(ids %in% rownames(K))) {
      missing <- setdiff(ids, rownames(K))
      stop(sprintf(
        "ginverse '%s' missing levels: %s",
        nm, paste(missing, collapse = ", ")
      ))
    }

    # Reorder
    K <- K[ids, ids, drop = FALSE]

    # Final safety check
    if (!identical(rownames(K), ids)) {
      stop(sprintf("Alignment failed for '%s'", nm))
    }

    # This is not an efficient conversion
    K <- as.matrix(K)
    K <- spam::as.spam(K)
    #K <- spam::as.spam.dgCMatrix(K)

    out[[nm]] <- K
  }

  return(out)
}



