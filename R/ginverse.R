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
#' @param levels Named list giving the levels corresponding to the
#'   rows and columns of each precision matrix.
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
#' K <- diag(1, 5)
#'
#' # Construct ginverse object
#' g <- as.ginverse(list(id = K),
#'                  levels = list(id = as.character(1:5)))
#' g
#'
#' @export

as.ginverse <- function(precisionMatrices,
                        levels,
                        tol = 1e-10) {

  if (!is.list(precisionMatrices)) {
    stop("precisionMatrices must be a list")
  }

  if (is.null(names(precisionMatrices))) {
    stop("precisionMatrices must be a named list")
  }

  if (!is.list(levels) || is.null(names(levels))) {
    stop("'levels' must be a named list")
  }

  if (!identical(names(precisionMatrices), names(levels))) {
    stop("'levels' must have the same names as 'precisionMatrices'")
  }

  for (nm in names(precisionMatrices)) {

    K <- precisionMatrices[[nm]]
    lev <- levels[[nm]]

    is_spam <- inherits(K, "spam")
    is_matrix <- inherits(K, c("matrix", "Matrix"))

    if (!is_spam && !is_matrix) {
      stop(sprintf(
        "'%s' must be a matrix, Matrix, or spam object",
        nm
      ))
    }

    ## Check dimensions
    if (nrow(K) != ncol(K)) {
      stop(sprintf(
        "'%s' must be square",
        nm
      ))
    }

    ## Check levels
    if (is.null(lev)) {
      stop(sprintf(
        "'levels[[%s]]' must not be NULL",
        nm
      ))
    }

    if (length(lev) != nrow(K)) {
      stop(sprintf(
        "'levels[[%s]]' must have length %d",
        nm, nrow(K)
      ))
    }

    if (anyDuplicated(lev)) {
      stop(sprintf(
        "'levels[[%s]]' contains duplicated levels",
        nm
      ))
    }

    ## Levels should be character, as they will be matched
    ## against factor levels in the data.
    if (!is.character(lev)) {
      stop(sprintf(
        "'levels[[%s]]' must be a character vector",
        nm
      ))
    }
  }

  structure(
    precisionMatrices,
    class = c("ginverse", "list"),
    levels = levels,
    tol = tol
  )
}

checkGinverseAgainstData <- function(ginverse, data, random_terms) {

  if (!inherits(ginverse, "ginverse")) {
    stop("ginverse must be created with as.ginverse()")
  }

  ginverse_levels <- attr(ginverse, "levels")

  if (is.null(ginverse_levels)) {
    stop("ginverse does not contain level information")
  }

  out <- list()

  for (nm in names(ginverse)) {

    if (!nm %in% random_terms) {
      stop(sprintf(
        "ginverse '%s' not present in random effects",
        nm
      ))
    }

    if (!nm %in% names(data)) {
      stop(sprintf(
        "Column '%s' not found in data",
        nm
      ))
    }

    ids <- levels(droplevels(as.factor(data[[nm]])))

    K <- ginverse[[nm]]
    K_levels <- ginverse_levels[[nm]]

    if (is.null(K_levels)) {
      stop(sprintf(
        "No level information available for ginverse '%s'",
        nm
      ))
    }

    if (length(K_levels) != nrow(K)) {
      stop(sprintf(
        "Number of levels for ginverse '%s' does not match matrix dimensions",
        nm
      ))
    }

    if (anyDuplicated(K_levels)) {
      stop(sprintf(
        "ginverse '%s' contains duplicated levels",
        nm
      ))
    }

    missing <- setdiff(ids, K_levels)

    if (length(missing) > 0) {
      stop(sprintf(
        "ginverse '%s' missing levels: %s",
        nm,
        paste(missing, collapse = ", ")
      ))
    }

    ind <- match(ids, K_levels)

    out[[nm]] <- K[ind, ind, drop = FALSE]
  }

  out
}



