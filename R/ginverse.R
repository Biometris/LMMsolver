#' Construct a ginverse Object from Precision Matrices
#'
#' Creates a \code{ginverse} object from a named list of precision
#' (inverse covariance) matrices. These matrices are typically used to
#' specify the inverse of covariance structures for random effects in
#' \code{LMMsolve}.
#'
#' Each precision matrix must be square. For \code{matrix} and \code{Matrix}
#' objects, the row and column names define the corresponding levels. For
#' \code{spam} objects, which do not use row and column names, the
#' corresponding levels must be supplied through \code{levels}.
#'
#' @param precisionMatrices A named list of square matrices. Each element
#'   must be a base \code{matrix}, an object inheriting from \code{Matrix},
#'   or a \code{spam} object. Each element represents a precision matrix
#'   corresponding to a random effect.
#' @param levels An optional named list giving the levels corresponding to
#'   the rows and columns of the precision matrices. This is required for
#'   \code{spam} objects, which do not have row and column names. For
#'   \code{matrix} and \code{Matrix} objects, levels are obtained from the
#'   row names; if supplied, they are checked for consistency with the row
#'   and column names.
#' @param tol A numeric tolerance used for numerical stability (e.g. during
#'   inversion or eigenvalue truncation). Stored as an attribute of the
#'   resulting object.
#'
#' @details
#' The function performs basic validation:
#' \itemize{
#'   \item \code{precisionMatrices} must be a named list.
#'   \item If supplied, \code{levels} must be a named list with matching
#'         names.
#'   \item Each precision matrix must be square.
#'   \item For \code{matrix} and \code{Matrix} objects, row and column names
#'         must be present and identical.
#'   \item For \code{spam} objects, \code{levels} must be supplied.
#'   \item Levels must be a character vector with length equal to the
#'         corresponding matrix dimension and contain no duplicates.
#'   \item If \code{levels} is supplied for a \code{matrix} or \code{Matrix}
#'         object, it must agree with its row and column names.
#' }
#' No reordering or alignment with the data is performed at this stage. This
#' is handled internally by \code{LMMsolve}.
#'
#' @return
#' An object of class \code{"ginverse"} (a named list) containing the supplied
#' precision matrices, with attributes \code{"levels"} and \code{"tol"}.
#'
#' @seealso \code{\link{LMMsolve}}
#'
#' @examples
#' K <- diag(1, 5)
#' dimnames(K) <- list(as.character(1:5), as.character(1:5))
#'
#' # Construct ginverse object from a matrix with names
#' g <- as.ginverse(list(id = K))
#' g
#'
#' # A spam matrix requires levels to be supplied
#' # Kspam <- spam::as.spam(K)
#' # g <- as.ginverse(list(id = Kspam),
#' #                  levels = list(id = as.character(1:5)))
#'
#' @export

as.ginverse <- function(precisionMatrices,
                        levels = NULL,
                        tol = 1e-10) {

  if (!is.list(precisionMatrices)) {
    stop("precisionMatrices must be a list")
  }

  if (is.null(names(precisionMatrices))) {
    stop("precisionMatrices must be a named list")
  }

  ## Check levels argument if supplied
  if (!is.null(levels)) {

    if (!is.list(levels) || is.null(names(levels))) {
      stop("'levels' must be a named list")
    }

    if (!identical(names(precisionMatrices), names(levels))) {
      stop("'levels' must have the same names as 'precisionMatrices'")
    }
  }

  ## Effective levels for each precision matrix
  ginverse_levels <- vector("list", length(precisionMatrices))
  names(ginverse_levels) <- names(precisionMatrices)

  for (nm in names(precisionMatrices)) {

    K <- precisionMatrices[[nm]]

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

    ## --------------------------------------------------------------
    ## Obtain and check levels
    ## --------------------------------------------------------------

    if (is_spam) {

      ## spam objects do not have row/column names
      if (is.null(levels) || is.null(levels[[nm]])) {
        stop(sprintf(
          "'levels[[%s]]' must be supplied for spam objects",
          nm
        ))
      }

      lev <- levels[[nm]]

    } else {

      rn <- rownames(K)
      cn <- colnames(K)

      if (is.null(rn) || is.null(cn)) {
        stop(sprintf(
          "'%s' must have row and column names",
          nm
        ))
      }

      if (!identical(rn, cn)) {
        stop(sprintf(
          "Row and column names of '%s' must be identical",
          nm
        ))
      }

      lev <- rn

      ## If levels were supplied, check consistency
      if (!is.null(levels)) {

        supplied <- levels[[nm]]

        if (!identical(supplied, lev)) {
          stop(sprintf(
            "'levels[[%s]]' does not match the row and column names of '%s'",
            nm, nm
          ))
        }
      }
    }

    ## Levels should be character, as they will be matched
    ## against factor levels in the data.
    if (!is.character(lev)) {
      stop(sprintf(
        "Levels for '%s' must be a character vector",
        nm
      ))
    }

    if (length(lev) != nrow(K)) {
      stop(sprintf(
        "Number of levels for '%s' must be %d",
        nm, nrow(K)
      ))
    }

    if (anyDuplicated(lev)) {
      stop(sprintf(
        "Levels for '%s' contain duplicated values",
        nm
      ))
    }

    ginverse_levels[[nm]] <- lev
  }

  structure(
    precisionMatrices,
    class = c("ginverse", "list"),
    levels = ginverse_levels,
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
