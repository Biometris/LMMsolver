#' Calculate standard errors for predictions
#'
#'
#' @keywords internal
calcStandardErrors <- function(C, U)
{
  ## !!! NOT CHANGE THE LINE OF CODE BELOW !!!
  ## It adds extra zeros ("fill-ins") to matrix C, needed
  ## to calculate the Partial Derivatives of Cholesky, not equal to zero.
  C = C + 0 * spam::crossprod.spam(U)
  cholC <- chol(C)

  ## calculate the partial derivatives of Cholesky
  cholC@entries <- partialDerivCholesky(cholC)

  ## convert to spam object and put in original order
  A <- spam::as.spam(cholC)
  A <- A[cholC@invpivot, cholC@invpivot]

  ## Equivalent to v <- diag(U %*% A %*% t(U))
  x <- spam::rowSums.spam((U %*% A) * U)
  se <- sqrt(x)
  return(se)
}

