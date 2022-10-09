#' Standard errors for predictions
#'
#' Calculates the standard errors for predictions \eqn{D \hat{u}},
#' see Welham et al. 2004 and Gilmour et al. 2004 for details.
#'
#' The prediction error variance is given by \eqn{D C^{-1} D'},
#' where \eqn{C} is the mixed model coefficient matrix, and \eqn{D} defines
#' linear combinations of fixed and random effects.
#' The standard errors are given by the the square root of
#' the diagonal. To calculate the standard errors in an efficient way we use that
#'
#' \deqn{\frac{\partial log|C + \xi_i d_i d_i'|}{\partial \xi_i} |_{\xi_i=0}
#'  = trace(C^{-1} d_i d_i') =
#' trace(d_i' C^{-1} d_i) = d_i' C^{-1} d_i, }
#' where \eqn{d_i} is row \eqn{i} of matrix \eqn{D}. The values of
#' \eqn{d_i' C^{-1} d_i} can be calculated more efficient, avoiding the
#' calculation of the inverse of \eqn{C}, by using Automated Differentiation
#' of the Choleksy algorithm, see section 2.3 in Smith (1995) for details.
#'
#' @param C a symmetric matrix of class spam
#' @param D a matrix of class spam
#'
#' @return a vector with standard errors for predictions \eqn{D \hat{u}}.
#'
#' @references
#' Welham, S., Cullis, B., Gogel, B., Gilmour, A., & Thompson, R. (2004).
#' Prediction in linear mixed models.
#' Australian & New Zealand Journal of Statistics, 46(3), 325-347.
#'
#' @references
#' Smith, S. P. (1995). Differentiation of the Cholesky algorithm.
#' Journal of Computational and Graphical Statistics, 4(2), 134-147.
#'
#' @references
#' Gilmour, A., Cullis, B., Welham, S., Gogel, B., & Thompson, R. (2004).
#' An efficient computing strategy for prediction in mixed linear models.
#' Computational statistics & data analysis, 44(4), 571-586.
#'
#' @keywords internal
calcStandardErrors <- function(C, D)
{
  ## !!! NOT CHANGE THE LINE OF CODE BELOW !!!
  ## It adds extra zeros ("fill-ins") to matrix C, needed
  ## to calculate the Partial Derivatives of Cholesky, not equal to zero.
  C = C + 0 * spam::crossprod.spam(D)
  cholC <- chol(C)

  p <- cholC@pivot
  tD <- t(D)
  tDp <- tD[p, ]
  x <- diagXCinvXt(cholC, tDp)
  se <- sqrt(x)
  return(se)
}

