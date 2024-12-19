#' Display the sparseness of the mixed model coefficient matrix
#'
#' @param object an object of class LMMsolve.
#' @param cholesky Should the cholesky decomposition of the coefficient matrix
#' be plotted?
#'
#' @returns A plot of the sparseness of the mixed model coefficient matrix.
#'
#' @examples
#' ## Fit model on john.alpha data from agridat package.
#' data(john.alpha, package = "agridat")
#'
#' ## Fit simple model with only fixed effects.
#' LMM1 <- LMMsolve(fixed = yield ~ rep + gen,
#'                 data = john.alpha)
#'
#' ## Obtain deviance.
#' displayMME(LMM1)
#'
#' @export
displayMME <- function(object,
                       cholesky = FALSE) {
  if (!inherits(object, "LMMsolve")) {
    stop("object should be an object of class LMMsolve.\n")
  }
  if (!cholesky) {
    spam::display(object$C)
  } else {
    cholC <- chol(object$C)
    L <- t(spam::as.spam(cholC))
    spam::display(L)
  }
}
