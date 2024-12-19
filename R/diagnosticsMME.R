#' Give diagnostics for mixed model coefficient matrix C and the cholesky
#' decomposition
#'
#' @param object an object of class LMMsolve.
#'
#' @returns A summary of the mixed model coefficient matrix and its choleski
#' decomposition.
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
#' diagnosticsMME(LMM1)
#'
#' @export
diagnosticsMME <- function(object) {
  if (!inherits(object, "LMMsolve")) {
    stop("object should be an object of class LMMsolve.\n")
  }
  cat("Summary of matrix C \n")
  print(spam::summary.spam(object$C))
  cat("\n Summary of cholesky decomposition of C \n")
  print(spam::summary.spam(chol(object$C)))
}
