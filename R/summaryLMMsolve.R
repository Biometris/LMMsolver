#' summary function
#'
#' @param object an object of class LMMsolve
#' @param which default dimensions (only option so far).
#' @param \dots some methods for this generic require additional arguments.
#' None are used in this method.
#'
#' @description Prints the effective dimensions of the model.
#'
#' @export
summary.LMMsolve <- function(object, which = "dimensions", ...) {

  # start and end of each variance component....
  Nres <- object$Nres
  e <- cumsum(object$varPar)
  s <- e - object$varPar + 1

  e <- e + Nres
  s <- s + Nres

  nVarComp <- length(object$term.labels.r)

  namesED <- names(object$ED)
  for (i in 1:nVarComp)
  {
    ndx <- c(s[i]:e[i])
    cat(object$term.labels.r[i], "with total effective dimension ",
        round(sum(object$ED[ndx]),2), "\n")
    if (length(ndx) > 1) {
      for (k in ndx) {
        cat("   ", namesED[k], "\t", round(object$ED[k],2 ), "\n")
      }
    }
  }
}
