chkInputLMMsolve <- function(fixed, random, data,
                             ginverse, residual, tolerance, maxit,
                             grpTheta, offset, family) {
  splTrms <- NULL
  if (!inherits(data, "data.frame")) {
    stop("data should be a data.frame.\n")
  }
  if (!inherits(fixed, "formula") || length(terms(fixed)) != 3) {
    stop("fixed should be a formula of the form \"resp ~ pred\".\n")
  }
  if (!is.null(random) &&
      (!inherits(random, "formula") || length(terms(random)) != 2)) {
    stop("random should be a formula of the form \" ~ pred\".\n")
  }
  if (!is.null(ginverse) &&
      (!is.list(ginverse) ||
       length(names(ginverse)) == 0 ||
       (!all(sapply(X = ginverse, FUN = function(x) {
         (is.matrix(x) || spam::is.spam(x)) && isSymmetric(x)}))))) {
    stop("ginverse should be a named list of symmetric matrices.\n")
  }
  if (!is.null(residual) &&
      (!inherits(residual, "formula") || length(terms(residual)) != 2)) {
    stop("residual should be a formula of the form \" ~ pred\".\n")
  }
  if (!is.numeric(tolerance) || length(tolerance) > 1 || tolerance < 0) {
    stop("tolerance should be a positive numerical value.")
  }
  if (!is.numeric(maxit) || length(maxit) > 1 || maxit < 0) {
    stop("maxit should be a positive numerical value.")
  }
  if (!is.null(grpTheta) &&
      (!is.numeric(grpTheta) || !isTRUE(all.equal(round(grpTheta),grpTheta)) ||
       max(grpTheta) != length(unique(grpTheta)))) {
    stop("grpTheta should be integers 1,2,...nGrp")
  }
  if (is.character(offset))
  {
    if (!hasName(data, offset)) {
      stop("offset ", offset, " not defined in the data.\n")
    }
    offset <- data[[offset]]
  }
  if (!is.numeric(offset)) {
    stop("offset should be numeric")
  }
  if (!(inherits(family, "family") || inherits(family, "familyLMMsolver"))) {
    stop("argument family not correct.\n")
  }
  return(offset)
}
