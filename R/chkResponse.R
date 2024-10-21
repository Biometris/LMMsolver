chkResponse <- function(y, residual, data)
if (is.null(residual)) {
  if (var(y) < .Machine$double.eps / 2) {
    stop("Variance response variable zero or almost zero.\n")
  }
} else {
  resVar <- all.vars(residual)
  varGrp <- tapply(X = y, INDEX = data[[resVar]], FUN = var)
  ndxVar <- which(varGrp < .Machine$double.eps / 2)
  if (length(ndxVar) > 0) {
    levels_f <- levels(data[[resVar]])
    levelsNoVar <- paste(levels_f[ndxVar], collapse = ", ")
    stop("Variance response variable zero or almost zero for levels:\n",
         levelsNoVar, "\n")
  }
}
