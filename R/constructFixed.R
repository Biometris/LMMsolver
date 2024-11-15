constructFixed <- function(fixed, data) {

  ## Make fixed part.
  mf <- model.frame(fixed, data, drop.unused.levels = TRUE)
  mt <- terms(mf)
  f.terms <- all.vars(mt)[attr(mt, "dataClasses") == "factor"]
  X <- Matrix::sparse.model.matrix(mt, data = mf,
                                   contrasts.arg = lapply(X = mf[, f.terms, drop = FALSE],
                                                          FUN = contrasts, contrasts = TRUE))
  term.labels.f <- attr(mt, "term.labels")

  q <- qr(as.matrix(X))
  if (q$rank != ncol(X)) {
    remCols <- q$pivot[-seq(q$rank)]
    ## Compare terms before and after removing extra columns.
    ## If a complete term is removed, it also has to be removed from the labels.
    f.terms.orig <- as.numeric(names(table(attr(X, "assign"))))

    dim.f.tab <- table(attr(X, "assign")[-remCols])
    dim.f <- as.numeric(dim.f.tab)
    X <- X[ , -remCols, drop = FALSE]
    f.terms.new <- as.numeric(names(dim.f.tab))
    if (!setequal(f.terms.orig, f.terms.new)) {
      term.labels.f <- term.labels.f[-setdiff(f.terms.orig, f.terms.new)]
    }
  } else {
    dim.f <- as.numeric(table(attr(X, "assign")))
  }
  ## Add intercept.
  if (attr(mt, "intercept") == 1) {
    term.labels.f <- c("(Intercept)", term.labels.f)
  }

  attr(X, "dim.f") <- dim.f
  attr(X, "term.labels.f") <- term.labels.f
  #attr(X, "mt") <- mt
  attr(X, "mf") <- mf
  return(X)
}
