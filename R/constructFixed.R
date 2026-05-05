
constructFixed_train <- function(fix, data) {

  ## Remove response
  fix[[2]] <- NULL

  mf <- model.frame(fix, data, drop.unused.levels = TRUE)

  ## Name check
  names_mf <- names(mf)
  IsValidName <- names_mf == make.names(names_mf)
  if (!all(IsValidName)) {
    stop("Syntactically invalid name(s): ",
         paste(names_mf[which(!IsValidName)], collapse=", "),
         call. = FALSE)
  }

  mt <- terms(mf)

  f.terms <- all.vars(mt)[attr(mt, "dataClasses") == "factor"]

  contrasts.arg <- lapply(
    mf[, f.terms, drop = FALSE],
    contrasts,
    contrasts = TRUE
  )

  X <- Matrix::sparse.model.matrix(
    mt,
    data = mf,
    contrasts.arg = contrasts.arg
  )

  term.labels.f <- attr(mt, "term.labels")

  ## ---- rank check ----
  q <- qr(as.matrix(X))

  if (q$rank != ncol(X)) {

    remCols <- q$pivot[-seq(q$rank)]

    f.terms.orig <- as.numeric(names(table(attr(X, "assign"))))

    dim.f.tab <- table(attr(X, "assign")[-remCols])
    dim.f <- as.numeric(dim.f.tab)

    X <- X[, -remCols, drop = FALSE]

    f.terms.new <- as.numeric(names(dim.f.tab))

    if (!setequal(f.terms.orig, f.terms.new)) {
      term.labels.f <- term.labels.f[-setdiff(f.terms.orig, f.terms.new)]
    }

  } else {
    dim.f <- as.numeric(table(attr(X, "assign")))
  }

  ## Add intercept
  if (attr(mt, "intercept") == 1) {
    term.labels.f <- c("(Intercept)", term.labels.f)
  }

  ## ---- BUILD SPEC ----
  spec <- list(
    fix            = fix,
    terms          = mt,
    xlevels        = .getXlevels(mt, mf),
    contrasts      = contrasts.arg,
    colnames       = colnames(X),

    assign         = attr(X, "assign"),
    dim.f          = dim.f,
    term.labels.f  = term.labels.f,
    intercept      = attr(mt, "intercept")
  )

  return(list(X = X, spec = spec))
}


constructFixed_pred <- function(spec, data) {
  mt <- spec$terms

  check_required_vars(mt, data)

  check_new_levels(data, spec$xlevels)

  mf <- model.frame(
    mt,
    data,
    xlev = spec$xlevels,
    drop.unused.levels = FALSE
  )

  f.terms <- all.vars(mt)[attr(mt, "dataClasses") == "factor"]

  contrasts.arg <- spec$contrasts

  X <- Matrix::sparse.model.matrix(
    mt,
    data = mf,
    contrasts.arg = contrasts.arg
  )

  # enforce identical column structure
  X <- X[, spec$colnames, drop = FALSE]

  return(X)
}







