
basis <- function(x, xmin, xmax, nseg, deg, pord) {

  var <- deparse(substitute(x))

  structure(
    list(
      var   = var,     # character, e.g. "x1"
      xmin  = xmin,
      xmax  = xmax,
      nseg  = nseg,
      deg   = deg,
      pord  = pord
    ),
    class = "basis_spec"
  )
}

eval_basis <- function(b, data) {

  x <- data[[b$var]]

  knots <- PsplinesKnots(
    xmin   = b$xmin,
    xmax   = b$xmax,
    degree = b$deg,
    nseg   = b$nseg
  )

  B <- Bsplines(knots, x)
  q <- ncol(B)

  M <- ortho_diff_matrix(p = b$deg, q = q)

  D   <- diff(diag(q), diff = b$pord)
  DtD <- crossprod(D)
  dx <- attr(knots, which="dx")
  sc <- (1 / dx)^(2 * b$pord - 1)
  P  <- sc * (t(M) %*% DtD %*% M)

  C <- spam::spam(0, nrow = q - 1, ncol = 1)
  C[1, 1] <- C[q - 1, 1] <- 1.0

  list(
    B     = B,
    M     = M,
    P     = P,
    C     = C,
    knots = knots
  )
}

orthoModel <- function(model, bases, data, trace=FALSE) {

  ## ---- 0. Extract response ----
  response_name <- all.vars(model)[attr(terms(model), "response")]

  y <- data[[response_name]]
  n <- length(y)

  ## ---- 1. Parse model ----
  tt <- terms(model)
  term_labels <- attr(tt, "term.labels")
  term_labels <- setdiff(term_labels, "1")

  main_terms <- term_labels[!grepl(":", term_labels)]
  interaction_terms <- term_labels[grepl(":", term_labels)]

  X <- spam::spam(x=1, nrow = n, ncol = 1)
  Z_list <- list()
  P_list <- list()
  C_list <- list()

  df_dim <- NULL
  for (nm in main_terms) {
    f <- eval_basis(bases[[nm]], data)

    Z_list[[nm]] <- f$B %*% f$M
    P_list[[nm]] <- f$P
    C_list[[nm]] <- f$C

    df_dim <- rbind(df_dim, data.frame(term=nm, dim=ncol(f$M)))
  }
  if (length(interaction_terms) >0 ) {
    for (it in interaction_terms) {
      vars <- strsplit(it, ":")[[1]]
      f1 <- eval_basis(bases[[vars[1]]], data)
      f2 <- eval_basis(bases[[vars[2]]], data)

      ## interaction design
      B12  <- RowKronecker(f1$B, f2$B)
      M12  <- f1$M %x% f2$M
      BM12 <- B12 %*% M12
      Z_list[[it]] <- BM12

      ## isotropic penalties
      P12 <- f1$P %x% (t(f2$M) %*% f2$M) + (t(f1$M) %*% f1$M) %x% f2$P
      P_list[[it]] <- P12

      ## pure interaction constraint
      C_list[[it]] <- f1$C %x% f2$C
      df_dim <- rbind(df_dim, data.frame(term=it, dim=ncol(f1$M)*ncol(f2$M)))
    }
  }
  e <- cumsum(df_dim$dim)
  s <- e - df_dim$dim + 1
  df_dim$s <- s
  df_dim$e <- e
  df_dim
  nTerms <- nrow(df_dim)
  TotDim <- sum(df_dim$dim)
  Z <- do.call(spam::cbind.spam, unname(Z_list))


  ## ---- build lP efficiently (NO sparse replacement) ----

  lP <- vector("list", nTerms)

  for (k in seq_len(nTerms)) {

    blocks <- vector("list", nTerms)

    for (j in seq_len(nTerms)) {
      if (j == k) {
        blocks[[j]] <- P_list[[k]]
      } else {
        qj <- df_dim$e[j] - df_dim$s[j] + 1
        blocks[[j]] <- spam::spam(0, qj, qj)
      }
    }

    lP[[k]] <- do.call(spam::bdiag.spam, blocks)
  }

  C_restrict <- spam::spam(x=0, ncol=length(C_list),nrow=TotDim)
  for (k in seq_len(nTerms)) {
    ndx <- c(df_dim$s[k] : df_dim$e[k])
    C_restrict[ndx, k] <- C_list[[k]]
  }
  C_restrict <- spam::cleanup(C_restrict)

  ## ---- 7. Residual precision ----
  lRinv <- list(spam::diag.spam(1, n))
  attr(lRinv, "cnt") <- n

  ## ---- 8. Fit ----
  fit <- sparseMixedModels(
    y          = y,
    X          = X,
    Z          = Z,
    lGinv      = lP,
    lRinv      = lRinv,
    C_restrict = C_restrict,
    trace      = trace
  )
  # add extra info
  df_dim$s <- df_dim$s + 1
  df_dim$e <- df_dim$e + 1
  df_dim <- rbind(data.frame(term = "Intercept", dim = 1, s=1, e=1), df_dim)

  fit$dim <- df_dim
  fit$response_name <- response_name
  fit$bases <- bases
  fit
}

predict_ortho <- function(object, newdat) {

  # helper to extract coefficients
  get_coef <- function(obj, term) {
    k <- match(term, obj$dim$term)
    ndx <- obj$dim$s[k]:obj$dim$e[k]
    obj$a[ndx]
  }

  f1_basis <- object$bases$f1
  f2_basis <- object$bases$f2

  # coefficients
  mu  <- get_coef(object, "Intercept")
  u1  <- get_coef(object, "f1")
  u2  <- get_coef(object, "f2")
  u12 <- get_coef(object, "f1:f2")

  # evaluate bases
  eval_f1 <- eval_basis(f1_basis, newdat)
  eval_f2 <- eval_basis(f2_basis, newdat)

  # marginals
  g1 <- as.vector(eval_f1$B %*% eval_f1$M %*% u1)
  g2 <- as.vector(eval_f2$B %*% eval_f2$M %*% u2)

  # interaction
  B12 <- RowKronecker(eval_f1$B, eval_f2$B)
  M12 <- eval_f1$M %x% eval_f2$M
  g12 <- as.vector((B12 %*% M12) %*% u12)

  # total fit
  fit_total <- as.vector(mu + g1 + g2 + g12)

  # return everything
  out <- cbind(newdat,
               g1 = g1,
               g2 = g2,
               g12 = g12,
               fit = fit_total)

  return(as.data.frame(out))
}

predict_marginal <- function(object, term, newdat) {

  # helper
  get_coef <- function(obj, term) {
    k <- match(term, obj$dim$term)
    ndx <- obj$dim$s[k]:obj$dim$e[k]
    obj$a[ndx]
  }

  # get basis spec
  basis_spec <- object$bases[[term]]

  # evaluate basis
  f_eval <- eval_basis(basis_spec, newdat)

  # coefficients
  u <- get_coef(object, term)

  # marginal effect
  g <- as.vector(f_eval$B %*% f_eval$M %*% u)

  # return
  out <- data.frame(newdat, value = g)
  names(out)[ncol(out)] <- term  # name column as f1 or f2

  return(out)
}

