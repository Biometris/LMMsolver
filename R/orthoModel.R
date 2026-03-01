
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
  
  knots <- LMMsolver:::PsplinesKnots(
    xmin   = b$xmin,
    xmax   = b$xmax,
    degree = b$deg,
    nseg   = b$nseg
  )
  
  B <- LMMsolver:::Bsplines(knots, x)
  q <- ncol(B)
  
  M <- ortho_diff_matrix(p = b$deg, q = q)
  
  D   <- diff(diag(q), diff = b$pord)
  DtD <- crossprod(D)
  P   <- t(M) %*% DtD %*% M
  
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

orthoModel <- function(model, bases, data) {
  
  ## ---- 0. Extract response ----
  response_name <- deparse(model[[2]])
  y <- data[[response_name]]
  n <- length(y)  
  
  ## ---- 1. Parse model ----
  tt <- terms(model)
  term_labels <- attr(tt, "term.labels")
  term_labels <- setdiff(term_labels, "1")  
  
  if (!all(term_labels %in% names(bases))) {
    stop("All model terms must be provided in 'bases'")
  }
  
  ## ---- 2. Evaluate bases and build Z ----
  Z_list <- list()
  P_list <- list()
  C_list <- list()
  
  for (nm in term_labels) {
    f_eval <- eval_basis(bases[[nm]], data)
    
    BM <- f_eval$B %*% f_eval$M
    Z_list[[nm]] <- BM
    P_list[[nm]] <- f_eval$P
    C_list[[nm]] <- f_eval$C
  }
  
  Z <- do.call(spam::cbind.spam, unname(Z_list))
  
  ## ---- 3. Fixed effects (intercept only) ----
  X <- spam::spam(1, nrow = n, ncol = 1)  
  
  ## ---- 4. Block-diagonal penalty list ----
  lP <- vector("list", length(P_list))
  for (i in seq_along(P_list)) {
    blocks <- vector("list", length(P_list))
    for (j in seq_along(P_list)) {
      blocks[[j]] <- if (i == j) P_list[[j]] else 0 * P_list[[j]]
    }
    lP[[i]] <- do.call(spam::bdiag.spam, blocks)
  }
  
  ## ---- 5. Constraint matrix (full coefficient space) ----
  q <- vapply(Z_list, ncol, 0L)
  C_restrict <- matrix(0, nrow = sum(q), ncol = length(C_list))
  
  r0 <- 0
  for (j in seq_along(C_list)) {
    idx <- r0 + seq_len(q[j])
    C_restrict[idx, j] <- as.vector(as.matrix(C_list[[j]]))
    r0 <- r0 + q[j]
  }
  C_restrict <- spam::as.spam(C_restrict)
  
  ## ---- 6. Residual precision ----
  lRinv <- list(spam::diag.spam(1, n))
  attr(lRinv, "cnt") <- n
  
  ## ---- 7. Fit ----
  fit <- LMMsolver:::sparseMixedModels(
    y          = y,
    X          = X,
    Z          = Z,
    lGinv      = lP,
    lRinv      = lRinv,
    C_restrict = C_restrict
  )
  
  fit
}

