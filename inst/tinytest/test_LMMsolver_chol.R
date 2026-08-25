expect_silent({
  library(spam)
  library(LMMsolver)

  n1 <- 2
  n2 <- 5
  n <- n1 * n2

  D1 <- diff.spam(diag.spam(n1), diff = 1)
  D2 <- diff.spam(diag.spam(n2), diff = 1)

  A <- diag.spam(1, n)
  B1 <- diag.spam(1, n1) + t(D1) %*% D1
  B2 <- diag.spam(1, n2) + t(D2) %*% D2
  B <- kronecker(B1, B2)

  f <- function(theta) theta[1] * A + theta[2] * B

  theta0 <- c(1, 1)
  C0 <- f(theta0)

  obj <- LMMsolver:::SparseCholesky(C0)

  expect_true(is(obj, "LMMsolver.chol"))
  expect_equal(length(obj@entries), length(obj@ADentries))

  ## Linear representation
  V <- LMMsolver:::vecList(obj, list(A, B))

  ## Individual derivatives
  g1 <- LMMsolver:::dlogdet(obj, A)
  g2 <- LMMsolver:::dlogdet(obj, B)
  g <- c(g1, g2)

  ## Homogeneity check:
  ## sum(theta * d log|C|/d theta) = n
  expect_equal(sum(theta0 * g), n, tolerance = 1e-10)

  ## Update to new theta
  theta <- c(1, 2)
  C <- f(theta)

  obj <- LMMsolver:::updateLinear(obj, V, theta)

  ## General update interface

  obj_update <- LMMsolver:::SparseCholesky(C0)
  obj_update <- LMMsolver:::update(obj_update, C)

  expect_equal(
    obj_update@entries,
    obj@entries,
    tolerance = 1e-12
  )

  expect_equal(
    obj_update@ADentries,
    obj@ADentries,
    tolerance = 1e-12
  )

  ## Individual derivatives after update
  g1_new <- LMMsolver:::dlogdet(obj, A)
  g2_new <- LMMsolver:::dlogdet(obj, B)
  g_new <- c(g1_new, g2_new)

  expect_equal(sum(theta * g_new), n, tolerance = 1e-10)

  ## Vectorized linear derivative
  g_linear <- LMMsolver:::dlogdetLinear(obj, V, theta)

  expect_equal(g_linear, g_new, tolerance = 1e-12)

  ## log determinant
  expect_equal(
    LMMsolver:::logdet(obj),
    as.numeric(determinant(C, logarithm = TRUE)$modulus),
    tolerance = 1e-10
  )

  ## Explicit normalization
  g_norm <- n * g_linear / sum(theta * g_linear)

  expect_equal(sum(theta * g_norm), n, tolerance = 1e-12)

  ## Solve
  set.seed(1234)
  b <- rnorm(n)

  x1 <- solve(C, b)
  x2 <- LMMsolver:::solve(obj, b)

  expect_equal(x2, x1, tolerance = 1e-10)

  #
  ## Backward-compatible ADchol interface
  #
  lP <- list(A, B)

  ADobj <- LMMsolver:::ADchol(lP)

  expect_true(inherits(ADobj, "ADchol"))
  expect_true(is(ADobj$chol, "LMMsolver.chol"))

  ## Same gradient as the new interface
  ## (used function name dlogdet before)
  g_ad <- LMMsolver:::dlogdet(ADobj, theta)

  expect_equal(
    g_ad,
    g_linear,
    tolerance = 1e-12
  )
})


