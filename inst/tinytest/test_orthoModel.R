
##
## Synthetic data
##
set.seed(1)

n <- 50

data <- data.frame(
  x1 = runif(n),
  x2 = runif(n)
)

data$y <- sin(2*pi*data$x1) + rnorm(n, sd = 0.1)


##
## Basis specification
##
b1 <- basis(x1, xmin = 0, xmax = 1, nseg = 6, deg = 3, pord = 2)
b2 <- basis(x2, xmin = 0, xmax = 1, nseg = 5, deg = 3, pord = 2)

bases <- list(
  x1 = b1,
  x2 = b2
)


##
## Test basis object
##
expect_true(inherits(b1, "basis_spec"))

expect_equal(b1$var, "x1")
expect_equal(b1$xmin, 0)
expect_equal(b1$xmax, 1)
expect_equal(b1$nseg, 6)


##
## Test eval_basis structure
##
f <- eval_basis(b1, data)

expect_true(is.list(f))

expect_true(all(c("B","M","P","C","knots") %in% names(f)))

q <- ncol(f$B)

expect_equal(nrow(f$B), n)
expect_equal(ncol(f$M), q-1)
expect_equal(nrow(f$M), q)

expect_equal(nrow(f$P), q-1)
expect_equal(ncol(f$P), q-1)

expect_equal(dim(f$C), c(q-1,1))


##
## Test penalty symmetry
##
P_dense <- as.matrix(f$P)

expect_true(isSymmetric(P_dense))


##
## Test orthogonality property of M
##
w <- ortho_int_condition(p = b1$deg, q = q)

expect_equal(drop(w %*% f$M), rep(0, q-1), tolerance = 1e-10)


##
## Test interaction construction
##
f1 <- eval_basis(b1, data)
f2 <- eval_basis(b2, data)

B12 <- RowKronecker(f1$B, f2$B)
M12 <- f1$M %x% f2$M

expect_equal(
  ncol(B12 %*% M12),
  ncol(f1$M) * ncol(f2$M)
)


##
## Test orthoModel (main effects only)
##
model <- y ~ x1 + x2

fit <- orthoModel(
  model = model,
  bases = bases,
  data  = data
)

expect_true(is.list(fit))

expect_true("dim" %in% names(fit))

expect_true("response_name" %in% names(fit))


##
## Dimension consistency
##
dimtab <- fit$dim

expect_true(all(c("term","dim","s","e") %in% names(dimtab)))

expect_equal(dimtab$term[1], "Intercept")


##
## Test model with interaction
##
model2 <- y ~ x1 + x2 + x1:x2

#fit2 <- orthoModel(
#  model = model2,
#  bases = bases,
#  data  = data
#)

#expect_true(is.list(fit2))
#expect_true(nrow(fit2$dim) >= 3)


##
## Check that constraints matrix matches dimension
##
TotDim <- max(fit$dim$e)

expect_true(TotDim > 0)


##
## Test prediction matrix dimensions
##
f <- eval_basis(b1, data)

Z <- f$B %*% f$M

expect_equal(nrow(Z), n)


##
## Penalty scaling check
##
dx <- attr(f$knots, "dx")

expect_true(dx > 0)
