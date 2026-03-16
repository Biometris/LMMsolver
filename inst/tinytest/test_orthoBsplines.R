
##
## Bernstein basis tests
##

expect_equal(bernstein(2, 0.5, 0), 0.25)
expect_equal(bernstein(2, 0.5, 1), 0.5)
expect_equal(bernstein(2, 0.5, 2), 0.25)

# Partition of unity
p <- 4
x <- runif(10)
for (xi in x) {
  s <- sum(sapply(0:p, function(k) bernstein(p, xi, k)))
  expect_equal(s, 1, tolerance = 1e-12)
}


##
## Bernstein matrix consistency
##

p <- 3
B <- bernstein_matrix(p, seq(0,1,length=p+1))

# rows sum to 1 (partition of unity)
expect_equal(rowSums(B), rep(1, p+1), tolerance = 1e-12)


##
## Bernstein Gram matrix symmetry
##

p <- 4
G <- bernstein_gram_coef(p)

expect_true(isSymmetric(G))

# positive definiteness
e <- eigen(G, symmetric = TRUE)$values
expect_true(all(e > 0))


##
## GramMatrix consistency with analytic cases
##

for (p in 0:3) {
  q <- p + 1
  G_num <- GramMatrix(p, q)
  G_ana <- Gram_analytic(p)

  expect_equal(G_num, G_ana, tolerance = 1e-10)
}


##
## B-spline integral of constant function
##

p <- 3
q <- 6

a <- rep(1, q)

int_val <- integrate_bspline(a, p)

expect_equal(int_val, 1, tolerance = 1e-10)


##
## Linear function integral
##

p <- 3
q <- 6

# approximate x using spline basis
a <- runif(q)

# integrate twice (should scale linearly)
I1 <- integrate_bspline(a, p)
I2 <- integrate_bspline(2*a, p)

expect_equal(I2, 2*I1, tolerance = 1e-10)


##
## First moment integral
##

p <- 3
q <- 6
a <- rep(1, q)

moment <- integrate_x_bspline(a, p)

# integral of x over [0,1]
expect_equal(moment, 0.5, tolerance = 1e-10)


##
## Orthogonality weights
##

p <- 3
q <- 8

w <- ortho_int_condition(p, q)

expect_equal(length(w), q)
expect_true(all(w > 0))


##
## Orthogonal difference matrix
##

p <- 3
q <- 8

M <- ortho_diff_matrix(p, q)

expect_equal(dim(M), c(q, q-1))


##
## Orthogonal x-difference matrix
##

M2 <- ortho_x_diff_matrix(p, q)

expect_equal(dim(M2), c(q, q-2))


##
## Cross product helper
##

u <- c(1,2,3)
v <- c(4,5,6)

w <- cross3(u,v)

expect_equal(w, c(-3,6,-3))


##
## Orthogonality check for constructed matrices
##

p <- 3
q <- 8

w0 <- ortho_int_condition(p,q)
w1 <- ortho_x_int_condition(p,q)

M <- ortho_x_diff_matrix(p,q)

# convert spam -> dense
M <- as.matrix(M)

expect_equal(drop(w0 %*% M), rep(0, q-2), tolerance = 1e-10)
expect_equal(drop(w1 %*% M), rep(0, q-2), tolerance = 1e-10)
