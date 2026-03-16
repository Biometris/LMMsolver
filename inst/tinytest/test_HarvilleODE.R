
##
## Mask used in tests
##
Mask <- rbind(
  c(1,1),
  c(2,2),
  c(2,1),
  c(3,3)
)

##
## Test matrix <-> vector conversion
##
C <- matrix(
  c(4,1,0,
    1,5,0,
    0,0,6),
  3,3, byrow = TRUE
)

theta <- LMMsolver:::mat2vec(C, Mask)
C2 <- LMMsolver:::vec2mat(theta, Mask)

expect_equal(C2, C)


##
## Test eta -> xi -> eta roundtrip
##
eta <- c(0.3, -0.4, 0.2, 0.5)

xi <- LMMsolver:::eta_to_xi(eta, Mask)
eta2 <- LMMsolver:::xi_to_eta(xi, Mask)

expect_equal(eta2, eta, tolerance = 1e-10)


##
## Test xi -> theta consistency
##
xi <- matrix(
  c(2,0,0,
    1,3,0,
    0,0,1),
  3,3, byrow = TRUE
)

theta <- LMMsolver:::xi_to_theta(xi, Mask)
thetaMat <- LMMsolver:::vec2mat(theta, Mask)

expect_equal(thetaMat, tcrossprod(xi))


##
## Test theta -> xi inversion
##
xi <- LMMsolver:::eta_to_xi(c(0.2,0.1,0.3,-0.2), Mask)

theta <- LMMsolver:::xi_to_theta(xi, Mask)
thetaMat <- LMMsolver:::vec2mat(theta, Mask)

xi2 <- LMMsolver:::theta_to_xi(thetaMat)

expect_equal(xi2, xi, tolerance = 1e-10)


##
## Test rotation helper
##
xi <- matrix(
  c(2,0,0,
    1,3,0,
    0,0,1),
  3,3, byrow = TRUE
)

M <- diag(c(1,2,3))

rot <- LMMsolver:::rotate_theta(M, xi)

expect_equal(rot, crossprod(xi, M %*% xi))


##
## Full transformation pipeline test
##
eta <- rnorm(4)

xi <- LMMsolver:::eta_to_xi(eta, Mask)
theta <- LMMsolver:::xi_to_theta(xi, Mask)

thetaMat <- LMMsolver:::vec2mat(theta, Mask)

xi2 <- LMMsolver:::theta_to_xi(thetaMat)
eta2 <- LMMsolver:::xi_to_eta(xi2, Mask)

expect_equal(eta2, eta, tolerance = 1e-10)


##
## Additional edge case: zero correlations
##
eta <- c(0, 0, 0, 0)

xi <- LMMsolver:::eta_to_xi(eta, Mask)
theta <- LMMsolver:::xi_to_theta(xi, Mask)

thetaMat <- LMMsolver:::vec2mat(theta, Mask)

expect_true(isSymmetric(thetaMat))
expect_true(all(diag(thetaMat) > 0))

# simulate some data
set.seed(1234)

n <- 100
x <- runif(n)
y <- exp(-2*x)*sin(x) + rnorm(n)

n <- length(y)

# as in Eilers and Marx 1996
xmin <- 0
xmax <- 1
nseg <- 20
deg <- 3

# construction of knots for P-splines amd Phi-splines
knots <- LMMsolver:::PsplinesKnots(xmin = xmin, xmax = xmax, degree=deg, nseg=nseg)
B <- LMMsolver:::Bsplines(knots, x)
G <- LMMsolver:::constructG(knots,scaleX=TRUE,pord=2)
G[, 1] <- 1

# desigm matrices
X <- B %*% G
Z <- B

# Construct Penalty, with built-in boundary constraints:
q <- ncol(B)
D <- spam::diff.spam(spam::diag.spam(q), diff=2)
DtD <- crossprod(D)
B_ref <- LMMsolver:::Bsplines(knots, x=c(xmin,xmax))
P <- DtD + crossprod(B_ref)
lGinv <- list(P)

lRinv <- list(spam::diag.spam(n))

# define the mask
Mask <- matrix(data=c(1,1, 2,2), nrow=2, ncol=2, byrow=TRUE)

obj_ODE <- LMMsolver:::HarvilleODE(y, X, Z, lGinv, lRinv, Mask, alpha=1.0, maxiter=100, thr=1.e-6)

dat <- data.frame(x=x,y=y)
obj <- LMMsolve(y~1, spline = ~spl1D(x, nseg=nseg, xlim=c(xmin,xmax), scale=TRUE), data=dat)
obj$logL
obj_ODE$logL
expect_equal(obj_ODE$logL, obj$logL)




