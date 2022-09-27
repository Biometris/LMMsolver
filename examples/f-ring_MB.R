# Smoothing ring image with array regression (Simulated data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(JOPS)
library(fields)
library(LMMsolver)
library(spam)

# Simulate the rings
nx = 1000
ny = 1000
x = seq(-1, 1, length = nx)
y = seq(-1, 1, length = ny)
ex = rep(1, nx)
ey = rep(1, ny)
X = outer(x, ey)
Y = outer(ex, y)
R1 = sqrt((X - 0.3)^2 + (Y - 0.3)^2)
R2 = sqrt((X + 0.2)^2 + (Y + 0.2)^2)
R3 = sqrt((X - 0.7)^2 + (Y + 0.7)^2)
R4 = sqrt((X + 0.7)^2 + (Y - 0.7)^2)
Z1 = exp(-50 * (R1 - 0.4)^2)
Z2 = exp(-50 * (R2 - 0.6)^2)
Z3 = exp(-50 * (R3 - 0.2)^2)
Z4 = exp(-50 * (R4 - 0.2)^2)
Z = pmax(pmax(pmax(Z1, Z2), Z3), Z4)

# Add noise
set.seed(2019)
Z = Z + matrix(rnorm(nx * ny), nx, ny)

nseg <- c(80, 70)

s1 <- proc.time()[3]

# Prepare bases
Bx = bbase(x, nseg = nseg[1])
By = bbase(y, nseg = nseg[2])
nbx = ncol(Bx)
nby = ncol(By)

# Prpare the penalty matrices
Dx = diff(diag(nbx), diff = 2)
Dy = diff(diag(nby), diff = 2)
lambdax = lambday = 1
Px = lambdax * t(Dx) %*% Dx
Py = lambday * t(Dy) %*% Dy
P = kronecker(Py, diag(nbx)) + kronecker(diag(nby), Px)

# Do the smoothing, using the array algorithm
W = matrix(runif(nx*ny), nx, ny)
#W = 0 * Z + 1
Tx = rowtens(Bx)
Ty = rowtens(By)
Q = t(Tx) %*% W %*% Ty
dim(Q) = c(nbx, nbx, nby, nby)
Q = aperm(Q, c(1, 3, 2, 4))
dim(Q) = c(nbx * nby, nbx * nby)
r = t(Bx) %*% (Z * W) %*% By
dim(r) = c(nbx * nby, 1)
A = solve(Q + P, r)
dim(A) = c(nbx, nby)
Zhat = Bx %*% A %*% t(By)

e1 <- proc.time()[3]

s2 <- proc.time()[3]

z <- as.vector(Z)
w <- as.vector(W)
knots1 <- LMMsolver:::PsplinesKnots(min(x),max(x),degree = 3, nseg = nseg[1])
Bx1 <- LMMsolver:::Bsplines(knots1, x)
knots2 <- LMMsolver:::PsplinesKnots(min(y),max(y),degree = 3, nseg = nseg[2])
Bx2 <- LMMsolver:::Bsplines(knots2, y)

sparseGLAM <- LMMsolver:::SparseGLAM(Bx2, Bx1)
BtWB <- LMMsolver:::calcBtWB(sparseGLAM, w)

lambdax <- lambday <- 1
Dx <- diff(diag.spam(nbx), diff = 2)
Dy <- diff(diag.spam(nby), diff = 2)

Px <- lambdax * t(Dx) %*% Dx
Py <- lambday * t(Dy) %*% Dy
P <- Py %x% diag.spam(nbx) + diag.spam(nby) %x% Px
#BtB <- crossprod(Bx2) %x% crossprod(Bx1)
C <- BtWB + P
BtZ <- LMMsolver:::KronProd2(t(Bx2), t(Bx1), w*z)
a <- solve(C, BtZ)
zhat <- LMMsolver:::KronProd2(Bx2, Bx1, a)

Zhat2 <- matrix(data=zhat, nx, ny)

range(Zhat-Zhat2)

e2 <- proc.time()[3]

cat("book:       ", e1-s1, " seconds\n")
cat("LMMsolver:  ", e2-s2, " seconds\n")



