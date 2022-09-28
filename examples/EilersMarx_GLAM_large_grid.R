#
# comparison of original code with LMMsolver using sparse GLAM.
# Martin Boer, Biometris, WUR, Wageningen.

# Smoothing ring image with array regression (Simulated data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(JOPS)
library(fields)

# Simulate the rings
nx = 500
ny = 500
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

nseg <- c(20, 20)

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

W = matrix(runif(nx*ny), nx, ny)
#W = 0 * Z + 1

s1 <- proc.time()[3]

# Do the smoothing, using the array algorithm
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

#
# part 2, using LMMsolver functions
#

library(LMMsolver)
library(spam)

# Make the B-splines and penalty matrix sparse, Z and W as vectors
Bx <- as.spam(Bx)
By <- as.spam(By)
P <- as.spam(P)
z <- as.vector(Z)
w <- as.vector(W)

s2 <- proc.time()[3]

sGLAMobj <- LMMsolver:::sparseGLAM(By, Bx)
BtWB <- LMMsolver:::calcBtWB(sGLAMobj, w)
BtZ <- LMMsolver:::calcBtY(sGLAMobj, w*z)
a <- solve(BtWB + P, BtZ)
zhat <- LMMsolver:::calcBa(sGLAMobj, a)

e2 <- proc.time()[3]

# compare original code with LMMsolver.
Zhat2 <- matrix(data=zhat, nx, ny)
all.equal(Zhat, Zhat2)

t1 <- e1-s1
t2 <- e2-s2

cat("book:       ", t1, " seconds\n")
cat("LMMsolver:  ", t2, " seconds\n")
cat("Factor:     ", round(t1/t2, 2), "times faster\n")



