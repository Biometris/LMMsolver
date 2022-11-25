#
# comparison of original code with LMMsolver using sparse GLAM.
# Martin Boer, Biometris, WUR, Wageningen.

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
Z0 = Z

# Add noise
set.seed(2019)
Z = Z + matrix(rnorm(nx * ny), nx, ny)
z <- as.vector(Z)

nseg <- c(100, 100)

dat <- data.frame(x=rep(x,times=ny),y=rep(y,each=nx),z=z)

s <- proc.time()[3]
obj1 <- LMMsolve(fixed=z~1,spline=~spl2D(x1=x,x2=y,nseg=nseg),data=dat,trace=TRUE)
e <- proc.time()[3]

cat("time: ", e-s, " seconds\n")

obj1$logL

# make predictions
pred <- obtainSmoothTrend(obj1, grid=c(nx, ny), includeIntercept = TRUE)

h1 <- (max(dat$x) - min(dat$x)) / nseg[1]
h2 <- (max(dat$y) - min(dat$y)) / nseg[2]
lambda1 <- as.numeric(obj1$sigma2e / obj1$tau2e[1]/ h1 ^ 3)
lambda2 <- as.numeric(obj1$sigma2e / obj1$tau2e[2]/ h2 ^ 3)
lambda1
lambda2
lambda2/lambda1
alpha <- lambda1/(lambda1+lambda2)
alpha

cols = gray(seq(0, 1, by = 0.01))
par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))

image(x, y, Z, col = cols, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
title("Data", cex.main = 1.5)

Zhat <- matrix(data=pred$ypred, nrow=nx,ncol=ny, byrow=TRUE)
image(x, y, Zhat, col = cols, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
title("Smoothed", cex.main = 1.5)

#Zhat <- matrix(data=pred$ypred, nrow=nx,ncol=ny, byrow=TRUE)
#image(x, y, Zhat-Z0, col = cols, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
#title("Difference", cex.main = 1.5)

summary(obj1, which='variances')


