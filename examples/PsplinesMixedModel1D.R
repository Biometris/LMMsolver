#' ---
#' title: P-splines in one dimension using mixed models
#' author: Martin Boer, Biometris, WUR
#' ---

rm(list=ls())

library(splines)
library(LMMsolver)

set.seed(1234)

# residual error
sigma2 <- 1

Nobs <- 100
xmin = 2
xmax = 3

x0 <- seq(xmin,xmax,by=0.001)
fsim <- function(x) {5*sin(2*pi*x)}

x = sort(runif(Nobs, min=2, max=3))
y = fsim(x) +  rnorm(Nobs, sd=sqrt(sigma2) )

sim.df <- data.frame(x=x, y=y)

degree = 3
k = degree+1
nseg = 10

dx = (xmax - xmin) / nseg

x = sim.df$x

knots = seq(xmin - degree * dx, xmax + degree * dx, by = dx)
B = splineDesign(knots, x, derivs=rep(0,length(x)), ord = degree+1)
q = ncol(B)
q

D <- diff(diag(q), diff=2)
DtD <- crossprod(D)
U <- eigen(DtD)$vectors[,-c(q-1,q)]
d <- eigen(DtD)$values[-c(q-1,q)]
Usc <- U %*% diag(1/sqrt(d))

lZ <- list()
lZ[[1]] = B %*% Usc
Z <- do.call("cbind",lZ)

df_ext = cbind(sim.df,Z)
df_ext$InterceptMB <- 1.0

dim <- sapply(lZ,ncol)
e <- cumsum(dim) + ncol(sim.df)
s <- e - dim + 1

# model without interactions:
lM <- list(B=c(s[1]:e[1]))
obj0 = LMMsolve(fixed=y~-1+InterceptMB+x, random=NULL,randomMatrices=lM, data=df_ext)
obj0$logL
obj0$ED

xgrid = seq(xmin, xmax, by=0.001)
Bx = splineDesign(knots, xgrid, derivs=rep(0,length(xgrid)), ord = degree+1,outer.ok=TRUE)

mu <- coef(obj0)$InterceptMB
beta  <- coef(obj0)$x
theta <- Usc %*% coef(obj0)$B
yfit <- mu + beta * xgrid + Bx %*% theta
plot(x=x, y=y)
lines(x=xgrid, y=fsim(xgrid), col='blue',lwd=2.5)
lines(x=xgrid, y=yfit, col='red' ,lwd=2.5)

#
#
# some extra code, to show that the non-linear part is orthogonal to the linear space:
#
# we can also show that the integrals over x are equal to zero,
# if we use an extended B-spline basis:
# int(f(x) g(x), x = -infty...infty) = 0
# int(f(x) h(x), x = -infty...infty) = 0
# where f(x) = sum(B_i(x), i=1...q), g(x) = 1, h(x) = x.
# So <f,g> = <f,h> = 0, which means thtat f is orthogonal to g and h.
#

# see sparse_mixed_model_splines.pdf:
# we can construct a sequence t_mn such that B(x) t_mn = x for [xmin,xmax]
# and therefore we can write X = B %*% A, with A orthogonal to Usc:
t_mn <- seq(xmin-(k/2)*dx+dx,by=dx,length=q)
A <- cbind(1, t_mn)
X <- cbind(1, x)
all.equal(X, B%*%A, check.attributes=FALSE)
range(t(A) %*% Usc)

# use extended B-spline basis:
xmin_ext <- xmin - degree*dx
xmax_ext <- xmax + degree*dx
knots_ext <- seq(xmin_ext - degree * dx, xmax_ext + degree * dx, by = dx)
xgrid_ext = seq(xmin_ext, xmax_ext, length = 5000)
Bx_ext = splineDesign(knots_ext, xgrid_ext, derivs=rep(0,length(xgrid_ext)), ord = degree+1)
q_ext = ncol(Bx_ext)

# we use extended B-spline basis, but extra B-splines equal to zero:
theta_ext <- c(rep(0,degree), theta, rep(0, degree))
fx = Bx_ext %*% theta_ext

# function g(x) = 1 and h(x) = x derived from extended B-spline basis:
gx <- Bx_ext %*% rep(1,q_ext)
t_mn <- seq(xmin_ext-(k/2)*dx+dx,by=dx,length=q_ext)
hx <- Bx_ext %*% t_mn

# function f(x)*g(x)
plot(x=xgrid_ext,y=fx*gx,type='l', main='f(x)*g(x)')
abline(h=0,col='red')
abline(v=xmin,lt=2)
abline(v=xmax,lt=2)
abline(v=xmin_ext,lt=2)
abline(v=xmax_ext,lt=2)
# numerical integration gives 0:
sum(fx*gx)
# integration per B-spline:
int_Bx <- rep(0, q_ext)
for (i in 1:q_ext) {
  int_Bx[i] = sum(Bx_ext[,i] * gx)
}
sum(int_Bx * theta_ext)

# function f(x)*h(x)
plot(x=xgrid_ext,y=fx*hx,type='l', main='f(x)*h(x)')
abline(h=0,col='red')
abline(v=xmin,lt=2)
abline(v=xmax,lt=2)
abline(v=xmin_ext,lt=2)
abline(v=xmax_ext,lt=2)
# numerical integration gives 0:
sum(fx*hx)

# integration per B-spline:
int_B2x <- rep(0, q_ext)
for (i in 1:q_ext) {
  int_B2x[i] = sum(Bx_ext[,i]^2 * gx * hx)
}
sum(int_B2x * theta_ext)

