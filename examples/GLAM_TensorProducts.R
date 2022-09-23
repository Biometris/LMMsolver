library(LMMsolver)
suppressPackageStartupMessages(library(spam))
set.seed(1234)

x1min <- 0
x1max <- 7
x2min <- 1
x2max <- 6
x3min <- 7
x3max <- 12
knots1 <- LMMsolver:::PsplinesKnots(x1min, x1max, degree=3, nseg=2)
knots2 <- LMMsolver:::PsplinesKnots(x2min, x2max, degree=3, nseg=4)
knots3 <- LMMsolver:::PsplinesKnots(x3min, x3max, degree=3, nseg=3)

x1 <- seq(x1min, x1max, by=1)
x2 <- seq(x2min, x2max, by=1)
x3 <- seq(x3min, x3max, by=1)
Bx1 <- LMMsolver:::Bsplines(knots1, x1)
Bx2 <- LMMsolver:::Bsplines(knots2, x2)
Bx3 <- LMMsolver:::Bsplines(knots3, x3)

q1 <- ncol(Bx1)
q2 <- ncol(Bx2)
q3 <- ncol(Bx3)
B <- Bx1 %x% Bx2 %x% Bx3
q <- q1*q2*q3
y <- rnorm(q)

yB1 <- as.vector(B %*% y)
yB2 <- LMMsolver:::prod3(Bx1, Bx2, Bx3, y)
yB3 <- LMMsolver:::prodList(list(Bx1,Bx2,Bx3), y)

all.equal(yB1, yB2)
all.equal(yB1, yB3)


