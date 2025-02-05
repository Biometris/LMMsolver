#' ---
#' title: Smoothing on torus, simulated data
#' author: Martin Boer, Biometris, WUR, Wageningen
#' ---

library(JOPS)
library(circular)
library(ggplot2)
library(LMMsolver)
suppressMessages(library(spam))

set.seed(1234)

cnt.data <- TRUE
if (cnt.data) {
  fam <- poisson()
} else {
  fam <- gaussian()
}

sim_fun <- function(x) {
  y <- exp(sin(2*pi*x[1])) + exp(cos(2*pi*x[2]))
  y
}

# simulate data on a grid:
n1 <- 100
n2 <- 100
n <- n1*n2
dat <- expand.grid(x1=seq(0,1,length=n1),x2=seq(0,1,length=n2))
dat$ytrue <- apply(dat, MARGIN=1, FUN = sim_fun)

# take subset for training
Ntraining <- 2500
k <- sample(x=c(1:n), size = Ntraining)
dat_train <- dat[k, ]
if (cnt.data) {
  dat_train$y <- sapply(dat_train$ytrue,
                        function(lambda) {rpois(n=1,lambda)})
} else {
  dat_train$y <- dat_train$ytrue + rnorm(n=Ntraining)
}

ggplot() +
  geom_tile(data = dat, aes(x=x1,y=x2,fill=ytrue)) +
  scale_fill_gradientn(colors = topo.colors(100)) +
  coord_fixed() + ggtitle("simulated true") + JOPS_theme()

ggplot() +
  geom_tile(data = dat_train, aes(x=x1,y=x2,fill=y)) +
  scale_fill_gradientn(colors = topo.colors(100)) +
  coord_fixed() + ggtitle("Training: simulated plus noise") +
  JOPS_theme()

nseg <- c(20,25)

obj1 <- LMMsolve(fixed = y~1,
                 spline = ~torus(x1, x2, nseg),
                 family = fam,
                 data = dat_train)
summary(obj1)

# prediction by foot, values should be in [0, 1]
grid <- expand.grid(x1 = seq(0,1,length=100),
                    x2 = seq(0,1,length=100))
x1 <- grid$x1
x2 <- grid$x2
knots <- obj1$splRes[[1]]$knots
B1 <- LMMsolver:::cBsplines(knots[[1]], x1)
B2 <- LMMsolver:::cBsplines(knots[[2]], x2)
B12 <- LMMsolver:::RowKronecker(B1, B2)

# define fixed part of model X:
z1 <- c(0:(nseg[1]-1))/nseg[1]
z2 <- c(0:(nseg[2]-1))/nseg[2]
G1 <- cbind(sin(2*pi*z1), cos(2*pi*z1))
G2 <- cbind(sin(2*pi*z2), cos(2*pi*z2))
X1 <- B1 %*% G1
X2 <- B2 %*% G2
G <- G1 %x% G2
X12 <- B12 %*% G
X <- cbind(1, X1, X2, X12)
U <- cbind(X, B12)

# make predictions:
a <- obj1$coefMME
yhat <- U %*% a
pred <- dat
pred$ypred <- fam$linkinv(yhat)

ggplot() +
  geom_tile(data = pred, aes(x=x1,y=x2, fill = ypred)) +
  scale_fill_gradientn(colors = topo.colors(100)) +
  coord_fixed() + ggtitle("fitted data LMMsolver torus") +
  JOPS_theme()
