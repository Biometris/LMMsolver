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

obj <- LMMsolve(fixed = y~1,
                 spline = ~spl2D(x1, x2, nseg, cyclic=c(TRUE, TRUE)),
                 family = fam,
                 data = dat_train)
summary(obj)

grid <- expand.grid(x1 = seq(0, 1, length = 100),
                    x2 = seq(0, 1, length = 100))
pred <- predict(obj, newdata = grid )

ggplot() +
  geom_tile(data = pred, aes(x=x1,y=x2, fill = ypred)) +
  scale_fill_gradientn(colors = topo.colors(100)) +
  coord_fixed() + ggtitle("fitted data LMMsolver torus") +
  JOPS_theme()



