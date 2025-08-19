library(LMMsolver)
library(ggplot2)
library(dplyr)
library(tictoc)

## simulate some data
alpha <- 0.45
f <- function(x) { 0.3 + 0.01*x + 0.01*x^2 + alpha*cos(pi*x)}

set.seed(12)
nTimePoints <- 25
nObsTimePoints <- 4
n <- nTimePoints*nObsTimePoints
timePoints <- c(1:nTimePoints)
z <- rep(timePoints, each=nObsTimePoints)
sigma2e <- 0.10
y <- f(z) + rnorm(n, sd = sqrt(sigma2e))
dat <- data.frame(z, y)

# some parameters for analysis:
xlim <- c(1,25)
nseg <- 50
tol <- 1.0e-8

## fit the model, using default starting values for theta:
obj1 <- LMMsolve(fixed = y ~ 1,
                spline = ~spl1D(x=z, nseg = nseg, xlim=xlim), data = dat, tolerance = tol)
summary(obj1)
obj1$nIter

newdat <- data.frame(z = seq(1, 25, length = 300))
pred1 <- predict(obj1, newdata = newdat, se.fit = TRUE)
pred1$y_true <- f(pred1$z)

ggplot(data = dat, aes(x = z, y = y)) +
  geom_point(col = "black", size = 1.5) +
  geom_line(data = pred1, aes(y = y_true), color = "red",
            linewidth = 1, linetype ="dashed") +
  geom_line(data = pred1, aes(y = ypred), color = "blue", linewidth = 1) +
  geom_ribbon(data= pred1, aes(x=z, ymin = ypred-2*se, ymax = ypred+2*se),
              alpha=0.2, inherit.aes = FALSE) +
  theme_bw() + ggtitle("Using default penalty theta")

## fit the model, with initial conditions theta_init:
theta_init <- c(100,1)
obj2 <- LMMsolve(fixed = y ~ 1,
                spline = ~spl1D(x=z, nseg = nseg,xlim=xlim), data = dat,
                theta=theta_init, tolerance = tol)
summary(obj2)
obj2$nIter

newdat <- data.frame(z = seq(1, 25, length = 300))
pred2 <- predict(obj2, newdata = newdat, se.fit = TRUE)
pred2$y_true <- f(pred2$z)

ggplot(data = dat, aes(x = z, y = y)) +
  geom_point(col = "black", size = 1.5) +
  geom_line(data = pred2, aes(y = y_true), color = "red",
            linewidth = 1, linetype ="dashed") +
  geom_line(data = pred2, aes(y = ypred), color = "blue", linewidth = 1) +
  geom_ribbon(data= pred2, aes(x=z, ymin = ypred-2*se, ymax = ypred+2*se),
              alpha=0.2, inherit.aes = FALSE) +
  ggtitle(paste("initial value theta smoothing: ", paste(theta_init, collapse=','))) +
  theme_bw()

# the REML estimates are different:
#  the 1st is the global maximum, the 2nd is a local maximum
logLik(obj1, includeConstant = FALSE)
logLik(obj2, includeConstant = FALSE)

# calculate logL on a grid, use log10 scale!
grid1 <- seq(-8, 8, length = 200)
grid2 <- seq(-8, 8, length = 200)
theta1 <- 10^grid1
theta2 <- 10^grid2

theta <- as.matrix(expand.grid(theta1=theta1, theta2=theta2))
dim(theta)

tic("start loglikelihood")
df_logL <- mLogLik(obj1, theta = theta)
toc()
# mv
sum(is.na(df_logL$logL))

# filter for values close to the (local) maxima:
thr <- -15
df2_logL <- df_logL %>% dplyr::filter(logL > thr)
ggplot(df2_logL) +
  geom_tile(aes(x=log10(theta1), y=log10(theta2), fill=logL)) +
  scale_fill_gradientn(colors = topo.colors(100)) +
  annotate("text", x=log10(obj1$theta[1]), y= log10(obj1$theta[2]), colour="black", size=4, label=c('Global \n Maximum')) +
  annotate("text", x=log10(obj2$theta[1]), y= log10(obj2$theta[2]), colour="black", size=4, label=c('Local \n Maximum')) +
  xlab("log10(penalty) smoothing") + ylab("log10(penalty) residual") +
  ggtitle("log-likelihood surface for simulated example")


