library(LMMsolver)
library(ggplot2)
library(dplyr)
library(tictoc)

## simulate some data
alpha <- 1.0
f <- function(x) { 0.3 + 0.01*x + 0.01*x^2 + alpha*cos(pi*x)}

set.seed(12)
nTimePoints <- 50
nObsTimePoints <- 3
n <- nTimePoints*nObsTimePoints
timePoints <- c(1:nTimePoints)
z <- rep(timePoints, each=nObsTimePoints)
sigma2e <- 0.15
y <- f(z) + rnorm(n, sd = sqrt(sigma2e))
dat <- data.frame(z, y)

# some parameters for analysis:
xlim <- c(1,50)
nseg <- 100

## fit the model, using default starting values for theta:
obj1 <- LMMsolve(fixed = y ~ 1,
                spline = ~spl1D(x=z, nseg = nseg, xlim=xlim), data = dat)
summary(obj1)

newdat <- data.frame(z = seq(1, 50, length = 300))
pred1 <- predict(obj1, newdata = newdat, se.fit = TRUE)
pred1$y_true <- f(pred1$z)

ggplot(data = dat, aes(x = z, y = y)) +
  geom_point(col = "black", size = 1.5) +
  geom_line(data = pred1, aes(y = y_true), color = "red",
            linewidth = 1, linetype ="dashed") +
  geom_line(data = pred1, aes(y = ypred), color = "blue", linewidth = 1) +
  geom_ribbon(data= pred1, aes(x=z, ymin = ypred-2*se, ymax = ypred+2*se),
              alpha=0.2, inherit.aes = FALSE) +
  theme_bw() + ggtitle("Using default theta")

## fit the model
theta_init <- c(100,1)
obj2 <- LMMsolve(fixed = y ~ 1,
                spline = ~spl1D(x=z, nseg = nseg,xlim=xlim), data = dat,
                theta=theta_init,trace=TRUE)
summary(obj2)

newdat <- data.frame(z = seq(1, 50, length = 300))
pred2 <- predict(obj2, newdata = newdat, se.fit = TRUE)
pred2$y_true <- f(pred2$z)

ggplot(data = dat, aes(x = z, y = y)) +
  geom_point(col = "black", size = 1.5) +
  geom_line(data = pred2, aes(y = y_true), color = "red",
            linewidth = 1, linetype ="dashed") +
  geom_line(data = pred2, aes(y = ypred), color = "blue", linewidth = 1) +
  geom_ribbon(data= pred2, aes(x=z, ymin = ypred-2*se, ymax = ypred+2*se),
              alpha=0.2, inherit.aes = FALSE) +
  ggtitle(paste("initial value theta smoothing: ", theta_init[1])) +
  theme_bw()

obj1$logL
obj2$logL

# calculate logL on a grid:
grid1 <- seq(-3,5, by = 0.05)
grid2 <- seq(-2,2, by = 0.05)
theta1 <- 10^grid1
theta2 <- 10^grid2

thetaM <- as.matrix(expand.grid(theta1=theta1,theta2=theta2))
grid <- expand.grid(grid1=grid1,grid2=grid2)
dim(thetaM)

tic("start loglikelihood")
object <- obj2
X <- object$X
Z <- object$Z
C <- object$C # needed, or rebuild anyway?
lGinv <- object$lGinv
lRinv <- object$lRinv
y <- object$y
logL <- LMMsolver:::logLikelihood(y,X,Z,lGinv,lRinv, theta = thetaM)
toc()
df <- cbind(grid, logL)

# define points for in the plot:
pnt1 <- log10(obj1$theta) # global max
pnt2 <- log10(obj2$theta) # local max

df2 <- df %>% dplyr::filter(logL > -130)
ggplot(df2) +
  geom_tile(aes(x=grid1, y=grid2, fill=logL)) +
  scale_fill_gradientn(colors = topo.colors(100)) +
  annotate("text", x=pnt1[1], y=pnt1[2], colour="black", size=4, label=c('Global \n Maximum')) +
  annotate("text", x=pnt2[1], y=pnt2[2], colour="black", size=4, label=c('Local \n Maximum')) +
  xlab("log10(penalty) smoothing") + ylab("log10(penalty) residual") +
  ggtitle("log-likelihood surface for simulated example")



