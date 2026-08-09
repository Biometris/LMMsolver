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

# calculate logL on a grid, use log10 scale!
grid1 <- seq(-2, 2, length = 5)
grid2 <- seq(-2, 2, length = 5)
theta1 <- 10^grid1
theta2 <- 10^grid2

theta <- as.matrix(expand.grid(theta1=theta1, theta2=theta2))
dim(theta)

df_logL <- mLogLik(obj1, theta = theta)
expect_equivalent_to_reference(df_logL, "mLogLik0")

expect_error(mLogLik(obj1, theta = theta[,-1,drop=FALSE]),
             "theta has wrong number of columns")
expect_error(mLogLik(obj1, theta = "x"),
             "theta should be matrix")

obj2 <- LMMsolve(fixed = y ~ z, data=dat)
theta <- matrix(data = obj2$theta, nrow=1,ncol=1)
r <- mLogLik(obj2, theta = theta)
expect_equal(obj2$theta, r$theta1)
expect_equal(0, r$dL_dtheta1)
expect_equal(obj2$logL, r$logL)
