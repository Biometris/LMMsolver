library(LMMsolver)

set.seed(1234)
x <- seq(-2, 2, by=0.1)
y <- x^3 - x + rnorm(n=length(x), sd=0.01)

plot(x=x,y=y)
abline(h=0)
dat <- data.frame(x=x, y=y)
dim(dat)

obj <- LMMsolver::LMMsolve(fixed = y~ 1,
                           spline = ~spl1D(x = x, nseg = 30,
                                           pord = 2, degree = 3,
                                           scaleX = FALSE),
                           data = dat, trace=TRUE)
summary(obj)

x0 <- seq(-2, 2, by=0.01)
newdat <- data.frame(x=x0)
pred <- obtainSmoothTrend(obj, newdata = newdat, includeIntercept = TRUE)
pred1 <- obtainSmoothTrend(obj, newdata = newdat, deriv = 1)
pred2 <- obtainSmoothTrend(obj, newdata = newdat, deriv = 2)

plot(dat$x, dat$y, main='x^3-x')
lines(pred$x, pred$ypred, type = "l", col ='red')
plot(pred1$x, pred1$ypred, type = "l",main='first derivative: 3x^2-1')
lines(x=x0,y=3*x0^2-1)
plot(pred2$x, pred2$ypred, type = "l",main='second derivative: 6x')
lines(x=x0,y=6*x0, col='blue')

