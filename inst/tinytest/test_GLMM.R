## Generate some random data following a poisson distribution.
set.seed(1234)

z1 <- rbinom(n=100, size=100, prob=0.5)
z2 <- rbinom(n=200, size=50,  prob=0.4)
z3 <- rbinom(n=150, size=200, prob=0.35)
z <- c(z1, z2, z3)
x <- c(1:100)
y <- sapply(x, FUN= function(x) {sum(x == z)})

dat = data.frame(x = x, y = y)

## Fit simple model.
obj1 <- LMMsolve(fixed = y ~ 1,
                 spline = ~spl1D(x = x, nseg = 10, degree = 3, pord=2),
                 data = dat,
                 family = poisson())

expect_equal_to_reference(obj1, "GLMMFull", check.attributes = FALSE)

## Check predictions.
pred1 <- obtainSmoothTrend(obj1, grid = 5)

expect_equal(pred1$ypred, c(1.735118429947, 3588.132243336660,
                            5337.660426348466, 3904.681375838120, 0.576329536209))
