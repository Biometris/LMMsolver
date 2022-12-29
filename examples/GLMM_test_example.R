library(LMMsolver)
library(JOPS)
set.seed(1234)

z1 <- rbinom(n=100, size=100, prob=0.5)
z2 <- rbinom(n=200, size=50,  prob=0.4)
z3 <- rbinom(n=150, size=200, prob=0.35)
z <- c(z1,z2,z3)
range(z)

x <- c(1:100)
y <- sapply(x, FUN= function(x) {sum(x==z)})

dat = data.frame(x = x, y = y)
dim(dat)

obj1 <- LMMsolve(fixed = y~1,
                 spline = ~spl1D(x = x, nseg = 50, degree = 3, pord=2),
                 data = dat,
                 family = poisson(),
                 trace = TRUE)
summary(obj1)

pred <- obtainSmoothTrend(obj1, includeIntercept = TRUE, grid=500)
Fit <- data.frame(x=pred$x,y=pred$ypred)

plt = ggplot(aes(x = x, y = y,  fill = I("wheat3")), data = dat) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0) +
  xlab("x") + ylab("Frequency") +
  ggtitle("histogram") +
  geom_line(data = Fit, col = I("steelblue"), size = 1.0) +
  JOPS_theme()
plt
