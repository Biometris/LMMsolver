library(JOPS)
library(ggplot2)
library(LMMsolver)

nseg <- 20

set.seed(1234)
n <- 500
z <- rnorm(n)
binwidth <- 0.10

# Make the histogram
xmin <- -4
xmax <- 4
brks <- seq(xmin, xmax, by = binwidth)
hst <- hist(z, breaks = brks, plot = FALSE)
x <- hst$mids
y <- hst$counts
dat_org <- data.frame(x=z)

dat <- data.frame(y=y,x=x)

pord <- 3
obj3 <- LMMsolve(fixed = y ~ 1,
                 spline = ~spl1D(x, nseg = 20, degree=3, pord=pord, xlim=c(xmin,xmax)),
                 family = poisson(),
                 data = dat, trace=TRUE)
summary(obj3)

newdat <- data.frame(x = seq(xmin, xmax, length = 300))
pred2 <- predict(obj3, newdata = newdat, se.fit = TRUE)
#pred2$y_true <- f2(pred2$x)

dat_normal <- data.frame(x=z, y= n*binwidth*dnorm(z))

ggplot(data = dat_org, aes(x = x)) +
  geom_histogram(fill = 'wheat3', binwidth=binwidth) +
  geom_line(data = pred2, aes(y = ypred), color = "blue", linewidth = 1) +
  geom_line(data = dat_normal, aes(y = y),color = 'red', linewidth = 1) +
  geom_ribbon(data= pred2, aes(x=x, ymin = ypred-2*se, ymax = ypred+2*se),
              alpha=0.2, inherit.aes = FALSE) + ggtitle(paste("pord=",pord)) +
  theme_bw()


