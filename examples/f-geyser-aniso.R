# Two-dimensional anisotropic density estimation (Old Faithful data)
# A graph from the book 'Practical Smoothing. The Joys of P-splines'
#Paul Eilers and Brian Marx, 2019

library(reshape2)
library(ucminf)
library(JOPS)
library(ggplot2)

# Get the data
data(faithful)
u = faithful[, 1]
v = faithful[, 2]

# Create 2-d histogram
h = hist2d(u, v, xlim = c(1, 6), ylim = c(40, 100))
Y = h$H

z <- as.vector(Y)
x <- rep(h$xgrid, times=100)
y <- rep(h$ygrid, each=100)
dat <- data.frame(x,y,z)
# Turn matrix into a 'long' data frame
#row.names(Y) = h$xgrid
#names(Y) = h$ygrid
#dat <- melt(Y)
#names(dat) = c("x", "y", "z")

#dat2 <- dat[dat$z!=0,]

# Here scaleX is FALSE in spl1D, to be consistent with model in JABES2020 paper.
obj1 <- LMMsolve(fixed = z~1,
                 spline = ~spl2D(x1 = y, x2=x, nseg = c(10,10), degree = 3, pord=2),
                 data = dat,
                 family = poisson(),
                 trace = FALSE)
summary(obj1)

pred <- obtainSmoothTrend(obj1, includeIntercept = TRUE, grid=c(100,100))

family <- poisson()
dense2 <- data.frame(x=pred$x,y=pred$y,z=sqrt(pred$ypred))

# Smooth histogram and return AIC
aicfun = function(lla) {
  lambda = 10^lla
  fit = hist2dsm(Y, lambdax = lambda[1], lambday = lambda[2], dx = 2)
  Mu = fit$Mu
  ok = Y > 0
  dev = 2 * sum(Y[ok] * log(Y[ok]/Mu[ok]))
  ed = fit$ed
  aic = dev + 2 * ed
}

# Search for best (log) lambdas
op = ucminf(c(0, 0), aicfun)
lambda = 10^op$par
cat("log10(lambdas), AIC:", op$par, op$value, "\n")
fit = hist2dsm(Y, lambdax = lambda[1], lambday = lambda[2], dx = 2)
Mu = fit$Mu

# Turn matrix into a 'long' data frame
row.names(Mu) = h$xgrid
names(Mu) = h$ygrid
dens <- melt(sqrt(Mu))
names(dens) = c("x", "y", "z")


# Plot density with contours
plt = ggplot(dense2, aes(x, y, fill = z)) +
  geom_raster(show.legend = F) +
  scale_fill_gradient(high = 'darkgreen', low = 'white') +
  xlab('Eruption length (min)') + ylab('Waiting time (min)') +
  ggtitle('Old Faithful, anisotropic smooth (square root of density)') +
  geom_contour(aes(z = z), color = "steelblue", show.legend = T) +
  geom_point(data=faithful,aes(x=eruptions,y=waiting), inherit.aes = FALSE) +
  JOPS_theme()

# Make and save plot
plot(plt)

plot(x=dens$z,y=dense2$z)
