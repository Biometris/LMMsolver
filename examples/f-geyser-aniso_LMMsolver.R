# Martin Boer, 2022, 2D density estimation using LMMsolver

# Two-dimensional anisotropic density estimation (Old Faithful data)
# A graph from the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(reshape2)
library(ucminf)
library(JOPS)
library(ggplot2)
library(LMMsolver)

# Get the data
data(faithful)
u = faithful[, 1]
v = faithful[, 2]

# Create 2-d histogram
nbins <- c(100, 100)
h <- hist2d(u, v, nb=nbins, xlim = c(1, 6), ylim = c(40, 100))
Y <- h$H

z <- as.vector(Y)
x <- rep(h$xgrid, times=nbins[2])
y <- rep(h$ygrid, each=nbins[1])
dat <- data.frame(x, y, z)

nseg <- c(10, 10)
obj1 <- LMMsolve(fixed = z~1,
                 spline = ~spl2D(x1 = x, x2=y, nseg = nseg),
                 data = dat,
                 family = poisson())
summary(obj1)

pred <- obtainSmoothTrend(obj1, includeIntercept = TRUE, grid=nbins)
densLMMsolver <- data.frame(x=pred$x,y=pred$y,z=sqrt(pred$ypred))

# Smooth histogram and return AIC
aicfun = function(lla) {
  lambda = 10^lla
  fit = hist2dsm(Y, nsegx=nseg[1],nsegy=nseg[2],
                 lambdax = lambda[1], lambday = lambda[2], dx = 2)
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
fit = hist2dsm(Y, nsegx=nseg[1], nsegy=nseg[2],
               lambdax = lambda[1], lambday = lambda[2], dx = 2)
Mu = fit$Mu

# Turn matrix into a 'long' data frame
z <- as.vector(t(Mu))
x <- rep(h$xgrid, each=nbins[2])
y <- rep(h$ygrid, times=nbins[1])
dens <- data.frame(x=x,y=y,z=sqrt(z))

# compare LMMsolver with code in JOPS book
plot(x=dens$z,y=densLMMsolver$z, main="compare methods",
     xlab="z (JOPS)",ylab='z (LMMsolver)')

# Plot density with contours using LMMsolver results:
plt = ggplot(densLMMsolver, aes(x, y, fill = z)) +
  geom_raster(show.legend = F) +
  scale_fill_gradient(high = 'darkgreen', low = 'white') +
  xlab('Eruption length (min)') + ylab('Waiting time (min)') +
  ggtitle('Old Faithful, anisotropic smooth (square root of density)') +
  geom_contour(aes(z = z), color = "steelblue", show.legend = T) +
  geom_point(data=faithful,aes(x=eruptions,y=waiting), inherit.aes = FALSE) +
  JOPS_theme()

# Make and save plot
plot(plt)

