library(ggplot2)
library(SemiPar)
library(gridExtra)
library(JOPS)
library(LMMsolver)

# Get the data
data(ethanol)

head(ethanol)

dat <- ethanol

obj <- LMMsolve(fixed = NOx~1,
                spline = ~spl1D(x = E, nseg = 100, degree = 3, pord = 2) +
                          spl1D(x = C, nseg =  50, degree = 3, pord = 2),
                 data = dat,
                 trace = FALSE)

summary(obj)
cf <- coef(obj)

displayMME(obj, cholesky=FALSE)
displayMME(obj, cholesky=TRUE)

predE <- obtainSmoothTrend(obj,grid=200,includeIntercept = TRUE, which=1)
predC <- obtainSmoothTrend(obj,grid=200,includeIntercept = TRUE, which=2)

# Add fitted components to data frame (for ggplot)
Fmod = ethanol
Fmod$f1 = obtainSmoothTrend(obj,newdata=ethanol,includeIntercept = FALSE, which=2)$ypred
Fmod$f2 = obtainSmoothTrend(obj,newdata=ethanol,includeIntercept = FALSE, which=1)$ypred

mu = cf$`(Intercept)`

# Create plots
plt1 = ggplot(aes(x = C, y = E), data = ethanol) +
  geom_point(color = "darkgrey") +
  xlab("Compression ratio (C)") +
  ylab("Equivalence ratio (E)") +
  ggtitle("Experiment design") +
  JOPS_theme()

plt2 = ggplot(aes(x = E, y = NOx - f1), data = Fmod) +
  geom_point(color = "darkgrey") +
  geom_line(aes(x = E, y = ypred), data = predE, size = 1, color = "blue") +
  xlab("Equivalence ratio") +
  ylab("Partial residuals") +
  ggtitle("Partial response") +
  JOPS_theme()

plt3 = ggplot(aes(x = C, y = NOx - f2), data = Fmod) +
  geom_point(color = "darkgrey") +
  geom_line(aes(x = C, y = ypred), data = predC, color = "blue", size = 1) +
  xlab("Compression ratio (C)") +
  ylab("Partial residuals") +
  ggtitle("Partial response") +
  JOPS_theme()

plt4 = ggplot(aes(x = f1 + f2 + mu, y = NOx), data = Fmod) +
  geom_point(color = "darkgrey") +
  geom_abline(slope = 1, intercept = 0, color = "blue", size = 1) +
  xlab("Fitted NOx") +
  ylab("Observed ") +
  ggtitle("Compare fit to data") +
  JOPS_theme()

grid.arrange(plt1, plt2, plt3, plt4, ncol = 2, nrow = 2)

