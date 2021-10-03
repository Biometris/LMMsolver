# Martin Boer, Biometris, Wageningen
# example sparse 1D, see Boer et al 2020, JABES for details
# of the example
rm(list=ls())

library(agridat)
library(LMMsolver)
library(ggplot2)

data(john.alpha)
dat <- john.alpha

N <- nrow(dat)

# Here scaleX is FALSE in spl1D, to be consistent with model in JABES2020 paper.
obj1 <- LMMsolve(fixed = yield~rep+gen,
                spline = ~spl1D(x = plot, nseg = N-1, degree = 1, pord= 1, scaleX=FALSE),
                data = dat,
                trace = FALSE,
                tolerance = 1.0e-10,
                omitConstant = FALSE)
summary(obj1)

round(obj1$dev, 2)

# residual variance, see JABES 2020 paper, table 1:
round(obj1$sigma2e, 5)

# genotype random, not in JABES paper
obj2 <- LMMsolve(fixed = yield~rep,
                 random = ~gen,
                spline = ~spl1D(x = plot, nseg = N-1, degree = 1, pord= 1, scaleX=FALSE),
                data = dat,
                trace = FALSE,
                tolerance = 1.0e-10)
summary(obj2)

obj2$EDmax
obj2$EDnominal

diagnosticsMME(obj2)
displayMME(obj2, cholesky=FALSE)
displayMME(obj2, cholesky=TRUE)

# obtain spatial trend with genotype fixed:
plotDat <- obtainSmoothTrend(obj1, grid=100)
head(plotDat)
ggplot(plotDat, aes(x = plot, y = ypred)) +
  geom_tile(show.legend = TRUE) +
  geom_line() +
  labs(title = "Spatial trend for the oats data", x = "plot", y = "plot effects") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())




