library(agridat)
library(LMMsolver)
library(ggplot2)
library(dplyr)
data("piepho.barley.uniformity")
dat <- piepho.barley.uniformity
dat$R <- as.factor(dat$row)
dat$C <- as.factor(dat$col)
dat <- filter(dat,!is.na(yield))

obj <- LMMsolve(fixed = yield~1,
                spline = ~spl1D(x = row, nseg = 50, degree = 3, pord = 2) +
                  spl1D(x = col, nseg =  50, degree = 3, pord = 2),
                data = dat,
                trace = FALSE)
summary(obj)

colors <- topo.colors(100)

plt1 <- ggplot(dat)+
  ggtitle(label="raw data") +
  geom_raster(mapping=aes(x=col,y=row,fill=yield))+
  scale_fill_gradientn(name="yield",colours= colors) +
  xlab("column") + ylab("row") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_equal()
plt1

predRow <- obtainSmoothTrend(obj,grid=200,includeIntercept = TRUE, which=1)
predCol <- obtainSmoothTrend(obj,grid=200,includeIntercept = TRUE, which=2)

plt2 = ggplot(data = predRow) + geom_line(aes(x = row, y = ypred))
plt3 = ggplot(data = predCol) + geom_line(aes(x = col, y = ypred))
plt2
plt3

