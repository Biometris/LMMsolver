library(LMMsolver)
library(ggplot2)
library(gridExtra)

set.seed(1234)

f1 <- function(x) {sin(2*pi*x) }

f2 <- function(x, y) { # Lee et al 2013
  y <- (cos(2*pi*sqrt((x - 0.5)^2 + (y - 0.5)^2)) + 0.2)/0.5
  y
}

N <- 500
x1 <- runif(N)
x2 <- runif(N)
x3 <- runif(N)
eps <- rnorm(N, sd=0.1)
y <- f1(x1) + f2(x2, x3) + eps

dat <- data.frame(x1,x2,x3,y)

obj <- LMMsolve(y~1,spline=~spl1D(x1, nseg=15) +
                            spl2D(x2, x3, nseg=c(15,15)),
                            data=dat)
summary(obj)

pred1 <- obtainSmoothTrend(obj, grid=c(200), which=1)
pred2 <- obtainSmoothTrend(obj, grid=c(200, 200), which=2)

p1 <- ggplot(pred1, aes(x = x1, y=ypred)) + geom_line() + ggtitle("f1")

p2 <- ggplot(pred2, aes(x = x2, y = x3, fill = ypred)) +
  geom_tile(show.legend = TRUE) +
  scale_fill_gradientn(colours = topo.colors(100))+
  coord_fixed() + ggtitle("f2") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

grid.arrange(p1,p2,nrow=1,ncol=2)

