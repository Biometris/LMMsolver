library(LMMsolver)
library(ggplot2)
library(gridExtra)

set.seed(1234)

f1 <- function(x) {sin(2*pi*x) }

f2 <- function(x, y) { # Lee et al 2013
  y <- (cos(2*pi*sqrt((x - 0.5)^2 + (y - 0.5)^2)) + 0.2)/0.5
  y
}

nTrain <- 500
nTest <- 100
nTot <- nTrain + nTest
x1 <- runif(nTot)
x2 <- runif(nTot)
x3 <- runif(nTot)
eps <- rnorm(nTot, sd=0.1)
y <- f1(x1) + f2(x2, x3) + eps

dat_all <- data.frame(x1,x2,x3,y)
dat <- dat_all[1:nTrain,]
newdat <- dat_all[-c(1:nTrain),]

obj <- LMMsolve(y~1,spline=~spl1D(x1, nseg=15) +
                            spl2D(x2, x3, nseg=c(15,15)),
                            data=dat)
summary(obj)

grid <- seq(0,1, length=200)
grid_f1 <- expand.grid(x1=grid, x2=0.5, x3=0.5)
grid_f2 <- expand.grid(x1=0.5, x2=grid, x3=grid)
pred1 <- predict(obj, newdata=grid_f1)
pred2 <- predict(obj, newdata=grid_f2)

p1 <- ggplot(pred1, aes(x = x1, y=ypred)) + geom_line() + ggtitle("f1")

p2 <- ggplot(pred2, aes(x = x2, y = x3, fill = ypred)) +
  geom_tile(show.legend = TRUE) +
  scale_fill_gradientn(colours = topo.colors(100))+
  coord_fixed() + ggtitle("f2") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

grid.arrange(p1,p2,nrow=1,ncol=2)

pred <- predict(obj, newdata=newdat)

ggplot(pred, aes(x=y, y=ypred)) + geom_point() +
         geom_abline(intercept=0, slope=1, col='red')



