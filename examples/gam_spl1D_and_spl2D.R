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
x1 <- runif(nTrain)
x2 <- runif(nTrain)
x3 <- runif(nTrain)
eps <- rnorm(nTrain, sd=0.1)
y <- f1(x1) + f2(x2, x3) + eps
dat <- data.frame(x1,x2,x3,y)

obj <- LMMsolve(y~1,spline=~spl1D(x1, nseg=17, xlim=c(0,1)) +
                            spl2D(x2, x3, nseg=c(12,17), x1lim=c(0,1), x2lim=c(0,1)),
                            data=dat)
summary(obj)

#
# make predictions on a grid for functions f1 and f2
#

grid_f1 <- expand.grid(x1=seq(0,1,length=200), x2=0.5, x3=0.5)
pred_f1 <- predict(obj, newdata=grid_f1)
p1 <- ggplot(pred_f1, aes(x = x1, y=ypred)) + geom_line() + ggtitle("f1")

grid_f2 <- expand.grid(x1=0.5, x2=seq(0,1,length=200), x3=seq(0,1,length=200))
pred_f2 <- predict(obj, newdata=grid_f2)
p2 <- ggplot(pred_f2, aes(x = x2, y = x3, fill = ypred)) +
  geom_tile(show.legend = TRUE) +
  scale_fill_gradientn(colours = topo.colors(100))+
  coord_fixed() + ggtitle("f2") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

grid.arrange(p1,p2,nrow=1,ncol=2)

#
# make predictions for an independent test set
#
nTest <- 100
x1 <- runif(nTest)
x2 <- runif(nTest)
x3 <- runif(nTest)
y_true <- f1(x1) + f2(x2, x3)
newdat <- data.frame(x1, x2, x3, y_true)
pred <- predict(obj, newdata=newdat)

ggplot(pred, aes(x=y_true, y=ypred)) + geom_point() +
  geom_abline(intercept=0, slope=1, col='blue')

#
# comparison with SOP
#
library(SOP)

# for comparison, don't set the range in LMMsolver...
obj0 <- LMMsolve(y~1,spline=~spl1D(x1, nseg=17) + spl2D(x2, x3, nseg=c(12,17)), data=dat)
obj1 <- sop(formula = y ~ f(x1, nseg = 17) + f(x2,x3, nseg=c(12,17)), data = dat)
pred <- predict(obj0, newdata=newdat, se.fit = TRUE)
pred_sop <- predict(obj1, newdata=newdat, se.fit=TRUE)

deviance(obj0) - deviance(obj1)

range(pred$ypred - pred_sop$fit)
range(pred$se - pred_sop$se.fit)

