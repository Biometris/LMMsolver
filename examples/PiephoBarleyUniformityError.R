library(LMMsolver)
library(agridat)
data("piepho.barley.uniformity")
dat <- piepho.barley.uniformity
head(dat)

# negative eff. dimension for s(col):
obj1 <- LMMsolve(yield~1,
                 spline=~spl2D(x1=col,x2=row,nseg=c(20, 20)), data=dat,
                 trace=TRUE)
summary(obj1)

# error...
obj2 <- LMMsolve(yield~1,
                 spline=~spl2D(x1=col,x2=row,nseg=c(15, 15),pord=2), data=dat,
                 trace=TRUE)
summary(obj2)


