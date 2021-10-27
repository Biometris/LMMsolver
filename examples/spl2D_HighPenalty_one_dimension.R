library(LMMsolver)
library(agridat)

data("piepho.barley.uniformity")
dat <- piepho.barley.uniformity
head(dat)

Nseg <- seq(10, 80, by=10)

for (i in Nseg){
  s <- proc.time()[3]
  nseg <- c(i,i)
  obj <- LMMsolve(yield~1,
                  spline=~spl2D(x1=col, x2=row, nseg=nseg),
                  data=dat)
  e <- proc.time()[3]
  tbl <- summary(obj)
  pen <- tbl$Penalty[3]
  cat(sprintf("%4d %4d %8.3g %8.1f  \n", i, obj$nIter, pen, e-s))
}


