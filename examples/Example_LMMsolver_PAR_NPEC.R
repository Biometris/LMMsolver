library(tidyr)
library(LMMsolver)

# Example data
input <- read.csv('./examples/PAR_NPEC.csv')
input

x1 <- input$row
x2 <- input$col
x3 <- input$time

# Fit spline
fitting2 <- LMMsolver:::spl3Dfast(y = input$par,
                                  x1 = input$row,
                                  x2 = input$col,
                                  x3 = input$time,
                                  nseg = c(2, 2, 50),
                                  tolerance = 1e-6)


#Prediction grid
g <- expand.grid(x1 = (c(1:35)*24-24),
                 x2 = (c(1:14)*58-58),
                 x3 = seq(min(x3), max(x3),
                          by = 0.02))

# Make predictions
predPAR <- predict(fitting2, g)

# error because different number of predicted vals
g <- cbind(g, pred[[2]])



#### Using new way of calling LMMsolve:

## Specify xlims to replicate results from spl3Dfast.
x1lim <- c(min(x1) - 0.01, max(x1) + 0.01)
x2lim <- c(min(x2) - 0.01, max(x2) + 0.01)
x3lim <- c(min(x3) - 0.01, max(x3) + 0.01)

fittingNew <- LMMsolve(fixed = par~1,
                       spline = ~spl3D(row, col, time, nseg = c(2, 2, 50),
                                       x1lim = x1lim, x2lim = x2lim,
                                       x3lim = x3lim),
                       data = input,
                       tolerance = 1e-6,
                       trace = TRUE)

fitting2$edf
fittingNew$ED

#Prediction grid
g2 <- expand.grid(row = (c(1:35)*24-24),
                  col = (c(1:14)*58-58),
                  time = seq(min(x3), max(x3),
                             by = 0.02))

pred2 <- obtainSmoothTrend(fittingNew, newdata = g2, includeIntercept = TRUE)


