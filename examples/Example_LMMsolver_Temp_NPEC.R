library(tidyr)
library(LMMsolver)
library(wesanderson)
library(tidyverse)
library(fields)
library(lubridate)

# Example data
input <- read.csv('./examples/Temp_NPEC.csv')
input

### Define sensor coord #####################
pos <- data.frame(sensor = c(60598, 19445, 4567,
                             49176, 62701, 14061,
                             37681, 35159, 56963),
                  coord = c('AA1', 'AA18', 'AA35',
                            'AH1', 'AH18', 'AH35',
                            'AN1', 'AN18', 'AN35'),
                  rownr = c(35,18,1,
                            35,18,1,
                            35,18,1),
                  colnr = c(14,14,14,
                            8,8,8,
                            1,1,1))

# Sensor positions in cm
rowdist <- 24
coldist <- 58
pos$colm <- (as.numeric(pos$colnr)) * coldist - coldist
pos$rowm <- as.numeric(pos$rownr) * rowdist - rowdist
head(input)
x1 <- input$row
x2 <- input$col
x3 <- input$time

#### Using new way of calling LMMsolve:

## Specify xlims to replicate results from spl3Dfast.
x1lim <- c(min(x1) - 0.01, max(x1) + 0.01)
x2lim <- c(min(x2) - 0.01, max(x2) + 0.01)
x3lim <- c(min(x3) - 0.01, max(x3) + 0.01)
length(unique(input$time))

fittingNew <- LMMsolve(fixed = par~1,
                       spline = ~spl3D(row, col, time, nseg = c(2, 2, 50),
                                       x1lim = x1lim, x2lim = x2lim,
                                       x3lim = x3lim),
                       data = input,
                       tolerance = 1e-6,
                       trace = TRUE)

fittingNew$ED

#Prediction grid
g2 <- expand.grid(row = (c(1:max(pos$rownr))*rowdist-rowdist),
                  col = (c(1:max(pos$colnr))*coldist-coldist),
                  time = seq(min(x3), max(x3),
                             by = 0.02))


pred2 <- obtainSmoothTrend(fittingNew, newdata = g2,
                           includeIntercept = TRUE)

# Plot predictions
for (j in unique(pred2$time)) {
  g2 <- subset(pred2, time == j)
  g2 <- g2 %>%
    arrange(col, row)

  # png with common scale within day
  X1 = matrix(data=g2$ypred, nrow= 35, ncol = 14)
  image.plot(t(X1), x = (unique(g2$col)), y = unique(g2$row),
             col = wes_palette("Zissou1", 100, type = "continuous"),
             xlab='column',ylab='row',
             main = paste0('Temp ', as_datetime(j*100000)),
             zlim = c(min(pred2$ypred),max(pred2$ypred)))
  text(y=pos$rowm, x= pos$colm, col = 'black', pos$coord)
}



