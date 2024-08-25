# note: run usethis::use_logo() to generate the generate correct size

library(hexSticker)
library(ggplot2)
library(LMMsolver)
library(spam)
library(fields)
library(maps)
library(dplyr)
library(gridExtra)

set.seed(12)
n <- 150
x <- seq(0, 1, length = n)
sigma2e <- 0.04
y <- f2(x) + rnorm(n, sd = sqrt(sigma2e))
dat2 <- data.frame(x, y)


obj2 <- LMMsolve(fixed = y ~ 1,
                 spline = ~spl1D(x, nseg = 50),
                 data = dat2)

newdat <- data.frame(x = seq(0, 1, length = 300))
pred2 <- predict(obj2, newdata = newdat, se.fit = TRUE)
pred2$y_true <- f2(pred2$x)

p <- ggplot(data = dat2, aes(x = x, y = y)) +
  geom_point(col='black', size = 0.2) +
  geom_line(data = pred2, aes(y = ypred), color = "red", linewidth = 0.5) +
  geom_ribbon(data= pred2, aes(x=x, ymin = ypred-2*se, ymax = ypred+2*se),
              alpha=0.50, inherit.aes = FALSE) +
  theme_bw()

p <- p + theme_void() + theme_transparent()
p

sticker(p, package="LMMsolver", p_size=20, s_x=1.0, s_y=.85, s_width=1.2, s_height=0.9,
        h_fill = "lightyellow", p_color="blue",h_color="green",
        filename="LMMsolver_hexSticker.png")


