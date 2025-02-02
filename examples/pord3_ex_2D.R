library(LMMsolver)
library(spam)
library(ggplot2)
library(tictoc)
## Get precipitation data from spam
data(USprecip, package = "spam")

## Only use observed data
USprecip <- as.data.frame(USprecip)
USprecip <- USprecip[USprecip$infill == 1, ]

pord <- 3

tic("LMMsolver")
obj5 <- LMMsolve(fixed = anomaly ~ 1,
                 spline = ~spl2D(x1 = lon, x2 = lat, nseg = c(41, 41), pord=pord),
                 data = USprecip)
toc()
summary(obj5)

lon_range <- range(USprecip$lon)
lat_range <- range(USprecip$lat)
newdat <- expand.grid(lon = seq(lon_range[1], lon_range[2], length = 200),
                      lat = seq(lat_range[1], lat_range[2], length = 300))
plotDat5 <- predict(obj5, newdata = newdat)

plotDat5 <- sf::st_as_sf(plotDat5, coords = c("lon", "lat"))
usa <- sf::st_as_sf(maps::map("usa", regions = "main", plot = FALSE))
sf::st_crs(usa) <- sf::st_crs(plotDat5)
intersection <- sf::st_intersects(plotDat5, usa)
plotDat5 <- plotDat5[!is.na(as.numeric(intersection)), ]

ggplot(usa) +
  geom_sf(color = NA) +
  geom_tile(data = plotDat5,
            mapping = aes(geometry = geometry, fill = ypred),
            linewidth = 0,
            stat = "sf_coordinates") +
  scale_fill_gradientn(colors = topo.colors(100))+
  labs(title = paste("Precipitation (anomaly) , pord = ", pord),
       x = "Longitude", y = "Latitude") +
  coord_sf() +
  theme(panel.grid = element_blank())
