## ------------------------------------------------------------------------
library(red)
data(red.records)
spRec = red.records

## ------------------------------------------------------------------------
eoo(spRec)
aoo(spRec)
countries(spRec)

## ------------------------------------------------------------------------
data(red.layers)
spRaster <- red.layers
altRaster <- spRaster[[3]]
elevation(spRec, altRaster)

## ----fig.width=8, fig.height = 6-----------------------------------------
raster::plot(altRaster)
points(spRec)

## ----fig.width=8, fig.height = 6-----------------------------------------
raster::plot(altRaster)
points(spRec, col="red")
spRec = move(spRec, spRaster)
points(spRec)

## ----fig.width=8, fig.height = 4-----------------------------------------
outliers(spRec, spRaster)

## ------------------------------------------------------------------------
distrRaw <- map.points(spRec, spRaster)

## ----fig.width=8, fig.height = 6-----------------------------------------
raster::plot(distrRaw[[1]])
distrRaw[[2]]

## ----fig.width=8, fig.height = 6-----------------------------------------
spHabitat <- spRaster[[4]]
spHabitat[spHabitat != 4] <- 0
spHabitat[spHabitat == 4] <- 1
raster::plot(spHabitat, legend = FALSE)
points(spRec)

## ----fig.width=8, fig.height = 6-----------------------------------------
distrHabitat <- map.habitat(spRec, spHabitat, move = FALSE)
raster::plot(distrHabitat[[1]], legend = FALSE)

## ------------------------------------------------------------------------
distrHabitat[[2]]

## ----fig.width=8, fig.height = 6-----------------------------------------
thinRec <- thin(spRec, 0.1)
plot(spRec, col = "red")
points(thinRec)

## ----fig.width=8, fig.height = 6-----------------------------------------
spRaster <- raster.reduce(red.layers[[1:3]], n = 2)
raster::plot(spRaster)

## ----fig.width=8, fig.height = 6-----------------------------------------
par(mfrow = c(1,1))
map.draw(spRec, distrHabitat[[1]])
kml(distrHabitat[[1]], filename = "spMap.kml")

## ----fig.width=8, fig.height = 6-----------------------------------------
rliData <- matrix(c("LC","LC","EN","EN","EX","EX","LC","CR","CR","EX"), ncol = 2, byrow = TRUE)
colnames(rliData) <- c("2000", "2010")
rliData
rli(rliData)
rliData <- cbind(c("Arthropods","Arthropods","Birds","Birds","Birds"), rliData)
rliData
rli.multi(rliData, boot = TRUE)

