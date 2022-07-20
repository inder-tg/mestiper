# -----------------------------------------------------------------------------
# Written on Aug 7, 2018
#
# En este script daremos una breve introducción al paquete raster y presentaremos
# algunas ideas para realizar análisis exploratorio de estructuras espaciales.
# Basado en Ch. 8 Remote sensing image analysis en http://rspatial.org/analysis/rst/9-remotesensing.html
#
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Preamble
# Load needed packages
library(raster)
# -----------------------------------------------------------------------------
# Read images 
rootDir <- getwd()
dirDataSample <- paste(rootDir, "/data/rspatial", sep = "")

# load a single Landsat band

# blue band

b2 <- raster(paste(dirDataSample, "/LC08_044034_20170614_B2.tif", sep = ""))

# green band
b3 <- raster(paste(dirDataSample, "/LC08_044034_20170614_B3.tif", sep = ""))

# red band
b4 <- raster(paste(dirDataSample, "/LC08_044034_20170614_B4.tif", sep = ""))

# near infrared (NIR) band
b5 <- raster(paste(dirDataSample, "/LC08_044034_20170614_B5.tif", sep = ""))

b2

# ---------------------------
# IMAGE info and statistics #
# ---------------------------

# coordinate reference system (CRS)
crs(b2)

# number of rows, cols, cells
dim(b2)

ncell(b2)

# spatial resolution
res(b2)

# number of bands
nlayers(b2)

# compare two rasters
compareRaster(b2, b5)

# ----------------------------------------------------
# CREATE a RasterStack/RasterBrick (multiband image) #
# ----------------------------------------------------

landsatRGB <- stack(b4, b3, b2)
landsatFCC <- stack(b5, b4, b3)

nlayers(landsatRGB)

# --------------------------------------------------------------
# CREATE a RasterStack/RasterBrick (multiband image) from disk #
# --------------------------------------------------------------

rasterList <- paste0(dirDataSample, "/LC08_044034_20170614_B", 1:11, ".tif")

landsat <- stack(rasterList)

nlayers(landsat)

# -----------------------------------------------------------------------------
# Layers: Ultra Blue, Blue, Green, Red, Near Infrared (NIR), Shortwave Infrared (SWIR) 1, 
# Shortwave Infrared (SWIR) 2, Panchromatic, Cirrus, Thermal Infrared (TIRS)1, 
# Thermal Infrared (TIRS)2.
# No se necesitan las últimas 4 bandas, cómo las borramos de "landsat"?
# -----------------------------------------------------------------------------

# ---------------
# VISUALIZATION #
# ---------------

par(mfrow = c(2,2))
plot(b2, main = "Landast Blue", col = gray(0:100 / 100) )
plot(b3, main = "Landast Green", col = gray(0:100 / 100))
plot(b4, main = "Landast Red", col = gray(0:100 / 100))
plot(b5, main = "Landast NIR", col = gray(0:100 / 100))

# ------------------------------------
# COMBINE bands to get better images #
# ------------------------------------

# True color composite
par(mfrow = c(1,1))
# par(mar = c(4.5, 5, 1.5, 2))
plotRGB(landsatRGB, r = 1, g = 2, b = 3, axes = TRUE, stretch = "lin",  
        main = "Landsat True Color Composite")

# False color composite
par(mfrow = c(1, 2))
plotRGB(landsatRGB, r = 1, g = 2, b = 3, axes = TRUE, stretch = "lin",  
        main = "Landsat True Color Composite")
plotRGB(landsatFCC, r = 1, g = 2, b = 3, axes = TRUE, stretch = "lin",  
        main = "Landsat False Color Composite")

# ----------------------------------
# SUBSET and RENAME spectral bands #
# ----------------------------------

# select first 3 bands ONLY
landsatsub1 <- subset(landsat, 1:3)
# same 
landsatsub2 <- landsat[[1:3]]

compareRaster(landsatsub1, landsatsub2)

# Para remover las 4 últimas bandas de "landsat" que no necesitamos:
landsat <- subset(landsat, 1:7)

# Rename bands:
names(landsat)

names(landsat) <- c("ultra-blue", "blue", "green", "red", "NIR", "SWIR1", "SWIR2")

names(landsat)

# --------------------------
# Spatial subset or "crop" #
# --------------------------

extent(landsat)

newExtent <- extent(624387, 635752, 4200047, 4210939)

landsatCrop <- crop(landsat, newExtent)

# par(mfrow=c(1,1))
# plotRGB(landsatCrop, 1,2,3, stretch = "lin")

# saving results to disk
test <- getwd()
setwd(paste(rootDir, "/data/rspatial", sep = ""))
writeRaster(landsatCrop, filename = "cropped-landsat.tif", format = "GTiff", 
            overwrite = TRUE)
setwd(test)

# ------------------------
# RELATION between bands #
# ------------------------

# plot between ultra-blue and blue bands
pairs( landsatCrop[[1:2]], main = "Scatterplot between Ultra-blue and Blue bands" )

pairs( landsatCrop[[4:5]], main = "Scatterplot between Red and NIR bands" )
# -----------------------------------------------------------------------------

# -----------------------
# EXTRACT raster values #
# -----------------------

# LOAD polygons (with land use land cover info)
poligono <- readRDS(paste(dirDataSample, "/samples.rds", sep = ""))

poligono

# get 300 points samples from polygon
pointSample <- spsample(poligono, 300, type = "random")

# add class info to the point samples
pointSample$class <- over(pointSample, poligono)$class

# extract values with points
dataFrame <- extract(landsat, pointSample)

# print the reflectance values
head(dataFrame)

# -------------------
# SPECTRAL Profiles #
# -------------------

ms <- aggregate(dataFrame, list(pointSample$class), mean)

rownames(ms) <- ms[,1]
ms <- ms[,-1]
ms

# plot of spectral profiles

# vector of color for the land cover classes
mycolor <- c("darkred", "yellow", "burlywood", "cyan", "blue")

# convert ms from a data.frame to a matrix
ms <- as.matrix(ms)

# first create an empty plot
plot(0, ylim = c(0, 0.6), xlim = c(1,7), type = "n", xlab = "Bands", 
     ylab = "Reflectance")

# add the different classes
for (i in 1:nrow(ms)) {
  lines(ms[i,], lwd = 3, col = mycolor[i])
}

# TITLE
title(main = "Spectral Profile from Landsat (subset)", font.main = 2)

# LEGEND
legend("topleft", rownames(ms), cex = 0.8, col = mycolor, lwd = 3, bty = "n" )

# -----------------------------------------------------------------------------

# -----------------------------------------------------
# BASIC mathematical operations (with raster package) #
# -----------------------------------------------------

# Let's compute NDVI

# Esta función calcula el NDVI utilizando como parámetros "imagen" y el número
# de dos bandas, a saber "k" e "i"
vegetationIndex <- function(image, k, i){
  bandK <- image[[k]] # get k-th band from image
  bandi <- image[[i]] # get i-th band from image
  index <- (bandK - bandi)/(bandK + bandi)
index
}

# For "landsat" NIR = 5, red = 4
ndvi <- vegetationIndex(landsat, 5, 4)
plot(ndvi, col = rev(terrain.colors(10)), main = "Landsat-NDVI")

# Ejercicio: Adapta la función para calcular índices (vegetationIndex) que remarquen
# agua. Sugerencia: Utiliza el perfil espectral para encontrar las bandas que tiene
# la máxima y mínima reflectancia de agua.

library(RColorBrewer)
display.brewer.all()
myBlues <- brewer.pal(9, "Blues")

waterIndex <- vegetationIndex(landsat, 7, 1)
plot(waterIndex, col = rev(myBlues), main = "Landsat-waterIndex")
# -----------------------------------------------------------------------------

# --------------------------------
# HISTOGRAMS!!!!!!!!!!!!!!!!!!!! #
# --------------------------------

# histogram of ndvi
hist(ndvi)

hist(ndvi,
     main = "Distribution of NDVI values",
     xlab = "NDVI",
     ylab="Frequency",
     col = "wheat",
     xlim = c(-0.5, 1),
     breaks = 30,
     xaxt = 'n')
axis(side=1, at = seq(-0.5,1, 0.05), labels = seq(-0.5,1, 0.05))

# Ejercicio: Haz un histogram de otro índice de vegetación que hayas derivado

# --------------
# Thresholding #
# --------------

# Assuming that values of ndvi above 0.4 represent definitely vegetation.
# Plot the vegetation cover 

testVeg <- ndvi
testVeg[testVeg < 0.4] <- NA

plot(testVeg, main = "Vegetation cover")

par(mar = c(4.5, 5, 1.5, 2))
plot(testVeg, main = "Vegetation cover")

veg <- calc(ndvi, function(x){x[x < 0.4]<-NA; x})
plot(veg, main = "Vegetation cover")

compareRaster(veg, testVeg)

# The histogram of ndvi has a "peak" around the interval (0.25, 0.3)
# Let's try to pin point that region

peakVeg <- ndvi
peakVeg[peakVeg > -Inf & peakVeg < 0.25] <- NA
peakVeg[peakVeg > 0.3 & peakVeg < Inf] <- NA
peakVeg[peakVeg >= 0.25 & peakVeg <= 0.3] <- 1

plot(peakVeg)

land <- reclassify(ndvi, c(-Inf, 0.25, NA, 0.25, 0.3, 1, 0.3, Inf, NA))
plot(land, main = "Qué es esto?")

compareRaster(land, peakVeg)

# Overlapping land and landsatRGB
plotRGB(landsatRGB, r = 1, g = 2, b = 3, axes = TRUE, stretch = "lin", 
        main = "Landsat False Color Composite")
plot(land, add = TRUE, legend = FALSE)

# creating classes for different values of NDVI

vegc <- reclassify(veg, c(-Inf, 0.25, 1, 0.25, 0.3, 2, 0.3, 0.4, 3, 0.4, 0.5, 
                          4, 0.5, Inf, 5))
plot(vegc, col = rev(terrain.colors(4)), main = 'NDVI based thresholding')
# -----------------------------------------------------------------------------

