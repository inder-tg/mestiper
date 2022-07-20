
# Written by Inder Tecuapetla, April 19, 2018
# Some examples of interpolation, applied harmonic analysis, trend analysis,
# structural change estimation on MODIS and Landsat data cubes

# -----------------------------------------------------------------------------
library(fields)
library(raster)
library(RColorBrewer)
library(tools)
# -----------------------------------------------------------------------------

rotate <- function(x) t(apply(x, 2, rev))

matrixToRaster <- function(matrix, raster){
  rasterTable <- data.frame(rasterToPoints(raster)) 
  
  temp_matFinal <- c(matrix) 
  
  df <- data.frame(x = rasterTable$x, y = rasterTable$y, values = temp_matFinal)
  
  coordinates(df) <- ~ x + y
  
  gridded(df) <- TRUE
  
  raster_df <- raster(df)
  
  projection(raster_df) <- projection(raster)
  
  raster_df
}

# =============================================================================
# MODIS
# =============================================================================

# -----------------------------------------------------------------------------
# IMPORTACION y VISUALIZACION: El caso .bin
# -----------------------------------------------------------------------------

dirQRoo <- paste(rootDir, "/data/site_QRoo", sep = "" )

masterFile <- raster(paste(dirQRoo, "/master_QuintanaRoo.tif", sep = ""))

nRow <- nrow(masterFile) # 500
nCol <- ncol(masterFile) # 600

listFiles <- list.files(path = dirQRoo, pattern = ".bin")

listFiles <- listFiles[-1]

# ----------------------------------------------------
# Importando un archivo .bin y exportando una imagen #
# ----------------------------------------------------

x <- 1:nCol
y <- 1:nRow

k <- 300
to.read <- file(paste(dirQRoo, "/", listFiles[k], sep = ""), "rb")

imagen <- matrix(NA, ncol = nCol, nrow = nRow, byrow = TRUE)
for(r in 1:nRow){
  imagen[r,] <- readBin(to.read, integer(), size = 2, n = nCol, endian = "little") / 1e4
}
close(to.read)

imagenRotada <- rotate(imagen)

image.plot(x, y, imagenRotada, main = paste(listFiles[k], sep = ""), useRaster = TRUE,
           xlab = "", ylab = "")

colorPalette <- designer.colors(n = 100, col = c("wheat4", "wheat", "forestgreen", "gray20"))

image.plot(x, y, imagenRotada, main = paste(listFiles[k], sep = ""), useRaster = TRUE,
           col = colorPalette, xlab = "", ylab = "")

myBlues <- brewer.pal(4, "Blues")

image.plot(x, y, imagenRotada, main = paste(listFiles[k], sep = ""), useRaster = TRUE,
           col = myBlues, xlab = "", ylab = "")

myGreens <- brewer.pal(4, "Greens")

image.plot(x, y, imagenRotada, main = paste(listFiles[k], sep = ""), useRaster = TRUE,
           col = myGreens, xlab = "", ylab = "")

myPalette <- c("#edf8fb", "#ccece6", "#99d8c9", "#66c2a4", "#2ca25f", "#006d2c")
image.plot(x, y, imagenRotada, main = paste(listFiles[k], sep = ""), useRaster = TRUE,
           col = myPalette, xlab = "", ylab = "")

# EXPORTANDO imagen
png(paste(rootDir, "/png_site_QRoo/", listFiles[k], ".png", sep = ""))
image.plot(x, y, imagenRotada, main = paste(listFiles[k], sep = ""), useRaster = TRUE,
           col = colorPalette, xlab = "", ylab = "")
graphics.off()

# EXPORTANDO imagen en formato .tif

# test <- matrixToRaster(apply(imagen,1, rev), masterFile)

# test <- matrixToRaster(imagenRotada, masterFile)
# 
# test <- matrixToRaster(apply(t(imagenRotada), 2, rev), masterFile)
# 
# test <- matrixToRaster(t(imagenRotada), masterFile)

image.plot(x, y, imagenRotada, main = paste(listFiles[k], sep = ""), useRaster = TRUE,
           xlab = "", ylab = "")
 
test <- matrixToRaster(t(imagen), masterFile)

plot(test, main = paste(listFiles[k], sep = ""), useRaster = TRUE,
     xlab = "", ylab = "")

writeRaster(test, filename = paste(rootDir, "/data/site_QRoo_tif/", 
                                   file_path_sans_ext(listFiles[k]), ".tif", 
                                   sep = ""), format = "GTiff")

# -----------------------------------------------------------------------------

# ----------------------------------------------
# Importando un set de imagens en formato .bin #
# ----------------------------------------------

nLayers <- length(listFiles)
ndvi_site_QRoo <- array(NA, c(nRow, nCol, nLayers))

for(k in 1:nLayers) {
  to.read <- file(paste(dirQRoo, "/", listFiles[k], sep = ""), "rb")
  
  imagen <- matrix(NA, ncol = nCol, nrow = nRow, byrow = TRUE)
  for(r in 1:nRow){
    imagen[r,] <- readBin(to.read, integer(), size = 2, n = nCol, endian = "little") / 1e4
  }
  close(to.read)
  
  ndvi_site_QRoo[, , k] <- imagen
  
  imagenRotada <- rotate(imagen)
  # 
  png(paste(rootDir, "/png_site_QRoo/", listFiles[k], ".png", sep = ""))
  image.plot(x, y, imagenRotada, main = paste(listFiles[k], sep = ""), useRaster = TRUE,
             col = colorPalette, xlab = "", ylab = "")
  graphics.off()
  
  temp <- matrixToRaster(t(imagen), masterFile)
  
  writeRaster(temp, filename = paste(rootDir, "/data/site_QRoo_tif/", 
                                     file_path_sans_ext(listFiles[k]), ".tif", 
                                     sep = ""), format = "GTiff")
}

str(ndvi_site_QRoo)


# ------------------------
# HACIENDO una animación #
# ------------------------

dirData <- paste( rootDir, "/data/site_QRoo_tif", sep = "")

listFiles <- list.files(path = dirData, pattern = ".tif")

brik <- brick()
setwd(dirData)
for(i in 1:length(listFiles)){
  importedLayer <- raster(listFiles[i])
  brik <- addLayer(brik, importedLayer)
}

brik

?animate

animate(brik, n = 1)





