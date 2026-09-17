
# --- Actualizado por Inder Tecuapetla, Marzo 9, 2023 ---
# --- Primera versión introTiSeG.R basada en Colditz et al 2008.
# --- Segunda versión intro_interpolation.R
# --- Tercera versión missingData.R, Feb 24, 2022

library(imputeTS)
library(terra)
# library(raster)
# library(numbers)
library(geoTS)
library(foreach)
library(doParallel)

library(ggplot2)
library(tidyverse)
library(tidyquant)
library(here)
library(kableExtra)

source(paste0(getwd(), "/Rscripts/auxFUN.R"))

MSE <- function(true, estimate){ mean( (true - estimate)^2 ) }

# -----------------------------------------------------------------------------

# --- Aplication a PR ---

dirDATA <- here( "data", "interpol" )  # paste0( getwd(), "/mestiper/data/interpol" )

ndviList <- list.files(path=dirDATA, pattern="",
                       full.names = TRUE)

ndvi <- t(sapply(1:length(ndviList), function(s) readRDS(ndviList[s])))

yr <- range(ndvi)
yr[1] <- yr[1]-0.1
yr[2] <- yr[2]+0.1

par(mar=c(4,4,1,2))
plot(ndvi[1,], type="l", ylim=yr, ylab="NDVI", xlab="Obs.", lwd=3)
lines(ndvi[2,], col="red", lwd=3)
lines(ndvi[3,], col="blue", lwd=3)
lines(ndvi[4,], col="green", lwd=3)
lines(ndvi[5,], col="darkorange", lwd=3)

# --- simulacion de datos faltantes

plot(ndvi[1,], type="l", col="red", ylim=yr, 
     ylab="NDVI", xlab="Obs.", lwd=3)
points(ndvi[1,], pch=16, col="red")

mask_md <- c(3,10,17)

ndvi_mask <- ndvi[1,]
ndvi_mask[mask_md] <- NA

plot(ndvi_mask, type="l", col="red", 
     ylim=yr, ylab="NDVI", xlab="Obs.", lwd=3)
points(ndvi_mask, pch=16, col="red")

ndvi_interpol1 <- approx(x=ndvi_mask, method="constant", f=0, n=23)
ndvi_interpol2 <- approx(x=ndvi_mask, method="constant", f=1, n=23)
ndvi_interpol3 <- approx(x=ndvi_mask, method="constant", f=0.5, n=23)
ndvi_interpol4 <- na_interpolation(x=ndvi_mask)
ndvi_interpol5 <- na_interpolation(x=ndvi_mask, option="spline")

plot(ndvi_mask, type="l", col="red", lwd=3,
     ylim=yr, ylab="NDVI", xlab="Obs.",
     main="Previo")
points(ndvi_interpol1$y, pch=16, col="blue")

plot(ndvi_mask, type="l", col="red", lwd=3,
     ylim=yr, ylab="NDVI", xlab="Obs.",
     main="Siguiente")
points(ndvi_interpol2$y, pch=16, col="darkorange")

plot(ndvi_mask, type="l", col="red", lwd=3,
     ylim=yr, ylab="NDVI", xlab="Obs.",
     main="Punto medio")
points(ndvi_interpol3$y, pch=16, col="purple")

plot(ndvi_mask, type="l", col="red", ylim=yr, ylab="NDVI", xlab="Obs.", lwd=3,
     main="Lineal")
points(ndvi_interpol4, pch=16, col="darkgreen")

plot(ndvi_mask, type="l", col="red", ylim=yr, ylab="NDVI", xlab="Obs.", lwd=3,
     main="Spline")
points(ndvi_interpol5, pch=16, col="lightseagreen")

# --- Análisis cuantitativo

mask <- matrix(nrow = 5, ncol = 11)
mask[1,] <- c(2:3, 6:7, 10, 13:16, 21:22)
mask[2,] <- c(4:5, 8:12, 15, 17, 19:20)
mask[3,] <- c(4, 7:8, 14:17, 20, 21:23)
mask[4,] <- c(2,4, 7:8, 13:16, 19, 21:22)
mask[5,] <- c(2, 4:5, 7:8, 16:19, 21:22)

estimatePrevious <- matrix(nrow = 5, ncol = 23)
estimateNext <- matrix(nrow = 5, ncol = 23)
estimateMedio <- matrix(nrow = 5, ncol = 23)
estimateLine <- matrix(nrow = 5, ncol = 23)
estimateSpline <- matrix(nrow = 5, ncol = 23)

# variando el patrón de datos faltantes
# dataset is fixed

dataTemp <- ndvi[2,]

plot(dataTemp, type="l", col="red", 
     ylim=yr, ylab="NDVI", xlab="Obs.", lwd=3)
points(dataTemp, pch=16, col="red")

# ndvi_mask <- ndvi[3,]
# ndvi_mask[mask_md] <- NA

for(i in 1:5){
  
  yMask <- dataTemp
  missingData <- mask[i,]
  yMask[missingData] <- NA
  
  estimatePrevious[i,] <- approx(x=1:23, y=yMask, method="constant", f=0, n=23)$y
  
  estimateNext[i,] <- approx(x=1:23, y=yMask, method="constant", f=1, n=23)$y
  
  estimateMedio[i,] <- approx(x=1:23, y=yMask, method="constant", f=0.5, n=23)$y
  
  estimateLine[i,] <- na_interpolation(x=yMask)
  
  estimateSpline[i,] <- na_interpolation(x=yMask, option="spline")
  
}

# MSE_Previous
mean(sapply( 1:5, function(s) MSE( estimatePrevious[s,], dataTemp ) ))

mse_previous <- mean(sapply( 1:5, function(s) MSE( estimatePrevious[s,], dataTemp ) ))

# MSE_Next
mean(sapply( 1:5, function(s) MSE( estimateNext[s,], dataTemp ) ))

mse_next <- mean(sapply( 1:5, function(s) MSE( estimateNext[s,], dataTemp ) ))

# MSE_Medio
mean(sapply( 1:5, function(s) MSE( estimateMedio[s,], dataTemp ) ))

mse_medio <- mean(sapply( 1:5, function(s) MSE( estimateMedio[s,], dataTemp ) ))

# MSE_Lineal
mean(sapply( 1:5, function(s) MSE( estimateLine[s,], dataTemp ) ))

mse_lineal <- mean(sapply( 1:5, function(s) MSE( estimateLine[s,], dataTemp ) ))

# MSE_Spline
mean(sapply( 1:5, function(s) MSE( estimateSpline[s,], dataTemp ) ))

mse_spline <- mean(sapply( 1:5, function(s) MSE( estimateSpline[s,], dataTemp ) ))


# --- output as a table

stats <- data.frame(
  ndvi1 = c( mse_previous, mse_next, mse_medio, mse_lineal, mse_spline  )
)

stats <- data.frame( stats, 
                     ndvi2 = c( mse_previous, mse_next, mse_medio, mse_lineal, mse_spline  )
                     )

row.names(stats) <- c("Previous", "Next", "Mean", "Linear", "Spline")

kable(stats, caption = "MSE") %>%
  kable_styling(bootstrap_options = c("striped"),
                full_width = FALSE, position = "center")

# -----------------------------------------------------------------------------
# EJERCICIO: Ver archivo MSE_applied_linearInterpolation


# -----------------------------------------------------------------------------
# --- Evaluamos el desempeño de imputación vía curva de climatología.
# --- Usamos como "verdad" las observaciones de NDVI de Mohinora presentadas
# --- en mohinora_QA.R

# --- Si el objeto mohinora_ndvi_rTp no está en tu sesión de trabajo, descomenta 
# --- y ejecuta las siguientes líneas


mohinoraDIR <- list.dirs(here( "data", "mohinora_2026" )) #paste0( getwd(), "/TIF" )

ndviFILES <- list.files(path = mohinoraDIR[3],
                        pattern = ".tif",
                        full.names = TRUE)

mohinora_NDVI <- rast(ndviFILES)

mohinora_NDVI_rTp <- spRast_valuesCoords(mohinora_NDVI)

set.seed(102)
SAMPLE <- sample(1:nrow(mohinora_NDVI_rTp$values), 10)
set.seed(NULL)

PIXELS_coords <- mohinora_NDVI_rTp$coords[SAMPLE, 1:2]
PIXELS_values <- mohinora_NDVI_rTp$values[SAMPLE, 1:(20+23*25)]

sapply(1:10, function(s) sum( is.na(PIXELS_values[s,])  ) )

for(i in 1:10){
  TEMP_ts <- ts(PIXELS_values[i,],
                start = c(2001,1),
                end = c(2009,23),
                frequency = 23)
  
  plot(TEMP_ts, main=paste0("Coords:", round(PIXELS_coords[i,1],3), ",",
                            round(PIXELS_coords[i,2],3)),
       ylab="NDVI")
}

# ---

ndvi_mohi_sinNA <- PIXELS_values[c(4,6),]

mask <- matrix(nrow = 5, ncol = 11)
mask[1,] <- c(2:3, 6:7, 10, 13:16, 21:22)
mask[2,] <- c(4:5, 8:12, 15, 17, 19:20)
mask[3,] <- c(4, 7:8, 14:17, 20, 21:23)
mask[4,] <- c(2,4, 7:8, 13:16, 19, 21:22)
mask[5,] <- c(2, 4:5, 7:8, 16:19, 21:22)

estimatePrevious <- matrix(nrow = 5, ncol = 23)
estimateNext <- matrix(nrow = 5, ncol = 23)
estimateMedio <- matrix(nrow = 5, ncol = 23)
estimateLine <- matrix(nrow = 5, ncol = 23)
estimateSpline <- matrix(nrow = 5, ncol = 23)
estimateClimatologyLower <- matrix(nrow = 5, ncol = 23)
estimateClimatologyMedian <- matrix(nrow = 5, ncol = 23)
estimateClimatologyUpper <- matrix(nrow = 5, ncol = 23)

# variando el patrón de datos faltantes
# dataset is fixed

# y <- ndvi_mohi_sinNA[2,c((20+23*5+1):(20+23*6))] * 1e-4
# x <- ndvi_mohi_sinNA[2,-c((20+23*5+1):(20+23*6))] * 1e-4

y <- ndvi_mohi_sinNA[1,-c(1:20)] * 1e-4

z <- y[c((23*5+1):(23*6))] # * 1e-4 # 2005-2006

x <- y[-c((23*5+1):(23*6))] #* 1e-4 # 2001-2025 exceptuando 2005


# plot(y, type="l", col="red",
#      ylim=yr, ylab="NDVI", xlab="Obs.", lwd=3)
# points(y, pch=16, col="red")
# 
# plot(ndvi_mohi_sinNA[1,], type="l", col="red",
#       ylab="NDVI", xlab="Obs.", lwd=3)
# points(ndvi_mohi_sinNA[1,], pch=16, col="red")

# plot(z, type="l", col="red",
#      ylim=yr, ylab="NDVI", xlab="Obs.", lwd=3)
# points(y, pch=16, col="red")
# 
# plot(ndvi_mohi_sinNA[1,], type="l", col="red",
#      ylab="NDVI", xlab="Obs.", lwd=3)
# points(ndvi_mohi_sinNA[1,], pch=16, col="red")

dataTemp <- z

for(i in 1:5){
  yMask <- dataTemp
  missingData <- mask[i,]
  yMask[missingData] <- NA
  
  estimatePrevious[i,] <- approx(x=1:23, y=yMask, method="constant", f=0, n=23)$y
  
  estimateNext[i,] <- approx(x=1:23, y=yMask, method="constant", f=1, n=23)$y
  
  estimateMedio[i,] <- approx(x=1:23, y=yMask, method="constant", f=0.5, n=23)$y
  
  estimateLine[i,] <- na_interpolation(x=yMask)
  
  estimateSpline[i,] <- na_interpolation(x=yMask, option="spline")
  
  estimateClimatologyLower[i,] <- gapfill_climatology(y=yMask,x=x,box="lower")
  
  estimateClimatologyMedian[i,] <- gapfill_climatology(y=yMask,x=x,box="median")
  
  estimateClimatologyUpper[i,] <- gapfill_climatology(y=yMask,x=x,box="upper")
}

# MSE_Previous
( mse_previous <- mean(sapply( 1:5, function(s) MSE( estimatePrevious[s,], dataTemp ) )) )

# MSE_Next
( mse_next <- mean(sapply( 1:5, function(s) MSE( estimateNext[s,], dataTemp ) )) )

# MSE_Medio
( mse_mean <- mean(sapply( 1:5, function(s) MSE( estimateMedio[s,], dataTemp ) )) )

# MSE_Lineal
( mse_lineal <- mean(sapply( 1:5, function(s) MSE( estimateLine[s,], dataTemp ) )) )

# MSE_Spline
( mse_spline <- mean(sapply( 1:5, function(s) MSE( estimateSpline[s,], dataTemp ) )) )

# MSE_Clima_Lower
( mse_clima_lower <- mean(sapply( 1:5, function(s) MSE( estimateClimatologyLower[s,], dataTemp ) )) )

# MSE_Clima_Median
( mse_clima_median <- mean(sapply( 1:5, function(s) MSE( estimateClimatologyMedian[s,], dataTemp ) )) )

# MSE_Clima_Upper
( mse_clima_upper <- mean(sapply( 1:5, function(s) MSE( estimateClimatologyUpper[s,], dataTemp ) )) )

# ---

stats <- data.frame(
  "2001" = c( mse_previous, mse_next, mse_medio, mse_lineal, mse_spline,
             mse_clima_lower, mse_clima_median, mse_clima_upper)
)

stats <- data.frame( stats, 
                     ndvi2 = c( mse_previous, mse_next, mse_medio, mse_lineal, mse_spline  )
)

row.names(stats) <- c("Previous", "Next", "Mean", "Linear", "Spline",
                      "Climatology_lower", "Climatology_median", "Climatology_upper")

kable(stats, caption = "MSE") %>%
  kable_styling(bootstrap_options = c("striped"),
                full_width = FALSE, position = "center")

# -----------------------------------------------------------------------------

