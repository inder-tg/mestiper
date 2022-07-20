
# Written by Inder Tecuapetla on March 3, 2020
# Based on Colditz et al. 2008 
# Based on Primer Taller TATSSI (https://github.com/GerardoLopez/TATSSI) 

# =============================================================================

library(raster) 
library(imputeTS)
library(zoo)
library(reticulate)

# Linux
# use_python("/usr/bin/python3/")
# Windows
use_python(py_config()$pythonhome)

scipyInt <- import("scipy.interpolate")

MSE <- function(estimate, true){
  mean(( estimate - true )^2)
}

# =============================================================================

# define datasets directory

rootDir <- getwd()

dataDir <- paste0( getwd(), "/data/interpol" )

data1 <- readRDS( paste0(dataDir, "/ndvi_trace1_tiseg") )
data2 <- readRDS( paste0(dataDir, "/ndvi_trace2_tiseg") )
data3 <- readRDS( paste0(dataDir, "/ndvi_trace3_tiseg") )
data4 <- readRDS( paste0(dataDir, "/ndvi_trace4_tiseg") )
data5 <- readRDS( paste0(dataDir, "/ndvi_trace5_tiseg") )

yr <- range( data1, data2, data3, data4, data5 )
yr[1] <- yr[1] - 0.1
yr[2] <- yr[2] + 0.1

par(mar = c(4, 3, 1, 2))
plot(data1, type = "l", ylim = yr)
points(1:23, data1, pch = 16)
par(new=T)
plot(data2, type = "l", col = "red", ylim = yr)
points(1:23, data2, pch = 16, col = "red")
par(new=T)
plot(data3, type = "l", col = "blue", ylim = yr)
points(1:23, data3, pch = 16, col = "blue")
par(new=T)
plot(data4, type = "l", col = "green", ylim = yr)
points(1:23, data4, pch = 16, col = "green")
par(new=T)
plot(data5, type = "l", col = "darkorange", ylim = yr)
points(1:23, data5, pch = 16, col = "darkorange")

# =============================================================================

# --- toy example

plot(data2, type = "l", col = "red", ylim = yr)
points(1:23, data2, pch = 16, col = "red")

masc <- c(3, 10, 17) # indices donde nos har'an falta datos

yMasc <- data2
yMasc[masc] <- NA

plot(yMasc, type = "l",  col = "red", ylim = yr)
points(1:23, yMasc, pch = 16, col = "red")

interpol1 <- approx(x = 1:23, y = yMasc, method = "constant", f = 0, n = 23)
interpol2 <- approx(x = 1:23, y = yMasc, method = "constant", f = 1, n = 23)
interpol3 <- approx(x = 1:23, y = yMasc, method = "constant", f = 0.5, n = 23)
interpol4 <- approx(x = 1:23, y = yMasc, n = 23)

x <- 1:23
xClean <- x[!is.na(yMasc)]
yClean <- yMasc[!is.na(yMasc)]

TEMP <- scipyInt$BarycentricInterpolator(xi = xClean, yi = yClean)
interpol5 <- TEMP(x)

interpol6 <- na_interpolation(x = yMasc, option = "spline")

# =============================================================================

plot(yMasc, type = "l",  col = "red", ylim = yr,
     ylab = "NDVI", xlab = "No. obs", main = "Previous")
points(1:23, interpol1$y, pch = 16, col = "blue")

# ---

plot(yMasc, type = "l",  col = "red", ylim = yr,
     ylab = "NDVI", xlab = "No. obs", main = "Next")
points(1:23, interpol2$y, pch = 16, col = "darkorange")

# ---

plot(yMasc, type = "l",  col = "red", ylim = yr,
     ylab = "NDVI", xlab = "No. obs", main = "Punto Medio")
points(1:23, interpol3$y, pch = 16, col = "purple")

# ---

plot(yMasc, type = "l",  col = "red", ylim = yr,
     ylab = "NDVI", xlab = "No. obs", main = "Lineal")
points(1:23, interpol4$y, pch = 16, col = "darkgreen")

# ---

plot(yMasc, type = "l",  col = "red", ylim = yr,
     ylab = "NDVI", xlab = "No. obs", main = "Barycentric")
points(1:23, interpol5, pch = 16, col = "green")

# ---

plot(yMasc, type = "l",  col = "red", ylim = yr,
     ylab = "NDVI", xlab = "No. obs" , main = "Spline")
points(1:23, interpol6, pch = 16, col = "lightseagreen")


# =============================================================================
# --- quantitative analysis ---

masks <- matrix(0, nrow = 5, ncol = 11)
masks[1,] <- c(2:3, 6:7, 10, 13:16, 21:22)
masks[2,] <- c(4:5, 8:12, 15, 17, 19:20)
masks[3,] <- c(4, 7:8, 14:17, 20, 21:23)
masks[4,] <- c(2,4, 7:8, 13:16, 19, 21:22)
masks[5,] <- c(2, 4:5, 7:8, 16:19, 21:22)

estimatePrevious <- matrix(NA, nrow = 5, ncol = 23)
estimateNext <- matrix(NA, nrow = 5, ncol = 23)
estimateMedio <- matrix(NA, nrow = 5, ncol = 23)
estimateLine <- matrix(0, nrow = 5, ncol = 23)
estimateBary <- matrix(0, nrow = 5, ncol = 23)
estimateSpline <- matrix(0, nrow = 5, ncol = 23)

bigData <- matrix(0, nrow = 5, ncol = 23)
bigData[1,] <- data1
bigData[2,] <- data2
bigData[3,] <- data3
bigData[4,] <- data4
bigData[5,] <- data5

# =============================================================================
# varying masks
# dataset is fixed

dataTemp <- bigData[1,]

for(i in 1:5){
  yMask <- dataTemp
  mask <- masks[i,]
  yMask[mask] <- NA
  
  estimatePrevious[i,] <- approx(x = 1:23, y = yMask, method = "constant", 
                                 f = 0, n = 23)$y
  
  estimateNext[i,] <- approx(x = 1:23, y = yMask, method = "constant", 
                             f = 1, n = 23)$y
  
  estimateMedio[i,] <- approx(x = 1:23, y = yMask, method = "constant", 
                              f = 0.5, n = 23)$y
  
  estimateLine[i,] <- na_interpolation(yMask)
  
  x <- 1:23
  xClean <- x[!is.na(yMask)]
  yClean <- yMask[!is.na(yMask)]
  
  TEMP <- scipyInt$BarycentricInterpolator(xi = xClean, yi = yClean)
  estimateBary[i,] <- TEMP(x)
  
  estimateSpline[i,] <- na.spline(yMask)
}

# MSE_Previous
mean(sapply( 1:23, function(s) MSE( estimatePrevious[,s], dataTemp[s] ) ))

# MSE_Next
mean(sapply( 1:23, function(s) MSE( estimateNext[,s], dataTemp[s] ) ))

# MSE_Medio
mean(sapply( 1:23, function(s) MSE( estimateMedio[,s], dataTemp[s] ) ))

# MSE_Lineal
mean(sapply( 1:23, function(s) MSE( estimateLine[,s], dataTemp[s] ) ))

# MSE_Barycentric
mean(sapply( 1:23, function(s) MSE( estimateBary[,s], dataTemp[s] ) ))

# MSE_Spline
mean(sapply( 1:23, function(s) MSE( estimateSpline[,s], dataTemp[s] ) ))



