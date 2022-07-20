# 
# Written on Sep 27, 2018
# This scripts introduces harmonic regression and its applications on fitting
# complete time series, filling in missing values and modeling the so-called
# vegetation phenological curve; an application of harmonic regression on trend
# analysis of two real datasets is also included.
# 
# Based on Eastman et al. (2009)
# https://mabouali.wordpress.com/projects/harmonic-analysis-of-time-series-hants/
# 
# -----------------------------------------------------------------------------
source("myFunctions.R")
# -----------------------------------------------------------------------------
LoadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env) 
}
# -----------------------------------------------------------------------------
expFourier <- function(amplitude, phase, t, f0) { 
  sum ( amplitude[1:length(amplitude)] * cos ( 2 * pi * f0 * t - phase[1:length(phase)] / (180/pi) )  )
}

dirData <- paste0(getwd(), "/data/hants")

library(trend)
library(raster)
#install.package geoTS
# -----------------------------------------------------------------------------

# --- MOTIVATION: harmonic functions ---

# the simplest harmonic

A <- 1
w <- 1
phi <- 0

t <- seq(0, 2*pi, length=100)

y <- A * sin(w * t + phi)

plot(t, y, type = "l", lwd = 3)

phi <- pi/2

y1 <- A * sin(w * t + phi)

plot(t, y1, type = "l", lwd = 3, col = "red")

plot(t, y, type = "l", lwd = 3, lty = 2)
par(new=T)
plot(t, y1, type = "l", lwd = 3, col = "red")
par(new = T)
points(c(rep(NA,25), t[26:100]), c(rep(NA,25), y1[1:75]), 
       type = "l", lwd = 3, col = "red")

A <- 1
w <- 2
phi <- 0

t <- seq(0, 2*pi, length=100)
y <- A * sin(w * t + phi)

phi <- pi/2
y1 <- A * sin(w * t + phi)

plot(t, y, type = "l", lwd = 3, lty = 2)
par(new=T)
plot(t, y1, type = "l", lwd = 3, col = "red")
par(new = T)
# -----------------------------------------------------------------------------

# --- Basics on harmonic analysis ---

acq.freq <- 100                    # data acquisition frequency (Hz)
time     <- 6                      # measuring time interval (seconds)
ts       <- seq(0, time, 1 / acq.freq) # vector of sampling time-points (s)
f.0      <- 1 / time                 # fundamental frequency (Hz)

dc.component       <- 0
component.freqs    <- c(3)      # frequency of signal components (Hz)
component.delay    <- c(0)       # delay of signal components (radians)
component.strength <- c(.5)    # strength of signal components

f <- function(dc.component, component.strength, component.freqs, component.delay, t) { 
  dc.component +
    sum( component.strength * sin( 2 * pi * f.0 * component.freqs * t + component.delay))
}

trajectory <- sapply(ts, function(t) f(dc.component, component.strength, 
                                       component.freqs, component.delay, t = t))

yr <- range(trajectory)
yr[1] <- yr[1] - 0.1
yr[2] <- yr[2] + 0.1
plot(ts, trajectory, type="l", ylim = yr)

X.k <- fft(trajectory)                   # find all harmonics with fft()
plot.frequency.spectrum(X.k, xlimits = c(0,20))

yr <- range(trajectory)
yr[1] <- yr[1] - 0.1
yr[2] <- yr[2] + 0.1

firstHarmonic <- getHarmonic(X.k, 1, ts = ts, acq.freq = acq.freq) / ( acq.freq * get.componentStrength(trajectory, j = 1) )
secondHarmonic <- getHarmonic(X.k, 2, ts = ts, acq.freq = acq.freq) / ( acq.freq * get.componentStrength(trajectory, j = 2) )
thirdHarmonic <- getHarmonic(X.k, 3, ts = ts, acq.freq = acq.freq) / ( acq.freq * get.componentStrength(trajectory, j = 3) )

par(mfrow = c(3, 1))
plot(ts, firstHarmonic, col = "red")
plot(ts, secondHarmonic, col = "blue")
plot(ts, thirdHarmonic, col = "green")
# dev.off()

yr <- range(trajectory)
yr[1] <- yr[1] - 0.1
yr[2] <- yr[2] + 0.1

par(mfrow = c(1, 1))
plot(ts, trajectory, type="l", ylim = yr, lwd = 3)
points(ts, thirdHarmonic, col = "green", pch = 16, ylim = yr)
# -----------------------------------------------------------------------------

# Example1: noisy time series

set.seed(102)
trajectoryNoise <- trajectory +  rnorm(601, sd = 0.05)
set.seed(NULL)

plot(ts, trajectoryNoise, type = "l")

X.k <- fft(trajectoryNoise)                   # find all harmonics with fft()
plot.frequency.spectrum(X.k, xlimits = c(0,20))

bestFit <- getHarmonic(X.k, i = 3, ts = 1:length(trajectoryNoise), acq.freq = acq.freq) / 
  (acq.freq * get.componentStrength(trajectoryNoise, j = 3))

test <- harmonicReg(y = trajectoryNoise, lenBasePeriod = length(trajectoryNoise),
                    numFreq = 3, delta = 0.01)

yr <- range( range(trajectoryNoise), range(Re(bestFit)) )
yr[1] <- yr[1] - 0.1
yr[2] <- yr[2] + 0.1

plot(1:length(trajectoryNoise), trajectoryNoise, pch = 16, ylim = yr, col = "lightgray")
lines(1:length(trajectory), trajectory, lwd = 3)
lines(1:length(trajectoryNoise), bestFit, col = "red", lwd = 3)
lines(1:length(trajectoryNoise), test$fitted, col = "green", lwd = 3)
# -----------------------------------------------------------------------------

# --- APPLICATIONS: harmonicR test ---

data <- read.csv(paste0(dirData, "/HANTS_MATLAB_test")) 

data2 <- data[-c(19,21,23,25),]

plot(data2$x, data2$y, pch = 16)

X.k <- fft(data2$y)                   # find all harmonics with fft()
plot.frequency.spectrum(X.k, xlimits = c(0,20))

test1 <- harmonicReg(y = data2$y, lenBasePeriod = length(data2$y), numFreq = 5, delta = 0.01)

plot(data2$y, pch = 16, col = "lightgray") #type = "l", lwd = 3, ylim = yr)
lines(test1$fitted)

data3 <- data[-c(16:18,20,22,24,26),]

plot(data3$x, data3$y, pch = 16)

X.k <- fft(data3$y)                   # find all harmonics with fft()
plot.frequency.spectrum(X.k, xlimits = c(0,20))

test2 <- harmonicReg(y = data3$y, lenBasePeriod = length(data3$y), numFreq = 5, delta = 0.01)
test3 <- harmonicReg(y = data3$y, lenBasePeriod = length(data3$y), numFreq = 10, delta = 0.01)
test4 <- harmonicReg(y = data3$y, lenBasePeriod = length(data3$y), numFreq = 8, delta = 0.01)

plot(data3$y, pch = 16, col = "lightgray") #type = "l", lwd = 3, ylim = yr)
lines(test2$fitted)
lines(test3$fitted, lty = 4)
lines(test4$fitted, lty = 2)

# -----------------------------------------------------------------------------
# hants test
# -----------------------------------------------------------------------------

plot(data$x, data$y, pch = 16)

hantsLow <- hants(numImages = length(data$y), lenBasePeriod = length(data$y),
                  numFreq = 3, y = data$y, ts = 1:length(data$y),
                  HiLo = "Lo", low = 0, high = 255, fitErrorTol = 5, degreeOverDeter = 1,
                  delta = 0.1)

plot(data$x, data$y, pch = 16)
lines(hantsLow$fitted, lty = 2)

hantsHigh <- hants(numImages = length(data$y), lenBasePeriod = length(data$y),
                  numFreq = 3, y = data$y, ts = 1:length(data$y),
                  HiLo = "Hi", low = 0, high = 255, fitErrorTol = 5, degreeOverDeter = 1,
                  delta = 0.1)

plot(data$x, data$y, pch = 16)
lines(hantsHigh$fitted, lty = 2)

harmR <- harmonicReg(y = data$y, lenBasePeriod = length(data$y), numFreq = 3, delta = 0.1)

plot(data$x, data$y, pch = 16)
lines(harmR$fitted, lty = 2)

# -----------------------------------------------------------------------------
# --- missing values ---
# -----------------------------------------------------------------------------

dataMissing <- read.csv(paste0(dirData, "/HANTS_MATLAB_testMissing"))

plot(dataMissing$x, dataMissing$y, pch = 16)

hantsLow <- hants(numImages = length(dataMissing$y), lenBasePeriod = length(dataMissing$y),
                  numFreq = 3, y = dataMissing$y, ts = 1:length(dataMissing$y),
                  HiLo = "Lo", low = 0, high = 255, fitErrorTol = 5, degreeOverDeter = 1,
                  delta = 0.1)

hantsHigh <- hants(numImages = length(dataMissing$y), lenBasePeriod = length(dataMissing$y),
                   numFreq = 3, y = dataMissing$y, ts = 1:length(dataMissing$y),
                   HiLo = "Hi", low = 0, high = 255, fitErrorTol = 5, degreeOverDeter = 1,
                   delta = 0.1)

harmR <- harmonicReg(y = dataMissing$y, lenBasePeriod = length(dataMissing$y), numFreq = 3, delta = 0.1)

lines(hantsLow$fitted, lty = 4, col = "red")
lines(hantsHigh$fitted, lty = 2, col = "blue")
lines(harmR$fitted, lty = 4, col = "green")

# -----------------------------------------------------------------------------

# Example 2: NDVI time series

mat_NDVI <- LoadToEnvironment( paste0(dirData, "/mat_NDVI_2001.RData") )

ndvi <- mat_NDVI$v / 1e4

pixel <- c(1, 257) # c(2,22)
tsNDVI <- ndvi[pixel[1], pixel[2],]
l <- length(tsNDVI[tsNDVI == 3.2767])

tsNDVI_toPlot <- tsNDVI
tsNDVI_toPlot[ tsNDVI_toPlot == 3.2767 ] <- NA

PP_sd <- tsNDVI_toPlot
rangePP_sd <- range( PP_sd, na.rm = T )

yrange <- range( PP_sd, na.rm = T )
yrange[1] <- yrange[1] - .05
yrange[2] <- yrange[2] + .05

par(mar = c(5,5,2,3))

plot( PP_sd, type = "p", col = "blue", pch = 16,
      lwd = 4, axes = F,
      ylim = yrange,
      cex.axis = 0.25, cex.lab = 1.5,
      xlab = "Days",
      main = paste0("Pixel: (", pixel[1], ",", pixel[2], ")"),
      ylab = "NDVI")
axis(side = 1, at = seq(1, 365, by = 50), lwd = 4, cex.axis = 1.25,
     labels = c( as.character(seq(1, 365, by = 50)) ) )
ran <- range(PP_sd, na.rm = T)
axis(side = 2, at = seq(ran[1], ran[2], length = 10), col = "blue", lwd = 4,
     labels = round(seq(ran[1], ran[2], length = 10), digits = 2)) 

fit <- hants(numImages = length(tsNDVI), lenBasePeriod = length(tsNDVI),
              ts = 1:length(tsNDVI), HiLo = "Lo",
              low = -1, high = 1,
              numFreq = 3, y = tsNDVI, fitErrorTol = 0.05,
              degreeOverDeter = 5, delta = 0.1)

pixelSmooth <- fit$fitted

lai <- pixelSmooth  
rangeLai <- range(lai, na.rm = T)

par(new = T, xpd = T)
plot( lai, type = "l", col = "red",  lwd = 4, yaxt = "n", axes = F,
      ylim = yrange, xlab = "", ylab = "" )
axis(side = 4, at = seq(rangeLai[1], rangeLai[2], length = 10), col = "red", 
     lwd = 4, labels = round(seq( rangeLai[1], rangeLai[2], length = 10), digits = 2))

# *** Ejercicio *** 
# Aplica el procedimiento anterior a los pixeles c(8, 210) y c(11,509)


# -----------------------------------------------------------------------------

# --- Trend analysis of shape parameters (NDVI) ---

dataTemp <- LoadToEnvironment( paste0(dirData, "/siteJalisco_interpol.RData") )
jalisco <- dataTemp$aux

ndvi_ts <- ts(jalisco[250, 250, ], start = c(2000, 1), end = c(2015, 23),
              freq = 23)

plot(1:(16*23), ndvi_ts, type = "l", col = "red")
plot(1:(3*23), ndvi_ts[1:(3*23)], type = "l", col = "red")
lines(1:23, ndvi_ts[1:23], type = "l", col = "blue")
lines(1:23+23, ndvi_ts[1:23+23], type = "l", col = "green")
lines(1:23+2*23, ndvi_ts[1:23+2*23], type = "l", col = "black")

# hants_year1 <- hants(numImages = length(ndvi_ts[1:23]), 
#                         lenBasePeriod = length(ndvi_ts[1:23]), 
#                         numFreq = 2, y = as.numeric(ndvi_ts[1:23]), 
#                         ts = 1:length(ndvi_ts[1:23]), fitErrorTol = 0.05, 
#                         HiLo = "Lo", low = -1, high = 1,
#                         degreeOverDeter = 5, delta = 0.1)

fit_year1 <- harmonicReg(y = as.numeric(ndvi_ts[1:23]), lenBasePeriod = 23, 
                     numFreq = 2, delta = 0.1)

plot(1:23, ndvi_ts[1:23], type = "l", col = "blue")
# lines(1:23, hants_year1$fitted, lty = 3, col = "red")
lines(1:23, fit_year1$fitted, lty = 3, col = "green")

fit_year2 <- harmonicReg(lenBasePeriod = 23, 
                        numFreq = 2, y = as.numeric(ndvi_ts[1:23+23]), 
                        delta = 0.1)

plot(1:23+23, ndvi_ts[1:23+23], type = "l", col = "green")
lines(1:23+23, fit_year2$fitted, lty = 3, col = "red")

fit_year3 <- harmonicReg(lenBasePeriod = 23, 
                        numFreq = 2, y = as.numeric(ndvi_ts[1:23+2*23]), 
                        delta = 0.1)

plot(1:23+2*23, ndvi_ts[1:23+2*23], type = "l", col = "black")
lines(1:23+2*23, fit_year3$fitted, lty = 3, col = "red")

# shape parameters for first three years
fit_year1
fit_year2
fit_year3

mean <- c(fit_year1$amplitude[1], fit_year2$amplitude[1], fit_year3$amplitude[1])
annualAmp <- c(fit_year1$amplitude[2], fit_year2$amplitude[2], fit_year3$amplitude[2])
semiAnnualAmp <- c(fit_year1$amplitude[3], fit_year2$amplitude[3], fit_year3$amplitude[3])

temp <- sens.slope(mean)
temp$estimates
temp$p.value

temp1 <- sens.slope(annualAmp)
temp1$estimates
temp1$p.value

temp2 <- sens.slope(semiAnnualAmp)
temp2$estimates
temp2$p.value

# in case you wanted to apply this idea on your own
# arrayAmpZero <- array(NA, dim = c(nrow=nrow(jalisco), ncol = ncol(jalisco), nlayers = 16))
# arrayAnnualAmp <- array(NA, dim = c(nrow=nrow(jalisco), ncol = ncol(jalisco), nlayers = 16))
# arraySemiAnnualAmp <- array(NA, dim = c(nrow=nrow(jalisco), ncol = ncol(jalisco), nlayers = 16))
# arrayAnnualPhase <- array(NA, dim = c(nrow=nrow(jalisco), ncol = ncol(jalisco), nlayers = 16))
# arraySemiAnnualPhase <- array(NA, dim = c(nrow=nrow(jalisco), ncol = ncol(jalisco), nlayers = 16))

# for(i in 1:nrow(jalisco)){
#   cat("Row:", i, "\n")
#   for(j in 1:ncol(jalisco)){
#     ndvi_ts <- jalisco[i,j,]
#     
#     amplitudeZero <- numeric(16)
#     annualAmplitude <- numeric(16)
#     annualPhase <- numeric(16)
#     semiAnnualAmplitude <- numeric(16)
#     semiAnnualPhase <- numeric(16)
#     for(l in 1:16){
#       # harmonicReg(y = as.numeric(ndvi_ts[1:23]), lenBasePeriod = 23, 
#       #             numFreq = 2, delta = 0.1)
#       
#       year <- harmonicReg(y = as.numeric(ndvi_ts[1:23 + (l-1) * 23]),
#                           lenBasePeriod = 23, numFreq = 2, 
#                           delta = 0.1)
#       
#       amplitudeZero[l] <- year$amplitude[1]
#       annualAmplitude[l] <- year$amplitude[2]
#       annualPhase[l] <- year$phase[2]
#       semiAnnualAmplitude[l] <- year$amplitude[3]
#       semiAnnualPhase[l] <- year$phase[3]
#     }
#     arrayAmpZero[i,j,] <- amplitudeZero
#     arrayAnnualAmp[i,j,] <- annualAmplitude
#     arraySemiAnnualAmp[i,j,] <- semiAnnualAmplitude
#     arrayAnnualPhase[i,j,] <- annualPhase
#     arraySemiAnnualPhase[i,j,] <- semiAnnualPhase
#   }
# }
# must save the objects just obtained

# amplitudeZero[1:5,1:5,1]
# arrayAmpZero[1:5,1:5,1]

amplitudeZero <- LoadToEnvironment( paste0(dirData, "/meanCube_siteJalisco.RData"))$meanCube
annualAmplitude <- LoadToEnvironment(paste0(dirData, "/annualAmpCube_siteJalisco.RData"))$annualAmpCube
semiAnnualAmplitude <- LoadToEnvironment(paste0(dirData, "/semiAnnualAmpCube_siteJalisco.RData"))$semiAnnualAmpCube

# -----------------------------------------------------------------------------
# Theil-Sen test
# -----------------------------------------------------------------------------

sens.slope(amplitudeZero[1,1,])

matEstimateAmpZero <- matrix(0, nrow = nrow(amplitudeZero), ncol = ncol(amplitudeZero))
matPValueAmpZero <- matrix(0, nrow = nrow(amplitudeZero), ncol = ncol(amplitudeZero))

matEstimateAnnualAmp <- matrix(0, nrow = nrow(amplitudeZero), ncol = ncol(amplitudeZero))
matPValueAnnualAmp <- matrix(0, nrow = nrow(amplitudeZero), ncol = ncol(amplitudeZero))

matEstimateSemiAnnualAmp <- matrix(0, nrow = nrow(amplitudeZero), ncol = ncol(amplitudeZero))
matPValueSemiAnnualAmp <- matrix(0, nrow = nrow(amplitudeZero), ncol = ncol(amplitudeZero))

for(i in 1:nrow(amplitudeZero)){
  if(i %% 100 == 0){
    cat("Row:", i, "\n")
  }
  for(j in 1:ncol(amplitudeZero)){
    temp <- sens.slope(amplitudeZero[i,j,])
    matEstimateAmpZero[i,j] <- temp$estimates
    matPValueAmpZero[i,j] <- temp$p.value
    
    temp1 <- sens.slope(annualAmplitude[i,j,])
    matEstimateAnnualAmp[i,j] <- temp1$estimates
    matPValueAnnualAmp[i,j] <- temp1$p.value
    
    temp2 <- sens.slope(semiAnnualAmplitude[i,j,])
    matEstimateSemiAnnualAmp[i,j] <- temp2$estimates
    matPValueSemiAnnualAmp[i,j] <- temp2$p.value
  }
}

stack_Jalisco <- raster(paste(dirData, "/stack_Jalisco.tif", sep = ""))

test1 <- matrixToRaster(matrix = t(matEstimateAmpZero), raster = stack_Jalisco)

plot(test1, main = "Trends of Amplitude 0")

test2 <- matrixToRaster(matrix = t(matPValueAmpZero), raster = stack_Jalisco)

plot(test2, main = "Pvalues of Trends of Amplitude 0")

alpha <- 0.05
red <- matEstimateAmpZero
red[matPValueAmpZero >= alpha] <- -1
red <- red + 1
r <- red
# 

ampZeroSignificant <- matrixToRaster(matrix = t(r), raster = stack_Jalisco)
writeRaster(ampZeroSignificant, "redSignificativo", format = "GTiff")
# 
plot(ampZeroSignificant, main = "Significant Trends of Amplitude 0")

#

test11 <- matrixToRaster(matrix = t(matEstimateAnnualAmp), raster = stack_Jalisco)

plot(test11, main = "Trends of Amplitude 1")

test21 <- matrixToRaster(matrix = t(matPValueAnnualAmp), raster = stack_Jalisco)

plot(test21, main = "Pvalues of Trends of Amplitude 1")

green <- matEstimateAnnualAmp
green[matPValueAnnualAmp >= alpha] <- -1
green <- green + 1
g <- green
# 
annualAmpSignificant <- matrixToRaster(matrix = t(g), raster = stack_Jalisco)
writeRaster(annualAmpSignificant, "greenSignificativo", format = "GTiff")
# 
plot(annualAmpSignificant, main = "Significant Trends of Amplitude 1")

#

test12 <- matrixToRaster(matrix = t(matEstimateSemiAnnualAmp), raster = stack_Jalisco)

plot(test12, main = "Trends of Amplitude 2")

test22 <- matrixToRaster(matrix = t(matPValueSemiAnnualAmp), raster = stack_Jalisco)

plot(test22, main = "Pvalues of Trends of Amplitude 2")

blue <- matEstimateSemiAnnualAmp
blue[matPValueSemiAnnualAmp >= alpha] <- -1
blue <- blue + 1
b <- blue
# 
semiAnnualAmpSignificant <- matrixToRaster(matrix = t(b), raster = stack_Jalisco)
writeRaster(semiAnnualAmpSignificant, "blueSignificativo", format = "GTiff")
# 
plot(semiAnnualAmpSignificant, main = "Significant Trends of Amplitude 2")


landsatRGB <- brick(ampZeroSignificant * 1e2, 
                    annualAmpSignificant * 1e2, 
                    semiAnnualAmpSignificant * 1e2)

plotRGB(landsatRGB, r = 1, g = 2, b = 3, axes = TRUE, stretch = "lin",
        main = "Landsat True Color Composite")
# 
# rgb_jalisco <- stack("stack_jalisco_rgb.tif")

# plotRGB(rgb_jalisco, axes = TRUE, stretch = "lin",
#         main = "Landsat True Color Composite")

# -----------------------------------------------------------------------------








