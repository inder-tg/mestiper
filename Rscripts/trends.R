
library(fpp)
library(forecast)
library(Kendall)
library(trend)

senSlope <- function(x){
  n <- length(x)
  
  temp <- unlist(sapply(1:(n-2), function(s) (x[(s+1):n] - x[s]) / ((s+1):n - s) ))
  
  median(temp)
}

LoadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env) 
}

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


# -----------------------------------------------------------------------------

data("ausbeer")

?ausbeer

ausbeer

tsBeer <- tail(head(ausbeer, 17*4+2), 17*4-4)

# --- getting trend using moving average ---

getMA <- function(y, h, bandwidth){
  sapply( (1+h):(length(y)-h), function(s) sum(y[(s-h):(s+h)])/bandwidth )
}

ventana <- 4

orden <- (ventana-1)/2

# h <- 4
# k <- (h-1)/2
# k1 <- ceiling(k)
# k2 <- floor(k)
# tsBeer[(3-k1):(3+k1)]
# sum(tsBeer[(3-k1):(3+k1)])/h
# tsBeer[(3-k2):(3+k2)]
# sum(tsBeer[(3-k2):(3+k2)])/h
# mean(c(324.5, 206.25))

temp1 <- getMA(tsBeer, h = ceiling(orden), bandwidth = ventana)

temp2 <- getMA(tsBeer, h = floor(orden), bandwidth = ventana)

trendTempBeer <- (temp1+temp2[2:(length(temp2)-1)])/2

yr <- range(as.numeric(tsBeer))
yr[1] <- yr[1] - 25
yr[2] <- yr[2] + 25

# pdf("trendMA.pdf")
par(mar = c(3,3,2,2))
plot(as.ts(tsBeer), type = "l", pch = 16, col = "gray", ylim = yr)
# plot(3:62, trendTempBeer, type = "l") 
par(new=T)
plot(3:62, trendTempBeer, type = "l", col = "green", lwd = 5, yaxt ="n", 
     xaxt="n", ylim = yr) #, xlim=c(1,10)
par(new=T)
plot(3:62, as.ts(temp1), type = "l", col = "red", yaxt ="n", xaxt = "n",
     ylim = yr) #, xlim=c(1,10)
par(new=T)
plot(3:62, as.ts(temp2[2:(length(temp2)-1)]), type = "l", col = "blue", yaxt ="n", xaxt = "n",
     ylim = yr) #, xlim=c(1,10)
legend("topleft", lwd = c(1,1,5), col = c("red", "blue", "green"), legend = c("h=2", "h=1", "mean"))
# dev.off()
# --- trend via MA END ---

# --- trend via LOWESS ---

# pdf("trendLOWESS.pdf")
par(mar = c(3,3,2,2))
plot(as.ts(tsBeer), type = "l", pch = 16, col = "gray", ylim = yr) #
lines(lowess(time(tsBeer), tsBeer), lwd=5, col="orange")
par(new=T)
plot(trendTempBeer, type = "l", col = "green", lwd = 5, yaxt ="n",
     xaxt = "n", ylim = yr) #
legend("topleft", lwd = c(5, 5), col = c("orange", "green"), legend = c("lowess", "MA"))
# dev.off()

# pdf("trendLOWESS_MA.pdf")
par(mar = c(3,3,2,2))
plot(as.ts(tsBeer), type = "l", pch = 16, col = "gray", ylim = yr) #
lines(lowess(time(tsBeer), tsBeer), lwd=5, col="orange")
par(new=T)
plot(trendTempBeer, type = "l", col = "green", lwd = 5, yaxt ="n",
     xaxt = "n", ylim = yr) #
legend("topleft", lwd = c(5, 5), col = c("orange", "green"), legend = c("lowess", "MA"))
# dev.off()

# --- using forecast step-by-step ---

trendBeer <- ma(tsBeer, order = 4, centre = T)

pdf("stl_trend.pdf")
par(mar = c(3,3,2,2))
plot(as.ts(tsBeer), type = "l", pch = 16, col = "gray")
lines(trendBeer, type = "l", col = "green", lwd = 5)
dev.off()
# plot(as.ts(trendBeer))

detrendBeer <- tsBeer - trendBeer
# pdf("stl_detrended.pdf")
par(mar = c(3,3,2,2))
plot(as.ts(detrendBeer), type = "l", pch = 16, col = "darkgray")
# dev.off()

t(matrix(data = trendBeer, nrow = 4))

t(matrix(data = detrendBeer, nrow = 4))

m_beer <- t(matrix(data = detrendBeer, nrow = 4))
seasonal_beer <- colMeans(m_beer, na.rm = T)

# pdf("stl_seasonality.pdf")
par(mar = c(3,3,2,2))
plot(as.ts(rep(seasonal_beer, 16)))
# dev.off()

plot(as.ts((seasonal_beer)), type = "p")

noiseBeer <- tsBeer - trendBeer - rep(seasonal_beer, 16)

# pdf("stl_noise.pdf")
par(mar = c(3,3,2,2))
plot(as.ts((noiseBeer)))
# dev.off()

putBackTogether <- trendBeer + rep(seasonal_beer, 16) + noiseBeer

# pdf("stl_trend_season_noise.pdf")
par(mar = c(3,3,2,2))
plot(as.ts(putBackTogether))
par(new=T)
plot(tsBeer, col = "red")
# dev.off()

# --- using decompose() --- 

ts_beer <- ts(tsBeer, frequency = 4)
decompose_beer <- decompose(ts_beer, "additive")

str(decompose_beer)

plot(as.ts(decompose_beer$seasonal))
plot(as.ts(decompose_beer$trend))
plot(as.ts(decompose_beer$random))
plot(decompose_beer)

# --- using STL() --- 

ts_beer <- ts(tsBeer, frequency = 4)
stl_beer <- stl(ts_beer, "periodic")

str(stl_beer)

# stl_beer$time.series{}

seasonal_stl_beer <- stl_beer$time.series[,1]
trend_stl_beer <- stl_beer$time.series[,2]
random_stl_beer <- stl_beer$time.series[,3]

plot(ts_beer)
plot(as.ts(seasonal_stl_beer))
plot(trend_stl_beer)
plot(random_stl_beer)
plot(stl_beer)

# --- Mann-Kendall test --- 

x <- 1:10
y <- rnorm(10)

plot(y)

out <- Kendall(x,y)
summary(out)

mk.test(y)

data(maxau)

maxau

Q <- maxau[,"Q"]

plot(Q)

mk.test(Q)

# --- Nile example ---

data(Nile)

plot(Nile)

mk.test(Nile, continuity = TRUE)
##
# n <- length(Nile)
# cor.test(x = 1:n, y = Nile, meth = "kendall", continuity = T)

# library(robustbase)
# 
# zdsNile <- (Nile - median(Nile, na.rm = T)) / Qn(Nile)
# 
# plot(zdsNile)
# 
# mk.test(zdsNile, continuity = T)

# ---- Nile example END

#
# Annual precipitation entire Great Lakes
# The time series plot with lowess smooth suggests an upward trend
# The autocorrelation in this data does not appear significant.
# The Mann-Kendall trend test confirms the upward trend.
#
data(PrecipGL)
plot(PrecipGL)
lines(lowess(time(PrecipGL),PrecipGL), lwd = 3, col = 2)
acf(PrecipGL)

MannKendall(PrecipGL)

sens.slope(PrecipGL)

senSlope(PrecipGL)


# --- Sen test ---

dataX <- maxau[, "s"]

# pdf("maxau.pdf")
par(mar = c(3,3,2,2))
plot(dataX)
# dev.off()

mk.test(maxau[,"s"]) 

# pdf("maxau_loess.pdf")
par(mar = c(3,3,2,2))
plot(dataX)
lines(lowess(time(dataX), dataX), lwd = 5, col = "blue")
# dev.off()

n <- length(dataX)
temp <- unlist(sapply(1:(n-2), function(s) (dataX[(s+1):n] - dataX[s]) / ((s+1):n - s) ))

median(temp)

sens.slope(maxau[,"s"])

senSlope(maxau[,"s"])

# ---

rootDir <- getwd()

dirData <- paste0(rootDir, "/data")

setwd(dirData)

temp_NDVI <- LoadToEnvironment( "ori_coa_2013_a.RData" )

datos_Coa <- temp_NDVI$ndvi

coaRaster <- raster(listFiles[4])

tempR <- matrixToRaster(datos_Coa[,,12], coaRaster)

plot(tempR)

library(raster)

listFiles <- list.files(getwd(), pattern = ".tif")

array_interpol <- array( NA, dim = c(nrow = 390, ncol = 395, nlayers = 108) )

matrix_trend <- matrix(NA, nrow = 390, ncol = 395)

matrix_slope_pValue <- matrix(NA, nrow = 390, ncol = 395)

tsTest <- datos_Coa[50, 50,]

tsTest[108] <- NA

tsTest_interpol <- na.interpolation(tsTest)

plot(tsTest_interpol)

tsTest_interpol[108]

array_interpol <- datos_Coa

array_interpol[1,1,108]

array_interpol[array_interpol == 0] <- NA

for(i in 1:nrow(array_interpol)){
  array_interpol[i,1:ncol(array_interpol),] <- t(sapply(1:ncol(array_interpol), function(s) na.interpolation(array_interpol[i,s,])))
}

# save(array_interpol, file = "ori_coa_interpol.RData")

test <- mk.test(array_interpol[50,50,])

test$p.value

for(i in 1:nrow(array_interpol)){
  matrix_trend[i, 1:ncol(array_interpol)] <- sapply(1:ncol(array_interpol), function(s) mk.test(array_interpol[i,s,])$p.value )
}

# save(matrix_trend, file = "ori_coa_trend.RData")

trendRaster <- matrixToRaster(matrix_trend, coaRaster)

plot(trendRaster)

copyTrendRaster <- matrix_trend

copyTrendRaster[copyTrendRaster >= 0.05] <- NA

trendRaster0p5 <- matrixToRaster(copyTrendRaster, coaRaster)

plot(trendRaster0p5)

