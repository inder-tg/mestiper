
# --- THIS MUST BE RUN ONLY ONCE ----------------------------------------------

library(raster)
library(fields)
library(mapview)
library(RColorBrewer) # http://colorbrewer2.org/#type=diverging&scheme=RdYlGn&n=5

rotate <- function(x) t(apply(x, 2, rev))

source("myFunctions.R")

dataFrames <- LoadToEnvironment(paste0(getwd(), "/data/marismas/dfBricks.RData"))
df_amp0 <- LoadToEnvironment(paste0(getwd(), "/data/marismas/dfAmps.RData"))$df_amp0

df_original <- dataFrames$df_original

df_interpol <- dataFrames$df_interpol

ndmi <- brick("stak_ndmi.tif")

ndmiHants <- brick("stak_ndmi_interpol.tif")

master <- raster("master_file.tif")

# -----------------------------------------------------------------------------

layer <- 20 # put the number you like between 1 and 612

nf <- layout(matrix( c(rep(0,4),0,1,2,0,0,3,3,0, rep(0,4)), 4, 4, byrow = TRUE),
             c(0.5, 5, 5, 0.5), c(0.5, 4, 4, 0.5), T)
layout.show(nf)

image(ndmi[[layer]], col = rev(terrain.colors(5)))

image(ndmiHants[[layer]], col = rev(terrain.colors(5)))

toPlot <- locator() # click on a pixel in the right-hand side plot then press ESC and continue

points(toPlot, pch = 8, cex = 2.5)

ts <- get_timeSeries_byClicking(toPlot, df = df_original)

tsInterpol <- get_timeSeries_byClicking(toPlot, df = df_interpol)

plot(ts)

lines(tsInterpol, col = "red")

dev.off() # if you want to display another pixel you must run this line and start over from line 20

# -----------------------------------------------------------------------------

r <- LoadToEnvironment(paste0(getwd(), "/RData/statsAmp0_2000_2016.RData"))
g <- LoadToEnvironment(paste0(getwd(), "/RData/statsAmp1_2000_2016.RData"))
b <- LoadToEnvironment(paste0(getwd(), "/RData/statsAmp2_2000_2016.RData"))

estimates_amp0 <- r$matEstimateAmpZero
pvalue_amp0 <- r$matPValueAmpZero

estimates_amp1 <- g$matEstimateAnnualAmp
pvalue_amp1 <- g$matPValueAnnualAmp

estimates_amp2 <- b$matEstimateSemiAnnualAmp 
pvalue_amp2 <- b$matPValueSemiAnnualAmp

sign_amp0 <- estimates_amp0
sign_amp0[pvalue_amp0 >= 0.05] <- -1
r <- sign_amp0 + 1
sign_amp1 <- estimates_amp1
sign_amp1[pvalue_amp1 >= 0.05] <- -1
g <- sign_amp1 + 1
sign_amp2 <- estimates_amp2
sign_amp2[pvalue_amp2 >= 0.05] <- -1
b <- sign_amp2 + 1

red <- matrixToRaster(r, master)
green <- matrixToRaster(g, master)
blue <- matrixToRaster(b, master)

MODISRGB <- stack(red * 1e2,
                  green * 1e2,
                  blue * 1e2)

plotRGB(MODISRGB, r = 1, g = 2, b = 3, axes = TRUE, stretch = "lin",
        main = "MODIS True Color Composite")

layer <- 450 # put the number you like between 1 and 612

# nf <- layout(matrix( c(rep(0,4),0,1,0,0,0,2,2,0, rep(0,4)), 4, 4, byrow = TRUE),
#              c(0.5, 7, 5, 0.5), c(0.5, 6, 4, 0.5), T)
# layout.show(nf)

nf <- layout(matrix( c(rep(0,4),
                       0,1,0,0,
                       0,1,2,0,
                       0,1,0,0,
                       rep(0,4)), 5, 4,
                       byrow = TRUE),
             c(0.5, 7, 5, 0.5), c(0.5, 6, 4, 0.5), T)
layout.show(nf)


plotRGB(MODISRGB, r = 1, g = 2, b = 3, axes = TRUE, stretch = "lin",
        main = "Trend Analysis of Shape Parameters")

toPlot <- locator() # click on a pixel in the right-hand side plot then press ESC and continue

points(toPlot, pch = 8, cex = 2.5)

ts <- get_timeSeries_byClicking(toPlot, df = df_original)

tsInterpol <- get_timeSeries_byClicking(toPlot, df = df_interpol)

plot(ts)

lines(tsInterpol, col = "red")

dev.off() # if you want to display another pixel you must run this line and start over from line 84/88

# --- Layer by layer and ts ---

r <- LoadToEnvironment(paste0(getwd(), "/statsAmp0_2000_2016.RData"))
# g <- LoadToEnvironment(paste0(getwd(), "/RData/statsAmp1_2000_2016.RData"))
# b <- LoadToEnvironment(paste0(getwd(), "/RData/statsAmp2_2000_2016.RData"))

estimates_amp0 <- r$matEstimateAmpZero
pvalue_amp0 <- r$matPValueAmpZero

# estimates_amp1 <- g$matEstimateAnnualAmp
# pvalue_amp1 <- g$matPValueAnnualAmp

# estimates_amp2 <- b$matEstimateSemiAnnualAmp 
# pvalue_amp2 <- b$matPValueSemiAnnualAmp

sign_amp0 <- estimates_amp0
sign_amp0[pvalue_amp0 >= 0.05] <- NA
# sign_amp1 <- estimates_amp1
# sign_amp1[pvalue_amp1 >= 0.05] <- NA
# sign_amp2 <- estimates_amp2
# sign_amp2[pvalue_amp2 >= 0.05] <- NA

master <- raster("master_file.tif")
amp0 <- matrixToRaster(sign_amp0, master) # mean
# amp1 <- matrixToRaster(sign_amp1, master) # annual amp
# amp2 <- matrixToRaster(sign_amp2, master) # semi annual amp

miAzul <- c("#d7191c", "#f8c73c", "#1a9641")
# cts_paleta <- c("#d7191c", "#f8c73c", "#a6d96a")

t <- range(values(amp0), na.rm = T)
breaks <- seq(t[1], t[2], length = 4)
breaksToPlot <- round(breaks, digits = 3)
breaksLegend <- matrix(NA, nrow=3, ncol=2)
breaksLegend[1,] <- c(breaksToPlot[1], breaksToPlot[2])
breaksLegend[2,] <- c(breaksToPlot[2]-0.0001, breaksToPlot[3])
breaksLegend[3,] <- c(breaksToPlot[3]+0.0001, breaksToPlot[4])
redLegend  <- paste0("(", breaksLegend[1,1], ",", breaksLegend[1,2], ")")
orangeLegend  <- paste0("(", breaksLegend[2,1], ", ", breaksLegend[2,2], ")")
greenLegend  <- paste0("( ", breaksLegend[3,1], ", ", breaksLegend[3,2], ")")

mapview(amp0,
        na.color = "lightgray",
        label = T,
        query.type = "mousemove",
        query.digits = 4)

nf <- layout(matrix( c(rep(0,4),
                       0,1,2,0,
                       0,1,3,0,
                       0,1,4,0,
                       rep(0,4)),
                     5, 4, byrow = TRUE),
             c(0.5, 7, 5, 0.5), c(0.5, 3, 3, 3, 0.5), T)
layout.show(nf)

# image(amp0, col = cts_paleta)

image(amp0, col = miAzul, breaks = breaks, main = "Significant trend in mean (at 5%)", 
      cex.main = 2)
legend("bottomleft", fill = miAzul, bty = "n", cex = 1.5,
       title = "Estimated trend",
       legend = c(redLegend, orangeLegend, greenLegend))

toPlot <- locator() # click on three (different) pixels in the right-hand side plot then press ESC and continue

locator(1, type="p", par(pch=16, cex =0.015))

output <- get_timeSeries_byClicking(toPlot, df = df_original)
outputInterpol <- get_timeSeries_byClicking(toPlot, df = df_interpol)

ts <- output$ts
tsInterpol <- outputInterpol$ts

iD <- getId(coord = output$coord, df = df_amp0, breaks = breaksLegend)
points(toPlot, pch = iD, bg = c("black", "purple", "blue"), cex = 1.75)

yRan <- range(ts, na.rm = T)
yRan[1] <- yRan[1] - 0.01
yRan[2] <- yRan[2] + 0.01

for(i in 1:3){
  ts_test <- ts(ts[i,], start = c(2000, 1), end = c(2016, 36), frequency = 36)
  
  tsInterpol_test <- ts(tsInterpol[i,], start = c(2000, 1), end = c(2016, 36), frequency = 36)
  
  par(cex.lab = 1.5, cex.axis = 1.5) 
  plot(ts_test, type = "p", ylab = "NDMI", xlab = "Years", 
       ylim = yRan, pch = 16)
  legend("topleft", legend = "", pch = iD[i], cex = 1.75, bty = "n", horiz = T,
         pt.bg = ifelse(i==1, "black", ifelse(i == 2, "purple", "blue")))
  lines(tsInterpol_test, col = "red")
}

# plot(ts)
# lines(tsInterpol, col = "red")

dev.off()
# -----------------------------------------------------------------------------

