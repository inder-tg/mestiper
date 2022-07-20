
# Written by Inder Tecuapetla on Nov. 7, 2018
# Introduction to BFAST, in particular we will see some analysis made on a data 
# cube located in Coahuila, Mexico and use bfast01

# =============================================================================
# Landsat imagery time series
# =============================================================================

rm(list = ls())
library(raster)
library(fields)
library(imputeTS)
library(forecast)
library(Kendall)
library(trend)
library(bfast)
library(mapview)

source("myFunctions.R")

# getBreak <- function(start = 2002, end = 2010, frequency = 12, data){
#   output <- NA
#   breakType <- NA
#   significance <- NA
#   
#   dataTS <- ts(data, start = c(start, 1), end = c(end, frequency), frequency = frequency)
#   
#   getBFAST <- bfast01(data = dataTS)
#   
#   if(getBFAST$breaks == 1){
#     output <- getBFAST$breakpoints
#     temp <- bfast01classify(getBFAST)
#     breakType <- as.numeric(temp[1])
#     significance <- as.numeric(temp[2])
#   }
#   
#   list(bPs = output, type = breakType, significance = significance)
# }

# -----------------------------------------------------------------------------

rootDir <- getwd()

dirData <- paste0(getwd(), "/data/BFAST")

dataCoahuila <- LoadToEnvironment( paste0(dirData, '/ori_coa_2013_a.RData') ) # NDVI_stack_interpol

example <- dataCoahuila$ndvi

str(example)

plot(example[45, 50, ], type = "l")

path <- ts(example[45, 50, ], start = c(2005, 1), end = c(2013, 12), freq = 12)

plot(path, type = "l")

stack_Coahuila <- stack( paste0(dirData, "/ori_coa_2013_a.tif" ) )#stack(listFiles[1])

# exampleCopy <- example

temp <- LoadToEnvironment( paste0(dirData, "/stackInterpol_Coahuila.RData") )
dataCoahuilaInterpol <- temp$index

pathInterpol <- ts(dataCoahuilaInterpol[45, 50, ], start = c(2005, 1), end = c(2013, 12), freq = 12)
plot(pathInterpol, type = "l" )

# --- trend analysis ---

# temp <- mk.test(dataCoahuilaInterpol[45, 50, ])

mk.test(dataCoahuilaInterpol[45, 50, ])

plot(dataCoahuilaInterpol[45, 50, ], type = "l")
lines(lowess(time(dataCoahuilaInterpol[45, 50, ]), dataCoahuilaInterpol[45, 50, ]), 
      lwd = 5, col = "blue")

# ---- application of this idea to entire region ---

mapTrend <- readRDS(file = paste0(dirData, "/mapTrend"))

tempMapTrend <- matrix(0, nrow = nrow(mapTrend), ncol = ncol(mapTrend))

tempMapTrend[mapTrend <= 0.05] <- 1

rasterMapTrend <- matrixToRaster(matrix = tempMapTrend, raster = stack_Coahuila[[10]])

plot(rasterMapTrend)

mapview(rasterMapTrend)

# -----------------------------------------------------------------------------

# --- abrupt changes analysis ---

path <- ts(dataCoahuilaInterpol[45, 50, ], start = c(2005, 1), end = c(2013, 12), freq = 12)

output <- bfast01(data = path)

plot(output)

output$breaks

output$breakpoints

bfast01classify(output)

# ---

Row <- 349 
Col <- 253 

path2 <- ts(dataCoahuilaInterpol[Row, Col, ], start = c(2005, 1), end = c(2013, 12), freq = 12)

plot(path2, type = "l", col = "blue")

output2 <- bfast01(path2)

plot(output2, main = "Interruption (decrease with positive break)")

# --- 

Row <- 18 
Col <- 140

path3 <- ts(dataCoahuilaInterpol[Row, Col, ], start = c(2005, 1), end = c(2013, 12), freq = 12)

plot(path3, type = "l", col = "blue")

output3 <- bfast01(path3)

plot(output3, main = "Increase to decrease")

# ---

Row <- 337 
Col <- 159 

path4 <- ts(dataCoahuilaInterpol[Row, Col, ], start = c(2005, 1), end = c(2013, 12), freq = 12)

plot(path4, type = "l", col = "blue")

output4 <- bfast01(path4)

plot(output4, main = "Monotonic decrease")

# ---

# mapAbruptChange <- matrix(0, nrow = nrow(example), ncol = ncol(example))
# mapTimeAbruptChange <- matrix(0, nrow = nrow(example), ncol = ncol(example))
# mapTypeAbruptChange <- matrix(0, nrow = nrow(example), ncol = ncol(example))
# mapSignificanceAbruptChange <- matrix(0, nrow = nrow(example), ncol = ncol(example))
# yearChanges <- c() # numeric(nrow(example) * ncol(example))

# system.time({
#   for(i in 1:nrow(mapAbruptChange)){
#     getBreaks <- matrix(NA, nrow = 3, ncol = ncol(mapAbruptChange)) 
#     getYears <- numeric(ncol(mapAbruptChange))
#     
#     getBreaks <- sapply(1:ncol(mapAbruptChange), function(s) getBreak(data = exampleInterpol[i,s,]))
#     
#     getYears[1:ncol(mapAbruptChange)] <- sapply(1:ncol(mapAbruptChange), 
#                                                 function(s) getYear(breakPoint = as.numeric(getBreaks[1,s])) )
#     
#     mapTimeAbruptChange[i, (1:ncol(mapAbruptChange))[!is.na(getBreaks[1,])]] <- as.numeric(getBreaks[1,])[!is.na(getBreaks[1,])]
#     
#     mapTypeAbruptChange[i, (1:ncol(mapAbruptChange))[!is.na(getBreaks[1,])]] <- as.numeric(getBreaks[2,])[!is.na(getBreaks[1,])]
#     
#     mapSignificanceAbruptChange[i, (1:ncol(mapAbruptChange))[!is.na(getBreaks[1,])]] <- as.numeric(getBreaks[3,])[!is.na(getBreaks[1,])]
#     
#   }
# })

# saveRDS(mapTimeAbruptChange, file = "~/Desktop/posgradosUNAM/R/timeSeries/mapTime", ascii = T)
# saveRDS(mapTypeAbruptChange, file = "~/Desktop/posgradosUNAM/R/timeSeries/mapType", ascii = T)
# saveRDS(mapSignificanceAbruptChange, file = "~/Desktop/posgradosUNAM/R/timeSeries/mapSignificance", ascii = T)

mapTimeAbruptChange <- readRDS(file = paste0(dirData, "/mapTime"))
mapTypeAbruptChange <- readRDS(file = paste0(dirData, "/mapType"))
mapSignificanceAbruptChange <- readRDS(file = paste0(dirData, "/mapSignificance"))

mapTimeAbruptChangeCopy <- mapTimeAbruptChange
mapTimeAbruptChangeCopy[mapTimeAbruptChangeCopy==0] <- NA

NCOL <- ncol(mapTimeAbruptChangeCopy)

for(i in 1:nrow(mapTimeAbruptChangeCopy)){
  mapTimeAbruptChangeCopy[i, 1:NCOL] <- unlist(sapply(1:NCOL, function(s) getYear(breakPoint = mapTimeAbruptChangeCopy[i,s])))
}

rasterMapTimeAbrupt <- matrixToRaster(matrix = mapTimeAbruptChangeCopy, 
                                      raster = stack_Coahuila[[10]])

plot(rasterMapTimeAbrupt, main = "Changes by years")

# ---

mapClassification <- matrix(NA, nrow=nrow(mapTimeAbruptChange), 
                            ncol = ncol(mapTimeAbruptChange))

# hist(mapTypeAbruptChangeCopy)
# 
# hist(mapSignificanceAbruptChange)

mapClassification[ mapTypeAbruptChange == 6 & mapSignificanceAbruptChange == 0] <- 1
mapClassification[ mapTypeAbruptChange == 6 & mapSignificanceAbruptChange == 1] <- 2
mapClassification[ mapTypeAbruptChange == 6 & mapSignificanceAbruptChange == 2] <- 3
mapClassification[ mapTypeAbruptChange == 6 & mapSignificanceAbruptChange == 3] <- 4

mapClassification[ mapTypeAbruptChange == 7 & mapSignificanceAbruptChange == 0] <- 5
mapClassification[ mapTypeAbruptChange == 7 & mapSignificanceAbruptChange == 1] <- 6
mapClassification[ mapTypeAbruptChange == 7 & mapSignificanceAbruptChange == 2] <- 7
mapClassification[ mapTypeAbruptChange == 7 & mapSignificanceAbruptChange == 3] <- 8

mapClassification[ mapTypeAbruptChange == 4 & mapSignificanceAbruptChange == 0] <- 9
mapClassification[ mapTypeAbruptChange == 4 & mapSignificanceAbruptChange == 1] <- 10
mapClassification[ mapTypeAbruptChange == 4 & mapSignificanceAbruptChange == 2] <- 11
mapClassification[ mapTypeAbruptChange == 4 & mapSignificanceAbruptChange == 3] <- 12

rasterClassification <- matrixToRaster(matrix = mapClassification, 
                                       raster = stack_Coahuila[[10]])

plot(rasterClassification, main = "Classification map")

# ---

mapCombined_interruptionDecrease <- matrix(NA, nrow = nrow(mapClassification), 
                                           ncol = ncol(mapClassification))

mapCombined_interruptionDecrease[mapClassification == 1 & mapTimeAbruptChangeCopy == 2005 ] <- 1
mapCombined_interruptionDecrease[mapClassification == 1 & mapTimeAbruptChangeCopy == 2006 ] <- 2
mapCombined_interruptionDecrease[mapClassification == 1 & mapTimeAbruptChangeCopy == 2007 ] <- 3

rasterCombined_interruptionDecrease <- matrixToRaster(matrix = mapCombined_interruptionDecrease, 
                                                      raster = stack_Coahuila[[10]])

plot(rasterCombined_interruptionDecrease, 
     main = "Interruption (decrease with positive break) by years")

# --- 

mapCombined_increaseDecrease <- matrix(NA, nrow = nrow(mapClassification), 
                                       ncol = ncol(mapClassification))

mapCombined_increaseDecrease[mapClassification == 5 & mapTimeAbruptChangeCopy == 2005 ] <- 1
mapCombined_increaseDecrease[mapClassification == 5 & mapTimeAbruptChangeCopy == 2006 ] <- 2
mapCombined_increaseDecrease[mapClassification == 5 & mapTimeAbruptChangeCopy == 2007 ] <- 3

rasterCombined_increaseDecrease <- matrixToRaster(matrix = mapCombined_increaseDecrease, 
                                                  raster = stack_Coahuila[[10]])

plot(rasterCombined_increaseDecrease, main = "Increase to decrease (by years)")

mapview(rasterCombined_increaseDecrease)

# ---

mapCombined_decrease <- matrix(NA, nrow = nrow(mapClassification), 
                               ncol = ncol(mapClassification))

mapCombined_decrease[mapClassification == 9 & mapTimeAbruptChangeCopy == 2005 ] <- 1
mapCombined_decrease[mapClassification == 9 & mapTimeAbruptChangeCopy == 2006 ] <- 2
mapCombined_decrease[mapClassification == 9 & mapTimeAbruptChangeCopy == 2007 ] <- 3

rasterCombined_decrease <- matrixToRaster(matrix = mapCombined_decrease, 
                                          raster = stack_Coahuila[[10]])

plot(rasterCombined_decrease, main = "Monotonic decrease (by years)")
