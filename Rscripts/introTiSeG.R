
# Written on Sept 13, 2018
# Se presenta MSE (mean squared error) como medida integral de precisión y exactitud.
# También se estudia las principales aportaciones del paper ligado al desarrollo TiSeG

library(imputeTS)
library(zoo)

MSE <- function(true, estimate){ mean( (true - estimate)^2 ) }

###############
# --- MSE --- #
###############

x <- rnorm(1000, mean = 2, sd = 0.75)

mu1 <- mean(x)
mu2 <- median(x)

mean( (mu1 - 2)^2 )

mean( (mu2 - 2)^2 )

sigma1 <- sum( (x - mu1)^2 ) / (length(x) - 1)

sigma2 <- (length(x) - 1) / length(x) * sigma1

mean( ( 0.75 - sqrt(sigma1) )^2 )

mean( ( 0.75 - sqrt(sigma2) )^2 )

# ---

sim <- 1000
mu1 <- numeric(length(sim))
mu1 <- numeric(length(sim))
sigma1 <- numeric(length(sim))
sigma1 <- numeric(length(sim))

for(i in 1:length(sim)){
  x <- rnorm(1000, mean = 2, sd = 0.75)
  mu1[i] <- mean(x)
  mu2[i] <- median(x)
  sigma1[i] <- sum( (x - mean(x))^2 ) / (length(x) - 1)
  sigma2[i] <- mean((x - mean(x))^2) #(length(x) - 1) / length(x) * sigma1[i]
}

cat("MSE de mu1:", MSE(true = 2, estimate = mu1), "\n")

cat("MSE de mu2:", MSE(true = 2, estimate = mu2), "\n")

cat("MSE de sigma1:", MSE(true = 0.75, estimate = sqrt(sigma1)), "\n")

cat("MSE de sigma2:", MSE(true = 0.75, estimate = sqrt(sigma2)), "\n")


##############################
# --- TiSeG Simulaciones --- #
##############################


prof1 <- c(0.38, 0.335, 0.365, 0.3625, .3725, .4425, .6125, .6275, .78, .8, 
           .7825, .595, .545,
           .2815, .1985, .201, .190, .202, .2125, .28,
           .175, .165, .4)

prof2 <- c(.125,
           .120,
           .11,
           .135,
           .105,
           .075, .215, .2, .205, .2105, 
           .7825, .795, .9125,
           .585, 
           .4, .595, .425, .285, .295, .197,
           .1785, .1625, .175)

prof3 <- c(.225, 0.385, .245, .215, .2125, .275, .195, .3125, .645, .465,
           .7815,
           .6135, .685, .68, .642, .5925, .4565, .585, .4105, .395, .285, .265, .3)

prof4 <- c(.115, .135, .11, .2, .105, .2, .205, .285, .295, .4, .475, 
           .585, .540,
           .565,
           .56,
           .495, .4, .385, .395, .365, .285, 
           .170, 
           .175)

prof5 <- c(.125, .209, .2, .395, .285, .225, .305, .235, .585, .665, 
           .785, 
           .775, .635, .505, .535, .525, .515, .395, 
           .295, 
           .1975, .286, .168, .165)

plot(prof1, type = "l", col = "blue", ylim = c(0,1))
points(1:23, prof1, pch = 6, col = "blue")
par(new=T)
plot(prof4, type = "l", col = "lightblue", ylim = c(0,1))
points(1:23, prof4, pch = 8, col = "lightblue")
par(new=T)
plot(prof5, type = "l", col = "darkmagenta", ylim = c(0,1))
points(1:23, prof5, pch = 4, col = "darkmagenta")

par(new=T)
plot(prof2, type = "l", col = "deeppink", ylim = c(0,1))
points(1:23, prof2, pch = 10, col = "deeppink")
par(new=T)
plot(prof3, type = "l", col = "mediumseagreen", ylim = c(0,1))
points(1:23, prof3, pch = 12, col = "mediumseagreen")

masks <- matrix(NA, nrow = 5, ncol = 13)
temp1 <- c(2:3, 6:7, 10, 13:16, 19, 21:22)
masks[1, 1:length(temp1)] <- temp1
temp2 <- c(4:5, 8:12, 15, 17, 19:20)
masks[2, 1:length(temp2)] <- temp2
temp3 <- c(1:2, 4, 7:8, 13:16, 19, 21:23)
masks[3, 1:length(temp3)] <- temp3
temp4 <- c(2,4, 7:8, 13:16, 19, 21:22)
masks[4, 1:length(temp4)] <- temp4
temp5 <- c(2, 4:5, 7:8, 12:13, 16:19, 21:22)
masks[5, 1:length(temp5)] <- temp5

profiles <- matrix(NA, nrow = 5, ncol = 23)
profiles[1, ] <- prof1
profiles[2, ] <- prof2
profiles[3, ] <- prof3
profiles[4, ] <- prof4
profiles[5, ] <- prof5

mse_RowMask_ColMethod <- matrix(NA, nrow = 5, ncol = 4)

# --- tryout with prof1 vs all masks and using linear interpolation ---

linearInterpolByMask <- matrix(NA, nrow = 5, 23)

dataTest <- profiles[1,]

mask1 <- masks[1,]
mask1 <- mask1[!is.na(masks[1,])]
dataTestMask1 <- dataTest
dataTestMask1[mask1] <- NA

mask2 <- masks[2,]
mask2 <- mask2[!is.na(masks[2,])]
dataTestMask2 <- dataTest
dataTestMask2[mask2] <- NA

mask3 <- masks[3,]
mask3 <- mask3[!is.na(masks[3,])]
dataTestMask3 <- dataTest
dataTestMask3[mask3] <- NA

mask4 <- masks[4,]
mask4 <- mask4[!is.na(masks[4,])]
dataTestMask4 <- dataTest
dataTestMask4[mask4] <- NA

mask5 <- masks[5,]
mask5 <- mask5[!is.na(masks[5,])]
dataTestMask5 <- dataTest
dataTestMask5[mask5] <- NA

plot(dataTestMask1, type = "p", col = "blue", ylim = c(0,1))
plot(dataTestMask2, pch = 2, col = "blue", ylim = c(0,1))
plot(dataTestMask3, pch = 3, col = "blue", ylim = c(0,1))
plot(dataTestMask4, pch = 4, col = "blue", ylim = c(0,1))
plot(dataTestMask5, pch = 5, col = "blue", ylim = c(0,1))

linearInterpolByMask[1,] <- na.interpolation(dataTestMask1)
linearInterpolByMask[2,] <- na.interpolation(dataTestMask2)
linearInterpolByMask[3,] <- na.interpolation(dataTestMask3)
linearInterpolByMask[4,] <- na.interpolation(dataTestMask4)
linearInterpolByMask[5,] <- na.interpolation(dataTestMask5)

plot(dataTestMask1, type = "p", col = "blue", ylim = c(0,1))
lines(1:23, linearInterpolByMask[1,])

plot(dataTestMask2, type = "p", pch = 2, col = "blue", ylim = c(0,1))
lines(1:23, linearInterpolByMask[2,])

plot(dataTestMask3, type = "p", pch = 3, col = "blue", ylim = c(0,1))
lines(1:23, linearInterpolByMask[3,])

plot(dataTestMask4, type = "p", pch = 4, col = "blue", ylim = c(0,1))
lines(1:23, linearInterpolByMask[4,])

plot(dataTestMask5, type = "p", pch = 5, col = "blue", ylim = c(0,1))
lines(1:23, linearInterpolByMask[5,])

MSE(true = profiles[1,], estimate = linearInterpolByMask[1,])
MSE(true = profiles[1,], estimate = linearInterpolByMask[2,])
MSE(true = profiles[1,], estimate = linearInterpolByMask[3,])
MSE(true = profiles[1,], estimate = linearInterpolByMask[4,])
MSE(true = profiles[1,], estimate = linearInterpolByMask[5,])

mse_RowMask_ColMethod[1,1] <- MSE(true = profiles[1,], estimate = linearInterpolByMask[1,])
mse_RowMask_ColMethod[2,1] <- MSE(true = profiles[1,], estimate = linearInterpolByMask[2,])
mse_RowMask_ColMethod[3,1] <- MSE(true = profiles[1,], estimate = linearInterpolByMask[3,])
mse_RowMask_ColMethod[4,1] <- MSE(true = profiles[1,], estimate = linearInterpolByMask[4,])
mse_RowMask_ColMethod[5,1] <- MSE(true = profiles[1,], estimate = linearInterpolByMask[5,])

apply(mse_RowMask_ColMethod, 2, mean)

copy <- mse_RowMask_ColMethod

# ---------------------------------------------------------------------

# ---
interpol <- function(x, typeMethod){
  switch(typeMethod,
         linear = na.interpolation(x),
         mean = na.mean(x),
         random = na.random(x),
         spline = na.spline(x))
}

linearInterpolByMask[1,]
x <- dataTestMask1
interpol(x, "linear")
# ---

# =============================================================================

mse_RowMask_ColMethod <- matrix(NA, nrow = 5, ncol = 4)
method <- "linear"
dataTest <- profiles[1,]

for(i in 1:4){
  interpolByMask <- matrix(NA, nrow = 5, 23)
  
  if(i == 1){
    method <- "linear"
  }

  if(i == 2){
    method <- "mean"
  }
  
  if(i == 3){
    method <- "random"
  }
  
  if(i == 4){
    method <- "spline"
  }
  
  mask1 <- masks[1,]
  mask1 <- mask1[!is.na(masks[1,])]
  dataTestMask1 <- dataTest
  dataTestMask1[mask1] <- NA
  
  mask2 <- masks[2,]
  mask2 <- mask2[!is.na(masks[2,])]
  dataTestMask2 <- dataTest
  dataTestMask2[mask2] <- NA
  
  mask3 <- masks[3,]
  mask3 <- mask3[!is.na(masks[3,])]
  dataTestMask3 <- dataTest
  dataTestMask3[mask3] <- NA
  
  mask4 <- masks[4,]
  mask4 <- mask4[!is.na(masks[4,])]
  dataTestMask4 <- dataTest
  dataTestMask4[mask4] <- NA
  
  mask5 <- masks[5,]
  mask5 <- mask5[!is.na(masks[5,])]
  dataTestMask5 <- dataTest
  dataTestMask5[mask5] <- NA
  
  interpolByMask[1,] <- interpol(dataTestMask1, typeMethod = method)
  interpolByMask[2,] <- interpol(dataTestMask2, typeMethod = method)
  interpolByMask[3,] <- interpol(dataTestMask3, typeMethod = method)
  interpolByMask[4,] <- interpol(dataTestMask4, typeMethod = method)
  interpolByMask[5,] <- interpol(dataTestMask5, typeMethod = method)
  
  mse_RowMask_ColMethod[1,i] <- MSE(true = profiles[1,], estimate = interpolByMask[1,])
  mse_RowMask_ColMethod[2,i] <- MSE(true = profiles[1,], estimate = interpolByMask[2,])
  mse_RowMask_ColMethod[3,i] <- MSE(true = profiles[1,], estimate = interpolByMask[3,])
  mse_RowMask_ColMethod[4,i] <- MSE(true = profiles[1,], estimate = interpolByMask[4,])
  mse_RowMask_ColMethod[5,i] <- MSE(true = profiles[1,], estimate = interpolByMask[5,])
}

apply(copy, 2, mean)

apply(mse_RowMask_ColMethod, 2, mean)

copy2 <- mse_RowMask_ColMethod

# =============================================================================

mse_RowMask_ColMethod_LayerProfile <- array(NA, dim = c(nrow = 5, ncol = 4, nlayers = 5) )

for(i in 1:5){ # Layers, i.e. Profiles
  dataTest <- profiles[i,]
  
  for(j in 1:4){
    interpolByMask <- matrix(NA, nrow = 5, 23)
    
    if(j == 1){
      method <- "linear"
    }
    
    if(j == 2){
      method <- "mean"
    }
    
    if(j == 3){
      method <- "random"
    }
    
    if(j == 4){
      method <- "spline"
    }
    
    mask1 <- masks[1,]
    mask1 <- mask1[!is.na(masks[1,])]
    dataTestMask1 <- dataTest
    dataTestMask1[mask1] <- NA
    
    mask2 <- masks[2,]
    mask2 <- mask2[!is.na(masks[2,])]
    dataTestMask2 <- dataTest
    dataTestMask2[mask2] <- NA
    
    mask3 <- masks[3,]
    mask3 <- mask3[!is.na(masks[3,])]
    dataTestMask3 <- dataTest
    dataTestMask3[mask3] <- NA
    
    mask4 <- masks[4,]
    mask4 <- mask4[!is.na(masks[4,])]
    dataTestMask4 <- dataTest
    dataTestMask4[mask4] <- NA
    
    mask5 <- masks[5,]
    mask5 <- mask5[!is.na(masks[5,])]
    dataTestMask5 <- dataTest
    dataTestMask5[mask5] <- NA
    
    interpolByMask[1,] <- interpol(dataTestMask1, typeMethod = method)
    interpolByMask[2,] <- interpol(dataTestMask2, typeMethod = method)
    interpolByMask[3,] <- interpol(dataTestMask3, typeMethod = method)
    interpolByMask[4,] <- interpol(dataTestMask4, typeMethod = method)
    interpolByMask[5,] <- interpol(dataTestMask5, typeMethod = method)
    
    mse_RowMask_ColMethod_LayerProfile[1,j,i] <- MSE(true = dataTest, estimate = interpolByMask[1,])
    mse_RowMask_ColMethod_LayerProfile[2,j,i] <- MSE(true = dataTest, estimate = interpolByMask[2,])
    mse_RowMask_ColMethod_LayerProfile[3,j,i] <- MSE(true = dataTest, estimate = interpolByMask[3,])
    mse_RowMask_ColMethod_LayerProfile[4,j,i] <- MSE(true = dataTest, estimate = interpolByMask[4,])
    mse_RowMask_ColMethod_LayerProfile[5,j,i] <- MSE(true = dataTest, estimate = interpolByMask[5,])
  }
}

apply(copy2, 2, mean)

apply(mse_RowMask_ColMethod_LayerProfile[,,1], 2, mean)

copy2

mse_RowMask_ColMethod_LayerProfile[,,1]

# ---

toPlot <- c(apply(mse_RowMask_ColMethod_LayerProfile[,,1], 2, mean), NA,
                apply(mse_RowMask_ColMethod_LayerProfile[,,2], 2, mean), NA,
                apply(mse_RowMask_ColMethod_LayerProfile[,,3], 2, mean), NA,
                apply(mse_RowMask_ColMethod_LayerProfile[,,4], 2, mean), NA,
                apply(mse_RowMask_ColMethod_LayerProfile[,,5], 2, mean))

yr <- range(toPlot, na.rm = T)
yr[2] <- yr[2] + .01

plot(toPlot, type = "h", lwd = 5, ylim = yr, ylab = "MSE",
     col = c(rep("blue", 4), NA, rep("deeppink", 4), NA, 
             rep("mediumseagreen", 4), NA, rep("lightblue", 4), NA,
             rep("darkmagenta", 4), NA) )
points(1:length(toPlot), toPlot+.0025, pch = rep(c(2, 4, 6, 8, NA),4) )
legend("topright", pch = c(2, 4, 6, 8, NA), legend = c("Linear", "Mean", "Random", "Spline"))
legend("topleft", col = c("blue", "deeppink", "mediumseagreen", "lightblue", "darkmagenta"), 
       legend = c("Prof1", "Prof2", "Prof3", "Prof4", "Prof5"), lty = rep(1, 5), lwd = 5,
       cex = 0.75)

