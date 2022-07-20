
plot.frequency.spectrum <- function(X.k, xlimits = c(0,length(X.k)), ...) {
  plot.data  <- cbind(0:(length(X.k)-1), Mod(X.k))
  
  # TODO: why this scaling is necessary?
  plot.data[2:length(X.k),2] <- 2 * plot.data[2:length(X.k),2] 
  
  plot(plot.data, t = "h", lwd = 2, ...,
       xlab = "Frequency (Hz)", ylab="Strength", 
       xlim = xlimits, ylim = c(0, max(Mod(plot.data[,2]))))
}

get.componentStrength <- function(y, j){
  N <- length(y)
  ks <- 2:(N-1)
  
  a <- sum( y[1] + y[N] + 2 * sum( y[ks] * cos ( 2 * pi * j * (ks-1) / (N-1) )  ) ) / (N-1)
  
  b <- 2 * sum( y[ks] * sin ( 2 * pi * j * (ks-1) / (N-1) )  ) / (N-1)
  
  sqrt(a^2 + b^2)
}

getMagnitudePhase <- function(y, j){
  N <- length(y)
  ks <- 2:(N-1)
  
  a <- sum( y[1] + y[N] + 2 * sum( y[ks] * cos ( 2 * pi * j * (ks-1) / (N-1) )  ) ) / (N-1)
  
  b <- 2 * sum( y[ks] * sin ( 2 * pi * j * (ks-1) / (N-1) )  ) / (N-1)
  
  phase <- atan2(b, a) #ifelse( a > 0,  atan2(b, a), atan2(b, a) + pi )
  
# list( magnitude = sqrt(a^2 + b^2), phase = tan(b/a) )
  c( sqrt(a^2 + b^2), ifelse(phase !=0, (phase + 2 * pi), phase) )
}


getHarmonic <- function(Xk, i, ts, acq.freq){
  Xk.h <- rep(0, length(Xk))
  Xk.h[i+1] <- Xk[i+1] # i-th harmonic
  get.trajectory(Xk.h, ts, acq.freq = acq.freq)
}

#
hantsMATLAB <- function(numImages, lenBasePeriod, numFreqs, y, ts, HiLo = c("Hi", "Lo"), 
                        low, high, fitErrorTol, degreeOverDeter, delta){
  
  # numImages = length(trajectoryNoise)
  # lenBasePeriod = 600
  # numFreqs = 10
  # y = trajectoryNoise
  # ts = ts
  # HiLo = "Lo"
  # low = -.5
  # high = .75
  # fitErrorTol = 0.05
  # degreeOverDeter = 5
  # delta = 0.1
  
  # path <- as.numeric(path)
  # numImages = length(path)
  # lenBasePeriod = length(path)
  # numFreqs = 5
  # y = as.numeric(path)
  # ts = 1:length(path)
  # HiLo = "Lo"
  # low = 1
  # high = 180
  # fitErrorTol = 20
  # degreeOverDeter = 5
  # delta = 0.1
  
  # numImages = 36 
  # lenBasePeriod = 36 
  # numFreqs = 3
  # y = data$y 
  # ts = 1:36 
  # HiLo = "Lo" 
  # low = 0 
  # high = 255 
  # fitErrorTol = 5
  # degreeOverDeter = 1 
  # delta = 0.1
  
  mat <- matrix(0, nrow = min(2 * numFreqs + 1, numImages), ncol = numImages)
  amp <- numeric(numFreqs + 1)
  phi <- numeric(numFreqs + 1)
  ySmooth <- numeric(numImages)
  
  HiLo <- match.arg(HiLo)
  sHiLo <- 0
  
  if(HiLo == "Hi"){
    sHiLo <- -1
  } else {
    sHiLo <- 1
  }
  
  numRows <- min(2 * numFreqs + 1, numImages)
  numOutMax <- numImages - numRows - degreeOverDeter
  degree <- 180 / pi
  mat[1,] <- 1
  
  ang <- 2 * pi * (0:(lenBasePeriod-1)) / lenBasePeriod
  cs <- cos(ang)
  sn <- sin(ang)
  
  i <- 1:numFreqs
  for( j in 1:numImages ){
    index <- 1 + (i * (ts[j]-1)) %% lenBasePeriod
    mat[2*i, j] <- cs[index]
    mat[2*i+1, j] <- sn[index]
  }
  
  p <- rep(1, numImages)
  
  p[y < low | y > high] <- 0
  numOut <- sum( p==0 )
  
  if(numOut > numOutMax){
    cat("numOut:", numOut, "numOutMax:", numOutMax)
    stop("Not enough data points")
  }
  
  ready <- F
  nLoop <- 0
  nLoopMax <- numImages
  
  while( !ready & nLoop < nLoopMax ){
    nLoop <- nLoop + 1
    za <- mat %*% ( p * y )
    
    A <- mat %*% diag(p) %*% t(mat)
    A <- A + diag(rep(1, numRows)) * delta
    A[1,1] <- A[1,1] - delta
    zr <- solve(A) %*% za
    
    ySmooth <- t(mat) %*% zr
    weightedResiduals <- sHiLo * (ySmooth - y)
    error <- p * weightedResiduals
    
    # rankVec <- sort(error)
    rankVec <- order(error)
    
    if( floor(rankVec[numImages]) == 0  ){
      maxError <- 0
    } else {
      maxError <- weightedResiduals[rankVec[numImages]]
    }
    
    ready <- (maxError < fitErrorTol) | (numOut == numOutMax)
    
    if(!ready){
      i <- numImages
      j <- rankVec[i]
      while( (p[j] * weightedResiduals[j] > maxError * 0.5) & (numOut < numOutMax) ){
        p[j] <- 0
        numOut <- numOut + 1
        i <- i - 1
        j <- rankVec[i]
      }
    }
  }
  
  amp[1] <- zr[1]
  phi[1] <- 0
  
  zr[numImages+1] <- 0
  
  i <- seq(2, numRows, by = 2)
  iFreq <- (i+2)/2
  ra <- zr[i]
  rb <- zr[i+1]
  amp[iFreq] <- sqrt( ra * ra + rb * rb )
  phase <- atan2(rb, ra) * degree
  phase[phase < 0] <- phase[phase <0] + 360
  phi[iFreq] <- phase
  
  list(amplitude = amp, phase = phi, fitted = ySmooth)  
}
# ---

hantsR <- function(numImages, lenBasePeriod, numFreqs, y, ts, HiLo = c("Hi", "Lo"), 
                        low, high, fitErrorTol, degreeOverDeter, delta){
  
  mat <- matrix(0, nrow = min(2 * numFreqs + 1, numImages), ncol = numImages)
  amp <- numeric(numFreqs + 1)
  phi <- numeric(numFreqs + 1)
  ySmooth <- numeric(numImages)
  
  HiLo <- match.arg(HiLo)
  sHiLo <- 0
  
  if(HiLo == "Hi"){
    sHiLo <- -1
  } else {
    sHiLo <- 1
  }
  
  numRows <- min(2 * numFreqs + 1, numImages)
  numOutMax <- numImages - numRows - degreeOverDeter
  degree <- 180 / pi

  mat <- matrix(0, nrow = numImages, ncol = min(2 * numFreqs + 1, numImages))
  vecTs <- 2 * pi * (0:(lenBasePeriod-1)) / lenBasePeriod
  mat[,1] <- 1
  mat[, seq(2, 2 * numFreqs, 2)] <- sapply(1:numFreqs, function(s) cos( s * vecTs ))
  mat[, seq(2, 2 * numFreqs, 2)+1] <- sapply(1:numFreqs, function(s) sin( s * vecTs ))
  
  p <- rep(1, numImages)
  
  p[y < low | y > high] <- 0
  numOut <- sum( p==0 )
  
  if(numOut > numOutMax){
    cat("numOut:", numOut, "numOutMax:", numOutMax)
    stop("Not enough data points")
  }
  
  ready <- F
  nLoop <- 0
  nLoopMax <- numImages
  
  while( !ready & nLoop < nLoopMax ){
    nLoop <- nLoop + 1
    za <- t(mat) %*% ( p * y )
    
    A <- t(mat) %*% diag(p) %*% mat
    A <- A + diag(rep(1, numRows)) * delta
    A[1,1] <- A[1,1] - delta
    zr <- solve(A) %*% za
    
    ySmooth <- (mat) %*% zr
    weightedResiduals <- sHiLo * (ySmooth - y)
    error <- p * weightedResiduals
    
    rankVec <- order(error)
    
    if( floor(rankVec[numImages]) == 0  ){
      maxError <- 0
    } else {
      maxError <- weightedResiduals[rankVec[numImages]]
    }
    
    ready <- (maxError < fitErrorTol) | (numOut == numOutMax)
    
    if(!ready){
      i <- numImages
      j <- rankVec[i]
      while( (p[j] * weightedResiduals[j] > maxError * 0.5) & (numOut < numOutMax) ){
        p[j] <- 0
        numOut <- numOut + 1
        i <- i - 1
        j <- rankVec[i]
      }
    }
  }
  
  amp[1] <- zr[1]
  phi[1] <- 0
  
  zr[numImages+1] <- 0
  
  i <- seq(2, numRows, by = 2)
  iFreq <- (i+2)/2
  ra <- zr[i]
  rb <- zr[i+1]
  amp[iFreq] <- sqrt( ra * ra + rb * rb )
  phase <- atan2(rb, ra) * degree
  phase[phase < 0] <- phase[phase <0] + 360
  phi[iFreq] <- phase
  
list(amplitude = amp, phase = phi, fitted = as.numeric(ySmooth), iterations = nLoop)  
}
# ---

#------------------------------------------------------------------------------
#
get.outliers <- function(x, alpha = 0.05) {
  
  m <- median(x)
  iqr <- diff(quantile(x, c(0.25, 0.75)))
  n <- length(x)
  crit <- iqr/2/qnorm(0.75) * qnorm(log(1 - alpha/2)/n, log.p = T)
  indices <- 1:length(x)
  
  # indices[abs(x - m) > crit]
  x[abs(x - m) > crit] <- NA
  x
}
#
#------------------------------------------------------------------------------
#
getAnomaly <- function(x, y, type = c("standard", "robust"), quantile = qnorm(0.05) ){
  
  if(length(y) == 0){
    stop("length of x must be greater than 2")
  }

  type <- match.arg(type)
  
  if(type == "standard"){
    centerY <- mean(y, na.rm = T)
    dispersionY <- sd(y, na.rm = T)
  } else {
    temp <- s_Qn(x = y, mu.too = T)
    centerY<- temp[1]
    dispersionY <- temp[2]
  }
  
  indices <- 1:length(y)
  anomalies <- indices[ y > centerY - quantile * dispersionY | y < centerY + quantile * dispersionY ]
  
  # center_y <- center_temp
  if(length(anomalies) != 0){
    aux <- ifelse(type == "standard", mean(y[anomalies], na.rm = T), median(y[anomalies], na.rm = T))
    centerY_Anomaly <- sign( centerY - aux ) * abs( centerY - aux )
    # center_y_Anomaly <- ifelse(type == "standard", temp, 
  } else {
    centerY_Anomaly <- centerY
  }
  
# list(x = x, y = y, x_Anomaly = anomalies, y_Anomaly = y[anomalies], mu = center_y, mu_anomaly = center_y_Anomaly)
  list(x_Anomaly = anomalies, y_Anomaly = y[anomalies], mu = centerY, mu_anomaly = centerY_Anomaly)
}
#
# z <- rnorm(150, sd = 0.5) + rnorm(150, sd = 1.5)
# plot(z)
# 
# temp <- getAnomaly(x = 1:150, y = z, type = "standard")
# 
# plot(temp$x, temp$y)
# points(temp$x_Anomaly, temp$y_Anomaly, col = "red")
#
#------------------------------------------------------------------------------