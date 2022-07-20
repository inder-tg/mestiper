
getAnomaly <- function(x, type = c("standard", "robust"), alpha = 0.05){
  
  if(length(x) == 0){
    stop("length of x must be greater than 2")
  }
  
  type <- match.arg(type)
  
  if(type == "standard"){
    center <- mean(x, na.rm = T)
    dispersion <- sd(x, na.rm = T)
  } else {
    temp <- s_Qn(x = x, mu.too = T)
    center <- temp[1]
    dispersion <- temp[2]
  }
  
  indices <- 1:length(x)
  anomalies <- indices[ x < center - qnorm(1-alpha/2) * dispersion | x > center + qnorm(1-alpha/2) * dispersion ]
  
  list(xAnomaly = anomalies, yAnomaly = x[anomalies], mu = center, sigma = dispersion, 
       quant = qnorm(1-alpha/2), alpha = alpha)
}

# ---

getXDailyMonthlyMeans <- function(x, freq, start, end){

  if(freq == 23){
    
    N <- length(x)
    
    if( N %% freq != 0 ){
      stop("length of x must be a multiple of freq")
    }
    
    monthMeans <- numeric(12 * (N/23) )
    for(k in 1:(N/23)){
      time <- ((k-1)*23 + 1):(k*23)
      
      pixel <- x[time]
      mu <- numeric(12)
      a <- seq(1,23,2)
      b <- seq(2,23,2)
      for(i in 1:12){
        if(i == 12){
          int <- a[12]
        } else {
          int <- a[i]:b[i]
        }
        mu[i] <- mean(pixel[int], na.rm = T)
      }
      
      monthMeans[((k-1)*12 + 1):(k*12)] <- mu
    }
  }
  
  if(freq == 1){
    
    totalNumberDaysPerYear <- ifelse(start:end %% 4 == 0, 366, 365)
    l <- length(totalNumberDaysPerYear)
    tempNumDays <- numeric(1)
    monthMeans <- numeric(l)
    
    for(k in 1:l){
      
      lengthMonths <- readRDS("lengthMonthsLeapYear")
      
      if( totalNumberDaysPerYear[k] == 365 ){
        lengthMonths[2] <- 28
      }
      
      numDaysPerMonth <- cumsum(lengthMonths) 
      
      monthsLength <- c(0, numDaysPerMonth)
      
      pixel <- x[(tempNumDays+1):totalNumberDaysPerYear[k]]
      tempNumDays <- tempNumDays + length(pixel)
      
      mu <- numeric(12)
      for(m in 1:(length(monthsLength) - 1)){
        mu[m] <- mean(x[(monthsLength[m]+1):monthsLength[m+1]], na.rm = T)
      }
      
      monthMeans[((k-1)*12 + 1):(k*12)] <- mu
    }
  }
  
monthMeans
}

# Given a time series from a period with frequency "freq", this function
# gets the difference between the global mean of the time series 
# (calculated over the whole span of monthly means) and the local mean (calculated each month)

getMonthAnomaly <- function(x, freq = 1, start, end){
  
  if( freq != 1 | freq != 23 ){
    stop("Analisys for that value of freq is under development")
  }
  
  if( freq == 1 ){
    if(missing(start) | missing(end)){
      stop("start and end years must be provided")
    }
  }
  
  mu <- getXDailyMonthlyMeans(x, freq = freq, start = start, end = end)
  
  overallMean <- mean(mu)
  
  monthlyAnomalies <- mu - overallMean
  
list(anomalies = monthlyAnomalies, monthlyMeans = mu)
}

# getMonthAnomaly(x = c(rep(1, 366), rep(1, 365)), start = 2000, end=2001)

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


# returns the x.n time series for a given time sequence (ts) and
# a vector with the amount of frequencies k in the signal (X.k)
get.trajectory <- function(X.k, ts, acq.freq) {
  
  N   <- length(ts)
  i   <- complex(real = 0, imaginary = 1)
  x.n <- rep(0,N)           # create vector to keep the trajectory
  ks  <- 0:(length(X.k) - 1)
  
  for(n in 0:(N-1)) {       # compute each time point x_n based on freqs X.k
    x.n[n+1] <- sum(X.k * exp(i * 2 * pi * ks * n / N)) / N
  }
  
  x.n * acq.freq
}

plot.frequency.spectrum <- function(X.k, xlimits = c(0,length(X.k))) {
  plot.data  <- cbind(0:(length(X.k)-1), Mod(X.k))
  
  # TODO: why this scaling is necessary?
  plot.data[2:length(X.k),2] <- 2 * plot.data[2:length(X.k),2] 
  
  plot(plot.data, t = "h", lwd = 2, main="", 
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

getHarmonic <- function(Xk, i, ts, acq.freq){
  Xk.h <- rep(0, length(Xk))
  Xk.h[i+1] <- Xk[i+1] # i-th harmonic
  get.trajectory(Xk.h, ts, acq.freq = acq.freq)
}


LoadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env) 
}


getYear <- function(start = 2002, end = 2012, breakPoint, totalDays = c(0, 12 * 1:9) ){
  period <- start:end
  year <- period[sum( totalDays - breakPoint < 0 )]
  year  
}

