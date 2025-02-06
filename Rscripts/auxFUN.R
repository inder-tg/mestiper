
# --- funciones usadas en clase introStats

getClases <- function(x){
  out <- NA
  if(x <= 1){
    out <- 1
  }
  
  if( x > 1 & x <= 4 ){
    out <- 2
  }
  
  if( x > 4 & x <= 9 ){
    out <- 3
  }
  
  out
}

getTrial <- function(...){
  L <- 0
  class <- rep(0, 3)
  flag <- c(1,2,3)
  z <- c()
  
  while( sum(class==flag) != 3 ){ # class != flag
    L <- L + 1
    
    x <- sample(0:9, 1)
    
    y <- getClases(x)
    
    z <- c(z, x)
    if(y == 1){
      class[1] <- 1
    }
    
    if(y == 2){
      class[2] <- 2
    }
    
    if(y == 3){
      class[3] <- 3
    }
  }
  
  list(output = z, numComponents = L)
}

LLN2 <- function(n) {
  output <- vector("list", n)
  
  for(i in seq_len(n)){
    output[[i]] <- rbinom(i,1,0.5)
  }
  
  MU <- sapply(seq_len(n), function(s) mean(unlist(output[[s]])))
  
  SIGMA <- sapply(seq_len(n), function(s) sd(unlist(output[[s]])))
  
  list(media = MU, sd = SIGMA)
}


tiroDadoRojo <- function(...){sample(1:6, 1)}

tiroDadoNegro <- function(...){sample(c(1:6,6,6), 1)}

tomaCartaMazo1 <- function(...){
  output <- c()
  TEMP <- rbinom(n=1,1,0.5)
  
  if(TEMP==1){
    output <- "ROJA"
  } else {
    output <- "NEGRA"
  }
  list(carta=output, flag=TEMP)
}

tomaCartaMazo2 <- function(...){
  output <- c()
  TEMP <- rbinom(n=1,1,0.)
  
  if(TEMP==1){
    output <- "ROJA"
  } else {
    output <- "NEGRA"
  }
  
  list(carta=output, flag=TEMP)
}

# --- added on May 19, 2022

# --- interpolation hybrid method

get_pixel_matrix <- function(x,lenPeriod=23){
  output <- matrix(nrow=length(x)/lenPeriod, ncol=lenPeriod)
  
  for(i in seq_len(nrow(output))){
    output[i,] <- x[((i-1) * lenPeriod + 1):(i * lenPeriod)]
  }
  output
}


# NOTE: for this to work, length(x) must be a multiple of lenPeriod
climatology <- function(x, lenPeriod){
  MAT <- get_pixel_matrix(x=x, lenPeriod=lenPeriod)
  
  BOXPLOT <- boxplot(MAT, plot=FALSE)
  
  list(matrix=MAT, boxplot=BOXPLOT)
}

# --- March, 18, 2022

LoadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env)
}


# --- added on May 19, 2022

## created on April 29, 2022

getSigma <- function(m){
  apply(m, MARGIN = 2, FUN=sd)
}

## created on April 28, 2022

harmonicR_sigma <- function(y, sigma=NULL, lenBasePeriod=length(y),
                            numFreq, delta){
  numImages <- length(y)
  amp <- numeric(numFreq + 1)
  phi <- numeric(numFreq + 1)
  ySmooth <- numeric(numImages)
  
  numRows <- min(2 * numFreq + 1, numImages)
  degree <- 180 / pi
  
  mat <- matrix(0, nrow = numImages, ncol = min(2 * numFreq + 1, numImages))
  vecTs <- 2 * pi * (0:(lenBasePeriod-1)) / lenBasePeriod
  mat[,1] <- 1
  mat[, seq(2, 2 * numFreq, 2)] <- sapply(1:numFreq, function(s) cos( s * vecTs ))
  mat[, seq(2, 2 * numFreq, 2)+1] <- sapply(1:numFreq, function(s) sin( s * vecTs ))
  
  matCopy <- mat
  
  if(!is.null(sigma)){
    if(is.vector(sigma)){
      y <- diag(1/sqrt(sigma)) %*% y
      mat <- diag(1/sqrt(sigma)) %*% mat
    }
  }
  
  za <- t(mat) %*% ( y )
  
  A <- t(mat) %*% mat
  A <- A + diag(rep(1, numRows)) * delta
  A[1,1] <- A[1,1] - delta
  zr <- solve(A) %*% za
  
  # ySmooth <- (mat) %*% zr
  
  if(!is.null(sigma)){
    ySmooth <- matCopy %*% zr
  } else {
    ySmooth <- (mat) %*% zr
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
  phase[phase < 0] <- phase[phase < 0] + 360
  phi[iFreq] <- phase
  
  list(a.coef = ra, b.coef = rb,
       amplitude = amp, phase = phi, fitted = as.numeric(ySmooth))
}


get_smoothing_sigma <- function(m, numFreq=4, sigma, delta=0.1){
  
  smooth_m <- matrix(nrow = nrow(m), ncol = ncol(m))
  
  for(i in 1:nrow(m)){
    smooth_m[i,] <- harmonicR_sigma(y=m[i,], sigma=sigma,
                                    numFreq=numFreq, delta=delta)$fitted
  }
  
  smooth_m
}

# --- this function applies haRmonics to each row of m
get_smoothing <- function(m, method="harmR", numFreq=4, delta=0.1){
  smooth_m <- matrix(nrow = nrow(m), ncol = ncol(m))
  for(i in 1:nrow(m)){
    smooth_m[i,] <- haRmonics(y=m[i,], method=method,
                              numFreq=numFreq, delta=delta)$fit
  }
  
  smooth_m
}

# --- created on May 2, 2022

getAICtable <- function(observed, fitted, freqNum, sigma=NULL,
                        method=c("iid", "heter")){
  
  AICm <- numeric(nrow(fitted))
  
  if(method=="iid"){
    
    for(i in 1:length(AICm)){
      AICm[i] <- 2 * ( 2 * freqNum + 1) - 2 * sum( (observed[i,] - fitted[i,])^2 / (2 * sd(observed[i,])^2) )
    }
    
  } else {
    
    # AICm <- numeric(nrow(mat_smooth_sigma)) # matrix(nrow=2, ncol=22)
    
    if(is.null(sigma)){
      
      stop("sigma must be provided")
      
    } else {
      
      for(i in 1:length(AICm)){
        AICm[i] <- 2 * ( 2 * freqNum + 1) - 2 * sum( (observed[i,] - fitted[i,])^2 / sigma[i]^2 )/2
      }
      
    }
  }
  
  AICm
}

# --- Added on Feb 5, 2025
# --- Taken from index.Rmd

getMA <- function(y, h, weight){# y must be ts object
  out <- y
  out[c(1:h,(length(y)-h+1):length(y))] <- NA
  out[(1+h):(length(y)-h)] <- sapply( (1+h):(length(y)-h), 
                                      function(s) sum(y[(s-h):(s+h)])/weight )
  out
}


getSD <- function(m){
  v <- numeric(nrow(m))
  H <- solve( t(m) %*% m )
  
  for(i in seq_len(length(v))){
    v[i] <- m[i,] %*% H %*% m[i,]
  }
  
  v
}

getDesignMat <- function(numFreq, freq=23){
  
  # numFreq <- 23
  # freq <- 23*9
  
  mat <- matrix(0, nrow = freq, ncol = min(2 * numFreq + 1, freq))
  vecTs <- 2 * pi * (0:(freq-1)) / freq
  mat[,1] <- 1
  mat[, seq(2, 2 * numFreq, 2)] <- sapply(1:numFreq, function(s) cos( s * vecTs ))
  mat[, seq(2, 2 * numFreq, 2)+1] <- sapply(1:numFreq, function(s) sin( s * vecTs ))
  
  mat
}

vecToMatrix <- function (x, lenPeriod = 23) {
  if (length(x)%%lenPeriod != 0) {
    stop("Length of 'x' must be a multiple of 'lenPeriod'")
  }
  output <- matrix(nrow = length(x)/lenPeriod, ncol = lenPeriod)
  for (i in seq_len(nrow(output))) {
    output[i, ] <- x[((i - 1) * lenPeriod + 1):(i * lenPeriod)]
  }
  output
}
