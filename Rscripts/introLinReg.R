# -----------------------------------------------------------------------------
# Written on Aug 24, 2018
#
# En este script mostraremos cómo calcular estadísticos asociados al modelo de 
# regresión lineal
# -----------------------------------------------------------------------------

test <- landsatCrop[[1]]

testR <- rasterToPoints(test)

test2 <- landsatCrop[[2]]

test2R <- rasterToPoints(test2)

head(test2R)

test2R$ultra.blue

ultra.blue <- testR[,3]
blue <- test2R[,3]

plot(ultra.blue, blue)

out <- lm(blue~ultra.blue)
summary(out)

# ---
# z-scores
# ---

getZscores <- function(x){
  (x - mean(x))/sd(x)
}

zy <- getZscores(blue)

zx <- getZscores(ultra.blue)

plot(zx, zy)

# ---
# correlation
# ----
getCor <- function(x,y){
  sum(x*y)/(length(x)-1)
}

r <- getCor(zx, zy)

#slope
b1 <- r * sd(blue) / sd(ultra.blue)

#intercept
b0 <- mean(blue) - b1 * mean(ultra.blue)

#estimated
yHat <- b0 + b1 * ultra.blue

se <- sqrt( sum( (blue - yHat)^2 ) / (length(blue)-2) )

# t-test
t <- b1 / (se/ (sqrt(length(blue)-1) * sd(ultra.blue) ) )

summary(out)

b0
b1
t
1-pt(t, length(blue)-2)

plot(ultra.blue, blue, pch = 20)
lines(ultra.blue, yHat, col = "blue")



