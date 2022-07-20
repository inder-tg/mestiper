# -----------------------------------------------------------------------------
#
# Written by Inder Tecuapetla on Aug 8, 2018
#
# En este script daremos una breve introducción al lenguaje de programación
# estadística R
# Basado en el Cap. 2 de "Análisis espacial con R: Usa R como un Sistema de 
# Información Geográfica de Jean-Francois Mas (European Scientific Institute, Publishing)
#
# -----------------------------------------------------------------------------

# ----------------
# CONOCE RStudio #
# ----------------

# ------------------------------------------------
# Objetos en R: vector, dataframe, matrix, lista #
# ------------------------------------------------

# ---

Prec <- c(15, 40, 37, 0, 0, 0, 0, 0, 7, 3, 77)
meses <- c("Enero", "Febrero", "Marzo", "Abril")
numeros <- 1:4

# ---

# ------------------
# IMPORTANDO datos #
# ------------------

getwd()

setwd(paste(getwd(), "/data/recursos-mx", sep = ""))

getwd()

# ------------------------------------------------------------------------------
# NOTA: En Windows, se utiliza la diagonal invertida \ o una doble diagonal // #
# ------------------------------------------------------------------------------

tabEnsenada <- read.csv("ensenada.csv")

# clase de objeto tabEnsenada
class(tabEnsenada)

# estructura de objeto tabEnsenada
str(tabEnsenada)

# primeras filas
head(tabEnsenada)

# primeras 3 filas
head(tabEnsenada, 3)

# últimas 2 filas
tail(tabEnsenada, 2)

# nombre de las variables
names(tabEnsenada)

# accesando la variable P
tabEnsenada$P

# son vectores los elementos de tabEnsenada?
is.vector(tabEnsenada$P)

tabEnsenada$PCvid # periodo de crecimiento de la uva

# -------------------------
# OPERACIONES matemáticas #
# -------------------------

# calcular rango y agreagarlo como nueva variable a tabEnsenada
tabEnsenada$rango <- tabEnsenada$Tmax - tabEnsenada$Tmin

tabEnsenada

# seleccionando elementos del objeto

tabEnsenada[1,2]

tabEnsenada[, c(1, 2, 6)]

tabEnsenada[, c("mes", "P", "Tprom")]

tabEnsenada[-(10:12), -4]

# selección basad en condiciones

tabEnsenada[tabEnsenada$Tmax > 25,]

subset(tabEnsenada, Tmax > 25)

tabEnsenada[tabEnsenada$PCvid == "si", ]

subset(tabEnsenada, PCvid == "si")

subset(tabEnsenada, PCvid == "no", select = c("mes","P","Tprom"))

# -----------------------------------------------------------------------------
# guardando tabEnsenada en un archivo de texto

write.table(tabEnsenada, file = "tablaEnsenada.txt")

# -----------------------------------------------------------------------------
# coersión de un objeto a otro de diferente clase

class(tabEnsenada[,1:7])

matriz <- as.matrix(tabEnsenada[,1:7])

class(matriz)

matriz

tabEnsenada[2,4]

matriz[2,4]
# -----------------------------------------------------------------------------

# ---------------
# Gráficas en R #
# ---------------

plot(tabEnsenada$dias, tabEnsenada$P)

par(mar = c(4.5, 5, 1.5, 2))
plot(tabEnsenada$dias, tabEnsenada$P, main = "Días de lluvia / Precipitación",
     xlab = "Días de lluvia", ylab = "Precipitación mensual (mm)", col = "darkblue",
     bg = "blue", pch = 22)

par(mar = c(5, 5, 3, 2))
plot(tabEnsenada$dias, tabEnsenada$P, main = "Días de lluvia / Precipitación",
     xlab = "Días de lluvia", ylab = "Precipitación mensual (mm)", col = "darkblue",
     bg = "blue", pch = 22, cex = 1.2, cex.lab = 0.5, cex.axis = 2)

# ?par

# --------------------------
# RELACIÓN entre variables #
# --------------------------

cor(tabEnsenada$dias, tabEnsenada$P)

fit <- lm(tabEnsenada$P ~ tabEnsenada$dias)

fit

fit$coefficients

summary(fit)

resumen <- summary(fit)

class(resumen)

resumen

str(resumen)

resumen$coefficients

yr <- range(tabEnsenada$P)
yr[1] <- yr[1] - 1.25
yr[2] <- yr[2] + 1.25

par(mar = c(5, 5, 3, 2))
plot(tabEnsenada$dias, tabEnsenada$P, main = "Días de lluvia / Precipitación",
     ylim = yr,
     xlab = "Días de lluvia", ylab = "Precipitación mensual (mm)", col = "darkblue",
     bg = "blue", pch = 16, cex = 1.2, cex.lab = 1.2, cex.axis = 1.2)
lines(fit$model$`tabEnsenada$dias`, fit$fitted.values)

# -------------------------------------------
# LISTA: un objeto para contenerlos a todos #
# -------------------------------------------

dock <- list(Prec, fit, "lista")

dock

# --------------------------------------
# OPERACIONES marginales sobre objetos #
# --------------------------------------

matriz <- matrix(c(1,2,2,3,6,0,4,7,9), ncol = 3)

matriz <- matrix(c(1,2,2,3,6,0,4,7,9), ncol = 3, byrow = T)

# suma todos las entradas de cada columns
colSums(matriz)

# suma todos las entradas de cada renglón
rowSums(matriz)

colMeans(matriz)

# utilizando la función "apply"
?apply

apply(matriz, 2, sum)

apply(matriz, 1, sum)

apply(matriz, 1, mean)

apply(matriz, 2, sd)

?lapply
?sapply

# -----------------------
# creación de funciones #
# -----------------------

# SINTAXIS
# nameFunction <- function(parametros) {
#   comandos que dependen de parametros
# salida
# }

getNombre <- function(a, b){
  cat("Mi nombre es", a, "\n")
  cat("Y mi apellido", b, "\n")
}

getNombre(a = "Inder", b = "Tecuapetla")

binomioCuadrado <- function(a, b){
  (a + b)^2
}

binomioCuadrado(10, 2)

# Ejercicio: Escribe una función que calcule el NDVI de una imagen,
# los parámetros de esta función son "imagen" y "bandaX", "bandaY"
# el resultado debe ser el NDVI asociado a esta imagen

# -------------------------------------------
# CONDICIONAL IF: if (condicion) {haz esto} #
# -------------------------------------------

decision <- function(x){
  if(x < 0){
    cat("El número ingresado es:", x, "\n")
  }
  print("El condicional IF es fácil")
}


# ----------------------------------------------------------------------
# CONDICIONAL IF - ELSE: if(condicion) {haz esto} else {haz esto otro} #
# ----------------------------------------------------------------------

# Función para verificar si un número es par o impar

imPar <- function(x){
  cat("Usted ha ingresado el número:", x, "\n")
  if( x %% 2 == 0){
    cat(x, "es un número par", "\n")
  }else{
    cat(x, "es un número impar", "\n")
  }
}

# Ejercicio: crea una función que relacione dos números usando "=", ">", o "<"

# ---------------
# FOR LOOP en R #
# ---------------

# SINTAXIS: for(i in lista-valores) { comandos a ejecutar para cada valor de i }

for(i in 1:6){
  cat("En ", i,"-ésima iteración, i =",i, "\n")
}

fac <- 1
for(i in 1:4){
  fac <- fac * i
  cat("En ", i,"-ésima iteración, fac =", fac, "\n")
}
fac

factorial(4)

# Función para calcular la suma de los primeros N números naturales

sumN <- function(x){
  cat("Usted ha ingresado el número:", x, "\n")
  
  if(x < 0){
    stop("x debe ser mayor que cero")
  }
    
  if(x - floor(x) != 0){
    stop("x debe ser un número natural")
  }
  
  cont <- 0
  for(i in 1:x){
    cont <- cont + i
  }
  
  cat("Sum =", cont)
}

# Ejercicio: Escribe una función para leer cada archivos .csv en el 
# directorio de trabajo



