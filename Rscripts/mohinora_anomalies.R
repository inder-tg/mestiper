
# --- Elaborado Mar 22, 2024
# --- PRELIM version in mohinora_anomalies.R para 
# --- el Diplomado en Geomática 2023

# --- En este script presentamos un análisis de anomalías.
# --- Usamos código en paralelo para eficientar el cómputo 

# --- DATASET: NDVI MOD13Q1 v061 en Cerro Mohinora, Chihuahua, 2000-2023 

# --- ADDicionalmente, este script requiere archivos
# --- MOD13Q1_061_250m_16_days_NDVI_interpol.tif
# --- creado con el archivo mohinora_temporal_gapfilling.R

# --- ACTUALIZACION: DATASET NDVI y EVI MOD13Q1 v061 Cerro Mohinora, 2000-2025


# --- Preámbulo
library(raster)
# library(terra)
library(rasterVis)
library(mapview)
library(RColorBrewer)
library(gtools)
library(foreach)
library(doParallel)
library(geoTS)
library(sf)
library( patchwork )

source( paste0( getwd(), "/Rscripts/auxFUN.R" ) )

# ---

# # mestiperDIR <- "C:/Users/inder/OneDrive/Desktop/mestiper" # USUARIO: modificar
# 
# mohinora_NDVI_DIRS <- list.dirs(path = paste0( getwd(), "/TIF" ))[-1]
# 
# # mohinora_NDVI_QA <- list.dirs(path = paste0( getwd(), "/data" ))[-1]
# 
# FILES_NDVI_imputation <- list.files(path = mohinora_NDVI_DIRS[1],
#                                     pattern = ".tif",
#                                     full.names = TRUE)

# --- DATA loading

FILES_NDVI_imputation <- list.files( path = here( "data", "outputs", "mohinora_imputation" ),
                                     pattern = "NDVI",
                                     full.names = TRUE )

FILES_NDVI_interpolation <- list.files(path = here( "data", "mohinora_2026", "250m_16_days_NDVI_interpol" ),
                                       pattern = ".tif",
                                       full.names = TRUE)

NDVI_files <- c(FILES_NDVI_imputation, FILES_NDVI_interpolation)

mohinora_NDVI <- rast(NDVI_files)

mohinora_NDVI_rTp <- spRast_valuesCoords(mohinora_NDVI) # rasterToPoints(mohinora_DATA) #spRast_valuesCoords(mohinora_DATA)

# mohinora_NDVI_rTp_coords <- mohinora_DATA_rTp[,1:2]
# 
# mohinora_NDVI_rTp_values <- mohinora_DATA_rTp[,3:ncol(mohinora_DATA_rTp)]

# mohinoraSHP <- paste0( getwd(), "/mestiper/RData" )
# 
# RDatafiles <- list.files(path = mohinoraSHP,
#                          pattern = ".RData",
#                          full.names = TRUE)
# 
# mohinora_SHP <- LoadToEnvironment(RDatafiles[1])$mohinora_SHP_sinusoidal

# mohinora_SHP_st <- st_transform(x=mohinora_SHP, crs=crs(mohinora_DATA_interpol))


DIR <- list.dirs( here("data") )
lst_dirs <- setNames( as.list(DIR), basename(DIR) )

usvFILES <- list.files(path = lst_dirs$mohinora_usv7,
                       full.names = TRUE,
                       pattern = ".shp$")

mohinora_USV <- read_sf(usvFILES)

usv_COLORS <- c("#A1E5A5", "#E9D66B", "#00A877", 
                "#66B032", "#83A4F0", "#FC8FAB", "#F500A1")

usv_NAMES <- c("Pino", "Pastizal", "Pino-Encino",
               "Ayarin", "Agro", "Arbustiva", "Arborea")

# ---

# -----------------------------------------------
# --- Análisis de cambio abrupto: anomalías --- #
# -----------------------------------------------

plot(subset(mohinora_NDVI, 453))
# lines( mohinora_SHP, lwd=4 )

lines(mohinora_USV, col = usv_COLORS, lwd = 6)
# plot(mohinora_USV, col = usv_COLORS, lwd = 6, add=TRUE)

# Añadir la leyenda
legend("topleft", # posición en el gráfico
       legend = usv_NAMES,         # nombres de las categorías
       col = usv_COLORS,           # colores de borde
       lwd = 6,                    # grosor de línea en la leyenda
       pt.cex = 1,                  # tamaño del texto
       bty = "n",                  # sin caja alrededor
       inset = c(0.015, 0.05))

# --- Para analizar la serie de tiempo de cualquier píxel 
# --- en la imagen sigue estos pasos:
# --- 1. Ejecuta la línea de abajo
XY <- locator()
# --- 2. Haz click (SOLO UNA VEZ) en algún píxel en la imagen
# --- 3. Presiona la tecla ESC de tu teclado
# --- 4. Continúa con el script a partir de la línea 70

xy <- get_timeSeries_byClicking(c(XY$x, XY$y),
                                        df=mohinora_NDVI_rTp$coords)

pixel <- mohinora_NDVI_rTp$values[xy$coord,]

pixel_mat <- get_pixel_matrix(pixel * 1e-4 )

# pixel_mat <- get_pixel_matrix(mohinora_DATA_interpol_rTp$values[2074,])

# --- promedio por fecha de adquisicion
mu <- apply(pixel_mat, 2, mean)

# --- desviacion estándar por fecha de adquisicion
sigma <- apply(pixel_mat, 2, sd)

# --- anomalias estandarizadas
anomalia <- (pixel_mat - mu)/sigma

pixel_ts <- ts(c(t(pixel_mat)), start = c(2000,1), end = c(2025,23), 
               frequency = 23)

anomalia_ts <- ts(c(t(anomalia)), start = c(2000,1), end = c(2025,23), 
                  frequency = 23)

# plot(pixel_ts)
# plot(anomalia_ts)

# --- Alternativa 1

# df <- data.frame(
#   time = time(pixel_ts),
#   series1 = as.numeric(pixel_ts) * 1e-4,
#   series2 = as.numeric(anomalia_ts)
# )
# 
# # Plot interactivo con dos series de tiempo
# plot_ly(df, x = ~time) %>%
#   add_lines(y = ~series1, name = "NDVI", 
#             line = list(color = "darkgreen")) %>%
#   add_lines(y = ~series2, name = "Anomalia", 
#             line = list(color = "blue")) %>%
#   layout(title = "NDVI vs. Anomalía",
#          xaxis = list(title = "Años"),
#          yaxis = list(title = "NDVI"))

# --- Alternativa 2

# Convertir ts a data.frame
df_pixel <- data.frame(
  time = time(pixel_ts),
  value = as.numeric(pixel_ts)
)

df_anomalia <- data.frame(
  time = time(anomalia_ts),
  value = as.numeric(anomalia_ts)
)

# Primer plot
p1 <- ggplot(df_pixel, aes(x = time, y = value)) +
  geom_line(color = "darkgreen") +
  labs(title = "NDVI", x = "Años", y = "NDVI")

# Segundo plot
p2 <- ggplot(df_anomalia, aes(x = time, y = value)) +
  geom_line(color = "red") +
  labs(title = "Serie de anomalías", x = "Años", y = "Anomalía NDVI")

p1 + p2

# --- Recordando la densidad normal
# --- y el cálculo de algunas probabilidades

x <- seq(-5,5, by = 0.05)
y <- dnorm(x = x)

# --- APROX
lB <- c(-1,-2,-3,-4) # lB:= lower bound, límite inferior o límite izquierdo
uB <- c(1,2,3,4) # uB:= upper bouund, límite superior o límite derecho

# --- EXACT
# lB <- c(-qnorm(0.85),
#         -qnorm(0.975),
#         -qnorm(0.999),
#         -qnorm(0.999975)) # lB:= lower bound, límite inferior o límite izquierdo
# 
# uB <- c(qnorm(0.85),
#         qnorm(0.975),
#         qnorm(0.999),
#         qnorm(0.999975)) # uB:= upper bouund, límite superior o límite derecho

# --- 1sigma, aprox 85% (qnorm(0.85);) de 
# --- la probabilidad total
# --- está concentrada en esta región
xz1 <- x
xz1[x >= uB[1]] <- NA
xz1[x <= lB[1]] <- NA

yz1 <- y
yz1[x >= uB[1]] <- NA
yz1[x <= lB[1]] <- NA

# --- 2 sigma, aprox 97.5% (qnorm(0.975);) de 
# --- la probabilidad total
# --- está concentrada en esta región

xz2 <- x
xz2[x >= uB[2]] <- NA
xz2[x <= lB[2]] <- NA

yz2 <- y
yz2[x >= uB[2]] <- NA
yz2[x <= lB[2]] <- NA

# --- 3 sigma, aprox 99.9% (qnorm(0.999);) de 
# --- la probabilidad total
# --- está concentrada en esta región

xz3 <- x
xz3[x >= uB[3]] <- NA
xz3[x <= lB[3]] <- NA

yz3 <- y
yz3[x >= uB[3]] <- NA
yz3[x <= lB[3]] <- NA

# --- 4 sigma, aprox 99.9975% (qnorm(0.999975);) de 
# --- la probabilidad total
# --- está concentrada en esta región

xz4 <- x
xz4[x >= uB[4]] <- NA
xz4[x <= lB[4]] <- NA

yz4 <- y
yz4[x >= uB[4]] <- NA
yz4[x <= lB[4]] <- NA

# ----
# --- PRESTAR atención a los colores

myPal <- brewer.pal('RdYlGn', n=8)

yRan <- range(y,yz1,yz2,yz3,yz4,na.rm = T)

plot(x, y, type = "l", col = "red", lwd = 3,
     ylab = "", xlab = "", main = "",
     ylim=yRan)

# --- 1sigma region
a <- x[!is.na(yz1)]
b <- yz1[!is.na(yz1)]

# cbind(a,b)

polygon(x=c(a[21:39],rev(a[21:39])),
        y=c(b[21:39],rep(0,19)),
        col = myPal[5], border = NA)

polygon(x=c(a[1:21],rev(a[1:21])),
        y=c(b[1:21],rep(0,21)),
        col = myPal[4], border = NA)

abline(v=lB[1], col=myPal[4], lwd=2)
abline(v=uB[1], col=myPal[5], lwd=2)

# --- 2sigma region
a2 <- x[!is.na(yz2)]
b2 <- yz2[!is.na(yz2)]

cbind(a2,b2)

polygon(x=c(a2[60:79],rev(a2[60:79])),
        y=c(b2[60:79],rep(0,20)),
        col = myPal[6], border = NA)

polygon(x=c(a2[1:20],rev(a2[1:20])),
        y=c(b2[1:20],rep(0,20)),
        col = myPal[3], border = NA)

abline(v=lB[2], col=myPal[3], lwd=2)
abline(v=uB[2], col=myPal[6], lwd=2)

# --- 3sigma region
a3 <- x[!is.na(yz3)]
b3 <- yz3[!is.na(yz3)]

cbind(a3,b3)

polygon(x=c(a3[101:119],rev(a3[101:119])),
        y=c(b3[101:119],rep(0,19)),
        col = myPal[7], border = NA)

polygon(x=c(a3[1:22],rev(a3[1:22])),
        y=c(b3[1:22],rep(0,22)),
        col = myPal[2], border = NA)

abline(v=lB[3], col=myPal[2], lwd=2)
abline(v=uB[3], col=myPal[7], lwd=2)

# --- 4sigma region
# par(new = TRUE)
# plot(x,yz4, type = "h",
#      ylim = yRan,
#      xlab = "z", ylab = "")
abline(v=lB[4], col=myPal[1], lwd=2)
abline(v=uB[4], col=myPal[8], lwd=2)

# --- obj anomalia como un obj 'ts'
anomalia_ts <- ts(c(t(anomalia)), 
                  start = c(2000,1), 
                  end = c(2023,23),
                  frequency = 23)

# --- PRESTAR atención a los colores
# --- ligar las líneas de abajo con las
# --- regiones de probabilidades definidas a partir
# --- de la densidad normal

plot(anomalia_ts, ylab="")
abline(h=1, col=myPal[5], lwd=3) # una desviación estándar
abline(h=-1, col=myPal[4], lwd=3) 

abline(h=2, col=myPal[6], lwd=3) # 2 desviaciones estándar
abline(h=-2, col=myPal[3], lwd=3)

abline(h=3, col=myPal[7], lwd=3) # 3 desviaciones estándar
abline(h=-3, col=myPal[2], lwd=3)

abline(h=4, col=myPal[8], lwd=3) # 4 desviaciones estándar
abline(h=-4, col=myPal[1], lwd=3)

# ---

startYear <- 2004
endYear <- 2005
Title <- paste0("pixel de ", startYear, " a ", endYear)

plot(anomalia_ts, ylab="", 
     xlim=c(startYear, endYear),
     main=Title)
abline(h=1, col=myPal[5], lwd=3) # una desviación estándar
abline(h=-1, col=myPal[4], lwd=3) 

abline(h=2, col=myPal[6], lwd=3) # 2 desviaciones estándar
abline(h=-2, col=myPal[3], lwd=3)

abline(h=3, col=myPal[7], lwd=3) # 3 desviaciones estándar
abline(h=-3, col=myPal[2], lwd=3)

abline(h=4, col=myPal[8], lwd=3) # 4 desviaciones estándar
abline(h=-4, col=myPal[1], lwd=3)


# ------------------------------------------
# --- Anomalías: cómputo en paralelo --- #
# ------------------------------------------

# --- CODIGO EN PARALELO

df_anomalias2000 <- matrix(nrow=nrow(mohinora_NDVI_rTp$values), 
                           ncol=26) # primeras 2 columnas x, y, restantes 23 los valores de las anomalias
df_anomalias2000[,1:2] <- mohinora_NDVI_rTp$coords

df_anomalias2004 <- matrix(nrow=nrow(mohinora_NDVI_rTp$values), 
                           ncol=26)
df_anomalias2004[,1:2] <- mohinora_NDVI_rTp$coords

df_anomalias2008 <- matrix(nrow=nrow(mohinora_NDVI_rTp$values), 
                           ncol=26)
df_anomalias2008[,1:2] <- mohinora_NDVI_rTp$coords

numCores <- detectCores()

# --- Asegurarse de crear /RData/progressReports/mohinora
progressReportFile <- paste0(getwd(), "/RData/progressReports/mohinora_anomalies.txt" )
file.create(path=progressReportFile, showWarnings=FALSE)

write("===ANOMALIES computation began at===",
      file=progressReportFile, append=TRUE)
write(as.character(Sys.time()[1]), file=progressReportFile,
      append=TRUE)

kluster <- parallel::makeCluster(numCores-1, outfile="")
registerDoParallel(kluster)

# --- Hacer una prueba pequeña antes de ejecutar la línea
# --- de abajo

output <- foreach(i=1:nrow(mohinora_NDVI_rTp$values), .combine="rbind") %dopar% { 
  
  pixel <- mohinora_NDVI_rTp$values[i, ] #* 1e-4
  
  pixel_mat <- get_pixel_matrix(pixel)
  
  mu <- apply(pixel_mat, 2, mean, na.rm=TRUE)
  
  sigma <- apply(pixel_mat, 2, sd, na.rm=TRUE)
  
  anomalia <- (pixel_mat - mu) / sigma
  
  anomalia_2000 <- anomalia[1,]
  
  anomalia_2004 <- anomalia[5,]
  
  anomalia_2008 <- anomalia[9,]
  
  s <- c(anomalia_2000, anomalia_2004, anomalia_2008)
  
  if(i %% 100 ==0){
    texto <- paste0("Working on ROW: ", i)
    write(texto, file=progressReportFile, append=TRUE)
  }
  
  return(s)
}
stopCluster(kluster)

write( as.character(Sys.time()[1]), file=progressReportFile, append=TRUE)
write( "===ANOMALIES analysis ended here===", file=progressReportFile, append=TRUE)
# ---

# --- Guardando las anomalías como objetos matrix
df_anomalias2000[,3:25] <- output[,1:23]
df_anomalias2004[,3:25] <- output[,24:46]
df_anomalias2008[,3:25] <- output[,47:69]

# --- Asegurarse de haber creado el folder /RData/mohinora_anomalies

save(df_anomalias2000, file=paste0(getwd(),"/RData/mohinora_anomalies/2000.RData"))
save(df_anomalias2004, file=paste0(getwd(),"/RData/mohinora_anomalies/2004.RData"))
save(df_anomalias2008, file=paste0(getwd(),"/RData/mohinora_anomalies/2008.RData"))
# ---

# -----------------------
# --- RASTERIZACION --- #
# -----------------------

PROJECTION <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"

# --- descomentar las siguientes 3 lineas
# --- si no tienes creados los objetos df_anomalias2000,
# --- df_anomalias2004 y df_anomalias2008
# df_anomalias2000 <- LoadToEnvironment(paste0(getwd(),"/RData/mohinora_anomalies/2000.RData"))$df_anomalias2000
# df_anomalias2004 <- LoadToEnvironment(paste0(getwd(),"/RData/mohinora_anomalies/2004.RData"))$df_anomalias2004
# df_anomalias2008 <- LoadToEnvironment(paste0(getwd(),"/RData/mohinora_anomalies/2008.RData"))$df_anomalias2008

map_anomalies2000 <- brick()
map_anomalies2004 <- brick()
map_anomalies2008 <- brick()

for( i in 1:23 ){
  temp2000 <- matrixToRaster(matrix=df_anomalias2000[,c(1,2,i+2)],
                             projection=PROJECTION)
  
  temp2004 <- matrixToRaster(matrix=df_anomalias2004[,c(1,2,i+2)], 
                             projection=PROJECTION)
  
  temp2008 <- matrixToRaster(matrix=df_anomalias2008[,c(1,2,i+2)], 
                             projection=PROJECTION)
  
  map_anomalies2000 <- addLayer(map_anomalies2000, temp2000)
  
  map_anomalies2004 <- addLayer(map_anomalies2004, temp2004)
  
  map_anomalies2008 <- addLayer(map_anomalies2008, temp2008)
}

DoY <- seq(1, 365, by=16)
LABELS <- sapply(1:length(DoY), function(s) paste0("DoY-", DoY[s]) )

names(map_anomalies2000) <- LABELS

names(map_anomalies2004) <- LABELS

names(map_anomalies2008) <- LABELS


# --- Asegurarse de haber creado el folder /TIF/mohinora_anomalies

writeRaster(map_anomalies2000,
            filename = here("data", "outputs", "mohinora_anomalies", "2000" ),
            
            # paste0( getwd(), "/TIF/mohinora_anomalies/2000" ),
            
            format="GTiff", datatype="FLT4S", overwrite=TRUE)

writeRaster(map_anomalies2004,
            filename = here( "data", "outputs", "mohinora_anomalies", "2004" ),
              # paste0( getwd(), "/TIF/mohinora_anomalies/2004" ),
            format="GTiff", datatype="FLT4S", overwrite=TRUE)

writeRaster(map_anomalies2008,
            filename = here( "data", "outputs", "mohinora_anomalies", "2008" ),
              # paste0( getwd(), "/TIF/mohinora_anomalies/2008" ),
            format="GTiff", datatype="FLT4S", overwrite=TRUE)

# -----------------------
# --- VISUALIZACION --- #
# -----------------------
# --- PRESTAR atención a los colores
# --- interpretación de las áreas coloreadas
# --- es similar a la mostrada con las regiones
# --- de probabilidad normal

myPal <- brewer.pal('RdYlGn', n=7)
myTheme <- rasterTheme(region = myPal)

levelplot(map_anomalies2000, main="Anomalías 2000", 
          par.settings = myTheme,
          at = seq(-6,6,by=1))

levelplot(map_anomalies2004, main="Anomalías 2004", 
          par.settings = myTheme,
          at = seq(-6,6,by=1))

levelplot(map_anomalies2008, main="Anomalías 2008",
          par.settings = myTheme,
          at = seq(-6,6,by=1))
# ---


plot(subset(map_anomalies2008,7))
