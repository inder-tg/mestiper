
# --- Escrito por Inder Tecuapetla, Abril 29, 2023
# --- Análisis de clasificación de tendencias sobre el 
# --- DATASET: NDVI MOD13Q1 v061 en Cerro Mohinora, Chihuahua, 2000-2023 

# --- En este script presentamos un análisis de punto de cambio
# --- empleando bfast01 
# --- Usamos código en paralelo para eficientar el cómputo 

# --- ADicionalmente, este script requiere archivos
# --- MOD13Q1_061_250m_16_days_NDVI_interpol.tif
# --- creado con el archivo mohinora_temporal_gapfilling.R

library(terra)
library(foreach)
library(doParallel)
library(geoTS)
library(sf)
library(tmap)
library(bfast)
library(tidyverse)
library(kableExtra)
library(here)

source( "Rscripts/auxFUN.R" )

# --- DATA Loading

mohinora_NDVI_DIRS <- list.dirs(path = here("data", "mohinora_2026") )

FILES_NDVI_imputation <- list.files(path = here("data", "outputs", "mohinora_imputation"),
                                    pattern = "NDVI",
                                    full.names = TRUE)

FILES_NDVI_interpolation <- list.files(path = mohinora_NDVI_DIRS[5],
                                       pattern = ".tif",
                                       full.names = TRUE)

mohinora_NDVI_files <- c(FILES_NDVI_imputation, FILES_NDVI_interpolation[1:595])

mohinora_DATA <- rast(mohinora_NDVI_files)

mohinora_DATA_rTp <- spRast_valuesCoords(mohinora_DATA) #rasterToPoints(mohinora_DATA) #spRast_valuesCoords(mohinora_DATA)

# mohinora_DATA_rTp_coords <- mohinora_DATA_rTp[,1:2]
# 
# mohinora_DATA_rTp_values <- mohinora_DATA_rTp[,3:ncol(mohinora_DATA_rTp)]

# mohinoraSHP <- paste0( getwd(), "/mestiper/RData" )
# 
# RDatafiles <- list.files(path = mohinoraSHP,
#                          pattern = ".RData",
#                          full.names = TRUE)
# 
# mohinora_SHP <- LoadToEnvironment(RDatafiles[1])$mohinora_SHP_sinusoidal

shpFiles <- list.files(path = here("data", "mohinora_usv7"),
                       pattern = ".shp",
                       full.names = TRUE)

mohinora_USV <- read_sf(shpFiles[1])

usv_COLORS <- c("#A1E5A5", "#E9D66B", "#00A877", 
                "#66B032", "#83A4F0", "#FC8FAB", "#F500A1")

usv_NAMES <- c("Pino", "Pastizal", "Pino-Encino",
               "Ayarin", "Agro", "Arbustiva", "Arborea")

# -----------------------------------------------
# --- Clasificación de tendencias --- #
# -----------------------------------------------

plot(subset(mohinora_DATA,453))
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


XY <- locator() 

xy <- get_timeSeries_byClicking(c(XY$x, XY$y), 
                                df=mohinora_DATA_rTp$coords)

# 2564, 2298, buenos pixeles para cps
pixel <- mohinora_DATA_rTp$values[xy$coord, ] * 1e-4

# --- objeto ts
pixel_ts <- ts(pixel, 
               start = c(2000,1), 
               end = c(2025,23),
               frequency = 23)
# ---

pixel_bfast01 <- bfast01( data = pixel_ts )

plot(pixel_bfast01)

pixel_bfast01$breakpoints

bfast01classify(pixel_bfast01)

getYear(start=2000, end=2025, 
        bp=pixel_bfast01$breakpoints, freq=23)

# ---

# --- TYPES
# --- 1 - monotonic increase
# --- 2 - monotonic decrease
# --- 3 - monotonic increase with positive break
# --- 4 - monotonic decrease with negative break
# --- 5 - interruption: increase with negative break
# --- 6 - interruption: decrease with positive break
# --- 7 - reversal: increase to decrease
# --- 8 - reversal: decrease to increase

# ---

# --- Clase de tendencia: 7; Año: 2006

XY <- list(x=--10698418, y=2887497)

xy <- get_timeSeries_byClicking(c(XY$x, XY$y),
                                df=mohinora_DATA_rTp_coords)

pixel <- mohinora_DATA_rTp_values[xy$coord, ] * 1e-4

pixel_ts <- ts(as.numeric(pixel), 
               start = c(2000,1), end = c(2023, 23),
               frequency = 23)

pixel_bfast01 <- bfast01(pixel_ts, bandwidth = 0.2)

plot(pixel_bfast01)

pixel_bfast01$breakpoints

bfast01classify(pixel_bfast01)

getYear(start=2000, end=2023, 
        bp=pixel_bfast01$breakpoints, 
        freq=23)

# --- Clase de tendencia: 8

XY <- list(x=-10687435, y=2894415)

xy <- get_timeSeries_byClicking(c(XY$x, XY$y),
                                df=mohinora_DATA_rTp_coords)

pixel <- mohinora_DATA_rTp_values[xy$coord, ] * 1e-4

pixel_ts <- ts(as.numeric(pixel), 
               start = c(2000,1), end = c(2023, 23),
               frequency = 23)

pixel_bfast01 <- bfast01(pixel_ts, bandwidth = 0.2)

plot(pixel_bfast01)

pixel_bfast01$breakpoints

bfast01classify(pixel_bfast01)

getYear(start=2000, end=2024, bp=pixel_bfast01$breakpoints, freq=23)

# --- CODIGO EN PARALELO

# TYPE, SIGN, STABLE, YEARS, CP

TYPE <- matrix(nrow=nrow(mohinora_DATA_rTp_coords), ncol=3)
TYPE[,1:2] <- mohinora_DATA_rTp_coords[,1:2]

SIGN <- matrix(nrow=nrow(mohinora_DATA_rTp_coords), ncol=3)
SIGN[,1:2] <- mohinora_DATA_rTp_coords[,1:2]

STABLE <- matrix(nrow=nrow(mohinora_DATA_rTp_coords), ncol=3)
STABLE[,1:2] <- mohinora_DATA_rTp_coords[,1:2]

YEARS <- matrix(nrow=nrow(mohinora_DATA_rTp_coords), ncol=3)
YEARS[,1:2] <- mohinora_DATA_rTp_coords[,1:2]

CP <- matrix(nrow=nrow(mohinora_DATA_rTp_coords), ncol=3)
CP[,1:2] <- mohinora_DATA_rTp_coords[,1:2]

# ---  

numCores <- detectCores()

progressReportFile <- paste0(getwd(), "/RData/progressReports/mohinora/progress_cps.txt" )
file.create(path=progressReportFile, showWarnings=FALSE)

write("===BFAST01 analysis began at===",
      file=progressReportFile, append=TRUE)
write(as.character(Sys.time()[1]), file=progressReportFile,
      append=TRUE)

kluster <- parallel::makeCluster(numCores-1, outfile="")
registerDoParallel(kluster)

output <- foreach(i=1:nrow(mohinora_DATA_rTp_coords), .combine="rbind",
                  .packages=c("bfast") ) %dopar% { # nrow(sp_ndvi_rTp)
                    
                    if(i %% 100 ==0){
                      texto <- paste0("Working on ROW: ", i)
                      write(texto, file=progressReportFile, append=TRUE)
                    }
                    
                    pixel <- mohinora_DATA_rTp_values[i,] * 1e-4
                    
                    pixel_ts <- ts(pixel, 
                                   start = c(2000,1), 
                                   end = c(2023, 23),
                                   frequency = 23)
                    
                    pixel_bfast01 <- bfast01( data=pixel_ts )
                    
                    # TYPE, SIGN, STABLE, YEARS, CP
                    
                    TEMP <- bfast01classify(pixel_bfast01)
                    
                    YEAR <- getYear(start=2000, end=2023, 
                                    bp=pixel_bfast01$breakpoints, freq=23)
                    
                    s <- c(TEMP$flag_type, TEMP$flag_significance, 
                           TEMP$flag_pct_stable, YEAR, 
                           pixel_bfast01$breakpoints)
                    
                    return(s)
                  }
stopCluster(kluster)

write( as.character(Sys.time()[1]), file=progressReportFile, append=T)
write( "===BFAST01 analysis ended here===", file=progressReportFile, append=T)

# ---

TYPE[,3] <- output[,1]
SIGN[,3] <- output[,2]
STABLE[,3] <- output[,3]
YEARS[,3] <- output[,4]
CP[,3] <- output[,5]

# --- asegurarse de haber creado /RData/mohinora_cps_class

dir.create(paste0(getwd(), "/RData/mohinora_cps"), recursive = TRUE)

save(TYPE, file=paste0(getwd(),"/RData/mohinora_cps/TYPE.RData"))
save(SIGN, file=paste0(getwd(),"/RData/mohinora_cps/SIGN.RData"))
save(STABLE, file=paste0(getwd(),"/RData/mohinora_cps/STABLE.RData"))
save(YEARS, file=paste0(getwd(),"/RData/mohinora_cps/YEARS.RData"))
save(CP, file=paste0(getwd(),"/RData/mohinora_cps/CP.RData"))

# -----------------------
# --- RASTERIZACION --- #
# -----------------------

PROJECTION <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"

map_TYPE <- matrixToRaster(matrix=TYPE, projection=PROJECTION)
map_SIGN <- matrixToRaster(matrix=SIGN, projection=PROJECTION)
map_YEARS <- matrixToRaster(matrix=YEARS, projection=PROJECTION)

dir.create( paste0( getwd(), "/TIF/mohinora_cps"), recursive = TRUE )

raster::writeRaster(map_TYPE,
                    filename = paste0( getwd(), "/TIF/mohinora_cps/map_TYPE"),
                    format="GTiff", datatype="INT2S", overwrite=TRUE)

raster::writeRaster(map_YEARS,
                    filename = paste0( getwd(), "/TIF/mohinora_cps/map_YEARS"),
                    format="GTiff", datatype="INT2S", overwrite=TRUE)

raster::writeRaster(map_SIGN,
                    filename = paste0( getwd(), "/TIF/mohinora_cps/map_SIGN"),
                    format="GTiff", datatype="INT2S", overwrite=TRUE)

# --- just the ANP Cerro Mohinora

TYPEmap <- rast(paste0( getwd(), "/TIF/mohinora_cps/map_TYPE.tif" ))

SIGNmap <- rast(paste0( getwd(), "/TIF/mohinora_cps/map_SIGN.tif" ))

signTypeMap <- TYPEmap
signTypeMap[ SIGNmap != 0 ] <- NA

mohinora_cps_TYPE <- terra::crop(x=signTypeMap, 
                                 y=mohinora_SHP,
                                 mask=TRUE)

# writeRaster(mohinora_cps_TYPE,
#             filename = paste0( getwd(), "/TIF/mohinora_cps_class/cps_signType.tif"),
#             datatype="INT2S", overwrite=TRUE)

# -----------------------
# --- VISUALIZACION --- #
# -----------------------

# --- objetos auxiliares para el tmap

COLORES_update <- c("#50C878", "#C08F73",
                    "#006B3C", "#E32636",
                    "#87A96B", "#2F221E",
                    "#F2C185", "#66FF00")

usv_COLORS <- c("#A1E5A5", "#E9D66B", "#00A877", "#66B032", "#83A4F0", "#FC8FAB", "#F500A1")
usv_NAMES <- c("Pino", "Pastizal", "Pino-Encino",
               "Ayarin", "Agro", "Arbustiva", "Arborea")

COLOR_USV <- c(usv_COLORS[1], 
               rep(usv_COLORS[3],2), 
               usv_COLORS[4], 
               rep(usv_COLORS[2],2),
               rep(usv_COLORS[5], length(7:12)), 
               rep(usv_COLORS[6], 9),
               usv_COLORS[7])

# --- definiendo un bbox a usar en el tmap

bbox_new <- st_bbox(mohinora_USV) # current bounding box

xrange <- bbox_new$xmax - bbox_new$xmin # rango x
yrange <- bbox_new$ymax - bbox_new$ymin # rango y

bbox_new[1] <- bbox_new[1] - (0.25 * xrange) # xmin - left
bbox_new[3] <- bbox_new[3] + (0.25 * xrange) # xmax - right
bbox_new[2] <- bbox_new[2] - (0.3 * yrange) # ymin - bottom
bbox_new[4] <- bbox_new[4] + (0.25 * yrange) # ymax - top

bbox_new <- bbox_new %>%  
  st_as_sfc() 

# --- definiedo un st multistring para usarlo en tm_lines()

visual_mohinora <- mohinora_USV %>%
  sf::st_cast("MULTILINESTRING")

visual_mohinora$COLOR <- COLOR_USV

# ---
library(RColorBrewer)

type_map = tm_shape(mohinora_cps_TYPE, 
                    bbox = bbox_new) +
  tm_raster(col.scale = tm_scale(), 
            palette = COLORES_update, 
            legend.show = TRUE,
            col.legend = tm_legend("Trend types")) +
  tm_shape(visual_mohinora) + tm_lines(col="COLOR", lwd=3) + 
  tm_title("", frame = FALSE, bg.color = NA) +
  tm_compass(type = "8star", position = c("right", "bottom")) +
  tm_scalebar(text.size = 0.65,
              position = c("right", "bottom")) +
  tm_add_legend("symbol", 
                labels=usv_NAMES, 
                col=usv_COLORS,
                border.col = "grey40",
                size=1,
                shape=18,
                is.portrait = TRUE)

type_map

# ---

# -------------------------------
# --- RESUMEN dE RESULTADOS --- #
# -------------------------------

# nombres <- unique( mohinora_SHP_USV_st$DESCRIPCIO )

getPercent <- function(v){
  round(100 * v/sum(v), 3)
}

trunk <- mohinora_USV %>%
  dplyr::select( DESCRIPCIO ) %>%
  st_drop_geometry() 

# --- terra::mask

df <- trunk %>%
  apply(MARGIN = 1,
        FUN = function(x) mask(mohinora_cps_TYPE,
                               dplyr::filter(mohinora_USV,
                                             DESCRIPCIO == x))) %>%
  lapply(FUN = spRast_values, na_rm=TRUE) %>%
  lapply(FUN = as.data.frame) %>%
  lapply(FUN = table) %>%
  lapply(FUN=getPercent)

names(df) <- trunk$DESCRIPCIO

df <- df[!duplicated(df)]

# --- terra::crop

df_crop <- trunk %>%
  apply(MARGIN = 1,
        FUN = function(x) crop(mohinora_cps_TYPE,
                               dplyr::filter(mohinora_USV,
                                             DESCRIPCIO == x),
                               mask=TRUE)) %>%
  lapply(FUN = spRast_values, na_rm=TRUE) %>%
  lapply(FUN = as.data.frame) %>%
  lapply(FUN = table) %>%
  lapply(FUN=getPercent)

names(df_crop) <- trunk$DESCRIPCIO

df_crop <- df_crop[!duplicated(df_crop)]

# df_crop <- df_crop[!duplicated(df_crop)]
# 
cbind(df$`BOSQUE DE PINO`, df_crop$`BOSQUE DE PINO`)
# 
# cbind(df$`BOSQUE DE PINO-ENCINO`, df_crop$`BOSQUE DE PINO-ENCINO`)
# 
# cbind(df$`BOSQUE DE AYARÍN`, df_crop$`BOSQUE DE AYARÍN`)
# 
# cbind(df$`PASTIZAL INDUCIDO`, df_crop$`PASTIZAL INDUCIDO`)
# 
# cbind(df$`AGRICULTURA DE TEMPORAL ANUAL`, df_crop$`AGRICULTURA DE TEMPORAL ANUAL`)
# 
# cbind(df$`VEGETACIÓN SECUNDARIA ARBUSTIVA DE BOSQUE DE PINO`,
#       df_crop$`VEGETACIÓN SECUNDARIA ARBUSTIVA DE BOSQUE DE PINO`)
# 
cbind(df$`VEGETACIÓN SECUNDARIA ARBÓREA DE BOSQUE DE PINO`,
      df_crop$`VEGETACIÓN SECUNDARIA ARBÓREA DE BOSQUE DE PINO`)

# --- terra:exact

df_exact <- trunk %>%
  apply(MARGIN = 1,
        FUN = function(x) terra::extract(mohinora_cps_TYPE,
                                         dplyr::filter(mohinora_USV,
                                                       DESCRIPCIO == x),
                                         exact = TRUE)) #%>%

names(df_exact) <- trunk$DESCRIPCIO

test <- df_exact$`VEGETACIÓN SECUNDARIA ARBÓREA DE BOSQUE DE PINO`

test_na <- test %>%
  filter(!is.na(map_TYPE))

a <- sum(test_na$fraction[test_na$map_TYPE == 8])/sum(test_na$fraction)
b <- sum(test_na$fraction[test_na$map_TYPE == 5])/sum(test_na$fraction)
c <- sum(test_na$fraction[test_na$map_TYPE == 6])/sum(test_na$fraction)

a+b+c

# --- CORRECT version

getCorrectPercent <- function(x){
  # x <- df_exact$`VEGETACIÓN SECUNDARIA ARBÓREA DE BOSQUE DE PINO`
  x_na <- x %>%
    filter(!is.na(map_TYPE))
  
  types <- sort(unique(x_na$map_TYPE))
  
  y <- unlist(lapply(types, function(s) sum(x_na$fraction[x_na$map_TYPE == s]) ))
  
  cbind.data.frame(TYPE=types, PERCENT=round(y/sum(x_na$fraction) * 100, digits=3) )
}

trunk <- mohinora_USV %>%
  dplyr::select( DESCRIPCIO ) %>%
  st_drop_geometry() 

df_exact <- trunk %>%
  apply(MARGIN = 1,
        FUN = function(x) terra::extract(mohinora_cps_TYPE,
                                         dplyr::filter(mohinora_USV,
                                                       DESCRIPCIO == x),
                                         exact = TRUE)) 

names(df_exact) <- trunk$DESCRIPCIO

df_exact <- df_exact[!duplicated(df_exact)]

df_exact_percent <- lapply(df_exact, function(s) getCorrectPercent(s)) 

# cbind(df$`VEGETACIÓN SECUNDARIA ARBÓREA DE BOSQUE DE PINO`,
#       df_exact_percent$`VEGETACIÓN SECUNDARIA ARBÓREA DE BOSQUE DE PINO`)

# ---

percent_type <- matrix(nrow=7, ncol=8)

for(i in 1:nrow(percent_type)){
  percent_type[i,as.numeric(df_exact_percent[[i]]$TYPE)] <- as.numeric(df_exact_percent[[i]]$PERCENT) 
}

df_type <- data.frame("Pino" = percent_type[1,],
                      "PinoEncino"= percent_type[2,],
                      "Ayarin" = percent_type[3,],
                      "Pastizal" = percent_type[4,],
                      "Agro" = percent_type[5,],
                      "ArbustivaPino" = percent_type[6,],
                      "ArboreaPino" = percent_type[7,])

row.names(df_type) <- 1:8

df_type[is.na(df_type)] <- 0

t(df_type) %>%
  kbl(digits=2, 
      caption = "Porcentaje de área clasificada por tipo de tendencias .") %>%
  kable_minimal(full_width = FALSE, html_font = "Cambria",
                font_size=20)

# --- Hermoseando la tabla anterior

usv_COLORS_sorted <- c("#A1E5A5", "#00A877", "#66B032",
                       "#E9D66B", "#83A4F0", "#FC8FAB", 
                       "#F500A1")

usv_NAMES_sorted <- c("Pino", "Pino-Encino", "Ayarin",
                      "Pastizales", "Agro", "Arbustiva", "Arborea")

t(df_type) %>%
  kbl(digits=2, booktabs = TRUE,
      caption = "Porcentaje de distintos tipos de
      tendencia por tipo de uso de suelo y vegetación usando 'terra::exact'") %>%
  kable_minimal(full_width = FALSE, font_size=10) %>%
  # kable_styling(latex_options = "striped", full_width = FALSE,
  #               font_size = 10) %>%
  column_spec(1, color = usv_COLORS_sorted) %>%
  row_spec(row=0, bold = TRUE)
# color = COLORES_update)

# ----

