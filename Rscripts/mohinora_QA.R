
# --- Elaborado Mar 9, 2023
# --- Código para calcular % de dato faltante -a nivel pixel-
# --- y maxGapLength
# --- DATASET: NDVI MOD13Q1 en Cerro Mohinora, Chihuahua, 2000-2023

# --- Actualizado: Feb 24, 2024, Abril 5. 2025
# --- DATASET: NDVI MOD13Q1 v061 en Cerro Mohinora, Chihuahua, 2000-2024

# --- Actualizado: Ago 6, 2026
# --- Se han actualizado algunas líneas de código empleando funciones actuales 
# --- para alinearnos a las expectativas del Diplomado en Geomática Edición XIX
# --- NOTA: Crear subdirectorios necesarios (ej. /data/mohinora)

# --- Actualizado: Sep 10, 2026
# --- Se generó código más eficiente
# --- Ahora consideramos tmb el producto EVI del MOD13Q1

library(terra)
library(sf)
library(geoTS)
library(foreach)
library(doParallel)
library(here)

source("Rscripts/auxFUN.R")

# ---

# DIRS <-  list.dirs( path = paste0( getwd(), "/data"  ) )
DIRS <- list.dirs( here("data", "mohinora_2026") )

lst_DIRS <- setNames( as.list(DIRS), basename(DIRS) )

NDVIfiles <- list.files(path = lst_DIRS$`250m_16_days_NDVI`,
                        pattern = ".tif", 
                        full.names = TRUE)

EVIfiles <- list.files( path = lst_DIRS$`250m_16_days_EVI`,
                        pattern = ".tif", 
                        full.names = TRUE )

mohinora_NDVI_DATA <- rast(NDVIfiles)

mohinora_EVI_DATA <- rast(EVIfiles)

RELIABILITYfiles <- list.files(path = lst_DIRS$`250m_16_days_pixel_reliability`,
                               pattern = ".tif", 
                               full.names = TRUE)

mohinora_DATA_reliability <- rast(RELIABILITYfiles)

# mohinoraRDataDIR <- paste0( mestiperDIR, "/RData" )
DIR_outputs <- here( "data", "outputs" )
SHPfiles <- list.files(path = DIR_outputs,
                       pattern = ".shp$",
                       full.names = TRUE)

mohinora_shp <- read_sf(SHPfiles[1])

# --- Ejemplos

TEMP_NDVI <- subset(mohinora_NDVI_DATA, 159)
TEMP_EVI <- subset(mohinora_EVI_DATA, 159)
AUX <- subset(mohinora_DATA_reliability, 159)

par(mfrow=c(1,2))
plot(TEMP_NDVI)
lines(mohinora_shp)

plot(AUX)
lines(mohinora_shp, col = "cyan")

plot(TEMP_EVI)
lines(mohinora_shp)

plot(AUX)
lines(mohinora_shp, col = "cyan")

TEMP_NDVI[ AUX >= 2 ] <- NA
TEMP_EVI[ AUX >= 2 ] <- NA

plot(TEMP_NDVI)
lines(mohinora_shp)

plot(AUX)
lines(mohinora_shp, col = "cyan")

plot(TEMP_EVI)
lines(mohinora_shp)

plot(AUX)
lines(mohinora_shp, col = "cyan")

# ---

# whereToSave <- paste0(DIRS[3], "/250m_16_days_NDVI_QA") 
whereToSaveNDVI <- here( lst_DIRS$mohinora_2026, "/250m_16_days_NDVI_QA" )
whereToSaveEVI <- here( lst_DIRS$mohinora_2026, "/250m_16_days_EVI_QA" )
dir.create(whereToSaveNDVI, recursive = TRUE)
dir.create(whereToSaveEVI, recursive = TRUE)

TEMP_NDVI <- subset(mohinora_NDVI_DATA, 1)
AUX <- subset(mohinora_DATA_reliability, 1)
TEMP_NDVI[ AUX >= 2 ] <- NA

nameFILE <- basename( NDVIfiles[1]  )
nameFILE <- paste0(strsplit( nameFILE, ".tif" )[[1]][1], "_QA.tif")
writeRaster(TEMP_NDVI, 
            filename = here(whereToSaveNDVI, nameFILE), # paste0( whereToSave, "/", nameFILE ),
            datatype = datatype(mohinora_NDVI_DATA)[1],
            overwrite = TRUE)


TEMP_EVI <- subset(mohinora_EVI_DATA, 1)
# AUX <- subset(mohinora_DATA_reliability, 1)
TEMP_EVI[ AUX >= 2 ] <- NA

nameFILE <- basename( EVIfiles[1]  )
nameFILE <- paste0(strsplit( nameFILE, ".tif" )[[1]][1], "_QA.tif")
writeRaster(TEMP_EVI, 
            filename = here(whereToSaveEVI, nameFILE), # paste0( whereToSave, "/", nameFILE ),
            datatype = datatype(mohinora_EVI_DATA)[1],
            overwrite = TRUE)

# --- for-loop: Resolviendo una tarea vía iteración
for(i in 2:nlyr(mohinora_NDVI_DATA)){
  
  if( i %% 50 == 0 ){
    cat("Working on layer: ", i, "\n")
  }
  
  TEMP_NDVI <- subset(mohinora_NDVI_DATA, i)
  TEMP_EVI <- subset(mohinora_EVI_DATA, i)
  AUX <- subset(mohinora_DATA_reliability, i)
  TEMP_NDVI[ AUX >= 2 ] <- NA 
  TEMP_EVI[ AUX >= 2 ] <- NA 
  
  nameFILE <- basename( NDVIfiles[i]  )
  nameFILE <- paste0(strsplit( nameFILE, ".tif" )[[1]][1], "_QA.tif")
  writeRaster(TEMP_NDVI, 
              filename = here(whereToSaveNDVI, nameFILE), # paste0( whereToSave, "/", nameFILE ),
              datatype = datatype(mohinora_NDVI_DATA)[1],
              overwrite = TRUE)
  
  nameFILE <- basename( EVIfiles[i]  )
  nameFILE <- paste0(strsplit( nameFILE, ".tif" )[[1]][1], "_QA.tif")
  writeRaster(TEMP_EVI, 
              filename = here(whereToSaveEVI, nameFILE), # paste0( whereToSave, "/", nameFILE ),
              datatype = datatype(mohinora_EVI_DATA)[1],
              overwrite = TRUE)
  
  if( i %% 600 == 0 ){
    cat("Terminó con éxito: ", Sys.time(), "\n")
  }
  
}

# --- 

ndviQAFILES <- list.files( path = here( lst_DIRS$mohinora_2026, "250m_16_days_NDVI_QA" ), # paste0( getwd(), "/data/mohinora/250m_16_days_NDVI_QA" ),
                           pattern = ".tif$",
                           full.names = TRUE )

mohinora_NDVI_QA <- rast(ndviQAFILES)

mohinora_NDVI_QA_shp <- crop(mohinora_NDVI_QA, mohinora_shp,
                             mask=TRUE)

eviQAFILES <- list.files( path = here( lst_DIRS$mohinora_2026, "250m_16_days_EVI_QA" ), # paste0( getwd(), "/data/mohinora/250m_16_days_NDVI_QA" ),
                           pattern = ".tif$",
                           full.names = TRUE )

mohinora_EVI_QA <- rast(eviQAFILES)

mohinora_EVI_QA_shp <- crop(mohinora_EVI_QA, mohinora_shp,
                             mask=TRUE)

# --- Los numeritos!!!
mohinora_NDVI_QA_rTp <- spRast_valuesCoords(mohinora_NDVI_QA_shp)
mohinora_EVI_QA_rTp <- spRast_valuesCoords(mohinora_EVI_QA_shp)

# --- EJEMPLO sobre un pixel

# --- Eg. con NDVI --- REVISAR end = c(2026, X)
pixel <- mohinora_NDVI_QA_rTp$values[1700,]

(pixel_percentMiss <- sum(is.na(pixel)) / length(pixel)) * 100 # length de cualquier pixel es 483

(pixel_maxgap <- maxLagMissVal(x=pixel)$maxLag)

pixel_ts <- ts( pixel[1:(20+(25)*23)], start = c(2000,1), end = c(2025,23), frequency = 23 )

par(mfrow=c(1,1))
plot(pixel_ts, ylab="NDVI (integer format)")

# --- Eg. con EVI
pixel <- mohinora_EVI_QA_rTp$values[1700,]

(pixel_percentMiss <- sum(is.na(pixel)) / length(pixel)) * 100 # length de cualquier pixel es 483

(pixel_maxgap <- maxLagMissVal(x=pixel)$maxLag)

pixel_ts <- ts( pixel[1:(20+(25)*23)], start = c(2000,1), end = c(2024,23), frequency = 23 )

par(mfrow=c(1,1))
plot(pixel_ts, ylab="EVI (integer format)")

# --- COMPUTO en PARALELO

df_perc_miss <- matrix(nrow=nrow(mohinora_NDVI_QA_rTp$values), ncol=3)
df_perc_miss[,1:2] <- mohinora_NDVI_QA_rTp$coords[,1:2]

maxgap_df <- matrix(nrow=nrow(mohinora_NDVI_QA_rTp$values), ncol=3)
maxgap_df[,1:2] <- mohinora_NDVI_QA_rTp$coords[,1:2]

# --- progress report file (to check out on the process)
# --- antes de ejecutar, crear /RData/progressReports/mohinora (sólo en caso de que los directorios no existan)

dir.create( here( "RData", "progressReports" ),
            recursive = TRUE )

numCores <- detectCores()

progressReportFile <- here( "RData", "progressReports", "mohinora_QA.txt" )
file.create(path=progressReportFile, showWarnings=FALSE)

write("===QA analysis began at===",
      file=progressReportFile, append=TRUE)
write(as.character(Sys.time()[1]), file=progressReportFile,
      append=TRUE)

kluster <- parallel::makeCluster(numCores-1, outfile="")
registerDoParallel(kluster)

output <- foreach(i=1:nrow(mohinora_NDVI_QA_rTp$values), .combine="rbind",
                  .packages="geoTS") %dopar% { 
                    
                    pixel <- mohinora_NDVI_QA_rTp$values[i,1:(20+(25)*23)]
                    
                    pixel_percentMiss <- sum(is.na(pixel)) / length(pixel) * 100 # length de cualquier pixel es 483
                    
                    pixel_maxgap <- maxLagMissVal(x=pixel)$maxLag
                    
                    s <- c(as.numeric(pixel_percentMiss), as.numeric(pixel_maxgap))
                    
                    if(i %% 100 ==0){
                      texto <- paste0("Working on ROW: ", i)
                      write(texto, file=progressReportFile, append=TRUE)
                    }
                    
                    return(s)
                  }
stopCluster(kluster)

write( as.character(Sys.time()[1]), file=progressReportFile, append=TRUE )
write( "===QA analysis ended here===", file=progressReportFile, append=TRUE )

# ---

# --- Saving output

df_perc_miss[,3] <- output[,1]
maxgap_df[,3] <- output[,2]

# --- antes de ejecutar, crear /RData/mohinora_QA

dir.create( here( "RData", "mohinora_QA" ), # paste0(getwd(),"/RData/mohinora_QA"), 
            recursive = TRUE )

save(df_perc_miss, file=paste0(getwd(),"/RData/mohinora_QA/percent_missingValue.RData"))
save(maxgap_df, file=paste0(getwd(),"/RData/mohinora_QA/maxgap.RData"))

# --- Rasterization

PROJECTION <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs" #projection(mohinora_DATA) # crs(mohinora_mask) # raster::projection(TEMP) # crs(STACK_sp_ndvi_subset)

map_percentMissing <- matrixToRaster(matrix=df_perc_miss, 
                                     projection=PROJECTION)
map_maxgap <- matrixToRaster(matrix=maxgap_df, 
                             projection=PROJECTION)

map_percentMissing
map_maxgap

# --- antes de ejecutar, crear /outputs/mohinora_QA

dir.create( paste0( getwd(), "/data/outputs/mohinora_QA" ),
            recursive = TRUE )

writeRaster( map_percentMissing,
             filename = paste0( getwd(), "/data/outputs/mohinora_QA/missingValue" ),
             format="GTiff", datatype="FLT4S", overwrite=TRUE )

writeRaster( map_maxgap,
             filename = paste0( getwd(), "/data/outputs/mohinora_QA/maxGap" ),
             format="GTiff", datatype="INT2U", overwrite=TRUE )

# ---

QAfiles <- list.files(path = here( "data", "outputs", "mohinora_QA" ),
                      pattern=".tif",
                      full.names=TRUE)

maxGap <- rast(QAfiles[1])
percent <- rast(QAfiles[2])

par(mfrow=c(1,2))
plot(maxGap, main="max-gap length")
plot(percent, main="% missing values")

maxGap <- crop(maxGap, mohinora_shp, mask=TRUE)
percent <- crop(percent, mohinora_shp, mask=TRUE)

plot(maxGap, main="max-gap length")
plot(percent, main="% missing values")
par(mfrow=c(1,1))
# ---
