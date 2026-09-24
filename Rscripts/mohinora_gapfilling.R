
# --- Elaborado Mar 8, 2024
# --- Actualizado Mar 6, 2025
# --- DATASET: NDVI MOD13Q1 v061 en Cerro Mohinora, Chihuahua, 2000-2023 
# --- Se hace uso del DATASET que se encuentra en el proyecto mestiper

# --- Este script interpola TEMPORALMENTE todos los datos faltantes
# --- en DATASET después de aplicar capa de calidad, es decir, después de aplicar
# --- el código del script mohinora_QA.R

# --- ACTUALIZACION: Sep. 23, 2026
# --- DATASET: NDVI y EVI MOD13Q!, Cerro Mohinora, 2000-2025

library(terra)
library(sf)

library(geoTS)
library(foreach)
library(doParallel)

library(imputeTS)
library(here)

source("Rscripts/auxFUN.R")

# ---

mestiperDIR <- paste0( getwd(), "/mestiper" )  # "C:/Users/inder/OneDrive/Desktop/mestiper" # USUARIO: modificar

# dataDIR <- list.dirs( path = here( "data" ) )
# 
# mohinora_NDVI_imputation <- 
#   # list.dirs(path = paste0( getwd(), "/TIF" ))[-1]
# 
# mohinora_NDVI_QA <-  #list.dirs(path = paste0( getwd(), "/data" ))[-1]

FILES_NDVI_imputation <- list.files( path = here( "data", "outputs", "mohinora_imputation" ),
                                     pattern = "NDVI",
                                     full.names = TRUE )

FILES_NDVI_QA <- list.files(path = here( "data", "mohinora_2026", "250m_16_days_NDVI_QA" ),
                            pattern = ".tif",
                            full.names = TRUE)

mohinora_DATA <- rast(FILES_NDVI_imputation, raw=TRUE, drivers = "GTiff")
add(mohinora_DATA) <- rast(FILES_NDVI_QA, raw=TRUE, drivers = "GTiff")

FILES_NDVI <- c(FILES_NDVI_imputation, FILES_NDVI_QA)

mohinora_DATA <- raster::stack(FILES_NDVI)

# mohinora_DATA_rTp <- spRast_valuesCoords(mohinora_DATA)

mohinora_DATA_rTp <- raster::rasterToPoints(mohinora_DATA)

mohinora_DATA_rTp_coords <- mohinora_DATA_rTp[,1:2]
mohinora_DATA_rTp_values <- mohinora_DATA_rTp[,3:551]

mohinora_interpol_linear <- matrix(nrow=nrow(mohinora_DATA_rTp$values), 
                                   ncol=ncol(mohinora_DATA_rTp$values))


pixel <- mohinora_DATA_rTp$values[295,]

pixel_ts <- ts(pixel, start = c(2000,1), end = c(2023,23), frequency = 23 )

plot(pixel_ts)


# mohinora_interpol_climatology <- matrix(nrow=nrow(mohinora_DATA_rTp$values), 
#                                         ncol=ncol(mohinora_DATA_rTp$values))

# --- TESTING code in parallel

numCores <- detectCores()

kluster <- parallel::makeCluster(numCores-1, outfile="")
registerDoParallel(kluster)

output <- foreach(i=c(5,3601), .combine="rbind",
                  .packages="imputeTS") %dopar% { # nrow(sp_ndvi_rTp)
                    
                    pixel <- mohinora_DATA_rTp$values[i,]
                    
                    out_linear <- pixel
                    out_climatology <- pixel
                    
                    if(length( is.na(pixel) ) > 0){
                      out_linear <- na_interpolation(pixel) 
                      out_climatology <- na_climatology(pixel)
                    }
                    
                    s <- c(as.numeric(out_linear),
                           as.numeric(out_climatology))
                    
                    # s <- c(as.numeric(out_linear))
                    
                    return(s)
                  }
stopCluster(kluster)

str(output)


# --- progress report file (to check out the process)

progressReportFile <- paste0(getwd(), "/RData/progressReports/mohinora/progress_temporal_gapfilling.txt" )
file.create(path=progressReportFile, showWarnings=FALSE)

write("===TEMPORAL GAPFILLING began at===",
      file=progressReportFile, append=TRUE)
write(as.character(Sys.time()[1]), file=progressReportFile,
      append=TRUE)

numCores <- detectCores()

kluster <- parallel::makeCluster(numCores-1, outfile="")
registerDoParallel(kluster)

output <- foreach(i=1:nrow(mohinora_DATA_rTp$values), .combine="rbind",
                  .packages="imputeTS") %dopar% { # nrow(sp_ndvi_rTp)
                    
                    pixel <- mohinora_DATA_rTp$values[i,]
                    
                    out_linear <- pixel
                    
                    if(length( is.na(pixel) ) > 0){
                      out_linear <- na_interpolation(pixel) 
                    }
                    
                    s <- c(as.numeric(out_linear))
                    
                    if(i %% 100 ==0){
                      texto <- paste0("Working on ROW: ", i)
                      write(texto, file=progressReportFile, append=TRUE)
                    }
                    
                    return(s)
                  }
stopCluster(kluster)

write("===TEMPORAL GAPFILLING ended at===",
      file=progressReportFile, append=TRUE)
write(as.character(Sys.time()[1]), file=progressReportFile,
      append=TRUE)

str(output)

mohinora_interpol_linear <- output
# mohinora_interpol_climatology[SAMPLE,] <- output[,550:(549*2)]

# --- rasterization

PROJECTION <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"

dirTIFS_toGet_names <- paste0(mestiperDIR, "/data/mohinora/250m_16_days_NDVI")
listTIFnames <- list.files(path = dirTIFS_toGet_names,
                           pattern = ".tif$",
                           full.names = TRUE)

vectorNAMES <- character(549)
for(i in 1:549){
  temp <- listTIFnames[i]
  aux <- strsplit(temp, "/")
  basename <- aux[[1]][ length(aux[[1]]) ]
  nameBASE <- strsplit(basename, ".tif", fixed=TRUE)[[1]]
  vectorNAMES[i] <- paste0(nameBASE, 
                           "_interpol")
}


# --- asegurarse de crear /TIF/mohinora_interpolation
# --- las primeras 3 columnas contienen valores de NDVI imputados
# --- a través del procedimiento de climatología, por tanto, no es necesario
# --- guardar esas capas nuevamente
for(i in 4:ncol(mohinora_interpol_linear)){
  
  if( i %% 100 == 0){
    cat("Working on layer ", i, "\n")
  }
  
  mat <- cbind(mohinora_DATA_rTp$coords[,1:2], 
               mohinora_interpol_linear[,i])
  
  layer <- matrixToRaster(matrix=mat, 
                          projection=PROJECTION) 
  
  raster::writeRaster(x=layer,
                      filename = paste0(getwd(),
                                        "/TIF/mohinora_interpolation/",
                                        vectorNAMES[i-3]),
                      format="GTiff",
                      datatype="INT2S",
                      overwrite=TRUE)
  
}

TIFilescheck <- list.files(path = paste0(getwd(), "/TIF/mohinora_interpolation"),
                           pattern = ".tif$",
                           full.names = TRUE)

rTest <- rast(TIFilescheck)
plot(rTest)


# --- re run before mohinora_anomalies.R
# --- guardando 552 capas en un solo archivo

imputeTIFs <- list.files(path = paste0(getwd(), "/TIF/mohinora_imputation"),
                         pattern = ".tif$",
                         full.names = TRUE)

TEMP <- rast(imputeTIFs)

mohinora_DATA_interpol <- TEMP

AUX <- rast(TIFilescheck)

add(mohinora_DATA_interpol) <- AUX

mohinora_DATA_interpol

writeRaster(x=mohinora_DATA_interpol,
            filename = paste0(getwd(), "/TIF/MOD13Q1_061_250m_16_days_NDVI_interpol.tif"),
            datatype="INT2S", overwrite=TRUE)

# --- verificando q todo está OK


maskTIF <- list.files(path=paste0(getwd(), "/data/mohinora/250m_16_days_NDVI_QA"),
                      pattern = ".tif$",
                      full.names = TRUE)

interpolTIFS <- list.files(path=paste0(getwd(), "/TIF/mohinora_interpolation"),
                           pattern = ".tif$",
                           full.names = TRUE)


a <- rast(maskTIF) # mask

b <- rast(interpolTIFS) # interpol

x <- 450

par(mfrow=c(1,2))
plot(subset(a,x), main="Sin interpolación")
# lines( mohinora_SHP_st, lwd=4)
plot(subset(b,x), main="Con interpolación")
# lines( mohinora_SHP_st, lwd=4)


# mohinoraSHP <- paste0( getwd(), "/RData" )
# 
# RDatafiles <- list.files(path = mohinoraSHP,
#                          pattern = ".RData",
#                          full.names = TRUE)
# 
# mohinora_SHP <- LoadToEnvironment(RDatafiles[1])$mohinora_SHP_st
# 
# mohinoraDIR <- list.dirs(path=paste0( getwd(), "/mestiper/data/mohinora" ),
#                          full.names = TRUE)
# 
# mohinora_SHP_st <- st_transform(x=mohinora_SHP, crs=crs(mohinora_DATA_interpol))



