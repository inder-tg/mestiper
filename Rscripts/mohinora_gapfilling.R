
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
library(plotly)

source("Rscripts/auxFUN.R")

# --- DATA loading

FILES_NDVI_imputation <- list.files( path = here( "data", "outputs", "mohinora_imputation" ),
                                     pattern = "NDVI",
                                     full.names = TRUE )

FILES_NDVI_QA <- list.files(path = here( "data", "mohinora_2026", "250m_16_days_NDVI_QA" ),
                            pattern = ".tif",
                            full.names = TRUE)

FILES_NDVI <- c(FILES_NDVI_imputation, FILES_NDVI_QA)

mohinora_DATA <- rast( FILES_NDVI[1:598] ) # raster::stack(FILES_NDVI)


mohinora_DATA_rTp <- spRast_valuesCoords(mohinora_DATA) #raster::rasterToPoints(mohinora_DATA)

mohinora_interpol_linear <- matrix(nrow=nrow(mohinora_DATA_rTp$values), 
                                   ncol=ncol(mohinora_DATA_rTp$values))

# --- AN example

maskQA <- rast( here( "data", "outputs", "mohinora_QA", "missingValue.tif" ) )

DIR_outputs <- here( "data", "outputs" )
SHPfiles <- list.files(path = DIR_outputs,
                       pattern = ".shp$",
                       full.names = TRUE)

mohinora_shp <- read_sf(SHPfiles[1])

maskQA <- crop(maskQA, mohinora_shp, mask=TRUE)

plot( maskQA )

# --- 1. Ejecuta la línea de abajo
XY <- locator()
# --- 2. Haz click (SOLO UNA VEZ) en algún píxel en la imagen
# --- 3. Presiona la tecla ESC de tu teclado
# --- 4. Continúa con el script a partir de la línea 107

xy <- get_timeSeries_byClicking(c(XY$x, XY$y),
                                df=mohinora_DATA_rTp$coords)

pixel <- mohinora_DATA_rTp$values[295,]

pixel <- mohinora_DATA_rTp$values[xy$coord,]

pixel_ts <- ts(pixel, start = c(2000,1), end = c(2025,23), frequency = 23 )

plot(pixel_ts)

out_linear <- na_interpolation(pixel)

out_linear_ts <- ts(out_linear, start = c(2000,1), end = c(2025,23), frequency = 23 )

plot(out_linear_ts)

# --- alternativa

# Convertir a data.frame
df <- data.frame(
  time = time(pixel_ts),
  series1 = as.numeric(pixel_ts),
  series2 = as.numeric(out_linear_ts)
)

# Plot interactivo con dos series de tiempo
plot_ly(df, x = ~time) %>%
  add_lines(y = ~series1, name = "NDVI original", 
            line = list(color = "darkgreen")) %>%
  add_lines(y = ~series2, name = "NDVI interpol", 
            line = list(color = "blue")) %>%
  layout(title = "Comparación de dos series de tiempo",
         xaxis = list(title = "Años"),
         yaxis = list(title = "NDVI"))

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

progressReportFile <- paste0(getwd(), "/RData/progressReports/mohinora_gapfilling.txt" )
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
                    
                    if( sum( !is.na(pixel) ) >= 2 ){
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

# --- RASTERIZATION

PROJECTION <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"

vectorNAMES <- character(598)
for(i in 1:598){
  # i = 598
  temp <- FILES_NDVI[i] #listTIFnames[i]
  aux <- strsplit(temp, "/")
  basename <- aux[[1]][ length(aux[[1]]) ]
  nameBASE <- strsplit(basename, ".tif", fixed=TRUE)[[1]]
  vectorNAMES[i] <- paste0(nameBASE, "_interpol")
}

# --- las primeras 3 columnas contienen valores de NDVI imputados
# --- a través del procedimiento de climatología, por tanto, no es necesario
# --- guardar esas capas nuevamente

dir.create( here( "data", "mohinora_2026", "250m_16_days_NDVI_interpol" ),
            recursive = TRUE )

for(i in 4:ncol(mohinora_interpol_linear)){
  
  if( i %% 100 == 0){
    cat("Working on layer ", i, "\n")
  }
  
  mat <- cbind(mohinora_DATA_rTp$coords[,1:2], 
               mohinora_interpol_linear[,i])
  
  layer <- matrixToRaster(matrix=mat, 
                          projection=PROJECTION) 
  
  raster::writeRaster(x=layer,
                      filename = here( "data", "mohinora_2026", "250m_16_days_NDVI_interpol",
                                       vectorNAMES[i] ),
                      format="GTiff",
                      datatype="INT2S",
                      overwrite=TRUE)
  
}


# --- verificando q todo está OK


ndviQA <- list.files(path = here( "data", "mohinora_2026", "250m_16_days_NDVI_QA" ),
                     pattern = ".tif$",
                     full.names = TRUE)

ndviINTERPOL <- list.files(path = here( "data", "mohinora_2026", "250m_16_days_NDVI_interpol" ),
                           pattern = ".tif$",
                           full.names = TRUE)


a <- rast(ndviQA[1:595]) # mask

b <- rast(ndviINTERPOL) # interpol

x <- 450

par(mfrow=c(1,2))
plot(subset(a,x), main="Sin interpolación")
# lines( mohinora_SHP_st, lwd=4)
plot(subset(b,x), main="Con interpolación")
# lines( mohinora_SHP_st, lwd=4)



