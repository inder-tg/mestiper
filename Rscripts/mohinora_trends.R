
# --- Escrito por Inder Tecuapetla, Marzo 23, 2023
# --- Análisis de tendencias sobre el datacube de
# --- DATASET: NDVI MOD13Q1 en Cerro Mohinora, Chihuahua

# --- ACTUALIZADO Abril 5, 2024

# --- En este script presentamos un análisis de tendencias (globales).
# --- Usamos código en paralelo para eficientar el cómputo 

# --- DATASET: NDVI MOD13Q1 v061 en Cerro Mohinora, Chihuahua, 2000-2023 

# --- ADDicionalmente, este script requiere archivos
# --- MOD13Q1_061_250m_16_days_NDVI_interpol.tif
# --- creado con el archivo mohinora_temporal_gapfilling.R

# library(fpp)
# library(forecast)
# library(Kendall)
# # library(trend)
# library(RColorBrewer)
# library(gtools)
# # library(foreach)
# # library(doParallel)
# library(raster)
# # library(geoTS)
# # library(tmap)

library(terra)
library(trend)
library(foreach)
library(doParallel)
library(geoTS)
library(sf)
library(tmap)

source( paste0( getwd(), "/Rscripts/auxFUN.R" ) )

# --- DATA Loading

listFILES_mohinora <- list.files(path=paste0(getwd(), "/TIF"), 
                                 pattern=".tif$", 
                                 full.names=TRUE)

mohinora_DATA_interpol <- rast(listFILES_mohinora[1])

mohinora_DATA_interpol_rTp <- spRast_valueCoords(mohinora_DATA_interpol)

mohinoraSHP <- paste0( getwd(), "/RData" )

RDatafiles <- list.files(path = mohinoraSHP,
                         pattern = ".RData",
                         full.names = TRUE)

mohinora_SHP <- LoadToEnvironment(RDatafiles[1])$mohinora_SHP_st

mohinora_SHP_st <- st_transform(x=mohinora_SHP, crs=crs(mohinora_DATA_interpol))


# -----------------------------------------------
# --- Análisis de tendencias --- #
# -----------------------------------------------


plot(subset(mohinora_DATA_interpol,453))
lines( mohinora_SHP_st, lwd=4 )

XY <- locator() 

# XY <- list(x=-10698785, y=2893289)

xy <- get_timeSeries_byClicking(c(XY$x, XY$y), 
                                df=mohinora_DATA_interpol_rTp$coords)
# 491, buen pixel para cps; 1975 buen pixel para trends
pixel <- mohinora_DATA_interpol_rTp$values[xy$coord, ]

# --- objeto ts
pixel_ts <- ts(pixel * 1e-4, start = c(2000,1), 
               end = c(2023,23),
               frequency = 23)

# --- aprox a tendencia via lowess
plot(pixel_ts)
lines(lowess(time(pixel_ts), pixel_ts), lwd = 5, col="purple")

# --- let's go and remove seasonality
pixel_mat <- get_pixel_matrix(pixel)

mu <- apply(pixel_mat, 2, mean)
mu_ts <- ts(rep(mu,24) * 1e-4, 
            start = c(2000,1), 
            end = c(2023,23), 
            frequency = 23)

plot(pixel_ts, lwd=2)
lines(mu_ts, col="red", lwd=3)

# --- pixel sin estacionalidad como obj 'ts'
deSeasonal_pixel_ts <- pixel_ts - mu_ts

plot(deSeasonal_pixel_ts)
lines(lowess(time(deSeasonal_pixel_ts), deSeasonal_pixel_ts), 
      lwd = 5, col="purple")

# --- getting p-value and slope

# --- pixel sin estacionalidad
pixel_deSeason <- pixel - rep(mu,24)

pixel_MannKendall <- mk.test(pixel_deSeason)

pixel_MannKendall$p.value

pixel_SenTheil <- sens.slope(pixel_deSeason)

pixel_SenTheil$estimates

b_hat <- as.numeric(pixel_SenTheil$estimates)
a_hat <- median( pixel_deSeason - b_hat * 1:length(pixel) )


pixel_deSeason_ts <- ts(pixel_deSeason, 
                        start = c(2000, 1),
                        end = c(2023, 23), 
                        frequency = 23)

lineaTheilSen <- ts(a_hat + b_hat * 1:length(pixel), 
                    start = c(2000, 1),
                    end = c(2023, 23), 
                    frequency = 23)

par(mfrow=c(1,1), mar = c(2,2,1,2), adj=0)
plot(pixel_deSeason_ts, type="l", col = "gray", ylab = "NDVI")
lines(lineaTheilSen, lwd = 5, col = "lightcoral")
legend("topright", legend = c("raw data", "linear trend"),
       col = c("gray", "lightcoral"), lty = rep(1,2), lwd = c(1,5), 
       bty = "n")

# --- CODIGO EN PARALELO

# --- en df_pvalue vamos a guardar los p-values de la prueba MK
# --- en df_slope vamos a guardar los estimadores de pendiente prueba TS

df_pvalue <- matrix(nrow=nrow(mohinora_DATA_interpol_rTp$coords), ncol=3)
df_pvalue[,1:2] <- mohinora_DATA_interpol_rTp$coords[,1:2]

df_slope <- matrix(nrow=nrow(mohinora_DATA_interpol_rTp$coords), ncol=3)
df_slope[,1:2] <- mohinora_DATA_interpol_rTp$coords[,1:2]

# ---

numCores <- detectCores()

progressReportFile <- paste0(getwd(), "/RData/progressReports/mohinora/progress_trends.txt" )
file.create(path=progressReportFile, showWarnings=FALSE)

write("===TREND analysis began at===",
      file=progressReportFile, append=TRUE)
write(as.character(Sys.time()[1]), file=progressReportFile,
      append=TRUE)

kluster <- parallel::makeCluster(numCores-1, outfile="")
registerDoParallel(kluster)

output <- foreach(i=1:nrow(mohinora_DATA_interpol_rTp$coords), .combine="rbind",
                  .packages=c("trend") ) %dopar% { 
                    
                    pixel <- mohinora_DATA_interpol_rTp$values[i,]
                    
                    pixel_mat <- get_pixel_matrix(pixel)
                    
                    mu <- apply(pixel_mat, 2, mean, na.rm=TRUE)
                    
                    pixel_deSeason <- pixel - rep(mu,24)
                    
                    pixel_MannKendall <- mk.test(pixel_deSeason)$p.value
                    
                    pixel_SenTheil <- sens.slope(pixel_deSeason)$estimates
                    
                    s <- c(as.numeric(pixel_MannKendall), as.numeric(pixel_SenTheil))
                    
                    if(i %% 100 ==0){
                      texto <- paste0("Working on ROW: ", i)
                      write(texto, file=progressReportFile, append=TRUE)
                    }
                    
                    return(s)
                  }
stopCluster(kluster)

write( as.character(Sys.time()[1]), file=progressReportFile, append=TRUE)
write( "===TREND analysis ended here===", file=progressReportFile, append=TRUE)

# ---

df_pvalue[,3] <- output[,1]
df_slope[,3] <- output[,2]

# --- Asegurarse de crear /mohinora_trends

save(df_pvalue, file=paste0(getwd(),"/RData/mohinora_trends/pvalue.RData"))
save(df_slope, file=paste0(getwd(),"/RData/mohinora_trends/slope.RData"))

# --- Si en tu sesion no se encuentran los objetos df_pvalue,
# y df_slope, entonces descomenta y ejecuta las siguientes 3 líneas

# df_pvalue <- LoadToEnvironment(paste0(getwd(),"/RData/mohinora_trends/pvalue.RData"))$df_pvalue
# df_slope <- LoadToEnvironment(paste0(getwd(),"/RData/mohinora_trends/slope.RData"))$df_slope

# -----------------------
# --- RASTERIZACION --- #
# -----------------------

PROJECTION <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"

map_pvalue <- matrixToRaster(matrix=df_pvalue, projection=PROJECTION)
map_slope <- matrixToRaster(matrix=df_slope, projection=PROJECTION)

writeRaster(map_pvalue,
            filename = paste0( getwd(), "/TIF/mohinora_trends/pvalueMap"),
            format="GTiff", datatype="FLT4S", overwrite=TRUE)

writeRaster(map_slope,
            filename = paste0( getwd(), "/TIF/mohinora_trends/slopeMap"),
            format="GTiff", datatype="FLT4S", overwrite=TRUE)

# ---

pvalueMap <- map_pvalue
pvalueMap[pvalueMap > 0.05] <- NA

slopeMap <- map_slope
slopeMap[is.na(pvalueMap)] <- NA

writeRaster(slopeMap,
            filename = paste0( getwd(), "/TIF/mohinora_trends/significantSlopeMap"),
            format="GTiff", datatype="FLT4S", overwrite=TRUE)


plot(slopeMap)

hist(slopeMap)
boxplot(slopeMap)

quantile(slopeMap)

QUANT <- quantile(slopeMap, probs= seq(0,1,by=0.1) )

plot(slopeMap, col=c("red", "yellow", "green"),
     breaks=c(QUANT[1], 0, QUANT[5], QUANT[11]))

# --- just the ANP Cerro Mohinora

# SHP_anp <- list.files( path = paste0( getwd(), "/data/anp_2021" ),
#                        pattern = ".shp", 
#                        full.names = TRUE )
# 
# shp_anp <- shapefile( SHP_anp[1] )
# 
# shp_anp_sinu <- spTransform(shp_anp, crs(r))
# 
# mohinora_trendAnalysis <- raster_intersect_sp(x=slopeMap, 
#                                               y=shp_anp_sinu[165,])
# 

mapSlope <- rast(paste0( getwd(), "/TIF/mohinora_trends/significantSlopeMap.tif" ))

mohinora_trends <- crop(x=mapSlope, 
                        y=mohinora_SHP_st,
                        mask=TRUE)

# plot(mohinora_trendAnalysis, col=c("red", "yellow", "green"),
#      breaks=c(QUANT[1], 0, QUANT[2], QUANT[11]))
# 
# # --- A (kinda) pro map
# 
nz_height_map = tm_shape(mohinora_trends, 
                         bbox = mohinora_SHP_st) +
  tm_raster(style = "cont", palette = "YlGn", 
            legend.show = TRUE,
            title="Linear trend slope") +
  tm_shape(mohinora_SHP_st) + tm_borders() +
  tm_layout(title = "Cerro Mohinora", frame = FALSE, bg.color = NA) +
  tm_compass(type = "8star", position = c("right", "bottom")) +
  tm_scale_bar(text.size = 0.65,
               position = c("right", "bottom"))
# 
nz_height_map


