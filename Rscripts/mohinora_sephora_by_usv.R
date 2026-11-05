
# --- Escrito por Inder Tecuapetla, Nov 5, 2026
# --- SEPHORA para Mohinora
# --- DATASET: NDVI MOD13Q1 en Cerro Mohinora, Chihuahua

# --- En este script aplicamos phenopar_polygon para
# ---- estimar parámetros fenológicos sobre cada cobertura de suelo
# --- en Cerro Mohinora

# --- DATASET: NDVI MOD13Q1 v061 en Cerro Mohinora, Chihuahua, 2000-2025 

library(terra)
library(sf)
library(tmap)
library(raster)
library(here)
library(dplyr)

library(sephora)

LoadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env)
}

spRast_valuesCoords <- function(spRaster, na_rm=FALSE){
  spPoints <- as.points(spRaster, na.rm=na_rm)
  
  spValues <- terra::extract(spRaster, spPoints)
  
  DIM <- dim(spValues)
  
  spRasterToPoints <- as.matrix(spValues[1:DIM[1],2:DIM[2]])
  
  spCoords <- crds(spRaster, na.rm=na_rm)
  
  list(values=spRasterToPoints, coords=spCoords)  
}

# --- DATA Loading

mohinora_NDVI_DIRS <- list.dirs(path = here("data", "mohinora_2026") )

FILES_NDVI_imputation <- list.files(path = here("data", "outputs", "mohinora_imputation"),
                                    pattern = "NDVI",
                                    full.names = TRUE)

FILES_NDVI_interpolation <- list.files(path = mohinora_NDVI_DIRS[5],
                                       pattern = ".tif",
                                       full.names = TRUE)

mohinora_NDVI_files <- c(FILES_NDVI_imputation, FILES_NDVI_interpolation)

mohinora_DATA <- rast(mohinora_NDVI_files)

shpFiles <- list.files(path = here("data", "mohinora_usv7"),
                       pattern = ".shp",
                       full.names = TRUE)

mohinora_USV <- read_sf(shpFiles[1])

# ---

NAMES <- unique( mohinora_USV$DESCRIPCIO )

for( nombre in NAMES ){
  dir.create(here( "RData", "mohinora_sephora", nombre ))
  
  mohinora_DATA_shp <- crop( mohinora_DATA, mohinora_USV, mask=TRUE )
  
  mohinora_DATA_shp_rTp <- spRast_valuesCoords(mohinora_DATA_shp)
  
  save( mohinora_DATA_shp_rTp, 
        file = here( "RData", "mohinora_sephora", nombre, "DATA_rTp.RData" )  )
}

# ---

keepDataIndices <- function(data) {
  valid <- rowSums(!is.na(data)) >= 2
  list(
    indices = which(valid),
    total_rows = nrow(data),
    kept_rows = sum(valid),
    dropped_rows = sum(!valid)
  )
}

DIRS <- list.dirs( here( "RData", "mohinora_sephora" ) )[-1]

for( dir in DIRS ){
  
  mohinora_DATA_shp_rTp <- LoadToEnvironment( paste0( dir, "DATA_rTp.RData" ) )$mohinora_DATA_shp_rTp 
  
  idx <- keepDataIndices(data=mohinora_DATA_rTp$values)
  
  rTp <- list()
  rTp$values <- mohinora_DATA_rTp$values[idx$indices,]
  rTp$coords <- mohinora_DATA_rTp$coords[idx$indices,]
  
  cat("SAVING usv", basename(dir), "\n")
  
  save(rTp, file = file.path(dir, "rTp.RData"))
  
}

BASIS <- drbasis(n=75, q=2)

for( dir in DIRS ){
  
  rTp <- LoadToEnvironment( paste0( dir, "rTp.RData" ) )$rTp 
  
  mohinora_dtw2 <- phenopar_polygon(product = "independent",
                                    data = rTp$values,
                                    frequency = 23,
                                    clusterSize = 0.65,
                                    save_all = TRUE,
                                    distance = "dtw2",
                                    basis = BASIS,
                                    numCores = 20,
                                    k=3,
                                    confiBand = TRUE,
                                    dirToSave = dir,
                                    outputFileBaseName = "mohinora_dtw2_sample_75_clusterSize_0p65")
  
}
