
# --- Escrito por Inder Tecuapetla, Nov 5, 2026
# --- SEPHORA para Mohinora
# --- DATASET: NDVI MOD13Q1 en Cerro Mohinora, Chihuahua

# --- En este script calculamos estadísticas sobre los parámetros fenológicos 
# --- estimados en mohinora_sephora_by_usv.R en Cerro Mohinora

# --- DATASET: NDVI MOD13Q1 v061 en Cerro Mohinora, Chihuahua, 2000-2025 

library(terra)
library(sf)
library(tmap)
library(raster)
library(here)
library(dplyr)
library(kableExtra)

library(sephora)

library(ggnewscale)
library(spiralize)
library(vcd)

source("Rscripts/getDist_phenoParam.R")
source("Rscripts/getSpiralPlot.R")
source("Rscripts/extract_sephora.R")
source("Rscripts/sephora_auxFUN.R")

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

DIRS <- list.dirs( here( "RData", "mohinora_sephora" ) )[-1]

list_DIRS <- setNames( as.list(DIRS), basename(DIRS) )

RDataFILES <- list.files( path = list_DIRS$`BOSQUE DE PINO-ENCINO`, pattern = ".RData",
                          full.names = TRUE )

sephora_phenoparams <- LoadToEnvironment( RDataFILES[2] )$output

str( sephora_phenoparams[[1]] )

# RDataFILES <- list.files( path = list_DIRS$`BOSQUE DE PINO`, pattern = ".RData",
#                           full.names = TRUE )
# 
# sephora_phenoparams <- LoadToEnvironment( RDataFILES[2] )$output
# 
# str( sephora_phenoparams[[1]] )

# --- SPIRAL plots

# --- colors used in spiralPlot below
cgu <- rgb(173/255,221/255,142/255)
csos <- rgb(120/255,198/255,121/255)
cmat <- rgb(49/255, 163/255,84/255)
csen <- rgb(217/255, 95/255, 14/255)
ceos <- rgb(254/255, 153/255, 41/255)
cdor <- rgb(208/255, 209/255, 230/255)

colores <- c(cgu,csos,cmat,csen,ceos,cdor)

# --- how to get a SpiralPlot
mohinora_params <- extract_sephora(x = sephora_phenoparams, object = "phenoparams")

getSpiralPlot(MAT=mohinora_params, 
              LABELS=month.name,
              vp_param=list(width=0.5, height=0.7))
vcd::grid_legend(x=1.215, y=0.125, pch=18, col=colores,
                 frame=FALSE,
                 labels=c("GU","SoS","Mat","Sen","EoS","Dor"),
                 title="Params")

# --- sephora STATS


a <- phenoDatesToDoY(phenoDates=unlist(sephora_basicStats(mohinora_params))[c(1,3,5,7,9,11)])
# b <- phenoDatesToDoY(phenoDates=unlist(sephora_basicStats(ayarin))[c(1,3,5,7,9,11)])
# c <- phenoDatesToDoY(phenoDates=unlist(sephora_basicStats(pinoEncino))[c(1,3,5,7,9,11)])
# d <- phenoDatesToDoY(phenoDates=unlist(sephora_basicStats(agro))[c(1,3,5,7,9,11)])

var_a <- redondea(unlist(sephora_basicStats(mohinora_params))[c(1,3,5,7,9,11)+1])
# var_b <- redondea(unlist(sephora_basicStats(ayarin))[c(1,3,5,7,9,11)+1])
# var_c <- redondea(unlist(sephora_basicStats(pinoEncino))[c(1,3,5,7,9,11)+1])
# var_d <- redondea(unlist(sephora_basicStats(agro))[c(1,3,5,7,9,11)+1])

df <- data.frame(#Date_Grasslands = paste( a$month, a$day, "(+", var_a, "days)"),
                 # Date_Ayarin = paste( b$month, b$day, "(+", var_b, "days)"),
                 Date_PineOak = paste( a$month, a$day, "(+", var_a, "days)") )
                 #Date_Agro = paste( d$month, d$day, "(+", var_d, "days)"))

# colnames(df) <- c("Grasslands", "Douglas Fir", "Pine-Oak", "Agriculture")
colnames(df) <- c("Pine-Oak")
row.names(df) <- c("GU", "SoS", "Mat", "Sen", "EoS", "Dor")

kbl(df,
    longtable = T, booktabs = T, caption = "Median and median absolute deivations of estimated phenological dates for the polygons shown in Figure 1.") %>%
  # add_header_above(c(" ", "Grasslands" = 2, "Douglas Fri" = 2, "Pine-Oak"= 2,
  #                    "Rainfed agriculture" = 2)) %>%
  kable_styling(latex_options = c("hold_position"))


