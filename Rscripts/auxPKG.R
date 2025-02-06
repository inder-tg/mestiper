# -----------------------------------------------------------------------------
#
# Elaborado por Inder Tecuapetla, May 31, 2023
#
# Modificado Julio 28, 2023
#
# Instalación de paquetes/bibliotecas a utilizar en este módulo
# 
# Hecho para SELPER/CEOS Working Group Chapter D Training Group  
#
# -----------------------------------------------------------------------------

neededPackages <- c("raster", "sf", "mapview", "geoTS", "fpp",
                    "forecast", "Kendall", "trend", "RColorBrewer",
                    "gtools", "sta", "changepoint", "bfast")
                    

packagesToInstall <- neededPackages[!(neededPackages %in% installed.packages()[,"Package"])]

if( length(packagesToInstall) ){
  for( i in 1:length(packagesToInstall) ){
    message("Installing package", packagesToInstall[i], "\n")
    install.packages(packagesToInstall[i], dependencies = TRUE)
  }
} 

# library(rgdal)
library(raster)
# library(sp)
library(sf)
library(mapview)
library(geoTS)

library(fpp)
library(forecast)
library(Kendall)
library(trend)
library(RColorBrewer)
library(gtools)

library(sta)
library(changepoint)
library(bfast)
