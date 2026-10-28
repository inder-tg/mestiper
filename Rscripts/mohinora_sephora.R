# --- Escrito por Inder Tecuapetla, Mayo 23, 2024
# --- SEPHORA para Mohinora
# --- DATASET: NDVI MOD13Q1 en Cerro Mohinora, Chihuahua

# --- En este script estimamos los parámetros fenológicos de 
# --- una curva idealizada de NDVI
# --- Usamos código en paralelo para eficientar el cómputo 

# --- DATASET: NDVI MOD13Q1 v061 en Cerro Mohinora, Chihuahua, 2000-2023 

# --- ADDicionalmente, este script requiere archivos
# --- MOD13Q1_061_250m_16_days_NDVI_interpol.tif
# --- creado con el archivo mohinora_temporal_gapfilling.R

library(terra)
library(sf)
# library(tmap)
library(raster)
library(here)
library(dplyr)
library(tidyterra)

library(sephora)

source( "Rscripts/auxFUN.R" )

# (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# 
# BiocManager::install("ComplexHeatmap", force = TRUE)
# 
# library(ComplexHeatmap)

# --- DATA Loading

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

mohinora_DATA_rTp <- spRast_valuesCoords(mohinora_DATA) #rasterToPoints(mohinora_DATA) #spRast_valuesCoords(mohinora_DATA)

shpFiles <- list.files(path = here("data", "mohinora_usv7"),
                       pattern = ".shp",
                       full.names = TRUE)

mohinora_USV <- read_sf(shpFiles[1])


# ------------------
# ---- SEPHORA --- #
# - Preliminares - #
# ------------------

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

visual_mohinora <- mohinora_USV %>%
  sf::st_cast("MULTILINESTRING")

visual_mohinora$COLOR <- COLOR_USV

# --- tmap ---

# plot(subset(mohinora_DATA,453))

mohinora_tmap = tm_shape(subset(mohinora_DATA,453),
                         bbox = mohinora_USV) +
  tm_raster(col.scale = tm_scale_continuous(values=rev(terrain.colors(12))), 
            # palette = "YlGn",
            col.legend = tm_legend_hide()) +
            # legend.show = FALSE,
            # ) +
  tm_shape(visual_mohinora) + tm_lines(col="COLOR", lwd=5) + 
  # tm_layout(tm_title(""), frame = FALSE, bg.color = NA) #+
  # tm_compass(type = "8star", position = c("right", "bottom")) +
  # tm_scalebar(text.size = 0.65,
  #              position = c("right", "bottom")) #+
  tm_add_legend("symbol",
                labels=usv_NAMES,
                fill=usv_COLORS,
                border.col = "grey40",
                size=1,
                shape=18,
                is.portrait = TRUE)

mohinora_tmap

# ---

plot(subset(mohinora_DATA,453))
lines( visual_mohinora, lwd=4 )

XY <- locator() # 2616, 2441, 3455

xy <- get_timeSeries_byClicking(c(XY$x, XY$y), 
                                df=mohinora_DATA_rTp$coords)

pixel <- mohinora_DATA_rTp$values[xy$coord, ] * 1e-4

# pixel <- mohinora_DATA_rTp$values[2616, ]

# --- GRAPHIX

pixel_ts <- ts( pixel, start = c(2000,1), end = c(2025,23), frequency = 23)

plot(pixel_ts)

tsPlot(pixel_ts, yLab = "NDVI")

rmPlot(pixel_ts, yLab = "NDVI")

# --- MODELING

BASIS <- drbasis(n=75, q=2)

sephora_fit <- sephora::phenopar(x=pixel_ts,
                        distance = "dtw2", basis = BASIS,
                        k=3)

str( sephora_fit )

tsPlot(sephora_fit)

# --- 

custom_axis_theme <- theme(
  axis.text = element_text(size = 20, face = "bold"),
  axis.title = element_text(size = 16, face = "bold")
)

tsPlot(sephora_fit, yLab = "NDVI", xLab = "", position_legend="bottom") + 
  # axis.text.x = element_text(angle = 45, hjust = 1)) + 
  # theme_minimal(axis.text.x = element_text(angle = 45, hjust = 1))
  # labs(y = "NDVI", 
  #      title = "Time-Series for Point1") +
  labs(y = "NDVI", 
       title = "(B)") +
  theme(plot.title = element_text(size = 35, face = "bold", hjust = 0)
  ) +
  custom_axis_theme +
  theme(legend.text = element_text(size = 20, face = "bold"))

# ---

pixel_smooth <- ts( c(t(sephora_fit$xSmooth$fitOptAIC)),
                      start = c(2000,1),
                      end = c(2025,23),
                      frequency = 23 )

tsPlot( pixel_smooth )

rmPlot( pixel_smooth, position_legend = "bottom", yLab = "NDVI" )

# --- 

diff_years <- (2000:2025)[setdiff(1:26, sephora_fit$series)]

# --- AUX object
df <- data.frame(
  time = time(pixel_smooth),
  value = as.numeric(pixel_smooth)
) %>%
  mutate(
    year = floor(time),
    value = ifelse(year %in% diff_years, NA, value)
  )

ggplot(df, aes(time, value)) +
  geom_line(color = "darkgreen")

pixel_cluster <- ts( df$value, start = c(2000, 1),
                     end = c(2025, 23),
                     frequency = 23 )
                     
rmPlot(pixel_smooth, position_legend = "right", yLab = "NDVI")

# rmPlot(pixel_cluster, position_legend = "bottom", yLab = "NDVI")



# dim(sephora_fit$xSmooth$fitOptAIC)
# length(sephora_fit$series)
# sephora_fit$series

years_cluster <- (2000:2025)[sephora_fit$series]

df_cluster <- data.frame(
  DoY = factor(rep(1:23, times = length(years_cluster))),
  years = factor(rep(years_cluster, each = 23)),
  values = as.numeric(t(sephora_fit$xSmooth$fitOptAIC[
    sephora_fit$series, ,
    drop = FALSE
  ]))
)


ggplot(
  df_cluster,
  aes(
    x = DoY,
    y = values,
    group = years,
    colour = years
  )
) +
  geom_line() +
  geom_point()


outlier_idx <- setdiff(1:26, sephora_fit$series)

df_outliers <- data.frame(
  DoY = factor(rep(1:23, times = length(outlier_idx))),
  years = factor(rep(diff_years, each = 23)),
  values = as.numeric(
    t(
      sephora_fit$xSmooth$fitOptAIC[
        outlier_idx,
        ,
        drop = FALSE
      ]
    )
  )
)

ggplot(
  df_outliers,
  aes(
    x = DoY,
    y = values,
    group = years,
    colour = years
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2)


#  TODO: include in rmPlot function
# if( inherits(x, "sephora") ){
#   
#   years_cluster <- (x$startYear:x$endYear)[x$series]
#   
#   df <- data.frame(
#     DoY = factor(rep(1:x$freq, times = length(years_cluster))),
#     years = factor(rep(years_cluster, each = x$freq)),
#     values = as.numeric(
#       t(x$xSmooth$fitOptAIC[x$series,,drop=FALSE])
#     )
#   )
# }

# rmPlot(sephora_fit)

# --- sephora_fit sobre puesto

# x_axis <- get_metadata_years(x=pixel_cluster,
#                              startYear=2000, endYear=2025, frequency=23)
# 
# 
# x_axis_correct <- list()
# x_axis_correct$xDates <- x_axis$xDates[-c( (11*23+1):(12*23), (16*23+1):(17*23) )]
# x_axis_correct$xLabels <- x_axis$xLabels[-c(12,17)]
# 
# DoY <- seq(1,365, by = ceiling(365/23) ) #seq(1,365, length.out=23)
# DoY_aux <- seq(1,365, length.out=23)-1
# fpca_DoY <- sephora_fit$fpcaFun0Der(t=DoY_aux)
# 
# ggFull <- tsPlot(sephora_fit, yLab = "NDVI", xLab = "")
# 
# COLORS <- unique( ggplot_build(ggFull)$data[1][[1]]$colour )[-c(12,17)]
# number_observations <- length(years_cluster)*23
# 
# df <- data.frame(values=c(pixel_cluster[-c( (11*23+1):(12*23), (16*23+1):(17*23) )], 
#                           fpca_DoY),
#                  years=as.factor(rep(c(x_axis_correct$xLabels, 
#                                        "Idealized Phenology"), each=23)),
#                  DoY=factor(DoY, levels=DoY),
#                  # DoY=factor(DoY, levels=nameMONTH),
#                  class=c(rep(1,number_observations), rep(2,23)))
# 
# gg_fpca <- ggplot(data=df,
#                   aes(x=DoY, y=values, group=years, colour=years)) +
#   ggplot2::geom_line(linewidth = c(rep(1,number_observations), rep(4,23))) +
#   geom_point(size = 2) + #sizePoint=2
#   # ggplot2::labs(y="NDVI", x="", color="years+FPCA") +
#   ggplot2::labs(y="NDVI", x="", color="") +
#   ggplot2::scale_color_manual(values = c(COLORS, "#FF4500")) +
#   ggplot2::theme(legend.position = "bottom")# +
# 
# gg_fpca


# --- ALT

curves_cluster <-
  sephora_fit$xSmooth$fitOptAIC[
    sephora_fit$series,
    ,
    drop = FALSE
  ]

values_cluster <- as.numeric(t(curves_cluster))

years_cluster <- (2000:2025)[sephora_fit$series]

df <- data.frame(
  values = c(values_cluster, fpca_DoY),
  years = factor(
    rep(
      c(years_cluster, "Idealized Phenology"),
      each = 23
    )
  ),
  DoY = factor(DoY, levels = DoY),
  class = c(
    rep(1, length(values_cluster)),
    rep(2, 23)
  )
)

gg_fpca <- ggplot(data=df,
                  aes(x=DoY, y=values, group=years, colour=years)) +
  ggplot2::geom_line(linewidth = c(rep(1,number_observations), rep(4,23))) +
  geom_point(size = 2) + #sizePoint=2
  # ggplot2::labs(y="NDVI", x="", color="years+FPCA") +
  ggplot2::labs(y="NDVI", x="", color="") +
  ggplot2::scale_color_manual(values = c(COLORS, "#FF4500")) +
  ggplot2::theme(legend.position = "bottom")# +

gg_fpca


# ---

# --- PhenoParams

sephora_fit$phenoparams / 365

fpca_DoY

# Build dataframe for the curve

df <- data.frame(
  DOY = DoY_aux/364, #seq_along(fpca_DoY)/length(fpca_DoY),  # normalized DOY (0–1 scale)
  NDVI = fpca_DoY
)

# Extract only the non-NA phenology points
phenology_points <- data.frame(
  DOY = sephora_fit$phenoparams[c("GU","SoS","Mat","EoS")]/365,
  Stage = c("GU","SoS","Mat","EoS")
)

# Interpolate NDVI values at those DOYs
phenology_points$NDVI <- approx(df$DOY, df$NDVI, xout = phenology_points$DOY)$y

cgu <- rgb(173/255,221/255,142/255)
csos <- rgb(120/255,198/255,121/255)
cmat <- rgb(49/255, 163/255,84/255)
csen <- rgb(217/255, 95/255, 14/255)
ceos <- rgb(254/255, 153/255, 41/255)
cdor <- rgb(208/255, 209/255, 230/255)

colores <- c(cgu,csos,cmat,csen,ceos,cdor)

mycolors <- c(
  "GU" = colores[1],
  "SoS" = colores[2],
  "Mat" = colores[3],
  "EoS" = colores[5]
)

# ggplot(df, aes(x = DOY, y = NDVI)) +
#   geom_line(color = "#FF4500", size = 3) +  # solid red curve
#   geom_point(data = phenology_points, 
#              aes(x = DOY, y = NDVI, shape = Stage, color=Stage), 
#              size = 3) +
#   geom_text(data = phenology_points, 
#             aes(x = DOY, y = NDVI, label = Stage), # <-- map color
#             vjust = -1, size = 3, color="black") +  # , show.legend = FALSE hide legend for text
#   geom_vline(data = phenology_points, aes(xintercept = DOY, color = Stage),
#              linetype = "dashed", size = 1.5, show.legend = FALSE) +
#   scale_color_manual(values = mycolors) +
#   scale_x_continuous(
#     name = "",
#     breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)/365,
#     labels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
#   ) +
#   labs(y = "NDVI", 
#        # shape="TEST1",
#        color="TEST2",
#        title = "Pixel-level Idealized Phenology Curve with Key Dates") +
#   guides(shape="none")

# ---

ggplot(df, aes(x = DOY, y = NDVI)) +
  geom_line(color = "#FF4500", size = 3) +
  geom_point(data = phenology_points,
             aes(x = DOY, y = NDVI, 
                 shape = Stage, color = Stage,
                 fill=Stage),
             size = 3) +
  geom_text(data = phenology_points,
            aes(x = DOY, y = NDVI, label = Stage),
            vjust = -1, size = 6, color = "black") +
  geom_vline(data = phenology_points,
             aes(xintercept = DOY, color = Stage),
             linetype = "dashed", size = 1.5,
             show.legend = FALSE) +   # suppress legend for lines
  scale_color_manual(values = mycolors,
                     labels = c("EoS", "GU", "Mat", "SoS")) +
  scale_shape_manual(values = c(25, 18, 15, 17),#c(15:18),
                     labels = c("EoS", "GU", "Mat", "SoS")) +
  scale_fill_manual(values = mycolors,
                    labels = c("EoS", "GU", "Mat", "SoS")) +
  # scale_x_continuous(
  #   name = "",
  #   breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)/365,
  #   labels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  # ) +
  scale_x_continuous(
    name = "",
    breaks = DoY/365,
    labels = as.character(DoY)
  ) +
  labs(y = "NDVI",
       # title = "Pixel-level Idealized Phenology Curve with Key Dates",
       title = "(D)",
       fill = "Pheno. dates",
       shape = "Pheno. dates",   # legend title for shape
       color = "Pheno. dates") +   # legend title for color
  theme(plot.title = element_text(size = 35, face = "bold", hjust = 0)
  ) +
  guides(
    shape = guide_legend(
      override.aes = list(size=6)
    )
  ) +
  # custom_axis_theme +
  theme(
    legend.title = element_text(size=20, face = "bold"),
    legend.text = element_text(size = 17, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold"))

# --- parallel code 

dir.create(here( "RData", "mohinora_sephora" ))

mohinora_dtw2 <- phenopar_polygon(product = "independent",
                              data = mohinora_DATA_rTp$values,
                              frequency = 23,
                              save_all = TRUE,
                              distance = "dtw2",
                              basis = BASIS,
                              numCores = 6,
                              k=3,
                              confiBand = TRUE,
                              dirToSave = here( "RData", "mohinora_sephora" ),
                                # paste0(getwd(), "/RData/chunk_3_2/output"),
                              outputFileBaseName = "mohinora_dtw2_sample_75")

# --- 

source( "Rscripts/extract_sephora.R" )

output_chunk = LoadToEnvironment( paste0( getwd(), "/RData/mohinora_sephora/mohinora_dtw2_sample_100_phenoParams.RData" ))$output

# --- extract_sephora está en extract_sephora.R 
fpca = extract_sephora(x = output_chunk, object = "fpca" )

phenoparams = extract_sephora(x = output_chunk, object = "phenoparams" )

# --- otros objetos que se pueden extraer fácilmente son:
# "x" := serie de tiempo (ts) original
# "m_aug_smooth" := matriz con ts suavizadas, aumentadas a 75 puntos; 
# "clustering" := objeto tsclust con resultados de clustering (basado en dtw2); 
# "fpca_harmfit_params" := objeto con parámetros de ajuste vía regresión armónica; 
# "fpca_fun_0der" := objeto function para reproducir-remuestrear el ajuste vía fpca; 
# "fpca_fun_1der" := objeto function para reproducir-remuestrear la primera derivada del ajuste vía fpca
# "fpca_fun_2der" := objeto function para reproducir-remuestrear la segunda derivada del ajuste vía fpca
# "fpca_fun_3der" := objeto function para reproducir-remuestrear la tercera derivada del ajuste vía fpca
# "fpca_fun_4der" := objeto function para reproducir-remuestrear la cuarta derivada del ajuste vía fpca
# NOTA: para extraer otro objeto se puede adaptar extract_sephora, en caso de ser necesario
# yo puedo generar la adaptación sólo háganmelo saber.
# PARAMETROS disponibles: x, sigma, x_smooth, m_aug_smooth, usedTotal, series,
# clusterSize, dominatingClusterSize, clustering, fpca, fpca_harmfit_params, 
# fpca_fun_0der, fpca_fun_1der, fpca_fun_2der, fpca_fun_3der, fpca_fun_4der,
# phenoparams, mixedmodel_status, mixedmodel_inverted, optimal_num_harms

# --- Guardando fpca & phenoparams en archivos RDS -más económicos que RData 
saveRDS(fpca, file = file.path( getwd(), "RData/mohinora_sephora/fpca.rds" ) )

saveRDS(phenoparams, file = file.path( getwd(), "RData/NDVI_2016_2024_chunk4/output/chunk4_dtw2_phenoparams.rds" ) )


# --- output_extraction ---

source( "Rscripts/extract_sephora.R" )

output_chunk = LoadToEnvironment( paste0( getwd(), "/RData/mohinora_sephora/mohinora_dtw2_sample_100_phenoParams.RData" ))$output

# --- extract_sephora está en extract_sephora.R 
fpca = extract_sephora(x = output_chunk, object = "fpca" )

phenoparams = extract_sephora(x = output_chunk, object = "phenoparams" )

# --- otros objetos que se pueden extraer fácilmente son:
# "x" := serie de tiempo (ts) original
# "m_aug_smooth" := matriz con ts suavizadas, aumentadas a 75 puntos; 
# "clustering" := objeto tsclust con resultados de clustering (basado en dtw2); 
# "fpca_harmfit_params" := objeto con parámetros de ajuste vía regresión armónica; 
# "fpca_fun_0der" := objeto function para reproducir-remuestrear el ajuste vía fpca; 
# "fpca_fun_1der" := objeto function para reproducir-remuestrear la primera derivada del ajuste vía fpca
# "fpca_fun_2der" := objeto function para reproducir-remuestrear la segunda derivada del ajuste vía fpca
# "fpca_fun_3der" := objeto function para reproducir-remuestrear la tercera derivada del ajuste vía fpca
# "fpca_fun_4der" := objeto function para reproducir-remuestrear la cuarta derivada del ajuste vía fpca
# NOTA: para extraer otro objeto se puede adaptar extract_sephora, en caso de ser necesario
# yo puedo generar la adaptación sólo háganmelo saber.
# PARAMETROS disponibles: x, sigma, x_smooth, m_aug_smooth, usedTotal, series,
# clusterSize, dominatingClusterSize, clustering, fpca, fpca_harmfit_params, 
# fpca_fun_0der, fpca_fun_1der, fpca_fun_2der, fpca_fun_3der, fpca_fun_4der,
# phenoparams, mixedmodel_status, mixedmodel_inverted, optimal_num_harms

# --- Guardando fpca & phenoparams en archivos RDS -más económicos que RData 
saveRDS(fpca, file = file.path( getwd(), "RData/mohinora_sephora/fpca.rds" ) )

saveRDS(phenoparams, file = file.path( getwd(), "RData/mohinora_sephora/phenoparams.rds" ) )


# --- spliting by groups

na_pattern <- apply(is.na(phenoparams), 1, 
                    function(x) paste(which(x), collapse = "_"))
na_pattern[na_pattern == ""] <- "complete"

unique(na_pattern)

# "complete"

group_A <- phenoparams[na_pattern == "complete", ] # 35 pixels

# "4_6"
group_B <- phenoparams[na_pattern == "4_6", ] # 4679 pixels

# "4_5_6"
group_C <- phenoparams[na_pattern == "4_5_6", ] # 122 pixels


# --- Getting actual pixels

pxl_group <- as.numeric(unlist(lapply(rownames(group_A), 
                                      function(s) strsplit(s, "result.")[[1]][2] )))

fpca_group <- fpca[pxl_group]

mat_fpca_group <- matrix(nrow = length(fpca_group), ncol = 75)

for(i in 1:length(fpca_group)){
  mat_fpca_group[i,] <-  fpca_group[[i]]
}

mins_grupo <- apply(mat_fpca_group, 1, min) 

positive_fpca_group <- which(mins_grupo > 0) # length(positive_fpca_group)

pxl_group_positive <- pxl_group[positive_fpca_group] # (B,211580 pxls); (H, 52 pxls)

fpca_group_positive <- fpca_group[positive_fpca_group]

P <- sample(1:length(pxl_group_positive), 1)
plot(fpca_group_positive[[P]], type="l", 
     main= paste0("Pixel-", P, " Grupo A"),
     ylab = "")

# --- Visuals

positive_groupA_coords <- mohinora_DATA_rTp$coords[pxl_group_positive,]

positive_groupA_sf <- st_as_sf(
  data.frame(positive_groupA_coords),
  coords = c("x", "y"),  # adjust column names if needed
  crs = crs(mohinora_DATA)    # match raster CRS
)

# 3. Create the plot
pa <- 
  ggplot() +
  # Plot the background raster layer
  geom_spatraster(data = subset(mohinora_DATA, 453)) +
  
  # Overlay the simple feature points
  geom_sf(data = positive_groupA_sf, 
          color = "red", 
          size = 2, 
          alpha = 0.8) +
  
  # Choose a color palette for your raster values (e.g., viridis)
  scale_fill_viridis_c(na.value = "transparent", name = "Value") +
  
  # Apply a clean theme for spatial maps
  theme_minimal() +
  labs(
    title = "Spatial Points Overlaid on Raster",
    x = "",
    y = ""
  )

pa

# --- group B

pxl_group <- as.numeric(unlist(lapply(rownames(group_B), 
                                      function(s) strsplit(s, "result.")[[1]][2] )))

fpca_group <- fpca[pxl_group]

mat_fpca_group <- matrix(nrow = length(fpca_group), ncol = 75)

for(i in 1:length(fpca_group)){
  mat_fpca_group[i,] <-  fpca_group[[i]]
}

mins_grupo <- apply(mat_fpca_group, 1, min) 

positive_fpca_group <- which(mins_grupo > 0)

pxl_group_positive <- pxl_group[positive_fpca_group] #

fpca_group_positive <- fpca_group[positive_fpca_group]

P <- sample(1:length(pxl_group_positive), 1)
plot(fpca_group_positive[[P]], type="l", 
     main= paste0("Pixel-", P, " Grupo B"),
     ylab = "")

# --- Visuals

positive_groupB_coords <- mohinora_DATA_rTp$coords[pxl_group_positive,]

positive_groupB_sf <- st_as_sf(
  data.frame(positive_groupB_coords),
  coords = c("x", "y"),  # adjust column names if needed
  crs = crs(mohinora_DATA)    # match raster CRS
)


# 3. Create the plot
pb <- 
  ggplot() +
  # Plot the background raster layer
  geom_spatraster(data = subset(mohinora_DATA, 453)) +
  
  # Overlay the simple feature points
  geom_sf(data = positive_groupB_sf, 
          color = "red", 
          size = 2, 
          alpha = 0.8) +
  
  # Choose a color palette for your raster values (e.g., viridis)
  scale_fill_viridis_c(na.value = "transparent", name = "Value") +
  
  # Apply a clean theme for spatial maps
  theme_minimal() +
  labs(
    title = "Spatial Points Overlaid on Raster",
    x = "",
    y = ""
  )

pb

# --- group C

pxl_group <- as.numeric(unlist(lapply(rownames(group_C), 
                                      function(s) strsplit(s, "result.")[[1]][2] )))

fpca_group <- fpca[pxl_group]

mat_fpca_group <- matrix(nrow = length(fpca_group), ncol = 75)

for(i in 1:length(fpca_group)){
  mat_fpca_group[i,] <-  fpca_group[[i]]
}

mins_grupo <- apply(mat_fpca_group, 1, min) 

positive_fpca_group <- which(mins_grupo > 0)

pxl_group_positive <- pxl_group[positive_fpca_group] # 

fpca_group_positive <- fpca_group[positive_fpca_group]

P <- sample(1:length(pxl_group_positive), 1)
plot(fpca_group_positive[[P]], type="l", 
     main= paste0("Pixel-", P, " Grupo B"),
     ylab = "")

# --- Visuals

positive_groupC_coords <- mohinora_DATA_rTp$coords[pxl_group_positive,]

positive_groupC_sf <- st_as_sf(
  data.frame(positive_groupC_coords),
  coords = c("x", "y"),  # adjust column names if needed
  crs = crs(mohinora_DATA)    # match raster CRS
)


# 3. Create the plot
pc <- 
  ggplot() +
  # Plot the background raster layer
  geom_spatraster(data = subset(mohinora_DATA, 453)) +
  
  # Overlay the simple feature points
  geom_sf(data = positive_groupC_sf, 
          color = "red", 
          size = 2, 
          alpha = 0.8) +
  
  # Choose a color palette for your raster values (e.g., viridis)
  scale_fill_viridis_c(na.value = "transparent", name = "Value") +
  
  # Apply a clean theme for spatial maps
  theme_minimal() +
  labs(
    title = "Spatial Points Overlaid on Raster",
    x = "",
    y = ""
  )

pc

# --- Visuals combined

temp <- st_transform(mohinora_USV, st_crs(positive_groupA_sf))

# 2. Crop the points down to the shapefile boundary
cropped_points_A <- st_intersection(positive_groupA_sf, temp)
cropped_points_B <- st_intersection(positive_groupB_sf, temp)
cropped_points_C <- st_intersection(positive_groupC_sf, temp)

ggplot() +
  # Background raster layer
  geom_spatraster(data = subset( mohinora_DATA, 453 )) +
  
  # Layer 1: Cropped Group A Points
  geom_sf(data = cropped_points_A, aes(color = "Group A"), size = 2.5, alpha = 0.8) +
  
  # Layer 2: Cropped Group B Points
  geom_sf(data = cropped_points_B, aes(color = "Group B"), size = 2.5, alpha = 0.8) +
  
  # Layer 3: Cropped Group C Points
  geom_sf(data = cropped_points_C, aes(color = "Group C"), size = 2.5, alpha = 0.8) +
  
  # Optional: Overlay the shapefile border for visual reference
  geom_sf(data = mohinora_USV, fill = NA, color = "black", linewidth = 0.5) +
  
  # Color palette for your 3 point categories
  scale_color_manual(
    name = "Point Groups",
    values = c("Group A" = "#E41A1C", "Group B" = "#377EB8", "Group C" = "#4DAF4A")
  ) +
  
  # Color palette for the underlying raster values
  scale_fill_viridis_c(na.value = "transparent", name = "Raster Value") +
  
  # Styling for clean maps
  theme_minimal() +
  labs(
    title = "Cropped Point Groups Overlaid on Raster",
    x = "",
    y = ""
  )















