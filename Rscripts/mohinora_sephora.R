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
library(tmap)
library(raster)

# instalar clv, longitudinalData & TSclust from tar files

library(TSclust)
library(sephora)


# 1. Install 'clv' (Required by longitudinalData)
clv_url <- "https://r-project.org"
install.packages(clv_url, repos = NULL, type = "source")


# Install the remotes manager if you don't have it already
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install TSclust from the author's repository
remotes::install_github("pmontman/TSclust")


library(TSclust)


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

mohinora_NDVI_files <- c(FILES_NDVI_imputation, FILES_NDVI_interpolation[1:595])

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

# library(sephora)

pixel_ts <- ts( pixel, start = c(2000,1), end = c(2025,23), frequency = 23)

plot(pixel_ts)

tsPlot(pixel_ts, yLab = "NDVI")

rmPlot(pixel_ts, yLab = "NDVI")

# --- MODELING

BASIS <- drbasis(n=100, q=2)

sephora_fit <- sephora::phenopar(x=pixel_ts,
                        distance = "dtw2", basis = BASIS,
                        k=3)

tsPlot(sephora_fit)

# tsPlot(sephora_fit, type="ms")

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

rmPlot( pixel_smooth )

# --- 
(2000:2025)[setdiff(1:26, sephora_fit$series)]

# --- AUX object
df <- data.frame(
  time = time(pixel_smooth),
  value = as.numeric(pixel_smooth)
) %>%
  mutate(
    year = floor(time),
    value = ifelse(year %in% c(2011, 2016), NA, value)
  )

ggplot(df, aes(time, value)) +
  geom_line(color = "darkgreen")

pixel_cluster <- ts( df$value, start = c(2000,1),
                     end = c(2025,23),
                     frequency = 23  )
                     
tsPlot(pixel_cluster)                     

rmPlot(pixel_cluster)

rmPlot(sephora_fit)

# --- sephora_fit sobre puesto

x_axis <- get_metadata_years(x=pixel_cluster,
                             startYear=2000, endYear=2025, frequency=23)


x_axis_correct <- list()
x_axis_correct$xDates <- x_axis$xDates[-c( (11*23+1):(12*23), (16*23+1):(17*23) )]
x_axis_correct$xLabels <- x_axis$xLabels[-c(12,17)]

DoY <- seq(1,365, by = ceiling(365/23) ) #seq(1,365, length.out=23)
DoY_aux <- seq(1,365, length.out=23)-1
fpca_DoY <- sephora_fit$fpcaFun0Der(t=DoY_aux)

# sephora_fit$fpcaFun0Der(t=DoY)
# sephora_fit$fpcaFun0Der(t=DoY_aux)
# 
# plot(DoY, sephora_fit$fpcaFun0Der(t=DoY))
# abline(v=0)
# plot(DoY_aux, sephora_fit$fpcaFun0Der(t=DoY_aux))
# abline(v=0)

ggFull <- tsPlot(sephora_fit, yLab = "NDVI", xLab = "")

COLORS <- unique( ggplot_build(ggFull)$data[1][[1]]$colour )[-c(12,17)]
number_observations <- 24*23

df <- data.frame(values=c(pixel_cluster[-c( (11*23+1):(12*23), (16*23+1):(17*23) )], 
                          fpca_DoY),
                 years=as.factor(rep(c(x_axis_correct$xLabels, 
                                       "Idealized Phenology"), each=23)),
                 DoY=factor(DoY, levels=DoY),
                 # DoY=factor(DoY, levels=nameMONTH),
                 class=c(rep(1,number_observations), rep(2,23)))

gg_fpca <- ggplot(data=df,
                  aes(x=DoY, y=values, group=years, colour=years)) +
  ggplot2::geom_line(linewidth = c(rep(1,number_observations), rep(4,23))) +
  geom_point(size = 2) + #sizePoint=2
  # ggplot2::labs(y="NDVI", x="", color="years+FPCA") +
  ggplot2::labs(y="NDVI", x="", color="") +
  ggplot2::scale_color_manual(values = c(COLORS, "#FF4500")) +
  ggplot2::theme(legend.position = "bottom")# +

gg_fpca


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
  Stage = c("Green-up","Start of Season","Maturity","End of Season")
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
  "Green-up" = colores[1],
  "Start of Season" = colores[2],
  "Maturity" = colores[3],
  "End of Season" = colores[5]
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
                              k=3,
                              confiBand = TRUE,
                              dirToSave = here( "RData", "mohinora_sephora" ),
                                # paste0(getwd(), "/RData/chunk_3_2/output"),
                              outputFileBaseName = "mohinora_dtw2_sample_100")






sephora_fit <- sephora::phenopar(x=pixel_ts,
                                 distance = "dtw2", 
                                 basis = BASIS,
                                 k=3)



















