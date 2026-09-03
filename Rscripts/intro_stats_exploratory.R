############################################################
# Curso CADI - Tecnológico de Monterrey
# ============================================================
# Módulo 3: Relación entre variables y pensamiento sistémico
# ============================================================
# Dataset: CalCOFI (datos_calcofi_corregido.csv)
# Autor: Inder Rafael Tecuapetla Gómez
############################################################

library(dplyr)
library(readr)
library(ggplot2)
library(GGally)   # para ggpairs
library(corrplot) # para matrices de correlación
library(ggcorrplot)
library(here)
library(patchwork)

# ----------------------------------------------------------
# 0. Lectura de datos
# ----------------------------------------------------------
# Leemos el dataset depurado de CalCOFI (corregido)
calcofi <- read_csv( here( "data", "datos_calcofi_corregido.csv" ) )

# ------------------------------------------------------------
# 1. Correlación y asociaciones entre variables
# ------------------------------------------------------------

# Seleccionamos todas las variables numéricas en el dataset
calcofi_subset <- calcofi %>%
  select( where(is.numeric) )

calcofi_subset <- calcofi_subset %>%
  rename(
    Chl = mean_ChlorA,
    NO2 = mean_NO2uM,
    NO3 = mean_NO3uM,
    O2 = mean_O2ml_L,
    O2Sat = mean_O2Sat,
    Phaeop = mean_Phaeop,
    PO4 = mean_PO4uM,
    Salinity = mean_Salnty,
    SiO3 = mean_SiO3uM,
    Temperature = mean_T_degC
  ) %>%
  select(-Distance, -Year, -Lat_Dec, -Lon_Dec, -Cst_Cnt, -Quarter, -Bottom_D)

cor_matrix <- cor(calcofi_subset, use = "complete.obs")

# Visualizamos la matriz de correlación
corrplot(cor_matrix, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45,
         title = "Matriz de correlación entre variables ambientales")

ggcorrplot(cor_matrix, hc.order = TRUE, type = "lower",
           lab = TRUE, lab_size = 3,
           colors = c("red", "white", "blue"),
           title = "Matriz de correlación entre variables ambientales")

# Las variables con mayor correlación positiva son PO4 y NO3
# Las variables con mayor correlación negativa son O2 y NO3

# ------------------------------------------------------------
# 2. Visualización de relaciones (pares)
# ------------------------------------------------------------

calcofi_mas_corr <- calcofi %>%
  select(mean_PO4uM, mean_O2ml_L, mean_NO3uM, Quarter, Dist_cat) %>%
  rename(
    PO4 = mean_PO4uM,
    O2 = mean_O2ml_L,
    NO3 = mean_NO3uM
  )

# Usamos GGally::ggpairs para ver relaciones bivariadas
calcofi_mas_corr_subset <- calcofi_mas_corr %>%
  select(-Quarter, -Dist_cat)

ggpairs(calcofi_mas_corr_subset,
        title = "Relaciones bivariadas entre variables ambientales")

# ------------------------------------------------------------
# 3. Interpretación contextual (no causalidad)
# ------------------------------------------------------------
# Ejemplo: scatterplot de Fosfatos vs nitratos
ggplot(calcofi_mas_corr, aes(x = PO4, y = NO3, color = Dist_cat)) +
  geom_point(alpha = 0.5) +
  # geom_smooth(method = "loess", se = FALSE) +
  labs(title = "Relación entre fosfato y nitrato por tipo de zona",
       x = "PO4 (mug/L)", y = "Nitrato (mumol/L)", color = "Zona") +
  theme_minimal()

# Ejemplo: scatterplot de Oxígeno vs nitratos
ggplot(calcofi_mas_corr, aes(x = O2, y = NO3, color = Dist_cat)) +
  geom_point(alpha = 0.5) +
  # geom_smooth(method = "loess", se = FALSE) +
  labs(title = "Relación entre oxígeno y nitrato por tipo de zona",
       x = "PO4 (ml/L)", y = "Nitrato (mumol/L)", color = "Zona") +
  theme_minimal()

# ------------------------------------------------------------
# 4. Introducción a variables ambientales interrelacionadas
# ------------------------------------------------------------

ggpairs(calcofi_mas_corr_subset,
        title = "Relaciones bivariadas entre variables ambientales (con alta correlación)")


# ------------------------------------------------------------
# 5. Planteamiento de preguntas estadísticas a partir de los datos
# ------------------------------------------------------------
# Ejemplo de preguntas que se pueden derivar:
# - ¿Qué variables muestran mayor correlación positiva o negativa?
# - ¿Cómo cambian las asociaciones según la categoría de distancia?
# - ¿Qué patrones sugieren estrés ambiental (ej. hipoxia)?
# - ¿Qué variables podrían ser más relevantes para biodiversidad o pesca?


# ------------------------------------------------------------
# 6. Aplicaciones: Productividad & Biodivesidad
# ------------------------------------------------------------

# De acuerdo al Diccionario-Calcofi.pdf, valores altos de PO4 (fosfatos) están relacionados
# con pérdida de Biodiversidad en los océanos.

PO4 <- calcofi %>%
  select(mean_PO4uM, mean_T_degC, mean_O2ml_L, mean_O2Sat,
         mean_NO3uM, mean_SiO3uM, Quarter, Dist_cat) %>%
  rename(
    PO4 = mean_PO4uM,
    Temperature = mean_T_degC,
    O2 = mean_O2ml_L,
    O2Sat = mean_O2Sat,
    NO3 = mean_NO3uM,
    SiO3 = mean_SiO3uM
  )

# Usamos GGally::ggpairs para ver relaciones bivariadas
PO4_subset <- PO4 %>%
  select(-Quarter, -Dist_cat)

cor_PO4 <- cor(PO4_subset, use = "complete.obs")

ggcorrplot(cor_PO4, hc.order = TRUE, type = "lower",
           lab = TRUE, lab_size = 3,
           colors = c("red", "white", "blue"),
           title = "Matriz de correlación (con interés en PO4)")

ggplot(PO4, aes(x = PO4, y = NO3, color = Dist_cat)) +
  geom_point(alpha = 0.6) +
  # geom_smooth(method = "loess", se = FALSE) +
  labs(title = "Relación entre fosfato y nitrato",
       x = "Fosfato (µg/L)", y = "Nitrato (µmol/L)",
       color = "Zona") +
  theme_minimal()

ggplot(PO4, aes(x = PO4, y = NO3, color = Quarter)) +
  geom_point(alpha = 0.6) +
  # geom_smooth(method = "loess", se = FALSE) +
  labs(title = "Relación entre fosfato y nitrato",
       x = "Fosfato (µg/L)", y = "Nitrato (µmol/L)",
       color = "Trimestre") +
  theme_minimal()

# ---

ggplot(PO4, aes(x = PO4, y = Temperature, color = Dist_cat)) +
  geom_point(alpha = 0.6) +
  # geom_smooth(method = "loess", se = FALSE) +
  labs(title = "Relación entre fosfato y temperatura",
       x = "Fosfato (µg/L)", y = "Temp (°C)",
       color = "Zona") +
  theme_minimal()

ggplot(PO4, aes(x = PO4, y = Temperature, color = Quarter)) +
  geom_point(alpha = 0.6) +
  # geom_smooth(method = "loess", se = FALSE) +
  labs(title = "Relación entre fosfato y temperatura",
       x = "Fosfato (µg/L)", y = "Temp (°C)",
       color = "Trimestre") +
  theme_minimal()

PO4_Q4 <- ggplot(filter(PO4, Quarter == 4), aes(x = PO4, y = Temperature, color = Quarter)) +
  geom_point(alpha = 0.6) +
  labs(title = "Relación entre fosfato y temperatura (Quarter 4)",
       x = "Fosfato (µg/L)", y = "Temp (°C)",
       color = "Trimestre") +
  theme_minimal()

PO4_Q1 <- ggplot(filter(PO4, Quarter == 1), aes(x = PO4, y = Temperature, color = Quarter)) +
  geom_point(alpha = 0.6) +
  labs(title = "Relación entre fosfato y temperatura (Quarter 1)",
       x = "Fosfato (µg/L)", y = "Temp (°C)",
       color = "Trimestre") +
  theme_minimal()

PO4_Q1 + PO4_Q4 

# ---

ggplot(PO4, aes(x = PO4, y = O2, color = Dist_cat)) +
  geom_point(alpha = 0.6) +
  labs(title = "Relación entre fosfato y oxígeno",
       x = "Fosfato (µg/L)", y = "Oxígeno (ml/L)",
       color = "Zona") +
  theme_minimal()

ggplot(PO4, aes(x = PO4, y = O2, color = Quarter)) +
  geom_point(alpha = 0.6) +
  labs(title = "Relación entre fosfato y oxígeno",
       x = "Fosfato (µg/L)", y = "Oxígeno (ml/L)",
       color = "Trimestre") +
  theme_minimal()

PO4_O2_Q4 <- ggplot(filter(PO4, Quarter == 4), aes(x = PO4, y = O2, color = Quarter)) +
  geom_point(alpha = 0.6) +
  labs(title = "Relación entre fosfato y oxígeno (Quarter 4)",
       x = "Fosfato (µg/L)", y = "Oxígeno (ml/L)",
       color = "Trimestre") +
  theme_minimal()

PO4_O2_Q3 <- ggplot(filter(PO4, Quarter == 3), aes(x = PO4, y = O2, color = Quarter)) +
  geom_point(alpha = 0.6) +
  labs(title = "Relación entre fosfato y oxígeno (Quarter 3)",
       x = "Fosfato (µg/L)", y = "Oxígeno (ml/L)",
       color = "Trimestre") +
  theme_minimal()

PO4_O2_Q2 <- ggplot(filter(PO4, Quarter == 2), aes(x = PO4, y = O2, color = Quarter)) +
  geom_point(alpha = 0.6) +
  labs(title = "Relación entre fosfato y oxígeno (Quarter 2)",
       x = "Fosfato (µg/L)", y = "Oxígeno (ml/L)",
       color = "Trimestre") +
  theme_minimal()


PO4_O2_Q1 <- ggplot(filter(PO4, Quarter == 1), aes(x = PO4, y = O2, color = Quarter)) +
  geom_point(alpha = 0.6) +
  labs(title = "Relación entre fosfato y oxígeno (Quarter 1)",
       x = "Fosfato (µg/L)", y = "Oxígeno (ml/L)",
       color = "Trimestre") +
  theme_minimal()

PO4_O2_Q1 + PO4_O2_Q2 + PO4_O2_Q3 + PO4_O2_Q4 

# ---

PO4_NO3_Q4 <- ggplot(filter(PO4, Quarter == 4), aes(x = PO4, y = NO3, color = Dist_cat)) +
  geom_point(alpha = 0.6) +
  labs(title = "Relación entre fosfato y nitrato (Quarter 4)",
       x = "Fosfato (µg/L)", y = "Nitrato (µmol/L)",
       color = "Categoría de distancia") +
  theme_minimal()

PO4_NO3_Q3 <- ggplot(filter(PO4, Quarter == 3), aes(x = PO4, y = NO3, color = Dist_cat)) +
  geom_point(alpha = 0.6) +
  labs(title = "Relación entre fosfato y nitrato (Quarter 3)",
       x = "Fosfato (µg/L)", y = "Nitrato (µmol/L)",
       color = "Categoría de distancia") +
  theme_minimal()

PO4_NO3_Q2 <- ggplot(filter(PO4, Quarter == 2), aes(x = PO4, y = NO3, color = Dist_cat)) +
  geom_point(alpha = 0.6) +
  labs(title = "Relación entre fosfato y nitrato (Quarter 2)",
       x = "Fosfato (µg/L)", y = "Nitrato (µmol/L)",
       color = "Categoría de distancia") +
  theme_minimal()

PO4_NO3_Q1 <- ggplot(filter(PO4, Quarter == 1), aes(x = PO4, y = NO3, color = Dist_cat)) +
  geom_point(alpha = 0.6) +
  labs(title = "Relación entre fosfato y nitrato (Quarter 1)",
       x = "Fosfato (µg/L)", y = "Nitrato (µmol/L)",
       color = "Categoría de distancia") +
  theme_minimal()

PO4_NO3_Q1 + PO4_NO3_Q2 + PO4_NO3_Q3 + PO4_NO3_Q4 

# ---




