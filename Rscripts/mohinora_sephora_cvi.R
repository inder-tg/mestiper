
# --- SUMMARY: All external cvi coincide in grouping the NDVI FPCA-idealized
# time series in 2 groups. All cvi found a dominating group, which basically coincides
# with the agro polygon; actually, the classification error is 0.006802721.
# This error would be impressive should we not have the results of the internal
# cvis from which we have observed that there are 2 different subgroups within
# the whole agro ROI. All cvis suggested the "single" method as the best clustering
# metric, but this is not surprising as this clustering method adopts a 
# ‘friends of friends’ clustering strategy.

# --- DEVELOPED by May 23, 2025
# --- FINDING Cluster Validity Indeces (cvi) for FPCA estimates in mohinora_agro
# --- GONNA use Silhouette, Davies-Bouldin, Calinski-Harabasz; if available I will
# --- try generalized Dunn, COP and SDbw
# --- PARTITION Similarity Measures used include: Variation Information (VI), Adjusted
# --- Rand and Jaccard

# --- UPDATED on Nov 19, 2026
# --- USED Agriculture polygon and applied cvi to find 2 clusters and 4 clusters according
# --- to Silhouette and DBStar cvi, respectively.

library(dtwclust)
library(terra)
library(here)
library(sephora)
library(sf)

LoadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env) 
}

output_mohinora_agro <- LoadToEnvironment( here("RData", "mohinora_sephora", 
                                                "AGRICULTURA DE TEMPORAL ANUAL",
                                                "mohinora_dtw2_sample_75_clusterSize_0p65_phenoParams.RData") )$output

mohinora_agro <- LoadToEnvironment( here("RData", "mohinora_sephora", 
                                         "AGRICULTURA DE TEMPORAL ANUAL",
                                         "rTp.RData") )$rTp
  
# --- PRELIMINARIES

# hclust methods:
# "ward.D"   "ward.D2"  "single"   "complete" "average"  "mcquitty"
# "median"   "centroid"

# --- ANALYSIS 

FPCA = vector("list", length(output_mohinora_agro))
for(i in 1:length(output_mohinora_agro)){
  FPCA[[i]] = output_mohinora_agro[[i]]$fpca
}

indices <- rep(FALSE, length(output_mohinora_agro))
for(i in 1:length(output_mohinora_agro)){
  indices[i] <- is.logical(FPCA[[i]])
}

which(indices == TRUE)

mohinora_agro$values[61,]
mohinora_agro$values[79,]

BASIS <- drbasis(n=75, q=2)

pixel_ts_61 <- ts(as.numeric(mohinora_agro$values[61,]),
                  start = c(2000,1), end = c(2025,23),
                  frequency = 23)

sephora_pixel_ts_61 <- sephora::phenopar(x=pixel_ts_61,
                                 distance = "dtw2", 
                                 basis = BASIS,
                                 clusterSize = 0.65,
                                 k=3)

tsPlot(sephora_pixel_ts_61)
rmPlot(sephora_pixel_ts_61)


# FPCA_trim <- vector("list", length(output_mohinora_agro)-sum(indices))
# str(mohinora_agro)
FPCA_trim <- FPCA
FPCA_trim <- FPCA_trim[!sapply(FPCA_trim, function(x) all(is.na(x)))]
length(FPCA_trim)

mohinora_agro_trim <- vector("list", length(output_mohinora_agro)-sum(indices))
mohinora_agro_trim$values <- mohinora_agro$values[-c(61,79),]
mohinora_agro_trim$coords <- mohinora_agro$coords[-c(61,79),]

# --- FIRST TRIAL

distanceMatrix <- proxy::dist(FPCA_trim, method = "dtw_basic", norm = "L2", 
                              step.pattern = symmetric2)

hc_out <- vector("list", length = sqrt(length(FPCA_trim)))

for(i in 2L:length(hc_out)){
  hc_out[[i]] <- tsclust(FPCA_trim, type = "hierarchical",
                         k = i, trace = TRUE,
                         control = hierarchical_control(method = "all",
                                                        distmat = distanceMatrix))
}

METHODS <- c("ward.D", "ward.D2", "single", "complete", 
             "average", "mcquitty", "median", "centroid")

# --- GETTING a results table for INTERNAL CVIs

cvi_eval_hc_out_k <- vector("list", length = length(hc_out))

for(i in 1:length(hc_out)){
  cvi_eval_hc_out_k[[i]] <- sapply(hc_out[[i]], cvi, b=rep(1,443))
}

save(hc_out, cvi_eval_hc_out_k, 
     file = here( "RData", "mohinora_sephora", 
                  "AGRICULTURA DE TEMPORAL ANUAL",
                  "hc_cvi.RData" ))


hc_out <- LoadToEnvironment( here("RData", "mohinora_sephora", 
                                             "AGRICULTURA DE TEMPORAL ANUAL",
                                             "hc_cvi.RData") )$hc_out

cvi_eval_hc_out_k <- LoadToEnvironment( here("RData", "mohinora_sephora", 
                                                     "AGRICULTURA DE TEMPORAL ANUAL",
                                                     "hc_cvi.RData") )$cvi_eval_hc_out_k


# --- Slihoutte

Sil_matrix <- matrix(nrow=length(hc_out), ncol=8)
for(i in 2:length(hc_out)){
  Sil_matrix[i,] <- cvi_eval_hc_out_k[[i]][6,]
}
colnames(Sil_matrix) <- METHODS

Sil_matrix[-1,]

max(Sil_matrix[-1,]) # ward.D2 [2] & complete [4]

which(Sil_matrix[-1,] == max(Sil_matrix[-1,]), arr.ind=TRUE) # 2 cluster

# --- Dunn

D_matrix <- matrix(nrow=length(hc_out), ncol=8)
for(i in 2:length(hc_out)){
  D_matrix[i,] <- cvi_eval_hc_out_k[[i]][11,]
}
colnames(D_matrix) <- METHODS

D_matrix[-1,]

max(D_matrix[-1,]) # single [3] & average [5] & centroid [8]

which(D_matrix[-1,] == max(D_matrix[-1,]), arr.ind=TRUE) # 2 clusters

# --- COP

COP_matrix <- matrix(nrow=length(hc_out), ncol=8)
for(i in 2:length(hc_out)){
  COP_matrix[i,] <- cvi_eval_hc_out_k[[i]][12,]
}
colnames(COP_matrix) <- METHODS

COP_matrix[-1,]

min(COP_matrix[-1,]) # ward.D [20]

which(COP_matrix[-1,] == min(COP_matrix[-1,]), arr.ind=TRUE)

# --- DB*

DBstar_matrix <- matrix(nrow=length(hc_out), ncol=8)
for(i in 2:length(hc_out)){
  DBstar_matrix[i,] <- cvi_eval_hc_out_k[[i]][10,]
}
colnames(DBstar_matrix) <- METHODS

DBstar_matrix[-1,]

min(DBstar_matrix[-1,]) # mcquitty [6]

which(DBstar_matrix[-1,] == min(DBstar_matrix[-1,]), arr.ind=TRUE) # 4 clusters

# --- CH

CH_matrix <- matrix(nrow=length(hc_out), ncol=8)
for(i in 2:length(hc_out)){
  CH_matrix[i,] <- cvi_eval_hc_out_k[[i]][8,]
}
colnames(CH_matrix) <- METHODS

CH_matrix[-1,]

max(CH_matrix[-1,]) # mcquitty [1]

which(CH_matrix[-1,] == max(CH_matrix[-1,]), arr.ind=TRUE) # 2 clusters


# ---

SF_matrix <- matrix(nrow=length(hc_out), ncol=8)
for(i in 2:length(hc_out)){
  SF_matrix[i,] <- cvi_eval_hc_out_k[[i]][7,]
}
colnames(SF_matrix) <- METHODS

SF_matrix[-1,]

max(SF_matrix[-1,]) # mcquitty [1]

# which(SF_matrix[-1,] == max(CH_matrix[-1,]), arr.ind=TRUE)

# --- Espacialización de los resultados ---

shpFiles <- list.files(path = here("data", 
                                   "mohinora_usv7"),
                       pattern = ".shp",
                       full.names = TRUE)

mohinora_USV <- read_sf(shpFiles[1])

# --- CLusters based on cvi Silhouettes

plot( st_geometry( mohinora_USV ), main = "Sil" )

Sil_indCluster1 = which(hc_out[[2]][[2]]@cluster == 1)
Sil_CLUSTER1 = vector("list", sum( hc_out[[2]][[2]]@cluster == 1 ))
for(i in Sil_indCluster1){
  Sil_CLUSTER1[[i]] = FPCA_trim[[i]]
}

Sil_indCluster2 = which(hc_out[[2]][[2]]@cluster == 2)
Sil_CLUSTER2 = vector("list", sum( hc_out[[2]][[2]]@cluster == 2 ))
for(i in Sil_indCluster2){
  Sil_CLUSTER2[[i]] = FPCA_trim[[i]]
}

for(i in Sil_indCluster1){
  points(x=mohinora_agro$coords[i,1], y=mohinora_agro$coords[i,2], 
         pch=18, col="red")
}
# ---
for(i in Sil_indCluster2){
  points(x=mohinora_agro$coords[i,1], y=mohinora_agro$coords[i,2], 
         pch=18, col="blue")
}

plot(Sil_CLUSTER1[[Sil_indCluster1[1]]], type="l", col="red", ylim=c(0,1e4),
     main = "Idealized NDVI via FPCA Sil Cluster1 vs. Cluster2", ylab="NDVI")
for(i in Sil_indCluster1){
  lines(Sil_CLUSTER1[[i]], col= "red")
}

for(i in Sil_indCluster2){
  lines(Sil_CLUSTER2[[i]], col= "blue")
}

# --- Clusters based on cvi DBstar

plot( st_geometry( mohinora_USV ), main = "DBstar" )

DBstar_indCluster1 = which(hc_out[[4]][[6]]@cluster == 1)
DBstar_CLUSTER1 = vector("list", sum( hc_out[[4]][[6]]@cluster == 1 ))
for(i in DBstar_indCluster1){
  DBstar_CLUSTER1[[i]] = FPCA_trim[[i]]
}

DBstar_indCluster2 = which(hc_out[[4]][[6]]@cluster == 2)
DBstar_CLUSTER2 = vector("list", sum( hc_out[[4]][[6]]@cluster == 2 ))
for(i in DBstar_indCluster2){
  DBstar_CLUSTER2[[i]] = FPCA_trim[[i]]
}

DBstar_indCluster3 = which(hc_out[[4]][[6]]@cluster == 3)
DBstar_CLUSTER3 = vector("list", sum( hc_out[[4]][[6]]@cluster == 3 ))
for(i in DBstar_indCluster3){
  DBstar_CLUSTER3[[i]] = FPCA_trim[[i]]
}

DBstar_indCluster4 = which(hc_out[[4]][[6]]@cluster == 4)
DBstar_CLUSTER4 = vector("list", sum( hc_out[[4]][[6]]@cluster == 4 ))
for(i in DBstar_indCluster4){
  DBstar_CLUSTER4[[i]] = FPCA_trim[[i]]
}


for(i in DBstar_indCluster1){
  points(x=mohinora_agro$coords[i,1], y=mohinora_agro$coords[i,2], 
         pch=18, col="red")
}
# ---
for(i in DBstar_indCluster2){
  points(x=mohinora_agro$coords[i,1], y=mohinora_agro$coords[i,2], 
         pch=18, col="blue")
}
# ---
for(i in DBstar_indCluster3){
  points(x=mohinora_agro$coords[i,1], y=mohinora_agro$coords[i,2], 
         pch=18, col="purple")
}
# ---
for(i in DBstar_indCluster4){
  points(x=mohinora_agro$coords[i,1], y=mohinora_agro$coords[i,2], 
         pch=18, col="yellow")
}



plot(DBstar_CLUSTER1[[DBstar_indCluster1[1]]], type="l", col="red", 
     ylim=c(0,1e4),
     main = "Idealized NDVI via FPCA DBstar Cluster1 vs. Cluster2 vs Cluster3", 
     ylab="NDVI")
for(i in DBstar_indCluster1){
  lines(DBstar_CLUSTER1[[i]], col= "red")
}
# ---
for(i in DBstar_indCluster2){
  lines(DBstar_CLUSTER2[[i]], col= "blue")
}
# ---
for(i in DBstar_indCluster3){
  lines(DBstar_CLUSTER3[[i]], col= "purple")
}
# ---
for(i in DBstar_indCluster4){
  lines(DBstar_CLUSTER4[[i]], col= "yellow")
}


# --- SPIRAL plot

library(ggnewscale)
library(vcd)
library(spiralize)
library(kableExtra)

source("Rscripts/getDist_phenoParam.R")
source("Rscripts/getSpiralPlot.R")
source("Rscripts/extract_sephora.R")
source("Rscripts/sephora_auxFUN.R")

# ---

# --- colors used in spiralPlot below
cgu <- rgb(173/255,221/255,142/255)
csos <- rgb(120/255,198/255,121/255)
cmat <- rgb(49/255, 163/255,84/255)
csen <- rgb(217/255, 95/255, 14/255)
ceos <- rgb(254/255, 153/255, 41/255)
cdor <- rgb(208/255, 209/255, 230/255)

colores <- c(cgu,csos,cmat,csen,ceos,cdor)

# --- how to get a SpiralPlot

DIRS <- list.dirs( here( "RData", "mohinora_sephora" ) )[-1]

list_DIRS <- setNames( as.list(DIRS), basename(DIRS) )

RDataFILES <- list.files( path = list_DIRS$`AGRICULTURA DE TEMPORAL ANUAL`, pattern = ".RData",
                          full.names = TRUE )

sephora_phenoparams_agro <- LoadToEnvironment( RDataFILES[3] )$output


mohinora_params_agro <- extract_sephora(x = sephora_phenoparams_agro, 
                                        object = "phenoparams")

# --- Clusters from Sil cvi

mohinora_params_agro_Sil_indCluster1 <- mohinora_params_agro[Sil_indCluster1,]
mohinora_params_agro_Sil_indCluster2 <- mohinora_params_agro[Sil_indCluster2,]

getSpiralPlot(MAT=mohinora_params_agro_Sil_indCluster1, 
              LABELS=month.name,
              vp_param=list(width=0.5, height=0.7))
vcd::grid_legend(x=1.215, y=0.125, pch=18, col=colores,
                 frame=FALSE,
                 labels=c("GU","SoS","Mat","Sen","EoS","Dor"),
                 title="Params")

getSpiralPlot(MAT=mohinora_params_agro_Sil_indCluster2, 
              LABELS=month.name,
              vp_param=list(width=0.5, height=0.7))
vcd::grid_legend(x=1.215, y=0.125, pch=18, col=colores,
                 frame=FALSE,
                 labels=c("GU","SoS","Mat","Sen","EoS","Dor"),
                 title="Params")

# --- sephora STATS


a <- phenoDatesToDoY(phenoDates=unlist(sephora_basicStats(mohinora_params_agro_Sil_indCluster1))[c(1,3,5,7,9,11)])
b <- phenoDatesToDoY(phenoDates=unlist(sephora_basicStats(mohinora_params_agro_Sil_indCluster2))[c(1,3,5,7,9,11)])
# b <- phenoDatesToDoY(phenoDates=unlist(sephora_basicStats(ayarin))[c(1,3,5,7,9,11)])
# c <- phenoDatesToDoY(phenoDates=unlist(sephora_basicStats(pinoEncino))[c(1,3,5,7,9,11)])
# d <- phenoDatesToDoY(phenoDates=unlist(sephora_basicStats(agro))[c(1,3,5,7,9,11)])

var_a <- redondea(unlist(sephora_basicStats(mohinora_params_agro_Sil_indCluster1))[c(1,3,5,7,9,11)+1])
var_b <- redondea(unlist(sephora_basicStats(mohinora_params_agro_Sil_indCluster2))[c(1,3,5,7,9,11)+1])
# var_b <- redondea(unlist(sephora_basicStats(ayarin))[c(1,3,5,7,9,11)+1])
# var_c <- redondea(unlist(sephora_basicStats(pinoEncino))[c(1,3,5,7,9,11)+1])
# var_d <- redondea(unlist(sephora_basicStats(agro))[c(1,3,5,7,9,11)+1])

df <- data.frame(#Date_Grasslands = paste( a$month, a$day, "(+", var_a, "days)"),
  # Date_Ayarin = paste( b$month, b$day, "(+", var_b, "days)"),
  # Date_PineOak = paste( a$month, a$day, "(+", var_a, "days)") )
Date_Agro_Sil1 = paste( a$month, a$day, "(+", var_a, "days)"),
Date_Agro_Sil2 = paste( b$month, b$day, "(+", var_b, "days)"))

# colnames(df) <- c("Grasslands", "Douglas Fir", "Pine-Oak", "Agriculture")
# colnames(df) <- c("Pine-Oak")
colnames(df) <- c("Agro_Sil1", "Agro_Sil2")
row.names(df) <- c("GU", "SoS", "Mat", "Sen", "EoS", "Dor")

kbl(df,
    longtable = T, booktabs = T, 
    caption = "Median and median absolute deivations of estimated phenological dates for Agriculture according to the 3 clusters suggested by Silhouette CVI.") %>%
  kable_styling(latex_options = c("hold_position"))

# --- Clusters from DBStar cvi

mohinora_params_agro_DBstar_indCluster1 <- mohinora_params_agro[DBstar_indCluster1,]
mohinora_params_agro_DBstar_indCluster2 <- mohinora_params_agro[DBstar_indCluster2,]
mohinora_params_agro_DBstar_indCluster3 <- mohinora_params_agro[DBstar_indCluster3,]
mohinora_params_agro_DBstar_indCluster4 <- mohinora_params_agro[DBstar_indCluster4,]


getSpiralPlot(MAT=mohinora_params_agro_DBstar_indCluster1, 
              LABELS=month.name,
              vp_param=list(width=0.5, height=0.7))
vcd::grid_legend(x=1.215, y=0.125, pch=18, col=colores,
                 frame=FALSE,
                 labels=c("GU","SoS","Mat","Sen","EoS","Dor"),
                 title="Params")

getSpiralPlot(MAT=mohinora_params_agro_DBstar_indCluster2, 
              LABELS=month.name,
              vp_param=list(width=0.5, height=0.7))
vcd::grid_legend(x=1.215, y=0.125, pch=18, col=colores,
                 frame=FALSE,
                 labels=c("GU","SoS","Mat","Sen","EoS","Dor"),
                 title="Params")

getSpiralPlot(MAT=mohinora_params_agro_DBstar_indCluster3, 
              LABELS=month.name,
              vp_param=list(width=0.5, height=0.7))
vcd::grid_legend(x=1.215, y=0.125, pch=18, col=colores,
                 frame=FALSE,
                 labels=c("GU","SoS","Mat","Sen","EoS","Dor"),
                 title="Params")

getSpiralPlot(MAT=mohinora_params_agro_DBstar_indCluster4, 
              LABELS=month.name,
              vp_param=list(width=0.5, height=0.7))
vcd::grid_legend(x=1.215, y=0.125, pch=18, col=colores,
                 frame=FALSE,
                 labels=c("GU","SoS","Mat","Sen","EoS","Dor"),
                 title="Params")

# --- sephora STATS

a <- phenoDatesToDoY(phenoDates=unlist(sephora_basicStats(mohinora_params_agro_DBstar_indCluster1))[c(1,3,5,7,9,11)])
b <- phenoDatesToDoY(phenoDates=unlist(sephora_basicStats(mohinora_params_agro_DBstar_indCluster2))[c(1,3,5,7,9,11)])
c <- phenoDatesToDoY(phenoDates=unlist(sephora_basicStats(mohinora_params_agro_DBstar_indCluster3))[c(1,3,5,7,9,11)])
d <- phenoDatesToDoY(phenoDates=unlist(sephora_basicStats(mohinora_params_agro_DBstar_indCluster4))[c(1,3,5,7,9,11)])

var_a <- redondea(unlist(sephora_basicStats(mohinora_params_agro_DBstar_indCluster1))[c(1,3,5,7,9,11)+1])
var_b <- redondea(unlist(sephora_basicStats(mohinora_params_agro_DBstar_indCluster2))[c(1,3,5,7,9,11)+1])
var_c <- redondea(unlist(sephora_basicStats(mohinora_params_agro_DBstar_indCluster3))[c(1,3,5,7,9,11)+1])
var_d <- redondea(unlist(sephora_basicStats(mohinora_params_agro_DBstar_indCluster4))[c(1,3,5,7,9,11)+1])

df <- data.frame(
  Date_Agro_DBstar1 = paste( a$month, a$day, "(+", var_a, "days)"),
  Date_Agro_DBstar2 = paste( b$month, b$day, "(+", var_b, "days)"),
  Date_Agro_DBstar3 = paste( c$month, c$day, "(+", var_c, "days)"),
  Date_Agro_DBstar4 = paste( d$month, d$day, "(+", var_d, "days)")
)

colnames(df) <- c("Agro_DBstar1", "Agro_DBstar2", "Agro_DBstar3", "Agro_DBstar4")
row.names(df) <- c("GU", "SoS", "Mat", "Sen", "EoS", "Dor")

kbl(df,
    longtable = T, booktabs = T, 
    caption = "Median and median absolute deivations of estimated phenological dates for Agriculture according to the 4 clusters suggested by DBstar CVI.") %>%
  # add_header_above(c(" ", "Grasslands" = 2, "Douglas Fri" = 2, "Pine-Oak"= 2,
  #                    "Rainfed agriculture" = 2)) %>%
  kable_styling(latex_options = c("hold_position"))



