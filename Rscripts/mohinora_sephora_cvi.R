
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

library(dtwclust)
library(terra)
library(here)
library(sephora)
library(sf)

LoadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env) 
}

output_mohinora_agro <- LoadToEnvironment( here("mestiper", "RData", "mohinora_sephora", 
                                                "AGRICULTURA DE TEMPORAL ANUAL",
                                                "mohinora_dtw2_sample_75_clusterSize_0p65_phenoParams.RData") )$output
  
  # paste0( getwd(), "/RData/mohinora_agro/mohinora_agro_phenoParams.RData" ))$output

mohinora_agro <- LoadToEnvironment( here("mestiper", "RData", "mohinora_sephora", 
                                         "AGRICULTURA DE TEMPORAL ANUAL",
                                         "rTp.RData") )$rTp
  
  # paste0( getwd(), "/RData/mohinora_agro/mohinora_agro.RData" ))$rTp

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
     file = here( "mestiper",  "RData", "mohinora_sephora", 
                  "AGRICULTURA DE TEMPORAL ANUAL",
                  "hc_cvi.RData" ))

# --- Slihoutte

Sil_matrix <- matrix(nrow=length(hc_out), ncol=8)
for(i in 2:length(hc_out)){
  Sil_matrix[i,] <- cvi_eval_hc_out_k[[i]][6,]
}
colnames(Sil_matrix) <- METHODS

# Sil_matrix[-1,]
# 
# max(Sil_matrix[-1,])
# 
# which(Sil_matrix[-1,] == max(Sil_matrix[-1,]), arr.ind=TRUE)

D_matrix <- matrix(nrow=length(hc_out), ncol=8)
for(i in 2:length(hc_out)){
  D_matrix[i,] <- cvi_eval_hc_out_k[[i]][11,]
}
colnames(D_matrix) <- METHODS

COP_matrix <- matrix(nrow=length(hc_out), ncol=8)
for(i in 2:length(hc_out)){
  COP_matrix[i,] <- cvi_eval_hc_out_k[[i]][12,]
}
colnames(COP_matrix) <- METHODS

DBstar_matrix <- matrix(nrow=length(hc_out), ncol=8)
for(i in 2:length(hc_out)){
  DBstar_matrix[i,] <- cvi_eval_hc_out_k[[i]][10,]
}
colnames(DBstar_matrix) <- METHODS

CH_matrix <- matrix(nrow=length(hc_out), ncol=8)
for(i in 2:length(hc_out)){
  CH_matrix[i,] <- cvi_eval_hc_out_k[[i]][8,]
}
colnames(CH_matrix) <- METHODS

SF_matrix <- matrix(nrow=length(hc_out), ncol=8)
for(i in 2:length(hc_out)){
  SF_matrix[i,] <- cvi_eval_hc_out_k[[i]][7,]
}
colnames(SF_matrix) <- METHODS

Sil_matrix[-1,]

which(Sil_matrix[-1,] == max(Sil_matrix[-1,]),
      arr.ind=TRUE) # ward.D2, complete

CH_matrix[-1,]

which(CH_matrix[-1,] == max(CH_matrix[-1,]),
      arr.ind=TRUE) # mcquitty

DBstar_matrix[-1,]

which(DBstar_matrix[-1,] == min(DBstar_matrix[-1,]),
      arr.ind=TRUE) # mcquitty

D_matrix[-1,]

which(D_matrix[-1,] == max(D_matrix[-1,]),
      arr.ind=TRUE) #single, average, centroid

COP_matrix[-1,]

which(COP_matrix[-1,] == min(COP_matrix[-1,]),
      arr.ind=TRUE) # ward.D

# SF_matrix[-1,]
# 
# which(SF_matrix[-1,] == max(SF_matrix[-1,]),
#       arr.ind=TRUE)

shpFiles <- list.files(path = here("mestiper", "data", "mohinora_usv7"),
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

DBstar_indCluster1 = which(hc_out[[3]][[6]]@cluster == 1)
DBstar_CLUSTER1 = vector("list", sum( hc_out[[3]][[6]]@cluster == 1 ))
for(i in DBstar_indCluster1){
  DBstar_CLUSTER1[[i]] = FPCA_trim[[i]]
}

DBstar_indCluster2 = which(hc_out[[3]][[6]]@cluster == 2)
DBstar_CLUSTER2 = vector("list", sum( hc_out[[3]][[6]]@cluster == 2 ))
for(i in DBstar_indCluster2){
  DBstar_CLUSTER2[[i]] = FPCA_trim[[i]]
}

DBstar_indCluster3 = which(hc_out[[3]][[6]]@cluster == 3)
DBstar_CLUSTER3 = vector("list", sum( hc_out[[3]][[6]]@cluster == 3 ))
for(i in DBstar_indCluster3){
  DBstar_CLUSTER3[[i]] = FPCA_trim[[i]]
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


