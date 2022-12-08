###########################################################################
####### Spatiotemporal decoupling of avian taxonomic and functional diversity, Jarzyna & Stagge
####### 11: Cluster analysis
###########################################################################

#############################################################################
## load required packages
#############################################################################
require(stats)
require(tidyverse)
require(here)
require(svglite)
require(viridis)
require(lubridate)
require(scales)
require(stars)
require(raster)
require(rnaturalearth)

require(grid)
require(maptools)
require(rgdal)
require(sf)
require(ggplot2)
require(dplyr)
require(cluster)
require(fastcluster)
require(ggpubr)

### Use the correct select
select <- dplyr::select

#############################################################################
## prep data for cluster analysis
#############################################################################
data_mat <- readRDS(file = "Data/PCA-output/data_analysis.RDS") #actual measurements for diversity
pca_fit <- readRDS(file = "Data/PCA-output/pca_fit_analysis.RDS") #results from the PCA
#Each site has it’s own long term mean (center) and standard deviation (scale) throughout the year for each variable. 
#PC score is a seasonal pattern shared by all sites/variables that can be turned up or down (or flipped)
#Each site and variable has a unique loading 
#Score * loading gives the number of standard deviations from the mean for a given site and variable


#############################################################################
## get the loadings and prep data for cluster analysis
#############################################################################
data_mat <- t(data_mat)
data_mat <- as.data.frame(data_mat)
data_mat$index_uniq_full <- rownames(data_mat)

loadings.full <- as.data.frame(data_mat$index_uniq_full)
colnames(loadings.full) <- "index_uniq_full"

for (i in 1:3){
  #i=1 ##PC 1
  loading_temp <- data.frame(index_uniq=rownames(pca_fit$rotation), load=pca_fit$rotation[,i], center = pca_fit$center, scale = pca_fit$scale)
  colnames(loading_temp) <- c("index_uniq_full",paste0("load_PC",i),"center","scale")
  loading_temp <- loading_temp %>% select(index_uniq_full, paste0("load_PC",i))
  
  loadings.full <- loadings.full %>% 
    full_join(loading_temp, by = "index_uniq_full")
}

loadings.full$load_PC2 <- loadings.full$load_PC2*(-1) #flip the loading for PC2 (it is easier to interpret, but the math does not change)

loadings.full <- loadings.full %>%
  separate(col = "index_uniq_full", into = c("a", "b", "c"), sep = "_") %>%
  mutate(index_uniq = paste0(b, "_", c)) %>%
  mutate(var = a)  %>%
  mutate(index_uniq_full = paste0(a,"_", b, "_", c)) %>%
  select(-a, -b, -c)

loadings.full <- loadings.full %>%
  select(index_uniq_full,index_uniq,var,load_PC1,load_PC2,load_PC3)

saveRDS(loadings.full, file=paste0("Data/PCA-postprocessing-data/loadings-indices-all.rds"))

#############################################################################
## separate loadings by index
#############################################################################
var_list <- c( "fdis", "feve", "fric", "fricres", "sr")

index.un <- as.data.frame(unique(loadings.full$index_uniq))
colnames(index.un) <- "index_uniq"

for (i in 1:length(var_list)){
  i.load <- loadings.full %>%
    filter(var == var_list[i]) %>%
    select(index_uniq,load_PC1,load_PC2,load_PC3)
  
  colnames(i.load) <- c("index_uniq",paste0("load_",var_list[i],"_PC1"),paste0("load_",var_list[i],"_PC2"),paste0("load_",var_list[i],"_PC3"))
  
  index.un <- index.un %>% 
    full_join(i.load, by = "index_uniq")
}

saveRDS(index.un, file=paste0("Data/PCA-postprocessing-data/loadings-indices-all-v2.rds"))


#############################################################################
## Run k-means clustering
#############################################################################

### We first determine what is the correct number of clusters using a sample
## (to computationally expensive to run this on the entire sample given ~1M grid cells for each variable)
for (j in 1:1){
  loadings.full <- readRDS(file=paste0("Data/PCA-postprocessing-data/loadings-indices-all-v2.rds"))
  #note: loadings for PC2 in this matrix are already flipped in terms of the sign
  
  ##################################################
  ##  Process distance matrix for the first three PCs, for cFRic, FEve, and FDis only for now
  ##################################################
  loadings.full <- loadings.full %>% select(index_uniq,load_sr_PC1,load_fric_PC1,load_fricres_PC1,load_feve_PC1,load_fdis_PC1,
                                            load_sr_PC2,load_fric_PC2,load_fricres_PC2,load_feve_PC2,load_fdis_PC2,
                                            load_sr_PC3,load_fric_PC3,load_fricres_PC3,load_feve_PC3,load_fdis_PC3)
  loadings.full <- na.omit(loadings.full) 
  loadings.full <- sample_n(loadings.full, 20000)   
  ids <- loadings.full[,1]
  saveRDS(ids, file=paste0("Data/PCA-postprocessing-data/kmeans-ids-random-sample",j,".rds"))
  
  loadings.full <- loadings.full %>% select(load_sr_PC1,load_fric_PC1,load_fricres_PC1,load_feve_PC1,load_fdis_PC1,
                                            load_sr_PC2,load_fric_PC2,load_fricres_PC2,load_feve_PC2,load_fdis_PC2,
                                            load_sr_PC3,load_fric_PC3,load_fricres_PC3,load_feve_PC3,load_fdis_PC3)
  
  
  ##################################################
  ##  Process Clustering - K-Means
  ##################################################
  cluster.list <- 2:15
  clusters <- matrix(NA,dim(loadings.full)[1],length(cluster.list))
  colnames(clusters) <- paste("k_",cluster.list, sep="")
  
  cluster_goodness <- matrix(NA,length(cluster.list), 5)
  rownames(cluster_goodness) <- paste("k_",cluster.list, sep="")
  colnames(cluster_goodness) <- c("k", "totss","betweenss","perc_betw","sil_width")
  
  ##  Calculate the k-means clustering
  for (k in cluster.list) {
    clust_result <- kmeansruns(loadings.full, krange=k, runs=10,criterion="asw")
    clusters[,k-1]=as.factor(clust_result$cluster)
    cluster_goodness[k-1,] <- c(k, clust_result$totss,clust_result$betweenss,100*(clust_result$betweenss/clust_result$totss),clust_result$crit[k])
    #cat(k)
  }
  
  here <- cbind(ids,loadings.full,clusters)
  saveRDS(here, file=paste0("Data/PCA-postprocessing-data/kmeans-clusters-random-5vars-sample",j,".rds"))
  saveRDS(cluster_goodness, file=paste0("Data/PCA-postprocessing-data/kmeans-clusters-goodness-random-5vars-sample",j,".rds"))
  
  cluster_goodness <- as.data.frame(cluster_goodness)
  p <- ggplot(cluster_goodness, aes(k,sil_width))
  #p <- p + theme_pub_majgrid()
  p <- p + geom_line()
  p <- p + geom_point()
  
  p <- p + xlab("Number of Clusters")  #insert the x-axis title
  p <- p + ylab("Avg. Silhouette Width")  #insert the y-axis title
  p <- p + scale_x_continuous(breaks=seq(0,15,2))
  
  ggsave(p, file=paste0("Data/Graphics/kmeans-clusters-goodness-5vars-sample",j,".png"), width=8, height=6, dpi=600)
}


### Based on the sample results for 20000 locations, the appropriate number of clusters is 7
### We now redo this knnowing 7 is the correct number of clusters
set.seed(1)
loadings.full <- readRDS(file=paste0("Data/PCA-postprocessing-data/loadings-indices-all-v2.rds"))
loadings.full <- loadings.full %>% select(index_uniq,load_sr_PC1,load_fric_PC1,load_fricres_PC1,load_feve_PC1,load_fdis_PC1,
                                          load_sr_PC2,load_fric_PC2,load_fricres_PC2,load_feve_PC2,load_fdis_PC2,
                                          load_sr_PC3,load_fric_PC3,load_fricres_PC3,load_feve_PC3,load_fdis_PC3)
loadings.full <- na.omit(loadings.full) 
#loadings.full <- sample_n(loadings.full, 1000)   

set.seed(1)
loadings.full <- as.data.frame(loadings.full)
ind <- sample(nrow(loadings.full), 20000)

loadings.full[["train"]] <- TRUE
loadings.full[["train"]][ind] <- FALSE

cl1 = kcca(loadings.full[loadings.full[["train"]]==TRUE, 2:16], k=7, kccaFamily("kmeans"))
saveRDS(cl1, file=paste0("Data/PCA-postprocessing-data/clustering7-5vars-kmeans-object.rds"))    

pred_train <- predict(cl1)
pred_test <- predict(cl1, newdata=loadings.full[loadings.full[["train"]]==FALSE, 2:16])

clusters.val <- c(pred_train,pred_test)
clusters.val <- as.data.frame(clusters.val)
clusters.val <- clusters.val[order(as.numeric(rownames(clusters.val))),,drop=FALSE]
loadings.full <- loadings.full %>%
  bind_cols(clusters.val)

saveRDS(loadings.full, file=paste0("Data/PCA-postprocessing-data/loadings-indices-all-v2-with7clusters-5vars.rds"))  



