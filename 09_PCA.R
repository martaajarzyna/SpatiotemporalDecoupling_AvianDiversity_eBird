###########################################################################
####### Spatiotemporal decoupling of avian taxonomic and functional diversity, Jarzyna & Stagge
####### 09: Principal Component Analysis (PCA) 
###########################################################################

#############################################################################
## load required packages
#############################################################################
require(tidyverse)
require(here)
require(svglite)
require(viridis)
require(lubridate)
require(scales)
require(stars)
require(raster)
require(rnaturalearth)
require(zoo)
require(grid)
require(maptools)
require(rgdal)
require(sf)
require(ggplot2)
require(dplyr)

### Use the correct select
select <- dplyr::select

### This is a frustrating bit about R. Two packages both want to use the command select. Make sure we are using the correct one
select <- dplyr::select

##################################################
##  PCA
##################################################
###########################################################################
## Set the Paths
###########################################################################
### Set here path
here_path <- here::here()

### Path for Data and Output
data_path <- file.path(here_path, "./data")
output_path <- file.path(here_path, "./output")

### Set up output folders
write_output_path <- file.path(output_path, "pca_analysis")
dir.create(write_output_path, recursive=TRUE, showWarnings = FALSE)

### Set up figure folder
write_figures_path <- file.path(output_path, "figures/pca_analysis")
dir.create(write_figures_path, recursive=TRUE, showWarnings = FALSE)

dir.create(file.path(write_figures_path, "png"), recursive=TRUE, showWarnings = FALSE)
dir.create(file.path(write_figures_path, "pdf"), recursive=TRUE, showWarnings = FALSE)
dir.create(file.path(write_figures_path, "svg"), recursive=TRUE, showWarnings = FALSE)

###########################################################################
## Set initial values
###########################################################################
### Set seed so everyone gets the same random numbers (reproducible example)
set.seed(7890)

###########################################################################
## Set projections
###########################################################################

### Original projection
crs_sinu <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

### CONUS Albers Equal Area
#new_proj = st_crs(5070)
### US National Atlas Equal Area
new_proj = st_crs(2163)

### Create a bounding box to crop
bbox <- st_bbox(c(xmin = -120, ymin = 17, xmax = -78, ymax = 41), crs = new_proj)

###########################################################################
## Read in data
###########################################################################

fdis_df <- readRDS(file.path(data_path, "index_matrix_FDis.rds"))
feve_df <- readRDS(file.path(data_path, "index_matrix_FEve.rds"))
fric_df <- readRDS(file.path(data_path, "index_matrix_FRic.rds"))
fricres_df <- readRDS(file.path(data_path, "index_matrix_FRic-res.rds"))
sr_df <- readRDS(file.path(data_path, "index_matrix_SR.rds"))

###########################################################################
## Check if there are any differences in index
###########################################################################

sum(feve_df$index_uniq != fdis_df$index_uniq)
sum(fric_df$index_uniq != fdis_df$index_uniq)
sum(feve_df$index_uniq != fric_df$index_uniq)
sum(sr_df$index_uniq != sr_df$index_uniq)

###########################################################################
## Reorganize data
###########################################################################
### Extract the index
index_c <- feve_df$index_uniq

### Extract and rotate
fdis_mat <- t(fdis_df[,6:57])
feve_mat <- t(feve_df[,6:57])
fric_mat <- t(fric_df[,6:57])
fricres_mat <- t(fricres_df[,6:57])
sr_mat <- t(sr_df[,6:57])

### Add column names to separate later
colnames(fdis_mat) <- paste0("fdis_",index_c)
colnames(feve_mat) <- paste0("feve_",index_c)
colnames(fric_mat) <- paste0("fric_",index_c)
colnames(fricres_mat) <- paste0("fricres_",index_c)
colnames(sr_mat) <- paste0("sr_",index_c)

### cbind to combine them
data_mat <- cbind(fdis_mat, feve_mat, fric_mat, fricres_mat, sr_mat)
dim(data_mat)

###########################################################################
## Clean up the loose data
###########################################################################
### Remove all the original matrices
rm(fdis_mat)
rm(feve_mat)
rm(fric_mat)
rm(fricres_mat)
rm(sr_mat)

### Create site_df
site_df <- feve_df %>%
  select(index_uniq, y_j, x_j, lon, lat)

### Save site_df
saveRDS(site_df, file = file.path(write_output_path, "site_df.RDS"))

### Remove the original data
rm(fdis_df)
rm(fric_df)
rm(fricres_df)
rm(feve_df)
rm(sr_df)


###########################################################################
## Prepare for PCA - Check for NAs
###########################################################################

### Find number of NAs
na_col_count <- apply(is.na(data_mat),2,sum)

### Create a dataframe
na_df <- data.frame(full_index = names(na_col_count), na_count = na_col_count)

### Plot to check NAs
plot_df <- site_df %>%
  mutate(var = "fdis") %>%
  bind_rows(site_df %>% mutate(var="fric")) %>%
  bind_rows(site_df %>% mutate(var="fricres")) %>%
  bind_rows(site_df %>% mutate(var="feve")) %>%
  bind_rows(site_df %>% mutate(var="sr")) %>%
  mutate(full_index = paste0(var, "_", index_uniq))

### Merge
plot_df <- plot_df %>%
  left_join(na_df , by = c("full_index" = "full_index"))

var_list <- unique(plot_df$var)

for (j in seq(1, length(var_list))){
  var_j <- var_list[[j]]
  
  ### Separate just this variable
  plot_temp <- plot_df %>%
    filter(var == var_j) %>%
    select(lon, lat, na_count)
  
  ### Create raster
  raster_temp <- rasterFromXYZ(plot_temp, crs= crs_sinu)
  
  ### Convert to stars and reproject
  raster_temp <- st_as_stars(raster_temp) %>%
    st_transform( new_proj)
  
  ### Crop the plot
  raster_temp <- raster_temp[bbox]
  
  ### Create plot
  p <- ggplot() %>%
    + geom_stars(data = raster_temp) %>%
    + scale_fill_viridis(name = paste0("NA Weeks\n", var_j), limits = c(1, 52))
  
  ### Save plot
  plot_name <- paste0(var_j, "_na_plot")
  
  ggsave(file.path(write_figures_path, paste0("png/", plot_name, ".png")), p,  width = 4, height = 6.5, dpi = 600)
  ggsave(file.path(write_figures_path, paste0("pdf/", plot_name, ".pdf")), p,  width = 4, height = 6.5)
  ggsave(file.path(write_figures_path, paste0("svg/", plot_name, ".svg")), p,  width = 4, height = 6.5)
  
}


###########################################################################
## Prepare for PCA - Remove NAs
###########################################################################
### Remove NA columns
na_col_test <- na_col_count > 0
data_mat <- data_mat[,!na_col_test]

### Remove constant values
const_col_test <- apply(data_mat,2,var) == 0
data_mat <- data_mat[,!const_col_test]

##################################################
##  Perform PCA
##################################################
##  Perform PCA with scaling using SVD
pca_fit <- prcomp(data_mat, retx=TRUE, scale=TRUE)

##  Show the loadings
pc_importance <- summary (pca_fit)
pc_importance

##################################################
##  Save objects
##################################################

saveRDS(data_mat, file = "Data/PCA-output/data_analysis.RDS")
saveRDS(pca_fit, file = "Data/PCA-output/pca_fit_analysis.RDS")



##################################################
##  Read PCA data and output in to plot diagnostics
##################################################
data_mat <- readRDS(file = "Data/PCA-output/data_analysis.RDS")
pca_fit <- readRDS(file = "Data/PCA-output/pca_fit_analysis.RDS")
pc_importance <- summary(pca_fit)
pc_importance


##################################################
##  Scree Plot
##################################################
#write_path <- write_figures_path
#file.path(write_figures_path, "fric")
#dir.create(write_path, recursive=TRUE, showWarnings = FALSE)

### Create Scree plot
Scree <- data.frame(cbind(Component=seq(1,dim(pc_importance$importance)[2]),t(pc_importance$importance)))
Scree$EigenVal <- Scree$Standard.deviation^2
Scree_portion <- Scree[1:52,]

###
p <- ggplot(Scree_portion, aes(x = Component, y  = Proportion.of.Variance*100)) %>%
	+ geom_line(size=2) %>%
  + geom_point(alpha=1,size=4, pch=16) %>%
	+ scale_x_continuous(name = "Principal Component", limits=c(1,12),breaks=c(seq(1,12,1))) %>%
	+ scale_y_continuous(name = "Proportion of variance (%)", limits=c(0,50),breaks=c(0,10,20,30,40,50)) %>%
  + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
        panel.border = element_rect(fill=NA, colour = "white", size=1),
        axis.line = element_line(color = 'black', size=1.5),
        plot.title = element_text(size=15, vjust=2, family="sans"),
        axis.text.x = element_text(colour='black',size=22),
        axis.text.y = element_text(colour='black',size=22),
        axis.title.x = element_text(colour='black',size=27),
        axis.title.y = element_text(colour='black',size=27),
        axis.ticks = element_line(color = 'black', size=1.5),
        axis.ticks.length=unit(0.3,"cm"),
        legend.position="none",
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p

ggsave(p, file="Graphics/prop_var_12pcs.png", width=8, height=6, dpi=600)


p <- ggplot(Scree_portion, aes(x = Component, y  = Proportion.of.Variance*100)) %>%
  + geom_line(size=2) %>%
  + geom_point(alpha=1,size=4, pch=16) %>%
  + scale_x_continuous(name = "Principal Component", limits=c(1,52),breaks=c(1,10,20,30,40,50)) %>%
  + scale_y_continuous(name = "Proportion of variance (%)", limits=c(0,50),breaks=c(0,10,20,30,40,50)) %>%
  + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
          panel.border = element_rect(fill=NA, colour = "white", size=1),
          axis.line = element_line(color = 'black', size=1.5),
          plot.title = element_text(size=15, vjust=2, family="sans"),
          axis.text.x = element_text(colour='black',size=22),
          axis.text.y = element_text(colour='black',size=22),
          axis.title.x = element_text(colour='black',size=27),
          axis.title.y = element_text(colour='black',size=27),
          axis.ticks = element_line(color = 'black', size=1.5),
          axis.ticks.length=unit(0.3,"cm"),
          legend.position="none",
          legend.text=element_text(size=20),
          legend.title=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
p

ggsave(p, file="Graphics/prop_var_allpcs.png", width=8, height=6, dpi=600)



###
p <- ggplot(Scree_portion, aes(x = Component, y  = Cumulative.Proportion*100)) %>%
  + geom_line(size=2) %>%
  + geom_point(alpha=1,size=4, pch=16) %>%
  + scale_x_continuous(name = "Principal Component", limits=c(1,12),breaks=seq(1,12,1)) %>%
  + scale_y_continuous(name = "Cumulative variance (%)", limits=c(0,100),breaks=seq(0,100,10)) %>%
  + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
          panel.border = element_rect(fill=NA, colour = "white", size=1),
          axis.line = element_line(color = 'black', size=1.5),
          plot.title = element_text(size=15, vjust=2, family="sans"),
          axis.text.x = element_text(colour='black',size=22),
          axis.text.y = element_text(colour='black',size=22),
          axis.title.x = element_text(colour='black',size=27),
          axis.title.y = element_text(colour='black',size=27),
          axis.ticks = element_line(color = 'black', size=1.5),
          axis.ticks.length=unit(0.3,"cm"),
          legend.position="none",
          legend.text=element_text(size=20),
          legend.title=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
p

ggsave(p, file="Graphics/cum_var_12pcs.png", width=8, height=6, dpi=600)


###
p <- ggplot(Scree_portion, aes(x = Component, y  = Cumulative.Proportion*100)) %>%
  + geom_line(size=2) %>%
  + geom_point(alpha=1,size=4, pch=16) %>%
  + scale_x_continuous(name = "Principal Component", limits=c(1,52),breaks=c(1,10,20,30,40,50)) %>%
  + scale_y_continuous(name = "Cumulative variance (%)", limits=c(0,100),breaks=seq(0,100,10)) %>%
  + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
          panel.border = element_rect(fill=NA, colour = "white", size=1),
          axis.line = element_line(color = 'black', size=1.5),
          plot.title = element_text(size=15, vjust=2, family="sans"),
          axis.text.x = element_text(colour='black',size=22),
          axis.text.y = element_text(colour='black',size=22),
          axis.title.x = element_text(colour='black',size=27),
          axis.title.y = element_text(colour='black',size=27),
          axis.ticks = element_line(color = 'black', size=1.5),
          axis.ticks.length=unit(0.3,"cm"),
          legend.position="none",
          legend.text=element_text(size=20),
          legend.title=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
p

ggsave(p, file="Graphics/cum_var_allpcs.png", width=8, height=6, dpi=600)




##################################################
##  Plot Scores Diagrams
##################################################
pc_list <- seq(1,3) #only for the first 3 PCs
var_list <- c( "fdis", "feve", "fric", "fricres", "sr")



### Loop through scores
for(i in seq(1, 52)){

	### Create loadings
	plot_df <- data.frame(week = seq(1,52), loading = pca_fit$x[,i])
	#plot_df$score <- plot_df$score*(-1)   ##do this only for PC2
	
	if (i == 2){
	  p <- ggplot(plot_df, aes(x=week, y=(loading*(-1)))) %>%
	    + geom_line(size=2) %>%
	    + geom_point(alpha=1,size=4, pch=16) %>%
	    + geom_hline(yintercept = 0, linetype = "dashed") %>%
	    + scale_x_continuous(name = "Time step (week)") %>%
	    + scale_y_continuous(name = paste0("PC",i," score")) %>%
	    + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
	            panel.border = element_rect(fill=NA, colour = "white", size=1),
	            axis.line = element_line(color = 'black', size=1.5),
	            plot.title = element_text(size=15, vjust=2, family="sans"),
	            axis.text.x = element_text(colour='black',size=22),
	            axis.text.y = element_text(colour='black',size=22),
	            axis.title.x = element_text(colour='black',size=27),
	            axis.title.y = element_text(colour='black',size=27),
	            axis.ticks = element_line(color = 'black', size=1.5),
	            axis.ticks.length=unit(0.3,"cm"),
	            legend.position="none",
	            legend.text=element_text(size=20),
	            legend.title=element_blank(),
	            panel.grid.major = element_blank(), 
	            panel.grid.minor = element_blank())
	} else {
	  p <- ggplot(plot_df, aes(x=week, y=loading)) %>%
	  + geom_line(size=2) %>%
	  + geom_point(alpha=1,size=4, pch=16) %>%
	  + geom_hline(yintercept = 0, linetype = "dashed") %>%
	  + scale_x_continuous(name = "Time step (week)") %>%
	  + scale_y_continuous(name = paste0("PC",i," score")) %>%
	  + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
	          panel.border = element_rect(fill=NA, colour = "white", size=1),
	          axis.line = element_line(color = 'black', size=1.5),
	          plot.title = element_text(size=15, vjust=2, family="sans"),
	          axis.text.x = element_text(colour='black',size=22),
	          axis.text.y = element_text(colour='black',size=22),
	          axis.title.x = element_text(colour='black',size=27),
	          axis.title.y = element_text(colour='black',size=27),
	          axis.ticks = element_line(color = 'black', size=1.5),
	          axis.ticks.length=unit(0.3,"cm"),
	          legend.position="none",
	          legend.text=element_text(size=20),
	          legend.title=element_blank(),
	          panel.grid.major = element_blank(), 
	          panel.grid.minor = element_blank())
	}
	
	### Save plot
	ggsave(p, file=paste0("Graphics/pc_",i,"_score-diagrams-figS2.png"), width=8, height=6, dpi=600)

}
