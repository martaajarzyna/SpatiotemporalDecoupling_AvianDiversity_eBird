###########################################################################
####### Spatiotemporal decoupling of avian taxonomic and functional diversity, Jarzyna & Stagge
####### 10: Plot PCA results: example diagnostics, plus scores diagrams and loading maps
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

require(grid)
require(maptools)
require(rgdal)
require(sf)
require(ggplot2)
require(dplyr)

### This is a frustrating bit about R. Two packages both want to use the command select. Make sure we are using the correct one
select <- dplyr::select

##################################################
##  Read PCA data and output in
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
        axis.title.x = element_text(colour='black',size=22),
        axis.title.y = element_text(colour='black',size=22),
        axis.ticks = element_line(color = 'black', size=1.5),
        axis.ticks.length=unit(0.3,"cm"),
        legend.position="none",
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p

ggsave(p, file="Graphics/prop_var.png", width=8, height=6, dpi=600)

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
          axis.title.x = element_text(colour='black',size=22),
          axis.title.y = element_text(colour='black',size=22),
          axis.ticks = element_line(color = 'black', size=1.5),
          axis.ticks.length=unit(0.3,"cm"),
          legend.position="none",
          legend.text=element_text(size=20),
          legend.title=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
p

ggsave(p, file="Graphics/cum_var.png", width=8, height=6, dpi=600)


###########################################################################
###  Run this part on OSC - too much memory required to create a map
###########################################################################

###########################################################################
###  Download background data
###########################################################################
#us <- readOGR(".", layer = "states")
#us <- fortify(us)
#index_sf = st_as_sf(index_df, coords = c("Xcoord", "Ycoord"), crs = 4326 )

usa_states <- ne_states(country = 'United States of America', returnclass = 'sf')
usa_states <- fortify(usa_states)
crs_sinu <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
#crs_sinu <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"
new_proj = st_crs(2163)
### Create a bounding box to crop
#bbox <- st_bbox(c(xmin = -126, ymin = 17, xmax = -108, ymax = 41), crs = new_proj)

iso_list <- usa_states$iso_3166_2
iso_list <- iso_list[! iso_list %in% c('US-AK', 'US-HI')]
region_test <- usa_states$iso_3166_2 %in% iso_list

region_states <- usa_states[region_test,]

### Merge into a single polygons
region_states <- region_states %>%
  st_transform(new_proj) 

#plot(region_states)
#plot(usa_states)

##################################################
##  Plot Loading Maps
##################################################
pc_list <- seq(1,3) #only for the first 3 PCs
var_list <- c( "fdis", "feve", "fric", "fricres", "sr")

### Loop through maps
for (i in seq(1, length(pc_list))){
  
  ### Create temporary loading dataframe
  loading_temp <- data.frame(index_uniq=rownames(pca_fit$rotation), load=pca_fit$rotation[,i])
  
  loading_temp <- loading_temp %>%
    separate( col = "index_uniq", into = c("a", "b", "c"), sep = "_") %>%
    mutate(index_uniq = paste0(b, "_", c)) %>%
    mutate(var = a) %>%
    select(var, index_uniq, load)
  
  ### Create plotting dataframe
  plot_df <- site_df %>%
    expand_grid(data.frame(var = var_list)) %>%
    mutate(full_index = paste0(var, "_", index_uniq))
  
  plot_df <- plot_df %>%
    full_join(loading_temp, by = c("var", "index_uniq"))
  
  ### Find the limits
  plot_lims <- c(-1,1)*max(abs(plot_df$load), na.rm=TRUE)
  
  for (j in seq(1, length(var_list))){
    var_j <- var_list[[j]]
    
    ### Create raster
    raster_temp <- rasterFromXYZ(plot_df %>% filter(var == var_j) %>% select(lon, lat, load), crs= crs_sinu)
    
    ### Convert to stars and reproject
    raster_temp <- st_as_stars(raster_temp) %>%
      st_transform( new_proj)
    
    ### Crop the plot
    #raster_temp_crop <- crop(raster_temp, region_states)
    #raster_temp <- raster_temp[bbox]
    raster_temp.c <- st_crop(raster_temp, region_states, crop=FALSE)
    
    ### Create plot
    p <- ggplot() +
      geom_stars(data = raster_temp) +
      scale_fill_distiller(name = paste0("Score\n", toupper(var_j), "\nPC ",i), palette = "RdBu" , limits = plot_lims) +
      geom_sf(data = region_states, fill = NA, alpha = 1, size=0.2) +
      xlab("Longitude") +
      ylab("Latitude") +
      theme_bw(8) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    ### Save plot
    ggsave(p, file=paste0("Data/Graphics/", var_j, "_pc_",i,"_loading_map.png"), width=6, height=4, dpi=600)
  }
  
}



##################################################
##  Plot PC scores
##################################################

### Loop through scores
for(i in seq(1, 3)){

	#pc_folder <- file.path(write_path, paste0("pc_", i))
	#dir.create(pc_folder, recursive=TRUE, showWarnings = FALSE)

	### Create loadings
	plot_df <- data.frame(week = seq(1,52), score = pca_fit$x[,i])
	#plot_df$score <- plot_df$score*(-1)   ##do this only for PC2
	
	### Create plot
	p <- ggplot(plot_df, aes(x=week, y=score)) %>%
	  + geom_line(size=2) %>%
	  + geom_point(alpha=1,size=4, pch=16) %>%
	  + geom_hline(yintercept = 0, linetype = "dashed") %>%
	  + scale_x_continuous(name = "Time step (week)") %>%
	  + scale_y_continuous(name = "Score") %>%
	  + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
	          panel.border = element_rect(fill=NA, colour = "white", size=1),
	          axis.line = element_line(color = 'black', size=1.5),
	          plot.title = element_text(size=15, vjust=2, family="sans"),
	          axis.text.x = element_text(colour='black',size=22),
	          axis.text.y = element_text(colour='black',size=22),
	          axis.title.x = element_text(colour='black',size=22),
	          axis.title.y = element_text(colour='black',size=22),
	          axis.ticks = element_line(color = 'black', size=1.5),
	          axis.ticks.length=unit(0.3,"cm"),
	          legend.position="none",
	          legend.text=element_text(size=20),
	          legend.title=element_blank(),
	          panel.grid.major = element_blank(), 
	          panel.grid.minor = element_blank())
	p
	
	### Save plot
	ggsave(p, file=paste0("Graphics/pc_",i,"_score.png"), width=8, height=6, dpi=600)

}

