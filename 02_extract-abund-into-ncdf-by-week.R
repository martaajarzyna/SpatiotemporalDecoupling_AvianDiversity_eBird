###########################################################################
####### Spatiotemporal decoupling of avian taxonomic and functional diversity, Jarzyna & Stagge
####### 02: Extract species abundance into netCDF files
###########################################################################

###########################################################################
## Load Packages
###########################################################################
require(ncdf4)
#install.packages("ebirdst")
require(ebirdst)
require(raster)
require(terra)
library(sf)
library(rnaturalearth)


###########################################################################
## Add values to each week's ncdf file (VERY COMPUTATIONALLY INTENSIVE)
###########################################################################
sp.nam <- readRDS(file="Data/sp-abun-paths-2019.rds")

system.time(
for (j in 1:807){
cat(j)
spj <- load_raster("abundance", path = sp.nam[j,5]) 

for (i in 1:52) {
nc_file <- file.path(paste0("Data/abundmean_allsp_week_",i,".nc"))
nc <- nc_open(nc_file, write=TRUE) 
# Extract week  and convert to a matrix
spj.i <- spj[[i]]
spj.i.mat <- as.matrix(spj.i)
# Put in values
ncvar_put(nc, "abundance", spj.i.mat, start = c(1, 1, j), count = c(-1,-1,1))
# Close out the NCDF file
nc_close(nc)
}}
)




