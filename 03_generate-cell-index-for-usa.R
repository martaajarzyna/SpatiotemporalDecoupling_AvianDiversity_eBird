###########################################################################
####### Spatiotemporal decoupling of avian taxonomic and functional diversity, Jarzyna & Stagge
####### 03: Generate cell index matrix for USA--to only analyze data for contiguous USA
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
require(rgdal)
require(tidyverse)

### Only for super hi res or states
require(rnaturalearthhires)
require(rnaturalearthdata)
require(stars)

###########################################################################
## Read in example netcdf file and overlay with us shp
###########################################################################
#i=week, i=26 for testing
i=26
nc_sp_file <- file.path(paste0("Data/abundmean_allsp_week_",i,".nc"))
#nc_sp <- read_ncdf(nc_sp_file, var = "abubdance", ncsub =cbind(start = c(1,1,1), count = c(-1,-1,1)))
nc_sp <- nc_open(nc_sp_file)

abun_mat <- ncvar_get(nc_sp, "abundance", start = c(1,1,3), count = c(-1,-1,1))
x_vec <- ncvar_get(nc_sp, "x")
y_vec <- ncvar_get(nc_sp, "y")

abun_raster <- raster(abun_mat, xmn = min(x_vec), xmx = max(x_vec), ymn = min(y_vec), ymx=max(y_vec))
projection(abun_raster) <- crs("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
plot(abun_raster)

### overlay shp to see if we have the same projection
### All countries in world as shapefiles. There are a bunch of ways you can filter this
world <- ne_countries(scale = "medium", returnclass = "sf")
head(world)

### Get a US shapefile
usa <- ne_countries(country = 'United States of America', scale = "medium", returnclass = "sf")

### You can do a custom download from the naturalearth website
### https://www.naturalearthdata.com/downloads/50m-cultural-vectors/50m-admin-0-details/
### There's a ton of stuff on this website. The subunit file has the lower 48

world_by_subunits <- ne_download(scale = 'medium', type = "map_subunits",  returnclass = "sf")
### You can make this higher resolution if you need it

### There are 3 subunits (48, alaska, hawaii)
world_by_subunits %>% filter(ADMIN == "United States of America")

### Get the lower 48
lower_48 <- world_by_subunits %>% filter(SUBUNIT == "United States")
plot(lower_48)
ggplot(lower_48) + geom_sf()
lower_48

### Get the outline
lower_48_outline <- as(st_geometry(lower_48), Class="Spatial")
crs(lower_48_outline)
us.shp <- spTransform(lower_48_outline, CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"))

plot(us.shp, add = TRUE)



###########################################################################
###  Create a row and column raster
###########################################################################
row_vec <- seq(1, length(y_vec))
col_vec <- seq(1,length(x_vec))
### Create a raster of row numbers
row_mat <- matrix(rep(row_vec, times = length(x_vec)), ncol = length(x_vec), nrow = length(y_vec))
row_raster <- raster(row_mat, xmn = min(x_vec), xmx = max(x_vec), ymn = min(y_vec), ymx = max(y_vec))
crs(row_raster) <- crs(abun_raster)

plot(row_raster)

### Create a raster of column numbers
col_mat <- matrix(rep(col_vec, each = length(y_vec)), ncol = length(x_vec), nrow = length(y_vec))
col_raster <- raster(col_mat, xmn = min(x_vec), xmx = max(x_vec), ymn = min(y_vec), ymx = max(y_vec))
crs(col_raster) <- crs(abun_raster)
plot(col_raster)


###########################################################################
###  Convert to stars object and crop
###########################################################################
### Project lower 48 to sinusoidal
### Then add a buffer. Its in meters
lower_48_buffer <- lower_48 %>%
  st_transform( st_crs(abun_raster)) %>%
  st_buffer( dist = 50000)

### Convert to stars object
require(stars)
row_stars <- row_raster %>%
  st_as_stars()
col_stars <- col_raster %>%
  st_as_stars()

### Subset to the buffered US
row_usa = row_stars[lower_48_buffer]
col_usa = col_stars[lower_48_buffer]

######################################
### Create an index dataframe
######################################
index_df <- data.frame(row = c(row_usa %>% pull(1)), col = c(col_usa %>% pull(1))) %>%
  drop_na()

index_df <- index_df %>%
  mutate(lon = x_vec[col]) %>%
  mutate(lat = y_vec[row]) %>%
  mutate(i = seq(1,dim(index_df)[1])) %>%
  dplyr::select(i, row, col, lon, lat)

### Check plot
ggplot(index_df, aes(x= lon, y=lat, fill = i)) + geom_raster()

### This is the object that you can loop over (I hope)
head(index_df)
dim(index_df)

saveRDS(index_df, file="Data/US_cell_index.rds")

nc_close(nc_sp) #don't forget to close the netCDF file


