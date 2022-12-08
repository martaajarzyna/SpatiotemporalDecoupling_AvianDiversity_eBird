###########################################################################
####### Spatiotemporal decoupling of avian taxonomic and functional diversity, Jarzyna & Stagge
####### 00: Download eBird data
###########################################################################


###########################################################################
## Load Packages
###########################################################################
require(ncdf4)
install.packages("ebirdst")
require(ebirdst)
require(raster)
require(terra)
library(viridis)
library(sf)
library(rnaturalearth)

set_ebirdst_access_key("qs1ap2m3uah") #this will change every time you download the new data, will need to make a separate request

###########################################################################
## Check the Ebird data
###########################################################################
sp_path <- ebirdst_download(species = "example_data")
abunds <- load_raster("abundance", path = sp_path)
dim(ebirdst_runs) #A full list of species and run names can be found in the ebirdst_runs dataframe
head(ebirdst_runs)

### Check the names of species and get the codes
spnames <- ebirdst_runs[,1:4]
spcodes <- spnames$species_code
  
### Set parameters to run through 
nweek <- 52
n_species <- length(spcodes)

###########################################################################
### Download the data first
###########################################################################
sp.nam <- as.data.frame(spnames)
sp.nam$path <- NA

for (j in 1:nrow(sp.nam)){ 
  cat(j)
  spj <- ebirdst_download(species=paste0(spnames[j,3]), path="Data/ind-species-abundance-2019/")
  sp.nam[j,5] <- spj
}

saveRDS(sp.nam,file="Data/sp-abun-paths-2019.rds")
