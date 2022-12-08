###########################################################################
####### Spatiotemporal decoupling of avian taxonomic and functional diversity, Jarzyna & Stagge
####### 06: Quantify avian diversity 
###########################################################################


###########################################################################
## Load Packages
###########################################################################
require(ncdf4)
require(FD)  #calculating multivariate and single trait distance-based FD indices
require(picante)  #calculating PD-equivalent indices and trait phylogenetic signal
require(dplyr)
require(funrar)
require(gawdis)

###########################################################################
## Get traits and species names
###########################################################################
## Read in eBird species list and EltonTraits species trait designations
sp.nam <- readRDS(file="Data/sp-abun-paths-2019.rds")
traits <- read.csv("Traits/Master_Funcdiv_ebirdnomenclature.csv", header=TRUE, row.names=1)
colnames(traits) <- c("WJSpecID","scientific_name","common_name",colnames(traits[,4:ncol(traits)]))

###########################################################################
## Remove pelagic species and calculate corrected gower distance for all species
###########################################################################
## Remove pelagic and nocturnal species, also remove one record not id'ed to species, that leaves us with 566 species
traits <- traits %>%
  filter(!(PelagicSpecialist == 1)) %>%
  filter(!(scientific_name == "Melanitta deglandi/stejnegeri")) #remove one remaining ambiguous id
rownames(traits) <- traits$scientific_name

###########################################################################
## Calculate distance matrix for all retained species
###########################################################################
## select traits: diet, foraging stratum, body mass, and nocturnality
traits <- traits[,c(4:18,20)] 

#designate trait groupings for gawdis (groups = different numbers)
trait.groups = c(rep(1, 7), rep(2, 7), 3, 4)
#identify the trait groups that represent fuzzy conding or dummy variables - numbers matching the group values above.
fuzzy.groups = c(1,2)  #diet, foraging stratum

#check trait groups
trait.names <- colnames(traits)
cbind.data.frame(trait.names, trait.groups)

#create functional dissimilarity matrix with (more) equal contribution of traits and trait groups (fuzzy coded categories and mass/head-body length (corr = 95%)) via gawdis
#must use optimized to avoid the negative weights of a few traits
#you only need to do this once, then just read in the trait space for all calculations
traits.dist <- gawdis(traits, w.type = "optimized", groups = trait.groups, fuzzy = fuzzy.groups, opti.maxiter = 300)
attr(traits.dist, "weights")
saveRDS(traits.dist,file="Data/traits-distance.rds")


###########################################################################
## Calculate SR & multivariate FD
###########################################################################
### Read US cell index: these will be the cell over which to loop and calculate FD (we're calculating FD only for the continental USA)
### this needs to be done only ONCE. once done, you can just read in the grids_subsets data frame
index_df <- readRDS(file="Data/US_cell_index.rds")
all_grids <- expand.grid(y_j = seq(1,5630), x_j = seq(1,7074)) #almost 40mln rows
### create unique columns in both data frames
all_grids <- all_grids %>% unite(index_uniq, y_j, x_j, sep = "_", remove = FALSE)
index_df <- index_df %>%unite(index_uniq, row, col, sep = "_", remove = FALSE)
grids_subset <- all_grids %>% filter(index_uniq %in% index_df$index_uniq)

#resave with new index
saveRDS(index_df, file="Data/US_cell_index.rds")
saveRDS(all_grids, file="Data/AllGrids_cell_index.rds")
saveRDS(grids_subset, file="Data/USGrids_cell_index.rds") #this is what we will be using for the calculation
#grids_subset <- readRDS(file="Data/USGrids_cell_index.rds")

### This is a very computatioonally intensive step
### We advise that the user splits the dataset (grids_subset) into chunks that can all run concurrently
### on a super computer
### in total, this took 11,000 hours
for (i in 1:52) {
    #open ncdf file to write from
    nc_sp_file <- file.path(paste0("Data/abundmean_allsp_week_",i,".nc"))
    nc_sp <- nc_open(nc_sp_file) 
    
    #open ncdf files to write to
    nc_fd_file <- file.path(paste0("Data/fd_week",i,".nc"))
    nc_fd <- nc_open(nc_fd_file, write=T) 
    
    ### grab the dimensions for future processing
    species_vec <- ncvar_get(nc_sp, "species")
    x_vec <- ncvar_get(nc_sp, "x")
    y_vec <- ncvar_get(nc_sp, "y")
    
    for (j in seq(1, dim(grids_subset)[1])){
      #for (j in c(100000:100010)){ ##test on a suubset of cells
      x_j <- grids_subset$x_j[j]
      y_j <- grids_subset$y_j[j]
      lon_j <- x_vec[x_j]
      lat_j <- y_vec[y_j]
      
      ### This is the command that extracts the time series
      abundance_j <- ncvar_get(nc_sp, "abundance", start = c(y_j,x_j,1), count = c(1,1,-1)) 
      #abundance_df <- data.frame(lon = lon_j, lat = lat_j, species = species_vec, abundance = abundance_j)
      #head(abundance_df)
      
      abundance_j <- t(as.matrix(abundance_j))
      colnames(abundance_j) <- sp.nam[,3]
      ### retain only select species
      abundance_j <- t(as.matrix(abundance_j[,colnames(abundance_j) %in% rownames(traits)]))
      
      ### if all species NaN : then automatically assign NaN to all metrics
      ### if some species NaN, but others not: assign 0 to NaN, and proceed with analysis (those will be excluded)
      ### if all species >=0, proceed with the analysis 
      ### it is computationoally prohibitive to ruun the dbFD fucntion on the entire matrix (each has ~1M cells)
      ### so instead we run a calculation for each grid cell separately
      ### to make this work we're rbinding with a vector of 1s because FD package does not take matrix where some species are present at all (all species must be present at at least 1 ste)
      ### the only way to ensure that is to have a dummy row of 1s which will later be ignored
      ### I checked and the values are the same regrdaless of whether the dummy site is there
      abundance_j <- rbind(abundance_j, rep(1,ncol(abundance_j)))
      #abundance_j[,1:20]
      ab_j <- abundance_j[1,!is.na(abundance_j[1,])]

      if (sum(is.na(abundance_j[1,])) == nrow(traits)) { #when we have no data for any species
        next #NaN is already there
        ### here place NaN into the ncdf grid cell, and all_grids_metrics columns
        #ncvar_put(nc_fd, "SR", NaN, start = c(y_j,x_j), count = c(1,1))
        #ncvar_put(nc_fd, "FRic", NaN, start = c(y_j,x_j), count = c(1,1))
        #ncvar_put(nc_fd, "FEve", NaN, start = c(y_j,x_j), count = c(1,1))
        #ncvar_put(nc_fd, "FDis", NaN, start = c(y_j,x_j), count = c(1,1))
      } else if (sum(is.na(abundance_j[1,])) < nrow(traits) & sum(abundance_j[1,], na.rm=T) == 0) { #when there is data but there are no species, all functional indices will be 0
        ncvar_put(nc_fd, "SR", 0, start = c(y_j,x_j), count = c(1,1))
        ncvar_put(nc_fd, "FRic", 0, start = c(y_j,x_j), count = c(1,1))
        ncvar_put(nc_fd, "FEve", 0, start = c(y_j,x_j), count = c(1,1))
        ncvar_put(nc_fd, "FDis", 0, start = c(y_j,x_j), count = c(1,1))
      } else if (length(ab_j[ab_j>0]) < 4) { #when there are < 4 species, dbFD will not run (fewer species than traits)
        #ncvar_put(nc_fd, "SR", length(abundance_j[1,abundance_j[1,]>0]), start = c(y_j,x_j), count = c(1,1))
        ncvar_put(nc_fd, "SR", length(ab_j[ab_j>0]), start = c(y_j,x_j), count = c(1,1))
        ncvar_put(nc_fd, "FRic", 0, start = c(y_j,x_j), count = c(1,1))
        ncvar_put(nc_fd, "FEve", 0, start = c(y_j,x_j), count = c(1,1))
        ncvar_put(nc_fd, "FDis", 0, start = c(y_j,x_j), count = c(1,1))
      } else {
        ### scale abundance to be relative: making abundance relative does not change the values of FD metrics, but it makes the calculation faster. 
        abundance_j[is.nan(abundance_j)] <- 0 #replace all NaN values with 0, some species might have data, others do not
        abundance_j <- make_relative(abundance_j)
        FD.ind <- dbFD(traits.dist, abundance_j, w.abun = T, stand.x = T, calc.CWM = F, stand.FRic = T, print.pco = T, m=3, corr="sqrt")  #takes ~20 sec for 2 sites
        ### Put metric values into netCDF file 
        ncvar_put(nc_fd, "SR", FD.ind$nbsp[1], start = c(y_j,x_j), count = c(1,1))
        ncvar_put(nc_fd, "FRic", FD.ind$FRic[1], start = c(y_j,x_j), count = c(1,1))
        ncvar_put(nc_fd, "FEve", FD.ind$FEve[1], start = c(y_j,x_j), count = c(1,1))
        ncvar_put(nc_fd, "FDis", FD.ind$FDis[1], start = c(y_j,x_j), count = c(1,1))
      }
    }
    nc_close(nc_sp)
    nc_close(nc_fd)
}



###########################################################################
## Convert netCDF into matrix
###########################################################################
grids_subset <- readRDS(file="Data/USGrids_cell_index.rds")

for (i in 1:52){
  grids_out <- grids_subset
  grids_out$week <- NA
  grids_out$lon <- NA
  grids_out$lat <- NA
  grids_out$SR <- NA
  grids_out$FRic <- NA
  grids_out$FEve <- NA
  grids_out$FDis <- NA
  saveRDS(grids_out, file=paste0("Data/fd_matrix_week",i,".rds"))
}


for (i in 1:52) {
    #open ncdf files to write to
    nc_fd_file <- file.path(paste0("Data/fd_week",i,".nc"))
    nc_fd <- nc_open(nc_fd_file) 
    grids_out <- readRDS(file=paste0("Data/fd_matrix_week",i,".rds"))
    
    for (k in seq(1, dim(grids_subset)[1])){
      x_k <- grids_subset$x_j[k]
      y_k <- grids_subset$y_j[k]
      
      ### Extract the variable
      x_k_extract <- ncvar_get(nc_fd, "x", start = c(x_k), count = c(1)) 
      y_k_extract <- ncvar_get(nc_fd, "y", start = c(y_k), count = c(1)) 
      sr_k <- ncvar_get(nc_fd, "SR", start = c(y_k,x_k), count = c(1,1)) 
      fric_k <- ncvar_get(nc_fd, "FRic", start = c(y_k,x_k), count = c(1,1)) 
      feve_k <- ncvar_get(nc_fd, "FEve", start = c(y_k,x_k), count = c(1,1)) 
      fdis_k <- ncvar_get(nc_fd, "FDis", start = c(y_k,x_k), count = c(1,1)) 
      
      ### Write in variable
      grids_out[k,4] <- i
      grids_out[k,5] <- x_k_extract
      grids_out[k,6] <- y_k_extract
      grids_out[k,7] <- sr_k
      grids_out[k,8] <- fric_k
      grids_out[k,9] <- feve_k
      grids_out[k,10] <- fdis_k
    }
    nc_close(nc_fd)
    saveRDS(grids_out, file=paste0("Data/fd_matrix_week",i,".rds"))
}
    