###########################################################################
####### Spatiotemporal decoupling of avian taxonomic and functional diversity, Jarzyna & Stagge
####### 08: Convert to per nidex matrix
###########################################################################

###########################################################################
## Load Packages
###########################################################################
require(dplyr)
grids_subset <- readRDS(file="Data/USGrids_cell_index.rds")

vars <- c("SR","FRic","FEve","FDis","FRic-res")

for (i in 1:length(vars)){
grids_out <- grids_subset
grids_out$lon <- NA
grids_out$lat <- NA
grids_out$week1 <- NA
grids_out$week2 <- NA
grids_out$week3 <- NA
grids_out$week4 <- NA
grids_out$week5 <- NA
grids_out$week6 <- NA
grids_out$week7 <- NA
grids_out$week8 <- NA
grids_out$week9 <- NA
grids_out$week10 <- NA
grids_out$week11 <- NA
grids_out$week12 <- NA
grids_out$week13 <- NA
grids_out$week14 <- NA
grids_out$week15 <- NA
grids_out$week16 <- NA
grids_out$week17 <- NA
grids_out$week18 <- NA
grids_out$week19 <- NA
grids_out$week20 <- NA
grids_out$week21 <- NA
grids_out$week22 <- NA
grids_out$week23 <- NA
grids_out$week24 <- NA
grids_out$week25 <- NA
grids_out$week26 <- NA
grids_out$week27 <- NA
grids_out$week28 <- NA
grids_out$week29 <- NA
grids_out$week30 <- NA
grids_out$week31 <- NA
grids_out$week32 <- NA
grids_out$week33 <- NA
grids_out$week34 <- NA
grids_out$week35 <- NA
grids_out$week36 <- NA
grids_out$week37 <- NA
grids_out$week38 <- NA
grids_out$week39 <- NA
grids_out$week40 <- NA
grids_out$week41 <- NA
grids_out$week42 <- NA
grids_out$week43 <- NA
grids_out$week44 <- NA
grids_out$week45 <- NA
grids_out$week46 <- NA
grids_out$week47 <- NA
grids_out$week48 <- NA
grids_out$week49 <- NA
grids_out$week50 <- NA
grids_out$week51 <- NA
grids_out$week52 <- NA
saveRDS(grids_out, file=paste0("Data/index_matrix_",vars[i],".rds"))
}


vars <- c("SR","FRic","FEve","FDis","FRic-res")

###SR
i <- 1
grids_out <- readRDS(file=paste0("Data/index_matrix_",vars[i],".rds"))
week_coords <- readRDS(file=paste0("Data/fd_matrix_week",1,".rds"))
grids_out$lon <- week_coords$lon
grids_out$lat <- week_coords$lat

for (j in 1:52) {
  week_out <- readRDS(file=paste0("Data/fd_matrix_week",j,".rds"))
  grids_out[,5+j] <- week_out$SR
}

saveRDS(grids_out, file=paste0("Data/index_matrix_",vars[i],".rds"))


###FRic
i <- 2
grids_out <- readRDS(file=paste0("Data/index_matrix_",vars[i],".rds"))
week_coords <- readRDS(file=paste0("Data/fd_matrix_week",1,".rds"))
grids_out$lon <- week_coords$lon
grids_out$lat <- week_coords$lat

for (j in 1:52) {
  week_out <- readRDS(file=paste0("Data/fd_matrix_week",j,".rds"))
  grids_out[,5+j] <- week_out$FRic
}

saveRDS(grids_out, file=paste0("Data/index_matrix_",vars[i],".rds"))


###FEve
i <- 3
grids_out <- readRDS(file=paste0("Data/index_matrix_",vars[i],".rds"))
week_coords <- readRDS(file=paste0("Data/fd_matrix_week",1,".rds"))
grids_out$lon <- week_coords$lon
grids_out$lat <- week_coords$lat

for (j in 1:52) {
  week_out <- readRDS(file=paste0("Data/fd_matrix_week",j,".rds"))
  grids_out[,5+j] <- week_out$FEve
}

saveRDS(grids_out, file=paste0("Data/index_matrix_",vars[i],".rds"))


###FDis
i <- 4
grids_out <- readRDS(file=paste0("Data/index_matrix_",vars[i],".rds"))
week_coords <- readRDS(file=paste0("Data/fd_matrix_week",1,".rds"))
grids_out$lon <- week_coords$lon
grids_out$lat <- week_coords$lat

for (j in 1:52) {
  week_out <- readRDS(file=paste0("Data/fd_matrix_week",j,".rds"))
  grids_out[,5+j] <- week_out$FDis
}

saveRDS(grids_out, file=paste0("Data/index_matrix_",vars[i],".rds"))


###FRic-Res (cFRic)
i <- 5
grids_out <- readRDS(file=paste0("Data/index_matrix_",vars[i],".rds"))
week_coords <- readRDS(file=paste0("Data/fd_matrix_week",1,".rds"))
grids_out$lon <- week_coords$lon
grids_out$lat <- week_coords$lat

for (j in 1:52) {
  week_out <- readRDS(file=paste0("Data/fd_matrix_week",j,".rds"))
  grids_out[,5+j] <- week_out$FRicRes
}

saveRDS(grids_out, file=paste0("Data/index_matrix_",vars[i],".rds"))
