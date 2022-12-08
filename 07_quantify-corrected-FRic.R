###########################################################################
####### Spatiotemporal decoupling of avian taxonomic and functional diversity, Jarzyna & Stagge
####### 07: Quantify corrected functional richness
###########################################################################

###########################################################################
## Load Packages
###########################################################################
#require(ncdf4)
require(dplyr)
require(INLA)

#grids_subset <- readRDS(file="Data/USGrids_cell_index.rds")

for (j in 1:52) {
    week_out <- readRDS(file=paste0("Data/fd_matrix_week",j,".rds"))
    week_out$SRlog <- log(week_out$SR)
    week_out$FRiclog <- log(week_out$FRic)
    week_out_sub <- week_out[!(is.na(week_out$FRic)),]
    week_out_sub <- week_out_sub[!(week_out_sub$FRic == 0),]
    week_out_sub <- week_out_sub[!(week_out_sub$SR == 0),]
    
    # simple lm regression
    lm2 <- lm(formula=FRic~SR, data=week_out_sub) 
    res2 <- lm2$residuals
    week_out_sub$FRicRes <- res2
    
    resid <- week_out_sub[,c(1,ncol(week_out_sub))]
    week_out_join <- left_join(week_out, resid, by="index_uniq")
    week_out_join[is.na(week_out_join)] <- NaN
    saveRDS(week_out_join,file=paste0("Data/fd_matrix_week",j,".rds"))
}

