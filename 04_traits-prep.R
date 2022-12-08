###########################################################################
####### Spatiotemporal decoupling of avian taxonomic and functional diversity, Jarzyna & Stagge
####### 04: Prep avian traits
###########################################################################

###########################################################################
## Load Packages
###########################################################################
require(ncdf4)
require(FD)  #calculating multivariate and single trait distance-based FD indices
require(picante)  #calculating PD-equivalent indices and trait phylogenetic signal
require(dplyr)

###########################################################################
## Get traits and adjust names
###########################################################################
## Read in eBird species list and EltonTraits species trait designations
sp.nam <- readRDS(file="Data/sp-abun-paths-2019.rds")

traits <- read.csv("Traits/Master_Funcdiv.csv", header=TRUE) #Taken from Wilman et al. (2014)
colnames(traits) <- c("WJSpecID","scientific_name","common_name",colnames(traits[,4:ncol(traits)]))

## Reconcile taxonomy between EltonTraits and eBird
## Sage sparrow was split into two species : Belle's sparrow and Sagebrush sparrow
## the two new species are in eBird, but not EltonTraits
## we will create rows for those two species that have the same trait values as Sage sparrow
traits.sparrow.1 <- traits %>%
  filter(scientific_name  == "Amphispiza belli") %>%
  mutate(scientific_name = ifelse(common_name == "Sage Sparrow", "Artemisiospiza belli", scientific_name)) %>%
  mutate(common_name = ifelse(common_name == "Sage Sparrow", "Bell's Sparrow", common_name_name)) 

traits.sparrow.2 <- traits %>%
  filter(scientific_name  == "Amphispiza belli") %>%
  mutate(scientific_name = ifelse(common_name == "Sage Sparrow", "Artemisiospiza nevadensis", scientific_name)) %>%
  mutate(common_name = ifelse(common_name == "Sage Sparrow", "Sagebrush Sparrow", common_name_name)) 

traits <- traits %>%
  bind_rows(traits.sparrow.1) %>%
  bind_rows(traits.sparrow.2)


## Whip-poor-whil was split into two species: Eastern and Mexian Whip poor will
## the two new species are in eBird, but not EltonTraits
## we will create rows for those two species that have the same trait values as whip-poor-will
#"Caprimulgus vociferus", scientific_name2))  # "Caprimulgus vociferus" is a synonym for "Antrostomus vociferus" 
traits.whip.1 <- traits %>%
  filter(scientific_name  == "Caprimulgus vociferus") %>%
  mutate(scientific_name = ifelse(common_name == "Whip-poor-will", "Caprimulgus vociferus", scientific_name)) %>%
  mutate(common_name = ifelse(common_name == "Whip-poor-will", "Eastern Whip-poor-will", common_name_name)) 
traits.whip.2 <- traits %>%
  filter(scientific_name  == "Caprimulgus vociferus") %>%
  mutate(scientific_name = ifelse(common_name == "Whip-poor-will", "Caprimulgus arizonae", scientific_name)) %>%
  mutate(common_name = ifelse(common_name == "Whip-poor-will", "Mexican Whip-poor-will", common_name_name)) 

traits <- traits %>%
  bind_rows(traits.whip.1) %>%
  bind_rows(traits.whip.2)

#Colibri thalassinus
traits.col.1 <- traits %>%
  filter(scientific_name  == "Colibri thalassinus") %>%
  mutate(scientific_name = ifelse(common_name == "Green Violet-ear", "Colibri thalassinus", scientific_name)) %>%
  mutate(common_name = ifelse(common_name == "Green Violet-ear", "Mexican Violet-ear", common_name_name)) 
traits.col.2 <- traits %>%
  filter(scientific_name  == "Colibri thalassinus") %>%
  mutate(scientific_name = ifelse(common_name == "Green Violet-ear", "Colibri cyanotus", scientific_name)) %>%
  mutate(common_name = ifelse(common_name == "Green Violet-ear", "Lesser Violet-ear", common_name_name)) 

traits <- traits %>%
  bind_rows(traits.col.1) %>%
  bind_rows(traits.col.2)


#Momotus
traits.mot.1 <- traits %>%
  filter(scientific_name  == "Momotus momota") %>%
  mutate(scientific_name = ifelse(common_name == "Blue-crowned Motmot", "Momotus coeruliceps", scientific_name)) %>%
  mutate(common_name = ifelse(common_name == "Blue-crowned Motmot", "Blue-capped Motmot", common_name_name)) 
traits.mot.2 <- traits %>%
  filter(scientific_name  == "Momotus momota") %>%
  mutate(scientific_name = ifelse(common_name == "Blue-crowned Motmot", "Momotus lessonii", scientific_name)) %>%
  mutate(common_name = ifelse(common_name == "Blue-crowned Motmot", "Lesson's Motmot", common_name_name)) 

traits <- traits %>%
  bind_rows(traits.mot.1) %>%
  bind_rows(traits.mot.2)


##iberian green woodpecker used to be a subspecies of eurasian woodpecker
traits.wood.1 <- traits %>%
  filter(scientific_name  == "Picus viridis") %>%
  mutate(scientific_name = ifelse(common_name == "Eurasian Green Woodpecker", "Picus sharpei", scientific_name)) %>%
  mutate(common_name = ifelse(common_name == "Eurasian Green Woodpecker", "Iberian Green Woodpecker", common_name_name)) 

traits <- traits %>%
  bind_rows(traits.wood.1)

##iberian winger magpie is the same as azure-winged magpie used to be a subspecies of eurasian woodpecker
traits.mag.1 <- traits %>%
  filter(scientific_name  == "Cyanopica cyanus") %>%
  mutate(scientific_name = ifelse(common_name == "Azure-winged Magpie", "Cyanopica cooki", scientific_name)) %>%
  mutate(common_name = ifelse(common_name == "Azure-winged Magpie", "Iberian Azure-winged Magpie", common_name_name)) 

traits <- traits %>%
  bind_rows(traits.mag.1)

##northern island robin is the same as azure-winged magpie used to be a subspecies of eurasian woodpecker
traits.rob.1 <- traits %>%
  filter(scientific_name  == "Petroica australis") %>%
  mutate(scientific_name = ifelse(common_name == "New Zealand Robin", "Petroica longipes", scientific_name)) %>%
  mutate(common_name = ifelse(common_name == "New Zealand Robin", "North Island Robin", common_name_name)) 

traits <- traits %>%
  bind_rows(traits.rob.1)

##white-eye
traits.eye.1 <- traits %>%
  filter(scientific_name  == "Zosterops pallidus") %>%
  mutate(scientific_name = ifelse(common_name == "Pale White-eye", "Zosterops virens", scientific_name)) %>%
  mutate(common_name = ifelse(common_name == "Pale White-eye", "Cape White-eye", common_name_name)) 
traits <- traits %>%
  bind_rows(traits.eye.1)

## european stonechat
traits.stone.1 <- traits %>%
  filter(scientific_name  == "Saxicola torquatus") %>%
  mutate(scientific_name = ifelse(common_name == "Common Stonechat", "Saxicola rubicola", scientific_name)) %>%
  mutate(common_name = ifelse(common_name == "Common Stonechat", "European Stonechat", common_name_name)) 
traits <- traits %>%
  bind_rows(traits.stone.1)

##wagtail
traits.wag.1 <- traits %>%
  filter(scientific_name  == "Motacilla flava") %>%
  mutate(scientific_name = ifelse(common_name == "Yellow Wagtail", "Motacilla tschutschensis", scientific_name)) %>%
  mutate(common_name = ifelse(common_name == "Yellow Wagtail", "Eastern Yellow Wagtail", common_name_name)) 
traits <- traits %>%
  bind_rows(traits.wag.1)


##crossbill
traits.cross.1 <- traits %>%
  filter(scientific_name  == "Loxia curvirostra") %>%
  mutate(scientific_name = ifelse(common_name == "Red Crossbill", "Loxia sinesciuris", scientific_name)) %>%
  mutate(common_name = ifelse(common_name == "Red Crossbill", "Cassia Crossbill", common_name_name)) 
traits <- traits %>%
  bind_rows(traits.cross.1)

#oriole
traits.orio.1 <- traits %>%
  filter(scientific_name  == "Icterus dominicensis")%>%
  mutate(scientific_name = ifelse(common_name == "Greater Antillean Oriole", "Icterus portoricensis", scientific_name)) %>%
  mutate(common_name = ifelse(common_name == "Greater Antillean Oriole", "Puerto Rican Oriole", common_name_name)) 
traits <- traits %>%
  bind_rows(traits.orio.1)

#baham warbler
traits.warl.1 <- traits %>%
  filter(scientific_name  == "Dendroica dominica") %>%
  mutate(scientific_name = ifelse(common_name == "Yellow-throated Warbler", "Setophaga flavescens", scientific_name)) %>%
  mutate(common_name = ifelse(common_name == "Yellow-throated Warbler", "Bahama Warbler", common_name_name)) 
traits <- traits %>%
  bind_rows(traits.warl.1)

#st lucia warbler
traits.warl.2 <- traits %>%
  filter(scientific_name  == "Dendroica adelaidae") %>%
  mutate(scientific_name = ifelse(common_name == "Adelaide's Warbler", "Setophaga delicata", scientific_name)) %>%
  mutate(common_name = ifelse(common_name == "Adelaide's Warbler", "St. Lucia Warbler", common_name_name))

traits <- traits %>%
  bind_rows(traits.warl.2)

#costa rica warbler
traits.warl.3 <- traits %>%
  filter(scientific_name  == "Basileuterus tristriatus") %>%
  mutate(scientific_name = ifelse(common_name == "Three-striped Warbler", "Basileuterus melanotis", scientific_name)) %>%
  mutate(common_name = ifelse(common_name == "Three-striped Warbler", "Costa Rican Warbler", common_name_name))

traits <- traits %>%
  bind_rows(traits.warl.3)


#river warbler
traits.warl.4 <- traits %>%
  filter(scientific_name  == "Phaeothlypis rivularis") %>%
  mutate(scientific_name = ifelse(common_name == "Neotropical River Warbler", "Myiothlypis rivularis", scientific_name)) %>%
  mutate(common_name = ifelse(common_name == "Neotropical River Warbler", "Riverbank Warbler", common_name_name))

traits <- traits %>%
  bind_rows(traits.warl.4)




## Select species surveyed by eBird that have equivalent in EltonTrait
traits.1 <- traits %>%
  filter(scientific_name %in% sp.nam$scientific_name) #630 species 
## Select species surveyed by eBird that do not have equivalent in EltonTrait (taxonomy issues)
sp.nam.absent <- sp.nam %>%
  filter(!(scientific_name %in% traits$scientific_name)) #there are 177 (127 for 2019) species that do not appear in EltonTraits by their scientific names
sp.nam.present <- sp.nam.absent %>%
  filter(common_name %in% traits$common_name) #of these 177, there are 143 that appear and 34 (16 for 2019) species that do not appear in EltonTraits by their common names

# find those 143 species in EltonTraits and bind with traits.1 data frame that contain species matched based on scientific name
sp.nam.present2 <- sp.nam.present %>%
  left_join(traits, by = c("common_name"))
traits.2 <- sp.nam.present2 %>%
  mutate(scientific_name = scientific_name.x) %>%
  select(c(WJSpecID,scientific_name,common_name,colnames(traits[,4:ncol(traits)])))
traits.3 <- traits.1 %>%
  bind_rows(traits.2) #this contains 773 species that have the same either sci or common names in eBird and eltontraits

# see what's still missing
sp.nam.absent <- sp.nam.absent %>%
  filter(!(common_name %in% traits$common_name)) #of these 177, there are 34 species that do not appear in EltonTraits by their common names or scientific names
sp.nam.absent <- sp.nam.absent %>%
  mutate(scientific_name2 = scientific_name)

#sp.nam.absent <- sp.nam.absent %>%
#  mutate(scientific_name_b = scientific_name)
#sp.nam.absent <- sp.nam.absent %>%
#  mutate(scientific_name2 = ifelse(common_name == "Eastern Whip-poor-will", "Caprimulgus vociferus", scientific_name2))  # "Caprimulgus vociferus" is a synonym for "Antrostomus vociferus" 
#sp.nam.absent <- sp.nam.absent %>%
#  mutate(scientific_name2 = ifelse(common_name == "Mexican Whip-poor-will", "Caprimulgus arizonae", scientific_name2))  # "Caprimulgus vociferus" is a synonym for "Antrostomus arizonae" 


sp.nam.absent <- sp.nam.absent %>% 
  mutate(scientific_name2 = ifelse(common_name == "Ridgway's Rail", "Rallus longirostris", scientific_name2)) #Ridgway's Rail used to be a subspecies of Clapper rail; they have been since split into "Rallus longirostris" and "Rallus crepitans"
sp.nam.absent <- sp.nam.absent %>% 
  mutate(scientific_name2 = ifelse(common_name == "Common Gallinule", "Gallinula chloropus", scientific_name2)) #Common Gallinule was split from Common moorhen 
sp.nam.absent <- sp.nam.absent %>% 
  mutate(scientific_name2 = ifelse(common_name == "Gray-headed Swamphen", "Porphyrio porphyrio", scientific_name2)) #Gray-headed Swamphen was split from Purple swamphen
sp.nam.absent <- sp.nam.absent %>% 
  mutate(scientific_name2 = ifelse(common_name == "Snowy Plover", "Charadrius alexandrinus", scientific_name2)) #Snowy Plover was split from Kentish plover
sp.nam.absent <- sp.nam.absent %>% 
  mutate(scientific_name2 = ifelse(common_name == "Wilson's Snipe", "Gallinago gallinago", scientific_name2)) #Wilson's Snipe was split from Common snipe
sp.nam.absent <- sp.nam.absent %>% 
  mutate(scientific_name2 = ifelse(common_name == "Gray Hawk", "Buteo nitidus", scientific_name2)) #Gray Hawk was split from Gray-lined Hawk
sp.nam.absent <- sp.nam.absent %>% 
  mutate(scientific_name2 = ifelse(common_name == "Gray Fantail", "Rhipidura fuliginosa", scientific_name2)) #synonym
sp.nam.absent <- sp.nam.absent %>% 
  mutate(scientific_name2 = ifelse(common_name == "Northern Shrike", "Lanius excubitor", scientific_name2)) #Northern Shrike was split from Great Grey Shrike
sp.nam.absent <- sp.nam.absent %>% 
  mutate(scientific_name2 = ifelse(common_name == "Woodhouse's Scrub-Jay", "Aphelocoma californica", scientific_name2)) #Woodhouse's Scrub-Jay was split from California Scrub Jay
sp.nam.absent <- sp.nam.absent %>% 
  mutate(scientific_name2 = ifelse(common_name == "Eurasian Blue Tit", "Parus caeruleus", scientific_name2)) #synonym of Cyanistes caeruleus
sp.nam.absent <- sp.nam.absent %>% 
  mutate(scientific_name2 = ifelse(common_name == "Pacific Wren", "Troglodytes troglodytes", scientific_name2)) #Pacific Wren was split from Winter wren
sp.nam.absent <- sp.nam.absent %>% 
  mutate(scientific_name2 = ifelse(common_name == "Eastern Yellow Wagtail", "Motacilla flava", scientific_name2)) #Eastern Yellow Wagtail was split from Western Yellow Wagtail
sp.nam.absent <- sp.nam.absent %>% 
  mutate(scientific_name2 = ifelse(common_name == "Black-throated Gray Warbler", "Dendroica nigrescens", scientific_name2)) #Black-throated Gray Warbler was spelled as Black-throated Grey Warbler
sp.nam.absent <- sp.nam.absent %>% 
  mutate(scientific_name2 = ifelse(common_name == "Gray-and-gold Warbler", "Basileuterus fraseri", scientific_name2)) #gray was spelled grey
sp.nam.absent <- sp.nam.absent %>% 
  mutate(scientific_name2 = ifelse(common_name == "Gray-throated Warbler", "Basileuterus cinereicollis", scientific_name2)) #gray was spelled grey
sp.nam.absent <- sp.nam.absent %>% 
  mutate(scientific_name2 = ifelse(common_name == "Morelet's Seedeater", "Sporophila torqueola", scientific_name2)) #Morelet's Seedeater was split from White-collared Seedeater
sp.nam.absent <- sp.nam.absent %>% 
  mutate(scientific_name2 = ifelse(common_name == "LeConte's Sparrow", "Ammodramus leconteii", scientific_name2)) #old nomenclature ofr LeConte's Sparrow was Ammodramus leconteii instead of "Ammospiza leconteii"
#sp.nam.absent <- sp.nam.absent %>% 
#  mutate(scientific_name2 = ifelse(scientific_name == "Melanitta deglandi/stejnegeri", "Melanitta fusca", scientific_name2)) #split from M. fusca

#not resolved: Mexican duck
#remove for now, add below
sp.nam.absent <- sp.nam.absent %>%
  filter(!(scientific_name2 == "Anas diazi"))

## join with main traits data
traits.x <- traits %>%
  mutate(scientific_name2 = scientific_name)
sp.nam.absent <- sp.nam.absent %>%
  left_join(traits.x, by = c("scientific_name2"))
traits.4 <- sp.nam.absent %>%
  mutate(common_name = common_name.x) %>%
  mutate(scientific_name = scientific_name.x) %>%
  select(c(WJSpecID,scientific_name,common_name,colnames(traits[,4:ncol(traits)])))
traits.5 <- traits.3 %>%
  bind_rows(traits.4)

## replicate Mallard row and rename to Anas Diazi--Mexican duck, Mexican duck was recently split from Mallard, so we will use Mallard's traits
require(stringr)
mallard <- traits.5 %>%
  dplyr::filter(grepl("Mallard", common_name)) 
mallard[,2] <- "Anas diazi"
mallard[,3] <- "Mexican Duck"
traits.6 <- traits.5 %>%
  bind_rows(mallard)

## check if all species are accounted for, should be 807
t1 <- traits.6 %>%
  filter(scientific_name %in% sp.nam$scientific_name)
dim(t1)

## order by sp.nam
traits.7 <- traits.6[match(sp.nam$scientific_name, traits.6$scientific_name),]

## save the final traits csv
write.csv(traits.7, file="Traits/Master_Funcdiv_ebirdnomenclature.csv")
