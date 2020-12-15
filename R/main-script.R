## load packages
packages_needed <- c("raster", "sp", "data.table", "truncdist", "gstat", "rgdal")
packages_to_install <- packages_needed[!( packages_needed %in% rownames(installed.packages()))]
if (length(packages_to_install) > 0)
  install.packages(packages_to_install)
lapply(packages_needed, require, character.only = TRUE)

#### load datasets #### 

# # volume dynamics data (TmFO)
# volDynData <- read.csv("volDynData.csv",  stringsAsFactors = F)
# # TmFO plots info (coordinates, year of logging, etc)
# plotInfo <- read.csv("plotInfo.csv",  stringsAsFactors = F)
# # RadamBrasil volume data
# volRadam <- read.csv("volRadam.csv",  stringsAsFactors = F)

#### load parameters posterior ####
paramsDam <- read.csv("data/paramsDam.csv")
paramsVolDyn <- read.csv("data/paramsVolDyn.csv")
paramsVolDynSite <- read.csv("data/paramsVolDynSite.csv")
paramsPropRecru <- read.csv("data/paramsPropRecru.csv")

#### load maps and spatially explicit variables ####

# spatially explicit variables from FORMIND
paramsFormind <- read.csv("data/paramsFormind.csv")

# open spatial variables
grd <- read.csv("data/spatVariables.csv", stringsAsFactors = F)
coordinates(grd) <- ~ longitude+latitude
gridded(grd) <- TRUE
proj4string(grd) <- CRS("+proj=longlat")

## concessions shapefiles 
### add here!
## area available for logging
concessions <- readOGR(dsn = "data/amazonia", layer = "test")
ext2 <- round(extent(raster_concessions2)) + c(-0.5, 0.5, -0.5, 0.5)
rst <- raster(resolution = 0.01, ext = ext2)
raster_concessions <- rasterize(concessions, rst)
raster_concessions2 <- !is.na(raster_concessions)
raster_concessions3 <- aggregate(raster_concessions2, fact = 100)

grd$pconcession <- raster::extract(raster_concessions3,grd)
grd$pconcession[is.na(grd$pconcession)] <- 0

## total area of pixels (in km2)
areaTot <- raster::extract(area(raster(grd)), grd)
## total area available for logging (in ha)
# grd$areaLogging <- areaTot*grd$pAreaLogging*100  ## km2 -> ha
grd$areaLogging <- areaTot*grd$pconcession*100  ## km2 -> ha


### inputs ####
## number of simulations
niter <- 100
## vext and trot combinations
logg_rules <- expand.grid(vext = c(10, 20, 30),
                          trot = c(15, 30, 65))
logg_rules <- subset(logg_rules, vext == 20 | trot == 30)
logg_rules$logid <- 1:nrow(logg_rules)
## simulation time: 300 years
tsim <- 300

#### source other R scripts #### 

source("R/functions.R")

source("R/predictionsFirstHarvest.R")

## Clear environment ## 
rm(list=setdiff(ls(), c("predictions", "tsim")))
source("R/functions.R")

source("R/predictionsTrajectory.R")
save(predicAm, file = "cache/predictions.Rdata")
