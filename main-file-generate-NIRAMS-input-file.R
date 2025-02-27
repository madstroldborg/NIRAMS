# ============================================================================ #
# Main file for loading and processing input data and saving the processed data
# in the required HDF5 format for NIRAMS.
# The processing of the climate data was developed by Adam Butler.
# Mads Troldborg 2024
# ============================================================================ #

# -------------------------------------------------------------------------- #
# LIBRARIES
# -------------------------------------------------------------------------- #

rm(list=objects())
rm(list=ls())

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("rhdf5")

library(sp)      # ok
library(ncdf4)   # for working with netCDF data
library(raster)  # for working with raster data
library(rhdf5)   # for saving input in hdf5 file format
library(stringr) # used for reading and loading agcensus files

# -------------------------------------------------------------------------- #
# PATHS
# -------------------------------------------------------------------------- #

# Path to input files
path_to_inputs <- "C:/Users/mt40375/Documents/NIRAMS/input-files/"

# Path to processed data folder
path_to_datafolder <- "Processed_RData/"

source("process-climate-and-PET-functions.R")
source("process-EMEP-data.R")
source("prepare-agcensus-data.R")
source("prepare-LCM-data.R")
source("process-land-use-inputs-for-NIRAMS-functions.R")  

meta <- read.csv("H5-naming-conventions.csv")
h5filename <- "Processed//nirams-inputs-2023.h5"


# -------------------------------------------------------------------------- #
# M A T C H   C E N S U S   &   N I R A M S   C A T E G O R I E S
# -------------------------------------------------------------------------- #
# POTENTIAL TO DO: Change names (reference to JAC and IACS no longer needed)

# Load files for matching the census categories to the nirams categories

# Animal categories (i.e., annual organic N input associated with livestock)
JAC_Livestock_N_Per_Animal_InputMatrix <- read.csv("JAC_Livestock_N_Per_Animal_InputMatrix.csv")

# Crop and land use categories (i.e., annual applied inorganic N and N uptake 
# associated with different land uses)
IACS_LU_CODE <- read.csv("InorgN_Uptake_PET_fact.csv")

# Load merged categories
Categories_AGCENSUS_NIRAMS <- read.csv("AGCENSUS_categories_fin.csv")
## replace '_' with '.' in category names
Categories_AGCENSUS_NIRAMS$AGCENSUS_cat <- gsub("-",".",Categories_AGCENSUS_NIRAMS$AGCENSUS_cat)
#Categories_AGCENSUS_NIRAMS[,c(1,2,3,4,6,7)]


# ============================================================================ #
# Subset of grid cells referring to Scotland/the NIRMAS grid
# ============================================================================ #
# -------------------------------------------------------------------------- #
# DEFINE THE NIRAMS GRID (BNG)
# -------------------------------------------------------------------------- #
# Below the NIRAMS grid is defined. This is used to extract only those data
# values from larger grids (UK or EU wide) that match the NIRAMS domain. 
# Specifically used for extracting data from EMEP and AGCensus.

# The grid used in NIRAMS is defined by the following BNG coordinates:
xmin = 0       # -8.17 (-8.2)
xmax = 485000  # -0.43 (-0.5)
ymin = 520000  # 54.42 (54.4)
ymax = 1235000 # 60.99 (61.0)

# The NIRAMS grid resolution is 1km x 1km. Hence, the grid spans 485 cells in 
# the x-direction/longitude (i.e. the number of columns is 485) and 715 cells 
# in the y-direction/latitude (i.e. the number of rows is 715)

# Create the grid used in NIRAMS
d.grid <- 1000  # grid resolution
x.grid <- seq(xmin+d.grid/2,xmax,by=d.grid)
nx     <- length(x.grid)
y.grid <- seq(ymin+d.grid/2,ymax,by=d.grid)
ny     <- length(y.grid)
loci   <- expand.grid(x=x.grid,y=y.grid)    # creating the grid center points
nirams.grid        	 <- SpatialPoints(loci) # convert to spatial points
gridded(nirams.grid) <- TRUE                # create the grid

# set the CRS to BNG
slot(nirams.grid, "proj4string")   <- CRS(SRS_string = "EPSG:27700")
#cat(wkt(nirams.grid))


# -------------------------------------------------------------------------- #
# Functions to get subset of climate data grid cells referring to Scotland
# -------------------------------------------------------------------------- #
getscot.elev <- function(z){ u <- t(z) ; u[,ncol(u):1] } 
## NOTE: this line isn't used in Python code - not sure why?
## MT 2023: not sure if the columns needs to be reversed here for Python

getscot.clim <- function(z, mt = FALSE){ 
  # The Met Office data are from a 5km grid (British National Grid projection)  
  # that covers all of the UK. This function extracts the values at the   
  # corresponding grid in NIRAMS (i.e. Scotland only).
  #
  # The Met Office grid is given my the coordinates:
  #    met.xmin = -200000; met.xmax = 700000; 
  #    met.ymin = -200000; met.ymax= 1250000
  #    At 5 km resolution it therefore consists of 180x290 cells
  # The NIRAMS grid is given by the coordinates: 
  #    xmin = 0; xmax = 485000; 
  #    ymin = 520000; ymax = 1235000
  #    At 1 km resolution this therefore consists of 485x715 cells.
  # xm and ym contain the indices of the Met Office grid cells covering Scotland.
  # This can be confirmed by loading the Met Office grid coordinates and extract
  # coordinates at indices - see check-spatial-grids-mar2023
  
  xm <- 41:137
  ym <- 145:287  # MT 2023: the NIRAMS code in Python reads the data upside-down
  # (I think due to the default way of plotting array); 
  # I therefore think this should be reversed to 287:145
  if(mt){ 
    z <- z[xm,ym,]
  }
  else{
    z <- z[xm,ym]
  }
  
  z
} 

scotdf2matrix <- function(z){z[,ncol(z):1] } 
# windows();image(scotdf2matrix(matrix(log(NPIG),nrow=485,ncol=715)))

## ########################################################################## ##
## G E T   I N P U T   D A T A   F O R   N I R A M S
## ########################################################################## ##

# ============================================================================ #
# STEP 1: Get daily CEDA climate data
# ============================================================================ #
# POTENTIAL TO DO: Change coverage so Northern Ireland is not included the data
varnames.day <- c("rainfall", "tasmin", "tasmax")

obj.clim.day <- get.ceda.data(dataset = "hadukgrid_uk_5km_day", varnames = varnames.day,
                              year.start = 2010, year.end = 2021, daily = TRUE, 
                              path.root = paste0(path_to_inputs,"Raw-data-climate/Daily/"), 
                              regfn = getscot.clim, bigval = 10000)
# ensure max temp > min temp (for some reason min temp is occasionally higher 
# than max temp)
for (i in 1:length(obj.clim.day$day)) {
  tmax <- obj.clim.day$tasmax[,,i]
  tmin <- obj.clim.day$tasmin[,,i]
  obj.clim.day$tasmax[,,i] <- pmax(tmax,tmin)
  obj.clim.day$tasmin[,,i] <- pmin(tmax,tmin)
}

#ndif0=0; for (i in 1:4383){ndif0 = ndif0+ length(which(obj.clim.day$tasmax[,,i]- obj.clim.day$tasmin[,,i]<0))}

# ============================================================================ #
# STEP 2: Get monthly CEDA climate data
# ============================================================================ #
# POTENTIAL TO DO: Change coverage so Northern Ireland is not included in the
#                  data. 
varnames.mon <- c("tas", "tasmin", "tasmax", "sun", "sfcWind", "rainfall", "pv", "hurs")

obj.clim.mon <- get.ceda.data(dataset = "hadukgrid_uk_5km_mon",
                              varnames = varnames.mon, 
                              year.start = 2010, year.end = 2021, daily = FALSE, 
                              path.root = paste0(path_to_inputs,"Raw-data-climate/Monthly/"), 
                              regfn = getscot.clim, bigval = 10000)

varnames.mon <- c(varnames.mon, c("pet.pm", "pet.thorn"))
#ndif0=0; for (i in 1:1:length(obj.clim.mon$month)){ndif0 = ndif0+ length(which(obj.clim.mon$tasmax[,,i]- obj.clim.mon$tasmin[,,i]<0))}

# ============================================================================ #
# STEP 3: Load and add non-climate data
# ============================================================================ #
# latitude (already loaded though??!) - needed for PET calculation
obj.clim.mon$lat  <- getscot.elev(as.matrix(read.table(paste0(path_to_inputs,"Raw-data-nonclimate/scot_latitude.txt"), skip = 6)))

# elevation - needed for PET calculation
obj.clim.mon$elev <- getscot.elev(as.matrix(read.table(paste0(path_to_inputs,"Raw-data-nonclimate/scot_elevation.txt"), skip = 6)))
obj.clim.mon$elev[obj.clim.mon$elev < -10] <- NA


# ============================================================================ #
# STEP 4: Get annual EMEP deposition data
# ============================================================================ #
# POTENTIAL TO DO: Change coverage so only Scotland included
obj.emep <- get.emep.data(path.to.files=paste0(path_to_inputs,"Raw-data-deposition/"), 
                          yr.start=2010, yr.end=2021, nirams.grid)

# Add emep deposition data to monthly climate data object - why??
obj.clim.mon$emep <- list(years = obj.emep$year, vals = obj.emep$n_deposition)

# ============================================================================ #
# STEP 5: Calculate PET using Penman-Monteith and Thornwaite
# ============================================================================ #
# Get mean temperature data for 2009 and 2022 (needed for calculating PET for 
# Jan 2010 and Dec 2021)
tas_2009 <- get.ceda.data(dataset = "hadukgrid_uk_5km_mon",
                              varnames = "tas", 
                              year.start = 2009, year.end = 2009, daily = FALSE, 
                              path.root = paste0(path_to_inputs,"Raw-data-climate/Monthly/"), 
                              regfn = getscot.clim, bigval = 10000)
tas_2022 <- get.ceda.data(dataset = "hadukgrid_uk_5km_mon",
                          varnames = "tas", 
                          year.start = 2022, year.end = 2022, daily = FALSE, 
                          path.root = paste0(path_to_inputs,"Raw-data-climate/Monthly/"), 
                          regfn = getscot.clim, bigval = 10000)
# PENMAN-MONTEITH PET
# Calculate monthly change in mean temperature
obj.clim.mon$tasdiff <- s2diff3d3(obj.clim.mon$tas)
obj.clim.mon$tasdiff[,,1] <- obj.clim.mon$tas[,,2]-tas_2009$tas[,,12]
rm(tas_2009)
obj.clim.mon$tasdiff[,,144] <- tas_2022$tas[,,1] - obj.clim.mon$tas[,,143]
rm(tas_2022)

# Parameters for PET calculations same as in 2017 calculations
pars.pet.pm <- data.frame(alb.coef = 0.23, # Albedo coefficients: for grassland
                          a.s = 0.25,      # Angstrom values: Fraction of et_rad reaching the earth on overcast days
                          b.s = 0.50)      # Angstrom values: (a_s + b_s) is fraction of et_rad reaching the earth on a clear day
# Calculate Penman-Monteith PET
obj.clim.mon$pet.pm <- calc.pet.pm(obj.clim.mon, pars.pet.pm = pars.pet.pm)

# THORNWAITE PET
# Calculate Thornwaite PET
obj.clim.mon$pet.thorn <- calc.pet.thorn(obj.clim.mon$tas, nuy = length(unique(obj.clim.mon$year)))



# ============================================================================ #
# STEP 6: Load and resample AGCensus raster files
# ============================================================================ #
agcensus.raster <- get.agcensus.data(dir.path=paste0(path_to_inputs,"Raw-data-AGCENSUS/Agcensus_raster_data_scotland_2003_2019/"), 
                                     st.yr=2010, end.yr=2019, nirams.grid) 

#save(agcensus.raster,file=paste0(path_to_datafolder,"agcensus_raster_2010_2019.Rdata"))
# ============================================================================ #
# STEP 7: Load LCM raster files
# ============================================================================ #
lcm_raster <- get.lcm.data(dir.path=paste0(path_to_inputs,"Raw-data-LCM/"), 
                           st.yr=2010, end.yr=2021, nirams.grid) 

#save(lcm_raster,file=paste0(path_to_datafolder,"lcm_raster_2010_2021.Rdata"))

# ============================================================================ #
# STEP 8: Assign NIRAMS land use categories
# ============================================================================ #
# The number of years for which agcensus and lcm data are available are not the 
# same. Here we are looking at the years 2010-2021, but we only have agcensus
# data from 2010-2019 and lcm data for 2015, 2017-2021. Hence, for assigning 
# NIRAMS land uses (for the subsequent N budget calculations) we assume:
#
# YEARS         AGCENSUS            LCM
# 2010-2015:    2010-2015           2015
# 2016-2017     2016-2017           2017
# 2018-2019     2018-2019           2018-2019
# 2020-2021     2019                2020-2021
obj.landuse <- list()
yrs <- seq(2010,2021, by=1)
nyr <- length(yrs)
lcm.yrs      <- as.character(c(2015,2015,2015,2015,2015,2015,2017,seq(2017,2021,by=1))) 
agcensus.yrs <- as.character(c(seq(2010,2019,by=1),2019,2019))

for (i in 1:nyr){
  obj.landuse[[i]] <- assign_nirams_landuse_categories(agcensus.raster[[agcensus.yrs[i]]], 
                                                          lcm_raster[[lcm.yrs[i]]], Categories_AGCENSUS_NIRAMS)
  names(obj.landuse)[[i]] <- yrs[i]
}

#test <- assign_nirams_landuse_categories(agcensus.raster$'2018', lcm_raster$'2018', Categories_AGCENSUS_NIRAMS)



# ============================================================================ #
# STEP 9: Calculate gridded organic N input from livestock 
# ============================================================================ #
# The organic N input is based on agcensus data on livestock. Because agcensus
# data are only available up until 2019, livestock numbers from 2019 are used 
# for 2020 and 2021 calculations. 

obj.livestock_N <- list()

for (i in 1:nyr){
  ifelse(i>1 && agcensus.yrs[i]==agcensus.yrs[i-1],
         {
           obj.livestock_N[[i]] <- obj.livestock_N[[i-1]]},
         {
           obj.livestock_N[[i]] <- Livestock_N(agcensus.raster[[agcensus.yrs[i]]],
                                      Categories_AGCENSUS_NIRAMS,
                                      JAC_Livestock_N_Per_Animal_InputMatrix)
           })
  names(obj.livestock_N)[[i]] <- yrs[i]
}

obj.orgN_kgha <- list()
for (i in 1:nyr){
  obj.orgN_kgha[[i]] <- org_n_distrib(obj.livestock_N[[i]], obj.landuse[[i]], 
                                      roughgraz_max = 10,impgraz_max = 250, 
                                      arable_max = 170)
  names(obj.orgN_kgha)[[i]] <- yrs[i]
  }

# save orgN in arrays and in similar format as the climate data
or_gr <- array(dim=c(nx,ny,nyr))
or_sp <- array(dim=c(nx,ny,nyr))
or_wi <- array(dim=c(nx,ny,nyr))
or_ot <- array(dim=c(nx,ny,nyr))

for (i in 1:nyr){
  rr <- rasterFromXYZ(obj.orgN_kgha[[i]])
  or_gr[,,i] <- getscot.elev(as.matrix(rr$kg_ha_grass))
  or_ot[,,i] <- getscot.elev(as.matrix(rr$kg_ha_othercrops))
  or_sp[,,i] <- getscot.elev(as.matrix(rr$kg_ha_springcrops))
  or_wi[,,i] <- getscot.elev(as.matrix(rr$kg_ha_wintercrops))
}

obj.orgN <- list()
obj.orgN$or_gr <- or_gr
obj.orgN$or_ot <- or_ot
obj.orgN$or_sp <- or_sp
obj.orgN$or_wi <- or_wi
obj.orgN$years <- yrs

# ============================================================================ #
# STEP 10: Calculate gridded inorganic N input and N uptake
# ============================================================================ #
# The inorganic N input and uptake are based on land cover.
obj.inorgN_uptake <- list()
for (i in 1:nyr){
  obj.inorgN_uptake[[i]] <- Inorg_N_Application_Uptake(obj.landuse[[i]],IACS_LU_CODE)
  names(obj.inorgN_uptake)[[i]] <- yrs[i]
}

# save inorganic N and uptake in arrays and in similar format as the climate data
in_gr <- array(dim=c(nx,ny,nyr))
in_ot <- array(dim=c(nx,ny,nyr))
in_sp <- array(dim=c(nx,ny,nyr))
in_wi <- array(dim=c(nx,ny,nyr))
up_gr <- array(dim=c(nx,ny,nyr))
up_ot <- array(dim=c(nx,ny,nyr))
up_sp <- array(dim=c(nx,ny,nyr))
up_wi <- array(dim=c(nx,ny,nyr))

for (i in 1:nyr){
  in_gr[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.inorgN_uptake[[i]]$in_gr)))
  in_ot[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.inorgN_uptake[[i]]$in_ot)))
  in_sp[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.inorgN_uptake[[i]]$in_sp)))
  in_wi[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.inorgN_uptake[[i]]$in_wi)))
  up_gr[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.inorgN_uptake[[i]]$up_gr)))
  up_ot[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.inorgN_uptake[[i]]$up_ot)))
  up_sp[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.inorgN_uptake[[i]]$up_sp)))
  up_wi[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.inorgN_uptake[[i]]$up_wi)))
}

obj.inorgN <- list()
obj.inorgN$in_gr <- in_gr
obj.inorgN$in_ot <- in_ot
obj.inorgN$in_sp <- in_sp
obj.inorgN$in_wi <- in_wi
obj.inorgN$years <- yrs

obj.uptake <- list()
obj.uptake$up_gr <- up_gr
obj.uptake$up_ot <- up_ot
obj.uptake$up_sp <- up_sp
obj.uptake$up_wi <- up_wi
obj.uptake$years <- yrs

# ============================================================================ #
# STEP 11: Calculate gridded PET_fact
# ============================================================================ #

obj.pet_fact <- list()
pet_fact     <- array(dim=c(nx,ny,nyr))

for (i in 1:nyr){
  obj.pet_fact[[i]] <- pet_facts_calc(obj.landuse[[i]], IACS_LU_CODE$PET_Fact)
  names(obj.pet_fact)[[i]] <- yrs[i]
  pet_fact[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.pet_fact[[i]])))
}

obj.pet_fact$pet_fact <- pet_fact
obj.pet_fact$years    <- yrs 

# ============================================================================ #
# STEP 12: Load data from previous HDF5 file
# ============================================================================ #
hdf5_path      <- "nirams_ii_input.h5"     # path to existing hdf5 file
#nirams_ii_2015 <- h5file("nirams_ii_input.h5" ,  mode = "a") 
# list groups and data sets in the hdf5 file
#dataset  <- list.datasets(nirams_ii_2015) # vs 1: list all datasets contained in the hdf5 file
#groups   <- list.groups(nirams_ii_2015)   #       list all groups in hdf5 file
#h5data <- h5ls(hdf5_path)
#h5groups <- unique(h5data[,1])

# Get soil properties from previous HDF5 file and add to list
soil.props    <- h5read(hdf5_path, 
                        "/one_km_grids/soil_properties/")

varnames.soil <- names(soil.props) 
obj.soil      <- list()
for (i in 1:length(varnames.soil)){
  obj.soil[[i]] <- scotdf2matrix(soil.props[[i]])
  names(obj.soil)[[i]] <- varnames.soil[i]
}

# Get land properties from previous HDF5 file and add to list
land.props  <- h5read("nirams_ii_input.h5", 
                      "/one_km_grids/old_land_properties/")
obj.land.props <- list()
obj.land.props$lcms88_pet_fact <- scotdf2matrix(land.props$lcms88_pet_fact)

# Get topographic data from previous HDF5 file and add to list
topo.data  <- h5read(hdf5_path, 
                     "/five_km_grids/topographical_data")

# Get idealized time series from previous HDF5 file
ts.data  <- h5read(hdf5_path, 
                   "/time_series")

ts_data <- ts.data$time_series_table[,c(1,2,6,7,8,9,14,3,4,5,14,10,12,13,11)]
names(ts_data) <- c("day","day_of_year","month","or_gr","or_sp","or_wi","or_ot","in_gr","in_sp","in_wi","in_ot","up_gr","up_sp", "up_wi","up_ot")
#fc  <- getscot.elev(as.matrix(read.table(paste0(path_to_inputs,"Raw-data-nonclimate/field_capacity.txt"), skip = 6)))
#sat <- getscot.elev(as.matrix(read.table(paste0(path_to_inputs,"Raw-data-nonclimate/saturation_capacity.txt"), skip = 6)))
# 
# raster.check= rasterFromXYZ(obj.orgN_kgha$"2010")
# image(matrix(NPIG,nrow=485,ncol=715))

# ============================================================================ #
# STEP 13: Create masking grid and mask non land-based emep values
# ============================================================================ #
# Currently basically using the lcm 2015 as the masking grid (which includes the 
# north of England). But wonder if I might as well just use the fc grid as the 
# masking grid. In the end, the NIRAMS calculations are only carried out for 
# cells for which all inputs are not NA.
# create masking grid
my.mask1 <- pet_fact[,,1]
my.mask1[!is.na(my.mask1)] <- 1
my.mask2 <- obj.soil$fc
my.mask2[my.mask2 >=0 ] <- 1

# mask non land-based emep values
for (i in 1:length(obj.emep$year)){
  emep.tmp <- obj.emep$n_deposition[,,i]
  emep.tmp[is.na(my.mask1)] <- NA # -9999
  obj.emep$n_deposition[,,i] <- emep.tmp
}

# Add emep deposition data to monthly climate data object - why??
obj.clim.mon$emep <- list(years = obj.emep$year, vals = obj.emep$n_deposition)




## ########################################################################## ##
## Single R object
## ########################################################################## ##
## Combine into a single R object

objmeta <- objdata <- as.list(NULL)

# -------------------------------------------------------------------------- #
# CREATE OBJMETA (basically a meta data object with the dates of all data sets)
# -------------------------------------------------------------------------- #
objmeta$day.year   <- obj.clim.day$year  # add the year for each of the daily climate data set
objmeta$day.julian <- obj.clim.day$day   # add the julian day for each of the daily climate data set

objmeta$mon.year   <- obj.clim.mon$year  # add the year for each of the monthly climate data set
objmeta$mon.month  <- obj.clim.mon$month # add the month for each of the monthly climate data set

objmeta$year       <- obj.clim.mon$emep$years  # add the year for each of the annual data sets

objmeta$fivekmgrid.longitude <- obj.clim.day$longitude
objmeta$fivekmgrid.latitude  <- obj.clim.day$latitude

## Add day and month to daily data
fn <- function(x,y){(y >= x[1]) & (y <= x[2])}

objmeta$day.month <- objmeta$day.day <- NA

leapy <- (objmeta$day.year == 2012 | objmeta$day.year == 2016 | objmeta$day.year == 2020)

for(k in 1:length(objmeta$day.year)){
  
  nim <- c(31, 28 + leapy[k], 31, 30, 31, 30, 31, 31, 30, 31, 30, 31) # days in months
  
  nim <- cbind(c(1, 1 + cumsum(nim)[-12]), cumsum(nim)) # start and end days in month
  
  objmeta$day.month[k] <- which(apply(nim, 1, fn, y = objmeta$day.julian[k]))
  
  objmeta$day.day[k]   <- objmeta$day.julian[k] + 1 - nim[objmeta$day.month[k],1]
}

# -------------------------------------------------------------------------- #
# CREATE OBJDATA (object with all the data sets)
# -------------------------------------------------------------------------- #
# First all non-climatic data sets are added to the list
objdata$onekmgrid.lcms88petfact  <- obj.land.props$lcms88_pet_fact # add default pet_fact to object
objdata$onekmgrid.fc       <- obj.soil$fc # add field capacity to object
objdata$onekmgrid.sat      <- obj.soil$sat # add soil saturation to object
objdata$onekmgrid.calibl   <- obj.soil$calibl # add calibration constant to object
objdata$onekmgrid.calibv   <- obj.soil$calibv # add calibration constant to object
objdata$onekmgrid.bfi      <- obj.soil$bfi # add bfi to object (needed?)
objdata$onekmgrid.spr      <- obj.soil$spr # add spr to object (needed?)

objdata$onekmgrid.petfact  <- pet_fact # add pet_fact to object
objdata$onekmgrid.ndep     <- obj.clim.mon$emep$vals # add N deposition data to object

objdata$onekmgrid.orgn.gr  <- obj.orgN$or_gr # add organic N data for grassland to object
objdata$onekmgrid.orgn.ot  <- obj.orgN$or_ot # add organic N data for other crops to object
objdata$onekmgrid.orgn.sp  <- obj.orgN$or_sp # add organic N data for spring crops to object
objdata$onekmgrid.orgn.wi  <- obj.orgN$or_wi # add organic N data for winter crops to object

objdata$onekmgrid.inorgn.gr  <- obj.inorgN$in_gr # add inorganic N data for grassland to object
objdata$onekmgrid.inorgn.ot  <- obj.inorgN$in_ot # add inorganic N data for other crops to object
objdata$onekmgrid.inorgn.sp  <- obj.inorgN$in_sp # add inorganic N data for spring crops to object
objdata$onekmgrid.inorgn.wi  <- obj.inorgN$in_wi # add inorganic N data for winter crops to object

objdata$onekmgrid.uptake.gr  <- obj.uptake$up_gr # add N uptake data for grassland to object
objdata$onekmgrid.uptake.ot  <- obj.uptake$up_ot # add N uptake data for other crops to object
objdata$onekmgrid.uptake.sp  <- obj.uptake$up_sp # add N uptake data for spring crops to object
objdata$onekmgrid.uptake.wi  <- obj.uptake$up_wi # add N uptake data for winter crops to object

objdata$fivekmgrid.lat  <- obj.clim.mon$lat # add latitude to object (needed?)
objdata$fivekmgrid.elev <- obj.clim.mon$elev # add elevation to object

nvday <- length(varnames.day) ; nvmon <- length(varnames.mon)

nd <- length(objdata)

# add daily climate data to object
for(k in 1:nvday){ 
  objdata[[nd+k]]      <- obj.clim.day[[varnames.day[k]]]
  names(objdata)[nd+k] <- paste("fivekmgrid.day", varnames.day[k], sep=".") 
}

# add daily average temperature to object
objdata$fivekmgrid.day.tas <- (objdata$fivekmgrid.day.tasmax + objdata$fivekmgrid.day.tasmin) / 2

# add monthly climate data to object
for(k in 1:nvmon){ 
  
  objdata[[nd+1+nvday+k]]      <- obj.clim.mon[[varnames.mon[k]]] 
  names(objdata)[nd+1+nvday+k] <- paste("fivekmgrid.mon", varnames.mon[k], sep=".") 
}

## ########################################################################## ##
## Prepare R objects to be saved in H5 format
## ########################################################################## ##
## MAKE SURE FILES ARE SAVED AS YYYY_MM_DD TO BE CONSISTENT WITH PREVIOUS H5 FILES
resolution <- meta$Resolution[match(names(objdata), as.character(meta$Variable))]

names(objdata) <- as.character(meta$H5naming[match(names(objdata), as.character(meta$Variable))])

nv <- length(resolution)

# create new list
new <- as.list(NULL)

m <- 1

for(k in 1:nv){
  
  if(resolution[k] == "Nontemporal"){
    
    new[[m]] <- objdata[[k]] 
    
    names(new)[m] <- names(objdata)[k]
    
    m <- m + 1
  }
  
  if(resolution[k] == "Annual"){
    
    for(j in 1:length(objmeta$year)){
      
      new[[m]] <- objdata[[k]][,,j]
      
      names(new)[m] <- gsub("YEAR", objmeta$year[j], names(objdata)[k])
      
      m <- m + 1
    }
  }
  
  if(resolution[k] == "Monthly"){
    
    for(j in 1:length(objmeta$mon.year)){
      
      new[[m]] <- objdata[[k]][,,j]
      
      names(new)[m] <- gsub("YEAR", objmeta$mon.year[j], names(objdata)[k])
      
      # slightly modified to ensure all months have 2 digits (i.e. add leading 
      # zero to 1 digit months)
      names(new)[m] <- gsub("MONTH", str_pad(objmeta$mon.month[j],2,pad="0"), names(new)[m])
      #names(new)[m] <- gsub("MONTH", objmeta$mon.month[j], names(new)[m])
      
      m <- m + 1
    }
  }
  
  if(resolution[k] == "Daily"){
    
    for(j in 1:length(objmeta$day.year)){
      
      new[[m]] <- objdata[[k]][,,j]
      
      names(new)[m] <- gsub("YEAR", objmeta$day.year[j], names(objdata)[k])
      
      
      names(new)[m] <- gsub("MONTH", str_pad(objmeta$day.month[j],2,pad="0"), names(new)[m])
      #names(new)[m] <- gsub("MONTH", formatC(objmeta$day.month[j], width=2, format="d",flag="0"), names(new)[m]) 
      #names(new)[m] <- gsub("MONTH", objmeta$day.month[j], names(new)[m])
      
      names(new)[m] <- gsub("DAY", str_pad(objmeta$day.day[j],2,pad="0"), names(new)[m])
      
      m <- m + 1
    }
  }
}

#image(new$`five-km-grids-meteorological-data-daily-max-temp-max-temp-2018-max-temp-2018-06-06`)
for(i in 1:length(new)){ new[[i]][is.na(new[[i]])] <- -9999 } # change NA to -9999

summary(factor(unlist(lapply(new, function(x){paste(dim(x),collapse="-")}))))

## ########################################################################## ##
## Save as H5 file
## ########################################################################## ##
# ------------------------------------------------------------------------- #
# Change the names of each data set to hdf5 format
# ------------------------------------------------------------------------- #
objnames <- names(new)

objnames <- gsub("-", "_", objnames)

objnames <- gsub("grids_", "grids/", objnames)
objnames <- gsub("data_", "data/", objnames)
objnames <- gsub("deposition_", "deposition/", objnames)
objnames <- gsub("n_uptake_", "n_uptake/", objnames)
objnames <- gsub("organic_n_", "organic_n/", objnames)
objnames <- gsub("inorganic_n_", "inorganic_n/", objnames)
objnames <- gsub("iacs_pet_facts_", "iacs_pet_facts/", objnames)
objnames <- gsub("soil_properties_", "soil_properties/", objnames)
objnames <- gsub("old_land_properties_", "old_land_properties/", objnames)

evar <- unique(meta$Evar)
evar_new <- gsub("-", "_", evar)

#for(k in 1:length(evar)){ objnames <- gsub(paste0("_", evar[k]), paste0("/", evar[k]), objnames) }
for(k in 1:length(evar_new)){ objnames <- gsub(paste0("_", evar_new[k]), paste0("/", evar_new[k]), objnames) }
## objnames <- paste0("/", objnames)

# ----------------------------------------------------------
# Create h5 file and add groups
# ----------------------------------------------------------
h5createFile(h5filename)

h5createGroup(h5filename, "lat_data")
h5createGroup(h5filename, "elev_data")

h5createGroup(h5filename, "one_km_grids")
h5createGroup(h5filename, "one_km_grids/n_deposition")
h5createGroup(h5filename, "one_km_grids/soil_properties")
h5createGroup(h5filename, "one_km_grids/old_land_properties")
h5createGroup(h5filename, "one_km_grids/iacs_pet_facts")
h5createGroup(h5filename, "one_km_grids/inorganic_n")
h5createGroup(h5filename, "one_km_grids/organic_n")
h5createGroup(h5filename, "one_km_grids/n_uptake")

h5createGroup(h5filename, "five_km_grids")
h5createGroup(h5filename, "five_km_grids/meteorological_data")
h5createGroup(h5filename, "five_km_grids/meteorological_data/daily")
h5createGroup(h5filename, "five_km_grids/meteorological_data/monthly")

h5createGroup(h5filename, "time_series")

#yrs <- 2010:2021

ed <- gsub("-", "_", as.character(meta$Evar[meta$Resolution == "Monthly"]))

for(i in 1:length(ed)){
  
  vi <- paste0("five_km_grids/meteorological_data/monthly/", ed[i])
  
  h5createGroup(h5filename, vi)
}

ed <- gsub("-", "_", as.character(meta$Evar[meta$Resolution == "Daily"]))

for(i in 1:length(ed)){
  
  vi <- paste0("five_km_grids/meteorological_data/daily/", ed[i])
  
  h5createGroup(h5filename, vi)
  
  for(j in 1:length(yrs)){
    
    vij <- paste0(vi, "/", ed[i], "_", yrs[j])
    
    h5createGroup(h5filename, vij)
  }
}


# ------------------------------------------------------------------------- #
# Add data sets to groups
# ------------------------------------------------------------------------- #
# Use scotdf2matrix function below to flip all data upside down; otherwise 
# won't be read correctly in Python model script
for(k in 1:length(new)){
  
  print(paste(k, date()))
  
  #h5write(new[[k]], file = h5filename, name = objnames[k])
  h5write(scotdf2matrix(new[[k]]), file = h5filename, name = objnames[k]) # flipping all data sets upside down to be consistent with NIRAMS
}

# add time series
ts.df <- ts.data$time_series_table

h5write(ts.data$time_series_table, file = h5filename, name = paste("time_series",names(ts.data),sep="/"),write.attributes=TRUE)
## ########################################################################## ##

# h5filename2 <- "Processed//nirams-inputs-2023-ts.h5"
# h5createFile(h5filename2)
# h5createGroup(h5filename2, "time_series")
# h5write(ts.df, file = h5filename2, name = paste("time_series",names(ts.data),sep="/"))
# #h5write(ts.df, file = h5filename2, name ="time_series/")


#h5close(h5file(h5filename))
#nirams_ii = h5file(h5filename)

# ndepo_check = h5read(h5filename, "/one_km_grids/n_deposition" ) 
# data5  <- h5read("../Met_Office_Data_w_weather_ip_1.h5", 
#                  "/met_data/rainfall")



