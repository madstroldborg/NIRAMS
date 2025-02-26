# ==============================================================================
# Main file for loading and processing future input data and saving the processed 
# data in the required HDF5 format for NIRAMS.
# Because the future climate data is in 1km resolution and hence a lot larger 
# files than the 5km resolution files used for the historical NIRAMS simulations,
# this creates memory issues in the way the input files are currently processed 
# and saved in the HDF5 file. This now therefore has to be done in steps.
# ==============================================================================

# ---------------------------------------------------------------------------- #
# LIBRARIES
# ---------------------------------------------------------------------------- #

rm(list=objects())
rm(list=ls())

library(sp)      # ok
library(ncdf4)   # for working with netCDF data
library(raster)  # for working with raster data
library(rhdf5)   # for saving input in hdf5 file format
library(stringr) # used for reading and loading agcensus files

# ---------------------------------------------------------------------------- #
# PATHS
# ---------------------------------------------------------------------------- #

# Path to input files
path_to_inputs <- "C:/Users/mt40375/Documents/NIRAMS/input-files/"

source("process-climate-and-PET-functions.R")
source("process-EMEP-data.R")
source("prepare-agcensus-data.R")
source("prepare-LCM-data.R")
source("process-future-land-use.R")
source("process-land-use-inputs-for-NIRAMS-functions.R")  

meta <- read.csv("H5-naming-conventions-future.csv")      # changed
#h5filename <- "Processed//nirams-future-calu-inputs-2060.h5"   # changed
path_to_datafolder <- "Processed_RData/" 

# ---------------------------------------------------------------------------
# M A T C H   C E N S U S   &   N I R A M S   C A T E G O R I E S
# ---------------------------------------------------------------------------
# Load files for matching the census categories to the nirams categories
# 
# Animal categories (i.e., annual organic N input associated with livestock)
JAC_Livestock_N_Per_Animal_InputMatrix <- read.csv("JAC_Livestock_N_Per_Animal_InputMatrix.csv")

# Crop and land use categories (i.e., annual applied inorganic N and N uptake
# associated with different land uses)
IACS_LU_CODE   <- read.csv("InorgN_Uptake_PET_fact.csv")
CRAFTY_LU_CODE <- read.csv("InorgN_Uptake_PET_fact_future.csv")

# Load merged categories
Categories_AGCENSUS_NIRAMS <- read.csv("AGCENSUS_categories_fin.csv")
## replace '_' with '.' in category names
Categories_AGCENSUS_NIRAMS$AGCENSUS_cat <- gsub("-",".",Categories_AGCENSUS_NIRAMS$AGCENSUS_cat)
# 

# ============================================================================= #
# Subset of grid cells referring to Scotland/the NIRMAS grid
# ============================================================================= #
# ---------------------------------------------------------------------------- #
# DEFINE THE NIRAMS GRID (BNG)
# ---------------------------------------------------------------------------- #
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


# ---------------------------------------------------------------------------- #
# Functions to get subset of the climate data grid cells referring to Scotland
# ---------------------------------------------------------------------------- #
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


getscot.clim.1k <- function(z, mt = FALSE){ 
  # The future climate data are from a 1km grid (British National Grid projection)  
  # that covers all of the UK, except Shetland! This function extracts the values   
  # at the corresponding grid in NIRAMS (i.e. Scotland only).
  #
  # The Chess-scape grid is given my the coordinates:
  #    met.xmin = 0; met.xmax = 656000; 
  #    met.ymin = 0; met.ymax= 1057000
  #    At 1 km resolution it therefore consists of 656x1057 cells
  # The NIRAMS grid is given by the coordinates: 
  #    xmin = 0;      xmax = 485000; 
  #    ymin = 520000; ymax = 1235000
  #    At 1 km resolution this therefore consists of 485x715 cells.
  # Note that the NIRAMS grid is not fully covered by the chess-scape grid!  
  # xm and ym contain the indices of the Chess-scape grid cells covering Scotland.
  
  xm <- 1:485
  ym <- 521:1057  # MT 2023: the NIRAMS code in Python reads the data upside-down
                  # (due to the default way of plotting array); should be reversed
  if(mt){ 
    z0 <- array(dim=c(485,715,dim(z)[3]))
    z0[1:485,1:537,] <- z[xm,ym,]
    
  }
  else{
    z0 <- array(dim=c(485,715))
    z0[1:485,1:537] <- z[xm,ym]
  }
  
  z0
} 


scotdf2matrix <- function(z){z[,ncol(z):1] } 


## ########################################################################## ##
## G E T   I N P U T   D A T A   F O R   N I R A M S
## ########################################################################## ##

# Chose years for which the future simulations will be carried out
year.name  = 2080 # name extension used for saving files # CHANGE HERE!
start.year = year.name-5
end.year   = year.name-1
nyrs       = length(start.year:end.year)
nms        = nyrs*12
nds        = nms*30   # all months have 30 days in Chess-Scape


# ============================================================================ #
# STEP 1: Get daily future climate data  ----
# ============================================================================ #
# NOTE: Future climate data saved in processed Rdata file! No need to run this!

varnames.day.future <- c("pr", "tasmin", "tasmax")

# Load future rainfall. Done in steps for memory reasons
rainfall.day.future <- get.future.ceda.data(dataset = "chess-scape_rcp60_bias-corrected_01", varnames = varnames.day.future[1],
                                            year.start = start.year, year.end = end.year, daily = TRUE, 
                                            path.root = paste0(path_to_inputs,"Raw-data-future-climate/Daily/"), 
                                            regfn = getscot.clim.1k, bigval = 10000)

# change precipitation flux (kg/m2/s) to rainfall accumulation (mm/day)
rainfall.day.future <- rainfall.day.future$pr*60*60*24
save(rainfall.day.future,file=paste0(path_to_datafolder,"rainfallDayFuture",year.name,".Rdata"))
rm(rainfall.day.future)
gc()
# Load future minimum temp. 
tasmin.day.future <- get.future.ceda.data(dataset = "chess-scape_rcp60_bias-corrected_01", varnames = varnames.day.future[2],
                                            year.start = start.year, year.end = end.year, daily = TRUE, 
                                            path.root = paste0(path_to_inputs,"Raw-data-future-climate/Daily/"), 
                                            regfn = getscot.clim.1k, bigval = 10000)

# change temperatures from kelvin to degree C
tasmin.day.future <- tasmin.day.future$tasmin - 273.15
save(tasmin.day.future,file=paste0(path_to_datafolder,"tasminDayFuture",year.name,".Rdata"))
rm(tasmin.day.future)
gc()
# Load future maximum temp. 
tasmax.day.future <- get.future.ceda.data(dataset = "chess-scape_rcp60_bias-corrected_01", varnames = varnames.day.future[3],
                                          year.start = start.year, year.end = end.year, daily = TRUE, 
                                          path.root = paste0(path_to_inputs,"Raw-data-future-climate/Daily/"), 
                                          regfn = getscot.clim.1k, bigval = 10000)

# change temperatures from kelvin to degree C
tasmax.day.future <- tasmax.day.future$tasmax - 273.15
save(tasmax.day.future,file=paste0(path_to_datafolder,"tasmaxDayFuture",year.name,".Rdata"))
rm(tasmax.day.future)

load(paste0(path_to_datafolder,"tasmaxDayFuture",year.name,".Rdata"))
load(paste0(path_to_datafolder,"tasminDayFuture",year.name,".Rdata"))
counter <- vector(length=nds)
# ensure max temp > min temp (for some reason min temp is occasionally higher than max temp)
for (i in 1:nds) {
  tmax <- tasmax.day.future[,,i]
  tmin <- tasmin.day.future[,,i]
  counter[i] <- sum(tmax<tmin,na.rm=TRUE)
  if(counter[i]>0) {
  tasmax.day.future[,,i] <- pmax(tmax,tmin)
  tasmin.day.future[,,i] <- pmin(tmax,tmin)
  }
  rm(tmax,tmin)
  gc()
}
#rm(tmax,tmin)

save(tasmax.day.future,file=paste0(path_to_datafolder,"tasmaxDayFuture",year.name,".Rdata"))
save(tasmin.day.future,file=paste0(path_to_datafolder,"tasminDayFuture",year.name,".Rdata"))
rm(tasmax.day.future,tasmin.day.future)
gc()
# # NOTE: The daily data are not run for storage/memory reasons
# obj.clim.day.future <- get.future.ceda.data(dataset = "chess-scape_rcp60_bias-corrected_01", varnames = varnames.day.future,
#                               year.start = start.year, year.end = end.year, daily = TRUE, 
#                               path.root = paste0(path_to_inputs,"Raw-data-future-climate/Daily/"), 
#                               regfn = getscot.clim.1k, bigval = 10000)
# # change temperatures from kelvin to degree C
# obj.clim.day.future$tasmax <- obj.clim.day.future$tasmax - 273.15
# obj.clim.day.future$tasmin <- obj.clim.day.future$tasmin - 273.15
# 
# # change precipitation flux (kg/m2/s) to rainfall accumulation (mm/day)
# obj.clim.day.future$rainfall <- obj.clim.day.future$pr*60*60*24
# 
# # ensure max temp > min temp (for some reason min temp is occasionally higher than max temp)
# for (i in 1:length(obj.clim.day.future$day)) {
#   tmax <- obj.clim.day.future$tasmax[,,i]
#   tmin <- obj.clim.day.future$tasmin[,,i]
#   obj.clim.day.future$tasmax[,,i] <- pmax(tmax,tmin)
#   obj.clim.day.future$tasmin[,,i] <- pmin(tmax,tmin)
# }
# rm(tmax,tmin)
# ---- End of STEP 1 ----

# ============================================================================ #
# STEP 2: Get monthly CEDA climate data  ----
# ============================================================================ #
# NOTE: Monthly future climate (incl evapotranspiration) and EMEP data (i.e. 
# Step 2-5) saved in processed Rdata file.
varnames.mon.future <- c("tas", "tasmin", "tasmax", "rsds", "sfcWind", "hurs")

obj.clim.mon.future <- get.future.ceda.data(dataset = "chess-scape_rcp60_bias-corrected_01",
                              varnames = varnames.mon.future, 
                              year.start = start.year, year.end = end.year, daily = FALSE, 
                              path.root = paste0(path_to_inputs,"Raw-data-future-climate/Monthly/"), 
                              regfn = getscot.clim.1k, bigval = 10000)

# change temperatures from kelvin to degree C
obj.clim.mon.future$tas    <- obj.clim.mon.future$tas - 273.15
obj.clim.mon.future$tasmax <- obj.clim.mon.future$tasmax - 273.15
obj.clim.mon.future$tasmin <- obj.clim.mon.future$tasmin - 273.15

varnames.mon.future <- c(varnames.mon.future, c("pet.pm", "pet.thorn"))

# ---- End of STEP 2 ----

# ============================================================================ #
# STEP 3: Load and add non-climate data ----
# ============================================================================ #
# # Option 1: load directly from 1km file (TO BE GENERATED!!)
# obj.clim.mon.future$lat  <- getscot.elev(as.matrix(read.table(paste0(path_to_inputs,"Raw-data-nonclimate/scot_latitude.txt"), skip = 6)))
obj.clim.mon.future$lat  <- getscot.elev(as.matrix(read.table(paste0(path_to_inputs,"Raw-data-nonclimate/scot_latitude_1km.asc"), skip = 6)))
# # elevation - needed for PET calculation
# obj.clim.mon.future$elev <- getscot.elev(as.matrix(read.table(paste0(path_to_inputs,"Raw-data-nonclimate/scot_elevation.txt"), skip = 6)))
# obj.clim.mon.future$elev[obj.clim.mon.future$elev < -10] <- NA
obj.clim.mon.future$elev <- getscot.elev(as.matrix(read.table(paste0(path_to_inputs,"Raw-data-nonclimate/scot_elevation_bil_1km.asc"), skip = 6)))
obj.clim.mon.future$elev[obj.clim.mon.future$elev < -10] <- NA

# # Option 2: resample 5km to 1km grid
# e.nirams  <- extent(nirams.grid) # extent(0, 485000, 520000, 1235000)
# my.crs    <- proj4string(nirams.grid) 
# raster1km <- raster(resolution=c(1000,1000), crs=my.crs, ext=e.nirams)
# 
# # latitude  - needed for PET calculation
# lat <- list()
# lat$x <- seq(2500,485000, by=5000)
# lat$y <- seq(522500,1235000,by=5000)
# lat$z <- getscot.elev(as.matrix(read.table(paste0(path_to_inputs,"Raw-data-nonclimate/scot_latitude.txt"), skip = 6)))
# lat$z[lat$z== -9999] <- NA
# r1      <- raster(lat)
# lat.rs  <- resample(r1, raster1km, method='ngb')
# obj.clim.mon.future$lat  <- getscot.elev(as.matrix(lat.rs))
# 
# # elevation - needed for PET calculation
# elev <- list()
# elev$x <- seq(2500,485000, by=5000)
# elev$y <- seq(522500,1235000,by=5000)
# elev$z <- getscot.elev(as.matrix(read.table(paste0(path_to_inputs,"Raw-data-nonclimate/scot_elevation.txt"), skip = 6)))
# elev$z[elev$z== -9999] <- NA
# r2 <- raster(elev)
# 
# elev.rs  <- resample(r2, raster1km, method='ngb')
# obj.clim.mon.future$elev <- getscot.elev(as.matrix(elev.rs))
# obj.clim.mon.future$elev[obj.clim.mon.future$elev < -10] <- NA
# 
# rm(lat,r1,lat.rs,elev,r2,elev.rs)
# ---- End of STEP 3 ----

# ============================================================================ #
# STEP 4: Get annual EMEP deposition data ----
# ============================================================================ #
# For the future simulations, the annual N deposition is assumed to be the 
# average N deposition from 2016-2021. 
# The N deposition can either be extracted for 2016-2021 from the input files or 
# can be loaded from the NIRAMS 2023 hdf5 input file.
hdf5_path_new      <- "Processed/nirams-inputs-2023.h5"     # path to existing hdf5 file

# Get emep data from previous HDF5 file and add to list
emep.tmp    <- h5read(hdf5_path_new,"/one_km_grids/n_deposition/")
emep.future <- Reduce("+",emep.tmp[7:12])/length(emep.tmp[7:12]) # average of 2016-2021 data
emep.future <- scotdf2matrix(emep.future)
emep.future[emep.future== -9999] <- NA
rm(emep.tmp)
# # option 2
# obj.emep <- get.emep.data(path.to.files=paste0(path_to_inputs,"Raw-data-deposition/"), 
#                           yr.start=2016, yr.end=2021, nirams.grid)
# get mean N deposition for 2016-2021
# emep.future <- apply(obj.emep$n_deposition,c(1,2),mean,na.rm=TRUE)

# Add emep deposition data to monthly climate data object - why??
obj.clim.mon.future$emep <- list(years = start.year:end.year, vals = emep.future)

# ---- End of STEP 4 ----

# ============================================================================ #
# STEP 5: Calculate PET using Penman-Monteith and Thornwaite ----
# ============================================================================ #
# Get mean temperature data for Dec before start year and Jan after end year 
# (needed for calculating PET )
myfilename = paste0(path_to_inputs,"Raw-data-future-climate/Monthly/",
                    "chess-scape_rcp60_bias-corrected_01_tas_uk_1km_monthly_19801201-20801130.nc")

tas_ini <- chessscape.read.specific.mon(filename = myfilename,varname = "tas", 
                              regfn = getscot.clim.1k, bigval = 10000,
                              yr = start.year-1, mm = 12)
tas_end <- chessscape.read.specific.mon(filename = myfilename,varname = "tas", 
                                        regfn = getscot.clim.1k, bigval = 10000,
                                        yr = end.year+1, mm = 1)
# PENMAN-MONTEITH PET
# Calculate monthly change in mean temperature
obj.clim.mon.future$tasdiff            <- s2diff3d3(obj.clim.mon.future$tas)
obj.clim.mon.future$tasdiff[,,1]       <- obj.clim.mon.future$tas[,,2]-tas_ini
obj.clim.mon.future$tasdiff[,,nms] <- tas_end - obj.clim.mon.future$tas[,,nms-1]
rm(tas_ini,tas_end)

# Parameters for PET calculations same as in 2017 calculations
pars.pet.pm <- data.frame(alb.coef = 0.23, # Albedo coefficients: for grassland
                          a.s = 0.25,      # Angstrom values: Fraction of et_rad reaching the earth on overcast days
                          b.s = 0.50)      # Angstrom values: (a_s + b_s) is fraction of et_rad reaching the earth on a clear day
# Calculate Penman-Monteith PET
obj.clim.mon.future$pet.pm <- calc.pet.pm.future(obj.clim.mon.future, pars.pet.pm = pars.pet.pm)

# THORNWAITE PET
# Calculate Thornwaite PET
obj.clim.mon.future$pet.thorn <- calc.pet.thorn(obj.clim.mon.future$tas, nuy = length(unique(obj.clim.mon.future$year)))

save(obj.clim.mon.future,file=paste0(path_to_datafolder,"objClimMonthly",year.name,".RData"))
rm(obj.clim.mon.future)


# load("monthlyClimate20112021.RData")   # ???

# ---- End of STEP 5 ----

# NOTE: FOR FUTURE CLIMATE SCENARIOS ONLY, STEP 6-11 CAN BE SKIPPED. FOR THESE
# SCENARIOS THER LAND USE IS ASSUMED FIXED

# ============================================================================ #
# STEP 6: Load and resample AGCensus raster files ----
# ============================================================================ #

agcensus.raster <- get.agcensus.data(dir.path=paste0(path_to_inputs,"Raw-data-AGCENSUS/Agcensus_raster_data_scotland_2003_2019/"), 
                                     st.yr=2015, end.yr=2019, nirams.grid) 
#save(agcensus.raster,file=paste0(path_to_datafolder,"agcensus_raster_2010_2019.Rdata"))
#load(file=paste0(path_to_datafolder,"agcensus_raster_2010_2019.Rdata"))
# ---- End of STEP 6 ----

# ============================================================================ #
# STEP 7: Load LCM raster files ----
# ============================================================================ #
# NOTE: LCM data saved in processed Rdata file.
lcm_raster <- get.lcm.data(dir.path=paste0(path_to_inputs,"Raw-data-LCM/"), 
                           st.yr=2015, end.yr=2021, nirams.grid) 

#save(lcm_raster,file=paste0(path_to_datafolder,"lcm_raster_2010_2021.Rdata"))
# ---- End of STEP 7 ----

# ============================================================================ #
# STEP 8: Assign NIRAMS land use categories ----
# ============================================================================ #
# THIS STEP IS NOT NEEDED FOR FUTURE SIMULATIONS
# The number of years for which agcensus and lcm data are available are not the 
# same. Here we are looking at the years 2010-2021, but we only have agcensus
# data from 2010-2019 and lcm data for 2015, 2017-2021. For the assignment of 
# NIRAMS land uses (for the subsequent N budget calculations) we therefore assume:
#
# YEARS         AGCENSUS            LCM
# 2010-2015:    2010-2015           2015
# 2016-2017     2016-2017           2017
# 2018-2019     2018-2019           2018-2019
# 2020-2021     2019                2020-2021
obj.landuse <- list()
yrs <- seq(2016,2019, by=1) # yrs <- seq(2010,2021, by=1)
nyr <- length(yrs)
lcm.yrs      <- as.character(c(2017,yrs[2:4]))# as.character(c(2015,2015,2015,2015,2015,2015,2017,seq(2017,2021,by=1))) 
agcensus.yrs <- as.character(yrs) # as.character(c(seq(2010,2019,by=1),2019,2019))

for (i in 1:nyr){
  obj.landuse[[i]] <- assign_nirams_landuse_categories(agcensus.raster[[agcensus.yrs[i]]], 
                                                       lcm_raster[[lcm.yrs[i]]], Categories_AGCENSUS_NIRAMS)
  names(obj.landuse)[[i]] <- yrs[i]
}

# ---- End of STEP 8 ----

# ============================================================================ #
# STEP 9: Calculate gridded organic N input from livestock ----
# ============================================================================ #
# The organic N input is based on agcensus data on livestock. Because agcensus
# data are only available up until 2019, livestock numbers from 2019 are used 
# for 2020 and 2021 calculations. 

# Calculate total org N input for each grid cell
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

# Distribute org N
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

# ---- End of STEP 9 ----

# ============================================================================ #
# STEP 10: Calculate gridded inorganic N input and N uptake ----
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

# ---- End of Step 10 ----

# ============================================================================ #
# STEP 11: Calculate gridded PET_fact ----
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

# ---- End of Step 11 ----

# ============================================================================ #
# STEP 7A: Load CRAFTY raster files ----
# ============================================================================ #
# Load future most dominant land use class as simulated with CRAFTY
crafty_raster <- get.crafty.data(dir.path=paste0(path_to_inputs,"Raw-data-future-landuse/Reproject/RCP6_0-SSP3/"), 
                                 nirams.grid) 
crafty.yrs <- names(crafty_raster)

# identify required crafty year
cyr_id          <- which(names(crafty_raster)==as.character(year.name))


# ---- End of STEP 7A ----
# ============================================================================ #
# STEP 9A: Calculate gridded organic N input from livestock CRAFTY ----
# ============================================================================ #
# The future organic N input is based on agcensus data on livestock from 2019. 
# The organic N is distributed based on CRAFTY future land use data. 
agcensus.raster <- get.agcensus.data(dir.path=paste0(path_to_inputs,"Raw-data-AGCENSUS/Agcensus_raster_data_scotland_2003_2019/"), 
                                     st.yr=2019, end.yr=2019, nirams.grid) 

# Calculate total org N input for each grid cell
obj.future_livestock_N <- list()

obj.future_livestock_N[[1]] <- Livestock_N(agcensus.raster[["2019"]],
                                           Categories_AGCENSUS_NIRAMS,
                                           JAC_Livestock_N_Per_Animal_InputMatrix)

names(obj.future_livestock_N)[[1]] <- 2019


# Distribute org N based on CRAFTY
obj.future_orgN_kgha <- list()

obj.future_orgN_kgha[[1]] <- future_org_n_distrib(obj.future_livestock_N[["2019"]], crafty_raster[[cyr_id]], 
                                                  roughgraz_max = 10,impgraz_max = 250, 
                                                  arable_max = 170)
names(obj.future_orgN_kgha)[[1]] <- crafty.yrs[cyr_id]

# save orgN in arrays and in similar format as the climate data
or_gr <- array(dim=c(nx,ny))
or_sp <- array(dim=c(nx,ny))
or_wi <- array(dim=c(nx,ny))
or_ot <- array(dim=c(nx,ny))


rr <- rasterFromXYZ(obj.future_orgN_kgha[[1]])
or_gr <- getscot.elev(as.matrix(rr$kg_ha_grass))
or_ot <- getscot.elev(as.matrix(rr$kg_ha_othercrops))
or_sp <- getscot.elev(as.matrix(rr$kg_ha_springcrops))
or_wi <- getscot.elev(as.matrix(rr$kg_ha_wintercrops))

# for (i in 1:length(crafty.yrs)){
#   obj.future_orgN_kgha[[i]] <- future_org_n_distrib(obj.future_livestock_N[["2019"]], crafty_raster[[i]],
#                                       roughgraz_max = 10,impgraz_max = 250,
#                                       arable_max = 170)
#   names(obj.future_orgN_kgha)[[i]] <- crafty.yrs[i]
# }

# # save orgN in arrays and in similar format as the climate data
# or_gr <- array(dim=c(nx,ny,length(crafty.yrs)))
# or_sp <- array(dim=c(nx,ny,length(crafty.yrs)))
# or_wi <- array(dim=c(nx,ny,length(crafty.yrs)))
# or_ot <- array(dim=c(nx,ny,length(crafty.yrs)))
# for (i in 1:length(crafty.yrs)){
#   rr <- rasterFromXYZ(obj.future_orgN_kgha[[i]])
#   or_gr[,,i] <- getscot.elev(as.matrix(rr$kg_ha_grass))
#   or_ot[,,i] <- getscot.elev(as.matrix(rr$kg_ha_othercrops))
#   or_sp[,,i] <- getscot.elev(as.matrix(rr$kg_ha_springcrops))
#   or_wi[,,i] <- getscot.elev(as.matrix(rr$kg_ha_wintercrops))
# }

obj.orgN.future <- list()
obj.orgN.future$or_gr <- or_gr
obj.orgN.future$or_ot <- or_ot
obj.orgN.future$or_sp <- or_sp
obj.orgN.future$or_wi <- or_wi
obj.orgN.future$years <- start.year:end.year # as.numeric(crafty.yrs)

rm(obj.future_orgN_kgha)
save(obj.orgN.future, file=paste0(path_to_datafolder,"orgNfuture_", year.name,".RData"))    # organic N
# ---- End of Step 9A ----
# ============================================================================ #
# STEP 10A: Calculate gridded inorganic N input and N uptake CRAFTY ----
# ============================================================================ #
# The inorganic N input and uptake are based on land cover.
obj.future_inorgN_uptake <- list()

obj.future_inorgN_uptake[[1]] <- Future_Inorg_N_Application_Uptake(crafty_raster[[cyr_id]],CRAFTY_LU_CODE)
names(obj.future_inorgN_uptake)[[1]] <- crafty.yrs[cyr_id]


# save inorganic N and uptake in arrays and in similar format as the climate data
in_gr <- array(dim=c(nx,ny))
in_ot <- array(dim=c(nx,ny))
in_sp <- array(dim=c(nx,ny))
in_wi <- array(dim=c(nx,ny))
up_gr <- array(dim=c(nx,ny))
up_ot <- array(dim=c(nx,ny))
up_sp <- array(dim=c(nx,ny))
up_wi <- array(dim=c(nx,ny))


in_gr <- getscot.elev(as.matrix(rasterFromXYZ(obj.future_inorgN_uptake[[1]]$in_gr)))
in_ot <- getscot.elev(as.matrix(rasterFromXYZ(obj.future_inorgN_uptake[[1]]$in_ot)))
in_sp <- getscot.elev(as.matrix(rasterFromXYZ(obj.future_inorgN_uptake[[1]]$in_sp)))
in_wi <- getscot.elev(as.matrix(rasterFromXYZ(obj.future_inorgN_uptake[[1]]$in_wi)))
up_gr <- getscot.elev(as.matrix(rasterFromXYZ(obj.future_inorgN_uptake[[1]]$up_gr)))
up_ot <- getscot.elev(as.matrix(rasterFromXYZ(obj.future_inorgN_uptake[[1]]$up_ot)))
up_sp <- getscot.elev(as.matrix(rasterFromXYZ(obj.future_inorgN_uptake[[1]]$up_sp)))
up_wi <- getscot.elev(as.matrix(rasterFromXYZ(obj.future_inorgN_uptake[[1]]$up_wi)))



# obj.future_inorgN_uptake <- list()
# for (i in 1:length(crafty.yrs)){
#   obj.future_inorgN_uptake[[i]] <- Future_Inorg_N_Application_Uptake(crafty_raster[[i]],CRAFTY_LU_CODE)
#   names(obj.future_inorgN_uptake)[[i]] <- crafty.yrs[i]
# }
# nyr_f <- length(crafty.yrs)
# # save inorganic N and uptake in arrays and in similar format as the climate data
# in_gr <- array(dim=c(nx,ny,nyr_f))
# in_ot <- array(dim=c(nx,ny,nyr_f))
# in_sp <- array(dim=c(nx,ny,nyr_f))
# in_wi <- array(dim=c(nx,ny,nyr_f))
# up_gr <- array(dim=c(nx,ny,nyr_f))
# up_ot <- array(dim=c(nx,ny,nyr_f))
# up_sp <- array(dim=c(nx,ny,nyr_f))
# up_wi <- array(dim=c(nx,ny,nyr_f))
# 
# for (i in 1:nyr_f){
#   in_gr[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.future_inorgN_uptake[[i]]$in_gr)))
#   in_ot[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.future_inorgN_uptake[[i]]$in_ot)))
#   in_sp[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.future_inorgN_uptake[[i]]$in_sp)))
#   in_wi[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.future_inorgN_uptake[[i]]$in_wi)))
#   up_gr[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.future_inorgN_uptake[[i]]$up_gr)))
#   up_ot[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.future_inorgN_uptake[[i]]$up_ot)))
#   up_sp[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.future_inorgN_uptake[[i]]$up_sp)))
#   up_wi[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.future_inorgN_uptake[[i]]$up_wi)))
# }

obj.inorgN.future <- list()
obj.inorgN.future$in_gr <- in_gr
obj.inorgN.future$in_ot <- in_ot
obj.inorgN.future$in_sp <- in_sp
obj.inorgN.future$in_wi <- in_wi
obj.inorgN.future$years <- start.year:end.year#as.numeric(crafty.yrs)

obj.uptake.future <- list()
obj.uptake.future$up_gr <- up_gr
obj.uptake.future$up_ot <- up_ot
obj.uptake.future$up_sp <- up_sp
obj.uptake.future$up_wi <- up_wi
obj.uptake.future$years <- start.year:end.year#as.numeric(crafty.yrs) # crafty years or simulation years????

rm(obj.future_inorgN_uptake)
save(obj.inorgN.future, file=paste0(path_to_datafolder,"inorgNfuture_", year.name,".RData"))  # inorganic N
save(obj.uptake.future, file=paste0(path_to_datafolder,"uptakefuture_", year.name,".RData"))  # N uptake


# ---- End of Step 10A ----


# ============================================================================ #
# STEP 11A: Calculate gridded PET_fact ----
# ============================================================================ #

obj.pet_fact.future <- list()
pet_fact     <- array(dim=c(nx,ny))

obj.pet_fact.future[[1]] <- future_pet_facts_calc(crafty_raster[[cyr_id]], CRAFTY_LU_CODE)
names(obj.pet_fact.future)[[1]] <- crafty.yrs[cyr_id]
pet_fact <- getscot.elev(as.matrix(rasterFromXYZ(obj.pet_fact.future[[1]])))

#pet_fact     <- array(dim=c(nx,ny,nyr_f))
# for (i in 1:nyr){
#   obj.pet_fact.future[[i]] <- future_pet_facts_calc(crafty_raster[[i]], CRAFTY_LU_CODE)
#   names(obj.pet_fact.future)[[i]] <- crafty.yrs[i]
#   pet_fact[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.pet_fact.future[[i]])))
# }

obj.pet_fact.future$pet_fact <- pet_fact
obj.pet_fact.future$years    <- start.year:end.year#crafty.yrs 

save(obj.pet_fact.future,file=paste0(path_to_datafolder,"petfactfuture_", year.name,".RData")) # petfact
# ---- End of Step 11A ----



# ==============================================================================
# STEP 12: Load data from previous HDF5 file
# ==============================================================================
hdf5_path      <- "nirams_ii_input.h5"     # path to existing hdf5 file

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

# # Get topographic data from previous HDF5 file and add to list
# topo.data  <- h5read(hdf5_path, 
#                      "/five_km_grids/topographical_data")

# Get idealized time series from previous HDF5 file
ts.data  <- h5read(hdf5_path, 
                   "/time_series")

ts_data <- ts.data$time_series_table[,c(1,2,6,7,8,9,14,3,4,5,14,10,12,13,11)]
names(ts_data) <- c("day","day_of_year","month","or_gr","or_sp","or_wi","or_ot","in_gr","in_sp","in_wi","in_ot","up_gr","up_sp", "up_wi","up_ot")

ts_data_future <- read.csv("ts_data_future.csv", header=TRUE)
ts_data_future <- ts_data_future[,c(1,2,6,7,8,9,14,3,4,5,14,10,12,13,11)]
names(ts_data_future) <- c("day","day_of_year","month","or_gr","or_sp","or_wi","or_ot","in_gr","in_sp","in_wi","in_ot","up_gr","up_sp", "up_wi","up_ot")



## ########################################################################## ##
## Load "fixed" data
## ########################################################################## ##
load(file=paste0(path_to_datafolder,"orgNfuture_", year.name,".RData"))    # organic N
load(file=paste0(path_to_datafolder,"inorgNfuture_", year.name,".RData"))  # inorganic N
load(file=paste0(path_to_datafolder,"uptakefuture_", year.name,".RData"))  # N uptake
load(file=paste0(path_to_datafolder,"petfactfuture_", year.name,".RData")) # petfact


## ########################################################################## ##
## Single R object
## ########################################################################## ##
## Combine into a single R object

objmeta <- objdata <- as.list(NULL)
load(file=paste0(path_to_datafolder,"objClimMonthly",year.name,".RData"))
# --------------------------------------------------------------------------
# CREATE OBJMETA (basically a meta data object with the dates of all data sets)
# --------------------------------------------------------------------------
objmeta$day.year   <- rep(x=start.year:end.year, each=360)# obj.clim.day$year  # add the year for each of the daily climate data set
objmeta$day.julian <- rep(x=1:360, nyrs)# obj.clim.day$day   # add the julian day for each of the daily climate data set

objmeta$mon.year   <- rep(x=start.year:end.year, each=12) #obj.clim.mon.future$year  # add the year for each of the monthly climate data set
objmeta$mon.month  <- rep(1:12, nyrs) #obj.clim.mon.future$month # add the month for each of the monthly climate data set

objmeta$year       <- start.year:end.year # obj.clim.mon.future$emep$years  # add the year for each of the annual data sets

objmeta$onekmgrid.longitude <- obj.clim.mon.future$longitude  # changed!
objmeta$onekmgrid.latitude  <- obj.clim.mon.future$latitude   # changed!

## Add day and month to daily data
objmeta$day.month <- rep( rep(x=1:12,each=30), nyrs)
objmeta$day.day   <- rep(1:30,nms)

#fn <- function(x,y){(y >= x[1]) & (y <= x[2])}

# objmeta$day.month <- objmeta$day.day <- NA

#leapy <- (objmeta$day.year == 2012 | objmeta$day.year == 2016 | objmeta$day.year == 2020)

# for(k in 1:length(objmeta$day.year)){
#   
#   nim <- c(31, 28 + leapy[k], 31, 30, 31, 30, 31, 31, 30, 31, 30, 31) # days in months
#   
#   nim <- cbind(c(1, 1 + cumsum(nim)[-12]), cumsum(nim)) # start and end days in month
#   
#   objmeta$day.month[k] <- which(apply(nim, 1, fn, y = objmeta$day.julian[k]))
#   
#   objmeta$day.day[k]   <- objmeta$day.julian[k] + 1 - nim[objmeta$day.month[k],1]
# }

# --------------------------------------------------------------------------
# CREATE OBJDATA (a list with all the data sets)
# --------------------------------------------------------------------------
# First all non-climatic data sets are added to the list
objdata$onekmgrid.lcms88petfact  <- obj.land.props$lcms88_pet_fact # add default pet_fact to object
objdata$onekmgrid.fc       <- obj.soil$fc # add field capacity to object
objdata$onekmgrid.sat      <- obj.soil$sat # add soil saturation to object
objdata$onekmgrid.calibl   <- obj.soil$calibl # add calibration constant to object
objdata$onekmgrid.calibv   <- obj.soil$calibv # add calibration constant to object
objdata$onekmgrid.bfi      <- obj.soil$bfi # add bfi to object (needed?)
objdata$onekmgrid.spr      <- obj.soil$spr # add spr to object (needed?)

objdata$onekmgrid.petfact  <- obj.pet_fact.future$pet_fact  # pet_fact # add pet_fact to object
objdata$onekmgrid.ndep     <- obj.clim.mon.future$emep$vals # add N deposition data to object

objdata$onekmgrid.orgn.gr  <- obj.orgN.future$or_gr # add organic N data for grassland to object
objdata$onekmgrid.orgn.ot  <- obj.orgN.future$or_ot # add organic N data for other crops to object
objdata$onekmgrid.orgn.sp  <- obj.orgN.future$or_sp # add organic N data for spring crops to object
objdata$onekmgrid.orgn.wi  <- obj.orgN.future$or_wi # add organic N data for winter crops to object

objdata$onekmgrid.inorgn.gr  <- obj.inorgN.future$in_gr # add inorganic N data for grassland to object
objdata$onekmgrid.inorgn.ot  <- obj.inorgN.future$in_ot # add inorganic N data for other crops to object
objdata$onekmgrid.inorgn.sp  <- obj.inorgN.future$in_sp # add inorganic N data for spring crops to object
objdata$onekmgrid.inorgn.wi  <- obj.inorgN.future$in_wi # add inorganic N data for winter crops to object

objdata$onekmgrid.uptake.gr  <- obj.uptake.future$up_gr # add N uptake data for grassland to object
objdata$onekmgrid.uptake.ot  <- obj.uptake.future$up_ot # add N uptake data for other crops to object
objdata$onekmgrid.uptake.sp  <- obj.uptake.future$up_sp # add N uptake data for spring crops to object
objdata$onekmgrid.uptake.wi  <- obj.uptake.future$up_wi # add N uptake data for winter crops to object

objdata$onekmgrid.lat  <- obj.clim.mon.future$lat # add latitude to object (needed?)
objdata$onekmgrid.elev <- obj.clim.mon.future$elev # add elevation to object

#save(objdata,file="objdata_2040_stage1.RData")
rm(obj.orgN.future,obj.inorgN.future,obj.uptake.future,obj.pet_fact.future,obj.soil)

# COPIED FROM STEP 1 AND 2
varnames.day.future <- c("pr", "tasmin", "tasmax")
varnames.mon.future <- c("tas", "tasmin", "tasmax", "rsds", "sfcWind", "hurs")
varnames.mon.future <- c(varnames.mon.future, c("pet.pm", "pet.thorn"))

nvday <- length(varnames.day.future) ; nvmon <- length(varnames.mon.future)

nd <- length(objdata)

varnames.day.new <- c("rainfall", "tasmin", "tasmax")

# add daily climate data to object
for(k in 1:nvday){ 

  objdata[[nd+k]]      <- 0# obj.clim.day[[varnames.day.new[k]]]
  names(objdata)[nd+k] <- paste("onekmgrid.day", varnames.day.new[k], sep=".") 
  #rm(temp.data.day)
}


# # add daily average temperature to object (THIS IS NOT NEEDED FOR NIRAMS SO REMOVED FOR STORAGE REASONS)
# objdata$onekmgrid.day.tas <- (objdata$onekmgrid.day.tasmax + objdata$onekmgrid.day.tasmin) / 2

nd <- length(objdata)
# add monthly climate data to object
for(k in 1:nvmon){ 
  
  objdata[[nd+k]]      <- obj.clim.mon.future[[varnames.mon.future[k]]] 
  names(objdata)[nd+k] <- paste("onekmgrid.mon", varnames.mon.future[k], sep=".") 
}

#save(objdata,file="objdata_2040_stage2.RData")
rm(obj.clim.mon.future)
## ##############################################################################
## Prepare R objects to be saved in H5 format

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
      
      new[[m]] <- objdata[[k]] #[,,j]   # changed! annual data are assumed the same for all years in the future
      
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
      
      new[[m]] <- 0#objdata[[k]][,,j]
      
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
#for(i in 1:length(new)){ new[[i]][is.na(new[[i]])] <- -9999 } # change NA to -9999

summary(factor(unlist(lapply(new, function(x){paste(dim(x),collapse="-")}))))

rm(objdata)
gc()
## ##############################################################################
## Save as H5 file, 26 Aug 2020
## ##############################################################################
# ----------------------------------------------------------
# Change the names of each data set to hdf5 format
# ----------------------------------------------------------
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
# is this working??????
#for(k in 1:length(evar)){ objnames <- gsub(paste0("_", evar[k]), paste0("/", evar[k]), objnames) }
for(k in 1:length(evar_new)){ objnames <- gsub(paste0("_", evar_new[k]), paste0("/", evar_new[k]), objnames) }
## objnames <- paste0("/", objnames)

# ----------------------------------------------------------
# Create h5 file and add groups
# ----------------------------------------------------------
h5filename <- paste0("Processed//nirams-future-calu-inputs-",year.name,".h5")   # changed
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

h5createGroup(h5filename, "one_km_grids/meteorological_data")
h5createGroup(h5filename, "one_km_grids/meteorological_data/daily")
h5createGroup(h5filename, "one_km_grids/meteorological_data/monthly")

h5createGroup(h5filename, "time_series")

yrs <- start.year:end.year

ed <- gsub("-", "_", as.character(meta$Evar[meta$Resolution == "Monthly"]))

for(i in 1:length(ed)){
  
  vi <- paste0("one_km_grids/meteorological_data/monthly/", ed[i])
  
  h5createGroup(h5filename, vi)
}

ed <- gsub("-", "_", as.character(meta$Evar[meta$Resolution == "Daily"]))
ed <- ed[2:4]  # changed - not including mean temperature
for(i in 1:length(ed)){
  
  vi <- paste0("one_km_grids/meteorological_data/daily/", ed[i])
  
  h5createGroup(h5filename, vi)
  
  for(j in 1:length(yrs)){
    
    vij <- paste0(vi, "/", ed[i], "_", yrs[j])
    
    h5createGroup(h5filename, vij)
  }
}

# USE scotdf2matrix BELOW TO FLIP ALL DATA UPSIDE DOWN!!!!!!!!!

# ----------------------------------------------------------
# Add data sets to groups
# ----------------------------------------------------------

# for(k in 1:length(new)){
#   
#   print(paste(k, date()))
#   
#   h5write(scotdf2matrix(new[[k]]), file = h5filename, name = objnames[k]) # flipping all data sets upside down to be consistent with NIRAMS
# }

# This is now done in steps - very clunky. Carefully check than objnames match
#  "entries" in new!
idrain <- grep("one_km_grids/meteorological_data/daily/rainfall",objnames) # id daily rainfall data   80 - 1879
idtmin <- grep("one_km_grids/meteorological_data/daily/min_temp",objnames) # id daily min temp data 1880 - 3679
idtmax <- grep("one_km_grids/meteorological_data/daily/max_temp",objnames) # id daily max temp data 3680 - 5479
idmon  <- grep("one_km_grids/meteorological_data/monthly/",objnames) # id monthly data              5480 - 5959

# STEP 1 - add nontemporal and annual data to h5files (1:79)
for(k in 1:79){
  
  print(paste(k, date()))
  
  tempdat <- new[[k]]
  tempdat[is.na(tempdat)] <- -9999 # change NA to -9999
  
  h5write(scotdf2matrix(tempdat), file = h5filename, name = objnames[k]) # flipping all data sets upside down to be consistent with NIRAMS
}
rm(tempdat)

# STEP 2 - add daily rainfall data
load(paste0(path_to_datafolder,"rainfallDayFuture",year.name,".RData"))
i=1
for(k in idrain){
  
  print(paste(k, date()))
  
  tempdat <- rainfall.day.future[,,i] #new[[k]]
  tempdat[is.na(tempdat)] <- -9999 # change NA to -9999
  
  h5write(scotdf2matrix(tempdat), file = h5filename, name = objnames[k]) # flipping all data sets upside down to be consistent with NIRAMS
  i=i+1
}
rm(rainfall.day.future, tempdat)
gc()

# STEP 3 - add daily min temp data
load(paste0(path_to_datafolder,"tasminDayFuture",year.name,".RData"))
i=1
for(k in idtmin){
  
  print(paste(k, date()))
  
  tempdat <- tasmin.day.future[,,i] #new[[k]]
  tempdat[is.na(tempdat)] <- -9999 # change NA to -9999
  
  h5write(scotdf2matrix(tempdat), file = h5filename, name = objnames[k]) # flipping all data sets upside down to be consistent with NIRAMS
  i=i+1
}
rm(tasmin.day.future, tempdat)
gc()

# STEP 4 - add daily max temp data
load(paste0(path_to_datafolder,"tasmaxDayFuture",year.name,".RData"))
i=1
for(k in idtmax){
  
  print(paste(k, date()))
  
  tempdat <- tasmax.day.future[,,i] #new[[k]]
  tempdat[is.na(tempdat)] <- -9999 # change NA to -9999
  
  h5write(scotdf2matrix(tempdat), file = h5filename, name = objnames[k]) # flipping all data sets upside down to be consistent with NIRAMS
  i=i+1
}
rm(tasmax.day.future, tempdat)
gc()
# STEP 5 - add monthly climate data

for(k in idmon){
  
  print(paste(k, date()))
  
  tempdat <- new[[k]]
  tempdat[is.na(tempdat)] <- -9999 # change NA to -9999
  
  h5write(scotdf2matrix(tempdat), file = h5filename, name = objnames[k]) # flipping all data sets upside down to be consistent with NIRAMS
  #print(max(tempdat))
}
rm(tempdat)



# add time series
#ts.df <- ts.data$time_series_table
# Get idealized time series from previous HDF5 file
ts.data  <- h5read(hdf5_path, 
                   "/time_series")
# # ts_data <- ts.data$time_series_table[,c(1,2,6,7,8,9,14,3,4,5,14,10,12,13,11)]
# # names(ts_data) <- c("day","day_of_year","month","or_gr","or_sp","or_wi","or_ot","in_gr","in_sp","in_wi","in_ot","up_gr","up_sp", "up_wi","up_ot")
# ts_data_future <- read.csv("ts_data_future.csv", header=TRUE)
# ts_data_future <- ts_data_future[,c(1,2,6,7,8,9,14,3,4,5,14,10,12,13,11)]
# names(ts_data_future) <- c("day","day_of_year","month","or_gr","or_sp","or_wi","or_ot","in_gr","in_sp","in_wi","in_ot","up_gr","up_sp", "up_wi","up_ot")
ts.data$time_series_table <- ts_data_future
h5write(ts.data$time_series_table, file = h5filename, name = paste("time_series",names(ts.data),sep="/"),write.attributes=TRUE)

## ##############################################################################




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


# # # ==============================================================================
# # # STEP 6: Load and resample AGCensus raster files
# # # ==============================================================================
# # agcensus.raster <- get.agcensus.data(dir.path=paste0(path_to_inputs,"Raw-data-AGCENSUS/Agcensus_raster_data_scotland_2003_2019/"), 
# #                                      st.yr=2015, end.yr=2019, nirams.grid) 
# # 
# # # ==============================================================================
# # # STEP 7: Load LCM raster files
# # # ==============================================================================
# # lcm_raster <- get.lcm.data(dir.path=paste0(path_to_inputs,"Raw-data-LCM/"), 
# #                            st.yr=2015, end.yr=2021, nirams.grid) 
# # 
# # # ==============================================================================
# # # STEP 8: Assign NIRAMS land use categories
# # # ==============================================================================
# # # The number of years for which agcensus and lcm data are available are not the 
# # # same. Here we are looking at the years 2010-2021, but we only have agcensus
# # # data from 2010-2019 and lcm data for 2015, 2017-2021. For the assignment of 
# # # NIRAMS land uses (for the subsequent N budget calculations) we therefore assume:
# # #
# # # YEARS         AGCENSUS            LCM
# # # 2010-2015:    2010-2015           2015
# # # 2016-2017     2016-2017           2017
# # # 2018-2019     2018-2019           2018-2019
# # # 2020-2021     2019                2020-2021
# # obj.landuse <- list()
# # yrs <- seq(2015,2021, by=1)
# # nyr <- length(yrs)
# # lcm.yrs      <- as.character(c(2015,2017,seq(2017,2021,by=1))) 
# # agcensus.yrs <- as.character(c(seq(2015,2019,by=1),2019,2019))
# # 
# # for (i in 1:nyr){
# #   obj.landuse[[i]] <- assign_nirams_landuse_categories(agcensus.raster[[agcensus.yrs[i]]], 
# #                                                           lcm_raster[[lcm.yrs[i]]], Categories_AGCENSUS_NIRAMS)
# #   names(obj.landuse)[[i]] <- yrs[i]
# # }
# # 
# # 
# # # ==============================================================================
# # # STEP 9: Calculate gridded organic N input from livestock 
# # # ==============================================================================
# # # The organic N input is based on agcensus data on livestock. Because agcensus
# # # data are only available up until 2019, livestock numbers from 2019 are used 
# # # for 2020 and 2021 calculations. 
# # 
# # obj.livestock_N <- list()
# # 
# # for (i in 1:nyr){
# #   ifelse(i>1 && agcensus.yrs[i]==agcensus.yrs[i-1],
# #          {
# #            obj.livestock_N[[i]] <- obj.livestock_N[[i-1]]},
# #          {
# #            obj.livestock_N[[i]] <- Livestock_N(agcensus.raster[[agcensus.yrs[i]]],
# #                                       Categories_AGCENSUS_NIRAMS,
# #                                       JAC_Livestock_N_Per_Animal_InputMatrix)
# #            })
# #   names(obj.livestock_N)[[i]] <- yrs[i]
# # }
# # 
# # obj.orgN_kgha <- list()
# # for (i in 1:nyr){
# #   obj.orgN_kgha[[i]] <- org_n_distrib(obj.livestock_N[[i]], obj.landuse[[i]], 
# #                                       roughgraz_max = 10,impgraz_max = 250, 
# #                                       arable_max = 170)
# #   names(obj.orgN_kgha)[[i]] <- yrs[i]
# #   }
# # 
# # # save orgN in arrays and in similar format as the climate data
# # or_gr <- array(dim=c(nx,ny,nyr))
# # or_sp <- array(dim=c(nx,ny,nyr))
# # or_wi <- array(dim=c(nx,ny,nyr))
# # or_ot <- array(dim=c(nx,ny,nyr))
# # for (i in 1:nyr){
# #   rr <- rasterFromXYZ(obj.orgN_kgha[[i]])
# #   or_gr[,,i] <- getscot.elev(as.matrix(rr$kg_ha_grass))
# #   or_ot[,,i] <- getscot.elev(as.matrix(rr$kg_ha_othercrops))
# #   or_sp[,,i] <- getscot.elev(as.matrix(rr$kg_ha_springcrops))
# #   or_wi[,,i] <- getscot.elev(as.matrix(rr$kg_ha_wintercrops))
# # }
# # 
# # obj.orgN <- list()
# # obj.orgN$or_gr <- or_gr
# # obj.orgN$or_ot <- or_ot
# # obj.orgN$or_sp <- or_sp
# # obj.orgN$or_wi <- or_wi
# # obj.orgN$years <- yrs
# 
# # **********************************************
# ## Option 1
# # obj.orgN.future <- list()
# # obj.orgN.future$or_gr <- apply(or_gr[,,7:12],c(1,2),mean,na.rm=TRUE)
# # obj.orgN.future$or_ot <- apply(or_ot[,,7:12],c(1,2),mean,na.rm=TRUE)
# # obj.orgN.future$or_sp <- apply(or_sp[,,7:12],c(1,2),mean,na.rm=TRUE)
# # obj.orgN.future$or_wi <- apply(or_wi[,,7:12],c(1,2),mean,na.rm=TRUE)
# # obj.orgN.future$years <- start.year:end.year
# 
# # Option 2
# # Get orgN data from previous HDF5 file and add to list
# hdf5_path_new  <- "Processed/nirams-inputs-2023.h5"     # path to existing hdf5 file
# orgN.tmp <- h5read(hdf5_path_new,"/one_km_grids/organic_n/")
# orgN.tmp <- lapply(orgN.tmp,function(x) ifelse( x== -9999,NA,x)) # replace -9999 with NA
# 
# or_gr.future <- Reduce("+",orgN.tmp[7:12])/length(orgN.tmp[7:12]) 
# or_ot.future <- Reduce("+",orgN.tmp[12+7:12])/length(orgN.tmp[12+7:12]) 
# or_sp.future <- Reduce("+",orgN.tmp[2*12+7:12])/length(orgN.tmp[2*12+7:12]) 
# or_wi.future <- Reduce("+",orgN.tmp[3*12+7:12])/length(orgN.tmp[3*12+7:12]) 
# 
# obj.orgN.future <- list()
# obj.orgN.future$or_gr <- scotdf2matrix(or_gr.future)
# obj.orgN.future$or_ot <- scotdf2matrix(or_ot.future)
# obj.orgN.future$or_sp <- scotdf2matrix(or_sp.future)
# obj.orgN.future$or_wi <- scotdf2matrix(or_wi.future)
# obj.orgN.future$years <- start.year:end.year
# save(obj.orgN.future, file="orgNfuture.RData")
# rm(orgN.tmp, or_gr.future, or_ot.future, or_sp.future, or_wi.future)
# # ==============================================================================
# # STEP 10: Calculate gridded inorganic N input and N uptake
# # ==============================================================================
# 
# # # The inorganic N input and uptake are based on land cover.
# # obj.inorgN_uptake <- list()
# # for (i in 1:nyr){
# #   obj.inorgN_uptake[[i]] <- Inorg_N_Application_Uptake(obj.landuse[[i]],IACS_LU_CODE)
# #   names(obj.inorgN_uptake)[[i]] <- yrs[i]
# # }
# # 
# # # save inorganic N and uptake in arrays and in similar format as the climate data
# # in_gr <- array(dim=c(nx,ny,nyr))
# # in_ot <- array(dim=c(nx,ny,nyr))
# # in_sp <- array(dim=c(nx,ny,nyr))
# # in_wi <- array(dim=c(nx,ny,nyr))
# # up_gr <- array(dim=c(nx,ny,nyr))
# # up_ot <- array(dim=c(nx,ny,nyr))
# # up_sp <- array(dim=c(nx,ny,nyr))
# # up_wi <- array(dim=c(nx,ny,nyr))
# # 
# # for (i in 1:nyr){
# #   in_gr[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.inorgN_uptake[[i]]$in_gr)))
# #   in_ot[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.inorgN_uptake[[i]]$in_ot)))
# #   in_sp[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.inorgN_uptake[[i]]$in_sp)))
# #   in_wi[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.inorgN_uptake[[i]]$in_wi)))
# #   up_gr[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.inorgN_uptake[[i]]$up_gr)))
# #   up_ot[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.inorgN_uptake[[i]]$up_ot)))
# #   up_sp[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.inorgN_uptake[[i]]$up_sp)))
# #   up_wi[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.inorgN_uptake[[i]]$up_wi)))
# # }
# # 
# # obj.inorgN <- list()
# # obj.inorgN$in_gr <- in_gr
# # obj.inorgN$in_ot <- in_ot
# # obj.inorgN$in_sp <- in_sp
# # obj.inorgN$in_wi <- in_wi
# # obj.inorgN$years <- yrs
# # 
# # obj.uptake <- list()
# # obj.uptake$up_gr <- up_gr
# # obj.uptake$up_ot <- up_ot
# # obj.uptake$up_sp <- up_sp
# # obj.uptake$up_wi <- up_wi
# # obj.uptake$years <- yrs
# 
# # ******************************
# # # Option 1
# # obj.inorgN.future <- list()
# # obj.inorgN.future$in_gr <- apply(in_gr[,,7:12],c(1,2),mean,na.rm=TRUE) 
# # obj.inorgN.future$in_ot <- apply(in_ot[,,7:12],c(1,2),mean,na.rm=TRUE) 
# # obj.inorgN.future$in_sp <- apply(in_sp[,,7:12],c(1,2),mean,na.rm=TRUE) 
# # obj.inorgN.future$in_wi <- apply(in_wi[,,7:12],c(1,2),mean,na.rm=TRUE) 
# # obj.inorgN.future$years <- start.year:end.year
# # 
# # obj.uptake.future <- list()
# # obj.uptake.future$up_gr <- apply(up_gr[,,7:12],c(1,2),mean,na.rm=TRUE) 
# # obj.uptake.future$up_ot <- apply(up_ot[,,7:12],c(1,2),mean,na.rm=TRUE) 
# # obj.uptake.future$up_sp <- apply(up_sp[,,7:12],c(1,2),mean,na.rm=TRUE) 
# # obj.uptake.future$up_wi <- apply(up_wi[,,7:12],c(1,2),mean,na.rm=TRUE) 
# # obj.uptake.future$years <- start.year:end.year
# 
# 
# # Option 2
# # Get inorgN data from previous HDF5 file and add to list
# hdf5_path_new  <- "Processed/nirams-inputs-2023.h5"     # path to existing hdf5 file
# inorgN.tmp <- h5read(hdf5_path_new,"/one_km_grids/inorganic_n/")
# inorgN.tmp <- lapply(inorgN.tmp,function(x) ifelse( x== -9999,NA,x)) # replace -9999 with NA
# 
# in_gr.future <- Reduce("+",inorgN.tmp[7:12])/length(inorgN.tmp[7:12])
# in_ot.future <- Reduce("+",inorgN.tmp[12+7:12])/length(inorgN.tmp[12+7:12])
# in_sp.future <- Reduce("+",inorgN.tmp[2*12+7:12])/length(inorgN.tmp[2*12+7:12])
# in_wi.future <- Reduce("+",inorgN.tmp[3*12+7:12])/length(inorgN.tmp[3*12+7:12])
# 
# obj.inorgN.future <- list()
# obj.inorgN.future$in_gr <- scotdf2matrix(in_gr.future)
# obj.inorgN.future$in_ot <- scotdf2matrix(in_ot.future)
# obj.inorgN.future$in_sp <- scotdf2matrix(in_sp.future)
# obj.inorgN.future$in_wi <- scotdf2matrix(in_wi.future)
# obj.inorgN.future$years <- start.year:end.year
# save(obj.inorgN.future, file="inorgNfuture.RData")
# rm(inorgN.tmp, in_gr.future, in_ot.future, in_sp.future, in_wi.future)
# 
# 
# # Get N uptake data from previous HDF5 file and add to list
# hdf5_path_new  <- "Processed/nirams-inputs-2023.h5"     # path to existing hdf5 file
# uptake.tmp <- h5read(hdf5_path_new,"/one_km_grids/n_uptake/")
# uptake.tmp <- lapply(uptake.tmp,function(x) ifelse( x== -9999,NA,x)) # replace -9999 with NA
# 
# up_gr.future <- Reduce("+",uptake.tmp[7:12])/length(uptake.tmp[7:12])
# up_ot.future <- Reduce("+",uptake.tmp[12+7:12])/length(uptake.tmp[12+7:12])
# up_sp.future <- Reduce("+",uptake.tmp[2*12+7:12])/length(uptake.tmp[2*12+7:12])
# up_wi.future <- Reduce("+",uptake.tmp[3*12+7:12])/length(uptake.tmp[3*12+7:12])
# 
# obj.uptake.future <- list()
# obj.uptake.future$up_gr <- scotdf2matrix(up_gr.future)
# obj.uptake.future$up_ot <- scotdf2matrix(up_ot.future)
# obj.uptake.future$up_sp <- scotdf2matrix(up_sp.future)
# obj.uptake.future$up_wi <- scotdf2matrix(up_wi.future)
# obj.uptake.future$years <- start.year:end.year
# save(obj.uptake.future, file="uptakefuture.RData")
# rm(uptake.tmp, up_gr.future, up_ot.future, up_sp.future, up_wi.future)
# 
# # ==============================================================================
# # STEP 11: Calculate gridded PET_fact
# # ==============================================================================
# # 
# # obj.pet_fact <- list()
# # pet_fact     <- array(dim=c(nx,ny,nyr))
# # 
# # for (i in 1:nyr){
# #   obj.pet_fact[[i]] <- pet_facts_calc(obj.landuse[[i]], IACS_LU_CODE$PET_Fact)
# #   names(obj.pet_fact)[[i]] <- yrs[i]
# #   pet_fact[,,i] <- getscot.elev(as.matrix(rasterFromXYZ(obj.pet_fact[[i]])))
# # }
# # 
# # obj.pet_fact$pet_fact <- pet_fact
# # obj.pet_fact$years    <- yrs 
# # 
# # # ****************************************************
# # # Option 1
# # obj.pet_fact.future <- list()
# # obj.pet_fact.future$pet_fact <- apply(pet_fact[,,7:12],c(1,2),mean,na.rm=TRUE)
# # obj.pet_fact.future$years    <- start.year:end.year
# 
# # Option 2
# #hdf5_path_new  <- "Processed/nirams-inputs-2023.h5"     # path to existing hdf5 file
# petfact.tmp <- h5read(hdf5_path_new,"/one_km_grids/iacs_pet_facts/")
# petfact.tmp <- lapply(petfact.tmp,function(x) ifelse( x== -9999,NA,x)) # replace -9999 with NA
# petfact.future <- Reduce("+",petfact.tmp[7:12])/length(petfact.tmp[7:12])
# 
# obj.pet_fact.future <- list()
# obj.pet_fact.future$pet_fact <- scotdf2matrix(petfact.future)
# obj.pet_fact.future$years <- start.year:end.year
# rm(petfact.tmp, petfact.future)
# save(obj.pet_fact.future, file="petfactfuture.RData")
# 



