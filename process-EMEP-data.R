# ==============================================================================
# Function to load the annual N deposition data from EMEP. The script has been 
# re-written from scratch compared to the Ina 2017 version, as the grid used 
# by EMEP has changed in coverage, resolution and projection. The deposition 
# data from EMEP is now available at 0.1x0.1 degree (long,lat) resolution. 
# The data can be downloaded in NetCDF format from the EMEP website 
# (https://emep.int/mscw/mscw_moddata.html#NCdata).
# 
# Previously, the data was read from WebDab, which is the EMEP emission database. 
# However, WebDab now only contains the national totals, while the gridded data
# can be downloaded in NetCDF format from the link above. 
#
# Here, the gridded annual N deposition data will be loaded and stored as an  
# array, which can later be saved in the hdf5 format required for NIRAMS. 
#
# Note that both monthly and daily deposition data are now also provided by EMEP.
# These data suggest that there is quite a lot of variation in the deposition 
# over the year so as a FUTURE MODIFICATION it is perhaps worthwhile to use 
# these rather than annual values
#
# MT July 2023
# ==============================================================================

# LOAD AND PROCESS ANNUAL EMEP DATA AND SAVE TO ARRAY

get.emep.data <- function(path.to.files, yr.start=2010, yr.end=2020, nirams.grid){
  # The NetCDF files from EMEP contains modelled data for various compounds, 
  # including gases and aerosols in air, ground level ozone and deposition for ALL
  # of Europe. The data is at a 0.1 degree resolution (long-lat). The whole of 
  # Europe is covered by 1200x520 grid cells.
  #
  # Here we are interested in the dry and wet deposition of oxidized and reduced
  # nitrogen data: 
  # DDEP_OXN_m2Grid    Dry deposition of oxidized nitrogen per m2 grid
  # WDEP_OXN           Wet deposition of oxidized nitrogen
  # DDEP_RDN_m2Grid    Dry deposition of reduced nitrogen per m2 grid
  # WDEP_RDN           Wet deposition of reduced nitrogen
  #
  # INPUT:
  # path.to.files     path to the emep ncdf files
  # yr.start          the first year for which emep data should be loaded
  # yr.end            the last year for which emep data should be loaded
  # nirams.grid       SpatialPixels object with the NIRAMS grid cell center 
  #                   coordinates (or for whatever gridded area that the emep 
  #                   data need to be extracted for)
  # OUTPUT
  # obj.emep          list containing the array (nx,ny,nyr) with the spatial 
  #                   deposition data for each of the years from yr.start to
  #                   yr.end.  
  
  # Open a NetCDF file with EMEP data; here the most recent one is opened. This
  # is only done to read the longitude and latitude for the data for initialization
  # note: the file name for the 2020 data file was modified to be consistent with previous years
  #path.to.files <- "C:/Users/mt40375/Documents/NIRAMS/input-files/Raw-data-deposition/"
  fileName      <- "EMEP01_rv4.45_year.2020met_2020emis_rep2022.nc"  
  ncdf_data     <- nc_open(paste(path.to.files,fileName,sep=""))
  
  # Get the longitude and latitude for this data set. Note this data set covers
  # not only land based grid cells and it covers all of Europe
  emep.lon <- round(ncvar_get(ncdf_data, "lon"),2)
  emep.lat <- round(ncvar_get(ncdf_data, "lat"),2)

  yrs <- seq(yr.start,yr.end, by=1)
  nyr <- length(yrs)
  
  # initialize
  tmp.all <- array(dim=c(nx,ny,nyr))
  
  # We want to extract the values for the NIRAMS grid only (i.e. Scotland). We do
  # this by first creating a raster with the geo-referenced values (in long-lat). 
  # This raster is then re-projected to BNG and the relevant values are extracted 
  # at the NIRAMS grid locations.
  
  # Get data
  for (i in 1:nyr){
    yr_i      <- yrs[i]
    fileName  <- paste0("EMEP01_rv4.45_year.",yr_i,"met_",yr_i,"emis_rep2022.nc")
    # if the file exist in the folder, then extract data. otherwise, reuse the 
    # extracted deposition data from the previous year. This assumes that a file 
    # exists for the specified starting year!
    if( fileName %in% list.files(path.to.files)) { 
      ncdf_data <- nc_open(paste(path.to.files,fileName,sep=""))
      
      # get all N deposition values for Europe
      n_depo    <- ncvar_get(ncdf_data, "DDEP_OXN_m2Grid") + ncvar_get(ncdf_data, "WDEP_OXN") + 
                   ncvar_get(ncdf_data, "DDEP_RDN_m2Grid") + ncvar_get(ncdf_data, "WDEP_RDN")
    
      # Save total N deposition data in a list. 
      # EMEP has data for all of Europe. Here, only those values that roughly 
      # correspond to covering Scotland are extracted.This appears to be a better
      # approach than extracting and re-projecting a raster for the whole of Europe. 
      emep.data   <- list()
      emep.data$x <- emep.lon[205:310]
      emep.data$y <- emep.lat[240:315]
      emep.data$z <- n_depo[205:310,240:315]
    
      # convert from mg/m2 to kg/ha
      emep.data$z <- emep.data$z/100 
    
      # create raster 
      emep.raster <- raster(emep.data)
    
      # re-project raster to BNG
      emep.rasterBNG <- projectRaster(emep.raster,crs = crs(nirams.grid))
    
      # extract N deposition values corresponding to NIRAMS grid 
      emep.ndep <- extract(emep.rasterBNG,nirams.grid)
    
      # save extracted deposition values in matrix
      tmp.ndep  <- matrix(nrow=485,ncol=715,emep.ndep) 
      }
  
    # save to array (UPSIDE DOWN)
    tmp.all[,,i] <- tmp.ndep # tmp.ndep[,ncol(tmp.ndep):1] # saving matrix upside-down to be consistent with format in Python
    }
  
  obj.emep <- as.list(NULL)
  obj.emep[[1]]   <- tmp.all
  names(obj.emep) <- "n_deposition"
  obj.emep$year   <- yrs
  
  return(obj.emep)
  }


# ## Assume deposition for 2021 is the average of previous years
# #mean_depo = apply(obj.emep$n_deposition,c(1,2),mean)
# mean_depo = rowMeans(obj.emep$n_deposition,dims=2)
# obj.emep$n_deposition[,,nyr+1] = mean_depo
# obj.emep$year[nyr+1] = 2021


