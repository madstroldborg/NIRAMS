## ########################################################################## ##
## Functions for loading, processing and saving the climate data needed for 
## NIRAMS. The script was developed by Adam Butler and only very minor changes 
## have been made to this here.
## 
#' @title Get gridded UK climate data from CEDA netCDF files
#' @description Import gridded UK climate data into R, using locally downloaded 
#' netCDF files from CEDA.
## Created 3 August 2020, last modified 4 August 2020
## ########################################################################## ##
#' @param dataset Name of dataset; a character string
#' @param varnames Names of variables to extract; a vector of character values. 
#' @param year.start First year of data to extract; an integer.
#' @param year.end Final year of data to extract; an integer.
#' @param daily Are the data daily or monthly? A logical value (TRUE/FALSE).
#' @param path.root Path to the directory containing the netCDF files; a character string.
#' @param regfn Function to extract required spatial region from either a 2 or 3 dimensional array
#' @param bigval Threshold above which values of the variable should be converted to be missing (NA).
## ########################################################################## ##
#' @return A list, containing climate variables
## ########################################################################## ##
#' @export

# ----------------------------- #
# F U N C T I O N   1 
# ----------------------------- #
get.ceda.data <- function(dataset, varnames, year.start, year.end, daily, path.root, regfn, bigval){
  
  nvars  <- length(varnames)
  years  <- year.start : year.end
  nyears <- length(years)

  ## ---- START INITIALIZATION ----- 
  #   basically to get latitude and longitude and x,y,t dimension of data array
  if(daily){
    
    nt    <- nyears * 365 + sum(is.leap(years))
    dyimv <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31) ## NOTE: no leap years
    dms   <- c(1, 1 + cumsum(dyimv[-12])) ## Julian days for start of each month
    dme   <- cumsum(dyimv) ## Julian days for end of each month
    
    filename.test <- ceda.filename.day(path.root, dataset = dataset, vn = varnames[1], yr = years[1], mn = 1, dyimv = dyimv)
  }
  else{
    
    nt     <- nyears * 12
    months <- rep(1:12, nyears)
    
    filename.test <- ceda.filename.mon(path.root, dataset = dataset, vn = varnames[1], yr = years[1])
  }
  
  ## ############################### ##
  
  tmp <- nc_open(filename.test)
  
  longitude <- regfn(ncvar_get(tmp, "longitude"), mt = FALSE)  # MT COMMENT: is longitude and latitude needed? 
  latitude  <- regfn(ncvar_get(tmp, "latitude"), mt = FALSE)
  
  nx <- dim(longitude)[1] ; ny <- dim(longitude)[2]
  
  ## ---- END OF INITIALIZATION ----
  
  
  ## ############################### ##
  
  obj.clim <- as.list(NULL)
  
  for(i in 1:nvars){  ## loop over variables
    
    tmp.all <- array(dim=c(nx,ny,nt))
    
    nprev <- 0
    
    for(j in 1:nyears){ ## loop over years
      
      if(daily){
        
        dyimv <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31) ## Days per month
        
        if(is.leap(years[j])){ dyimv[2] <- 29 } ## Correct if leap year
        
        dms <- c(1, 1 + cumsum(dyimv[-12])) ## Julian days for start of each month
        
        dme <- cumsum(dyimv) ## Julian days for end of each month
        
        for(k in 1:12){      ## loop over months
          
          filename <- ceda.filename.day(path.root, dataset = dataset, vn = varnames[i], yr = years[j], mn = k, dyimv)
          
          ox <- (dms[k]:dme[k]) + nprev
          
          tmp.each <- ceda.read(filename, varname = varnames[i], bigval = bigval) # read the NetCDF file
          
          tmp.all[,,ox] <- regfn(tmp.each, mt = TRUE) # extract values for Scotland and add to array
        }
        
        nprev <- nprev + sum(dyimv)
      } # end of daily values loop
      else{
        
        filename <- ceda.filename.mon(path.root, dataset = dataset, vn = varnames[i], yr = years[j])
        
        ox <- (1:12)+12*(j - 1)
        
        tmp.each <- ceda.read(filename, varname = varnames[i], bigval = bigval) # read the NetCDF file
        
        tmp.all[,,ox] <- regfn(tmp.each, mt = TRUE)  # extract values for Scotland and add to array
      } # end of monthly values loop
    } # end of loop over years
    
    obj.clim[[i]] <- tmp.all  # add array to output list
  }
  
  names(obj.clim) <- varnames # specify variable names
  
  obj.clim$latitude <- latitude # add latitude to output list
  
  obj.clim$longitude <- longitude # add longitude to output list
  
  if(daily){
    
    t.year <- t.day <- NULL
    
    for(i in 1:nyears){
      
      t.year <- c(t.year, rep(years[i], 365 + is.leap(years[i]))) 
      
      t.day <- c(t.day, 1:(365 + is.leap(years[i])))
    }
    
    obj.clim$year <- t.year # vector with the years of each data set
    
    obj.clim$day <- t.day   # vector with the days of each data set
  }
  else{
    
    obj.clim$year <- rep(years, each = 12) # vector with the years of each data set
    
    obj.clim$month <- months               # vector with the months of each data set
  }
  
  obj.clim
}

## ########################################################################## ##
## Functions to find relevant filenames for CEDA climate data
## Written from scratch 3 Aug 2020, as now using netCDF rather than ASCII files
## Written 4 Aug 2020 to make "dataset" a user-specific argument, for maximum
##   flexibility
## ########################################################################## ##

## Filenames for monthly data:

ceda.filename.mon <- function(path.root, dataset, vn, yr){ paste0(path.root, vn, "_", dataset, "_", yr, "01-", yr, "12.nc") }

## Filenames for daily data:

ceda.filename.day <- function(path.root, dataset, vn, yr, mn, dyimv){ 
  
  mnl <- paste0((c("", "0"))[1 + (substr(mn,2,2) == "")], mn)
  
  dyim <- (dyimv)[mn]

  if(is.leap(yr) & (mn == 2)){ dyim <- 29 }
    
  paste0(path.root, vn, "_", dataset, "_", yr, mnl, "01-", yr, mnl, dyim, ".nc") 
}

## ########################################################################## ##

is.leap <- function(yr){ ((yr - 4 * floor(yr/4)) < 0.1) }

## ########################################################################## ##
## Created 3 August 2020, last modified 4 August 2020
## Note: use to calculate differences of a monthly CEDA variable
## ########################################################################## ##
#' @title Calculate 2 step differences in the 3rd dimension of a 3d array
#' @param z A three-dimensional array
#' @return A three-dimensional array, of the size of `z`
## ########################################################################## ##
#' @export

s2diff3d3 <- function(z){ 
  
  dz <- array(dim=dim(z)) ; nm <- dim(z)[3] ; dz[,,-c(1,nm)] <- z[,,-(1:2)] - z[,,-((nm-1):nm)] ; dz
}

## ########################################################################## ##
## Created 3 August 2020, last modified 4 August 2020
## ########################################################################## ##
#' @title Create a H5 file containing CEDA climate data
#' @description Read in a CEDA netCDF file, and extract a variable 'varname'
#' @param filename Name of netCDF file to read; a character string
#' @param varname Name of variable to read; a character string. Must match to one of the variable names within the netCDF file.
#' @param bigval Threshold above which values of the variable should be converted to be missing (NA).
#' @return An R object containing the variable 'varname', as extracted from the netCDF file.
## ########################################################################## ##
#' @export

ceda.read <- function(filename, varname, bigval){
  
  tmp.each <- nc_open(filename)  # open NetCDF file
  
  tmp.each <- ncvar_get(tmp.each, varid = varname) #  get variable
  
  tmp.each[tmp.each > bigval] <- NA
  
  tmp.each
}




## ########################################################################## ##
#' @title Create a H5 file containing CEDA climate data
#' @description Create a H5 file containing daily or monthly CEDA climate data, in the format required by NIRAMS
#' @param obj.clim
#' @param h5filename The name of the H5 file to be created; a character string. The filename must have the extension ".h5".
#' @param varnames Names of variables to be included in the H5 file; a vector of character values. These are the names as they appear in the object 'obj.clim'. 
#' @param varalt Names of variables to be included in the H5 file; a vector of character values of the same length as "varnames". These are the names (in the same order as the entires to 'varnames') as they will appear in the H5 file.
#' @param is.tvv Which variables are time varying? A vector of logical value (TRUE/VALUE) of the same length as "varnames".
#' @param daily Are the data daily or monthly? A logical value (TRUE/FALSE).
#' @return Empty (NULL); the function is run for the side effect that it creates a H5 file.
## ########################################################################## ##
#' @export

# ## MT comment 2023: This function is not used in the end in the 
# ##                  process-climate-and-PET-run.R script. 
# ceda2h5 <- function(obj.clim, h5filename, varnames, varalt, is.tvv, daily){ 
#   
#   h5createFile(h5filename)
#   
#   nv <- length(varnames)
#   
#   tvn <- c("month", "day")[1 + daily]
#   
#   tvar <- obj.clim[[which(names(obj.clim) == tvn)]]
#   
#   nt <- length(tvar)
#   
#   for(i in 1:nv){
#     
#     print(paste(i,"/",nv,"     ",date()))
#     
#     grpname <- varalt[i]
#     
#     h5createGroup(h5filename, grpname)
#     
#     k <- which(names(obj.clim) == varnames[i])
#     
#     if(is.tvv[i]){
#       
#       for(j in 1:nt){
#       
#         dat <- (obj.clim[[k]])[,,j]
#       
#         objname <- paste0(grpname, obj.clim$year[j], "_", tvar[j], sep="")
#       
#         h5write(dat, h5filename, objname)
#       }
#     }
#     else{
#       
#       dat <- (obj.clim[[k]])
#       
#       objname <- paste0(grpname, "_all")
#       
#       h5write(dat, h5filename, objname)
#     }
#         
#   }
#   
#   NULL
# }


## ############################################################################
#' @title Calculate PET using the Penman-Monteith equation
#' @description Calculate PET using the Penman-Monteith equation
#' @details A conversion of the Python code created by James
#'   Sample. 
#'
## Steps in conversion from python code:
##   A. create this separate function for core PM calculations
##   B. additionally included steps 0.1, 0.2, 0.3
##   C. but dropped step 4.1, which is now done before inputs
##   D. calculation of 'month.dict' also now done before this
##   E. also tmean.diff (tmean.diff <- tmean.nxt-tmean.pre)
## 
## Key bits of R to Python conversion:
##   1. drop "np." from function names
##   2. change "arccos" to "acos"
##   3. change "**" to "^"
##   4. change "_" to "."
##   5. change "np.pi" to "pi"
##   6. renamed "c" to "cetr"
##
## Additional changes 3 Aug 2020: renamed variables to match
##   CEDA naming conventions:
##
##    tmean > tasmean
##    tmax > tasmax
##    tmin > tasmin
##    vap.p > vp
##    u.10 > sfcWind
##    rel.hum > hurs
##    sun.hrs > sun
## ########################################################################## ##
## Created 3 August 2020, last modified 10 August 2020
## ########################################################################## ##
#' @param obj.clim.mon A list containing monthly climate data.
#'    The list contains monthly spatio-temporal climate data for
#'    climate variables, as named in the Hadley centre gridded data. Each of these
#'    variables is a three-dimensional array: 'tasmax' (maximum temperature),
#'    'taxmin' (minimum temperature), 'tas' (mean temperature), 'pv' (vapour pressure),
#'    'fscWind' (mean windspeed at 10m), 'hurs' (relative humidity), 'sun' (sunshine hours).
#'    Also contains 'lat' (spatial grid, with latitude values), 'evel' (spatial grid with 
#'    elevation values) and 'month'.
#' @param pars.pet.pm A data frame, with a single row, containing parameters to use for P-M
#'    equation calculations; contains columns: 
#'    'a.s' (Angstrom values: Fraction of et_rad reaching the earth on overcast days), 
#'    'b.s' (Angstrom values: (a_s + b_s) is fraction of et_rad reaching the earth on a clear day) 
#'    and 'alb.coef' (Albedo coefficient)'
#' @return    
## ########################################################################## ##
#' @export

calc.pet.pm <- function(obj.clim.mon, pars.pet.pm){
  
  zoo <- function(z,k){ if(is.array(z)){ if(length(dim(z)) > 2){ z <- z[,,k] } } ; z }
  
  nm <- dim(obj.clim.mon$tas)[3] # total number of months
  
  pet.pm <- array(dim=dim(obj.clim.mon$tas))
  
  for(k in 1:nm){
    
    ## ######################################################## ##
    ## Extract a single month of data from monthly CEDA data
    ## Rewritten 10 Aug 2020
    
    datk <- lapply(obj.clim.mon, zoo, k = k) 
    
    datk$year <- datk$year[k] ; datk$month <- datk$month[k] 
    
    ## ######################################################## ##
    ## PM calculation of PET
        
    pet.pm[,,k] <- calc.singlemonth.pet.pm(datk, pars.pet.pm = pars.pet.pm)
    
    ## ######################################################## ##
  }
  
  pet.pm
}

calc.singlemonth.pet.pm <- function(dat, pars.pet.pm){
  
  ## ##################################################################
  ## Extract data (not in python - add because in R 'dat' is a data frame)
  ## 3 August 2020: changed so that it now takes CEDA names as inputs
  ## 4 August 2020: "month" now taken from "dat"
  
  a.s <- pars.pet.pm$a.s
  b.s <- pars.pet.pm$b.s
  alb.coef <- pars.pet.pm$alb.coef
  
  tmean.cur  <- dat$tas     # mean temp [degC]
  tmax.cur   <- dat$tasmax  # max temp [degC]
  tmin.cur   <- dat$tasmin  # min temp [degC]
  tmean.diff <- dat$tasdiff # temp difference [degC]
  vap.p   <- dat$pv         # partial pressure of water vapour [hPa] 
  u.10    <- dat$sfcWind    # wind speed at 10m [m/s]
  rel.hum <- dat$hurs       # relative humidity [%]
  sun.hrs <- dat$sun        # sunshine hours [hr]
  
  latitude <- dat$lat
  elev     <- dat$elev
  
  month <- dat$month
  
  ## ##################################################################
  ## Step 0 - preliminaries
  
  ## Step 0.1: month length (as in Python version)
  
  months.dict <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  
  month.dict <- months.dict[month]
  
  ## Step 0.2: calculate Julian day of year midpoint
  ## (NOTE: entirely rewritten for R)
  
  days.of.yr <- c(0, cumsum(months.dict)[-12]) + 15
  
  day.of.yr <- days.of.yr[month]
  
  # Step 0.3: Convert latitude decimal degrees to radians
  lat.rad <- pi*latitude / 180.0
  
  # Step 0.4: Atmos P at altitude (in kPa)
  atm.p <- 101.3*(((293-0.0065*elev)/293)^5.26)
  
  # Step 0.5: Psychrometric constant (kPa/degC)
  psy.c <- 6.65E-4 * atm.p
  
  ## ##################################################################
  ## Step 1 - Convert units
  
  ## 1.2: vap.p to kPa
  vap.p <- vap.p / 10
  
  ## 1.3: u.10 to m/s
  #u.10 <- 0.5144 * u.10   ## MT 2023: Step removed. Wind speed data is already in m/s, not in knots
  
  ## ##################################################################
  ## Step 2 - Estimate wind speed at 2m above ground from u.10
  
  ## 2.1
  u.2 <- (4.87*u.10)/(log(678-5.42))
  
  ## ##################################################################
  ## Step 3 - Estimate key atmospheric parameters
  
  # 3.3: Slope of saturation vapour pressure curve (kPa/degC)
  slp.sat.vp<-((4098*0.6108*exp(17.27*tmean.cur/(tmean.cur+237.3)))/
                 ((tmean.cur+237.3)^2))
  
  # 3.4: Saturation vapour pressure for min and max T
  e.min <- 0.6108*exp(17.27*tmin.cur/(tmin.cur+237.3))
  e.max <- 0.6108*exp(17.27*tmax.cur/(tmax.cur+237.3))
  
  # 3.5a: Modified 27/11/2014. Average saturation vapour pressure (kPa)
  e.s <- (e.min+e.max)/2.0
  
  # 3.5b: Added 27/11/2014. Calculate e.a (kPa)
  e.a <- e.s*rel.hum/100.
  
  # 3.6: Modified 27/11/2014. Vapour pressure deficit (kPa)
  vap.p.def <- e.s - e.a
  
  # ##################################################################
  # Step 4 - Estimate key radiative parameters
  
  # 4.2: Inverse relative Earth-Sun distance
  rel.dist <- 1+0.033*cos(2*pi*day.of.yr/365)
  
  # 4.3: Solar declination (radians)
  sol.dec <- 0.409*sin((2*pi*day.of.yr/365)-1.39)
  
  # 4.4: Sunset hour angle
  set.hr <- acos(-1.*tan(lat.rad)*tan(sol.dec))
  
  # 4.5: Extra-terrestrial radiation (MJ m-2 day-1)
  cetr <- 0.082 * 24 * 60 / pi
  et.rad <- cetr * rel.dist * (set.hr*sin(lat.rad)*sin(sol.dec) +
                                 cos(lat.rad)*cos(sol.dec)*sin(set.hr))
  
  # 4.6: Daylight hours
  day.hrs <- 24*set.hr/pi
  
  # 4.7: Solar radiation (MJ m-2 day-1)  
  #      MT 2023: CHANGED! sun.hrs is monthly sunshine hours so needs to be 
  #      divided by number of days in month
  # sol.rad <- (a.s+(b.s*sun.hrs/day.hrs))*et.rad
  sol.rad <- (a.s+(b.s*(sun.hrs/month.dict)/day.hrs))*et.rad
  
  # 4.8: Clear-sky solar radiation (MJ m-2 day-1) 
  #      MT 2023: Small change so fraction of Ra reaching earth on clear days 
  #      equals as+bs
  #clr.sky.rad <- (0.75+(elev*2E-5))*et.rad
  clr.sky.rad <- (a.s+b.s+(elev*2E-5))*et.rad
  
  # 4.9: Net shortwave radiation (MJ m-2 day-1) 
  sw.rad <- (1-alb.coef)*sol.rad
  
  # 4.10: Modified 27/11/2014. Net longwave radiation (MJ m-2 day-1) 
  # First convert Tmin and Tmax to kelvin
  tmin.K <- tmin.cur + 273.16
  tmax.K <- tmax.cur + 273.16
  lw.rad <- (4.903E-9 * ((tmin.K^4+tmax.K^4)/2) *
               (0.34 - 0.14*(e.a^0.5)) *
               ((1.35*sol.rad/clr.sky.rad)-0.35))
  
  # 4.11: Net radiation (MJ m-2 day-1) 
  net.rad <- sw.rad - lw.rad
  
  # ##################################################################
  # Step 5 - Estimate soil heat flux (MJ m-2 day-1)
  
  # 5.1: Soil heat flux. Assumes soil heat cap <- 2.1 MJ/m3/C
  #      and effective soil depth of 2m
  
  soil.heat <- 0.07*tmean.diff
  
  # ##################################################################
  # Step 6 - Estimate PET for grassland in mm/day
  
  # 6.1: The Penman-Monteith equation (av. PET over month in mm/day)
  # Currently includes a fudge to set any values less than zero to
  # zero. Need to investigate this!
  # MT 2023 notes on units: 
  #   0.408 (kg MJ-1). Inverse of the latent heat of vaporization (i.e. energy required to evaporate 1kg of water; 2.45 MJ/kg)
  #   slp.sat.vp (kPa degC-1)
  #   net.rad (MJ m-2 day-1)
  #   soil.heat (MJ m-2 day-1)
  #   psy.c (kPa degC-1)
  #   900 (s kg K day-1 kJ-1)
  #   tmean.cur (degC)
  #   u.2 (m s-1)
  #   vap.p.def (kPa)
  #   0.34 (s m-1). This is based on rs/ra = rs*u.2/208 and assuming rs=70 s m-1
   
  pet <- (((0.408*slp.sat.vp*(net.rad-soil.heat)) +
             (psy.c*(900/(tmean.cur+273.16))*u.2*vap.p.def)) /
            (slp.sat.vp + psy.c*(1+0.34*u.2)))
  
  pet[pet<0] <- 0
  
  # ##################################################################
  ## Step 7 - Convert to monthly totals i.e. mm/month
  
  pet.mon <- pet * month.dict
  
  # ##################################################################
  
  pet.mon
}

## End of Penman-Monteith #####################################################

## ############################################################################
#' @title Compute PET using Thorntwaite equation
#' @description Compute PET using Thornthwaite equation from monthly mean temperature data
## ########################################################################## ##
#' @details From James Sample: "The Penman-Monteith equation is the generally accepted standard for PET.
#' However, the Thornthwaite equation is simpler and will serve as a useful
#' check of the P-M estimates. The standard Thornthwaite equation can only be
#' applied to monthly data (not daily as I had hoped), but monthly could still
#' be useful. The Thornthwaite equation generally under-predicts PET compared to P-M."
#'
#' "The Thornthwaite equation is ET in mm = 16c(10T/I)^a, where I and a
#' are constants depending on monthly average temperatures in the year
#' of interest. c is the day length correction factor."
## ########################################################################## ##
## 3 August 2020: changed so that meantemp for single year calcs is array [nx,ny,12]
##  rather than [ngrid,12], and renamed "meantemp" to "tasmean"
## ########################################################################## ##
#' @param tasmean A three dimensional array containing gridded monthly mean 
#'   temperature values. The first two dimensions refer to the spatial grid, and
#'   the final dimension to months.
#' @param nuy Number of years; an integer.
#' @return A three dimensional array, of the same size as `tasmean`, containing estimated
#'   monthly PET values for each square on the grid.   
## ########################################################################## ##
#' @export

calc.pet.thorn <- function(tasmean, nuy){
  
  pet.thorn <- array(dim=dim(tasmean))
  
  for(j in 1:nuy){
    
    jk <- (1:12)+12*(j - 1)
    
    pet.thorn[,,jk] <- calc.singleyear.pet.thorn(tasmean = tasmean[,,jk])
  }
  
  pet.thorn
}

calc.singleyear.pet.thorn <- function(tasmean){
  
  ## ################################
  
  out <- array(dim=dim(tasmean))
  
  ## ################################
  ## Calculate I
  ## I = sum(i.month), month=1..12
  ## i.month = (T.month/5)^1.514
  ## ################################
  
  t.month <- tasmean
  
  t.month[t.month < 0] <- 0
  
  i.month <- (t.month/5)^1.514
  
  Io <- apply(i.month, 1:2, sum) ## 4 Aug 2020: changed from "1" to "1:2"
  
  ## ################################
  ## Calculate a
  ## a = (6.75x10^-7)I^3 - (7.71x10^-5)I^2 + (1.79x10^-2)I + 0.49
  ## ################################
  
  a <- (6.75E-7)*Io**3 - (7.71E-5)*Io**2 + (1.79E-2)*Io + 0.49
  
  ## ################################
  # Now estimate the PET for each month
  ## ################################
  
  # correction factor (where are these from?)
  day.dict <- c(0.65, 0.73, 0.99, 1.15, 1.37, 1.43, 1.43, 1.27, 1.03, 0.88, 0.68, 0.59)
  
  for(j in 1:12){
    
    tmean.cur <- tasmean[,,j] ## 3 Aug 2020: tasmean[,j] now tasmean[,,j]
    
    tmean.cur[tmean.cur < 0] <- 0
    
    # Calculate the Thornthwaite PET in mm (mm/month)
    out[,,j] <- day.dict[j] * 16 * (10 * tmean.cur/Io)^a ## 3 Aug 2020: out[j,] now out[,,j]
  }
  
  ## ################################
  
  if(all(is.na(out))){ browser() }
  
  out
}

## End of Thornwaite ########################################################## 

## ##############################################################################
## Created 12 August 2020: function to add in deposition data from EMEP
## ##############################################################################

add.emep.data <- function(obj, path.emep, pat, groupname, variablename){
  
  filenames <- list.files(path.emep, pattern = pat, full.names = TRUE)
  
  nf <- length(filenames)
  
  for (k in 1:nf){
    
    year <- str_extract(basename(filenames[k]), "[0-9]{4}")
    
    kn <- paste(groupname, variablename, year, sep="")    
    
    obj[kn] <- as.matrix(read.table(filenames[k], sep = " ", skip = 6))
  }
  
  obj
}


## ##############################################################################
## F U N C T I O N S   F O R   F U T U R E   D A T A
## ##############################################################################

# **********************************************************
# FUNCTION1: LOAD  FUTURE CLIMATE  ---- 
# **********************************************************
get.future.ceda.data <- function(dataset, varnames, year.start, year.end, 
                                 daily, path.root, regfn, bigval){
  # note that the future dataset simulate daily climate data where all months 
  # are 30 days long
  
  nvars  <- length(varnames)
  years  <- year.start : year.end
  nyears <- length(years)
  
  ## ---- START INITIALIZATION -----
  #  basically to get latitude and longitude and x,y,t dimension of data array
  if(daily){
    
    nt    <- nyears*360 #nyears * 365 + sum(is.leap(years))
    dyimv <- c(30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30) ## NOTE: no leap years
    dms   <- c(1, 1 + cumsum(dyimv[-12])) ## Julian days for start of each month
    dme   <- cumsum(dyimv) ## Julian days for end of each month
    
    filename.test <- future.ceda.filename.day(path.root, dataset = dataset, 
                                              vn = varnames[1], yr = years[1], 
                                              mn = 1, dyimv = dyimv)
  }
  else{
    
    nt     <- nyears * 12
    months <- rep(1:12, nyears)
    
    filename.test <- future.ceda.filename.mon(path.root, dataset = dataset, vn = varnames[1])
  }
  
  ## ************************************
  
  tmp <- nc_open(filename.test)
  
  longitude <- regfn(ncvar_get(tmp, "lon"), mt = FALSE)  # MT COMMENT: is longitude and latitude needed? 
  latitude  <- regfn(ncvar_get(tmp, "lat"), mt = FALSE)
  nx <- dim(longitude)[1] ; ny <- dim(longitude)[2]
  
  ## ---- END OF INITIALIZATION ----
  
  ## ************************************ 
  
  obj.clim <- as.list(NULL)
  
  for(i in 1:nvars){  ## loop over variables
    
    tmp.all <- array(dim=c(nx,ny,nt))
    
    nprev <- 0
    
    for(j in 1:nyears){ ## loop over years
      
      if(daily){
        
        dyimv <- c(30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30) 
        
        dms <- c(1, 1 + cumsum(dyimv[-12])) ## Julian days for start of each month
        
        dme <- cumsum(dyimv) ## Julian days for end of each month
        
        for(k in 1:12){      ## loop over months
          
          filename <- future.ceda.filename.day(path.root, dataset = dataset, 
                                               vn = varnames[i], yr = years[j], 
                                               mn = k, dyimv)
          
          ox <- (dms[k]:dme[k]) + nprev
          
          tmp.each <- ceda.read(filename, varname = varnames[i], bigval = bigval) # read the NetCDF file
          
          tmp.all[,,ox] <- regfn(tmp.each, mt = TRUE) # extract values for Scotland and add to array
          rm(tmp.each)
        }
        
        nprev <- nprev + sum(dyimv)
      } # end of daily values loop
      else{
        
        filename <- future.ceda.filename.mon(path.root, dataset = dataset, vn = varnames[i])
        
        ox <- (1:12)+12*(j - 1)
        
        tmp.each <- chessscape.mon.read(filename, varname = varnames[i], 
                                        bigval = bigval, years[j], years[j]) # read the NetCDF file
        
        tmp.all[,,ox] <- regfn(tmp.each, mt = TRUE)  # extract values for Scotland and add to array
        rm(tmp.each)
        
      } # end of monthly values loop
    } # end of loop over years
    
    obj.clim[[i]] <- tmp.all  # add array to output list
  }
  
  names(obj.clim) <- varnames # specify variable names
  
  obj.clim$latitude <- latitude # add latitude to output list
  
  obj.clim$longitude <- longitude # add longitude to output list
  
  if(daily){
    
    t.year <- t.day <- NULL
    
    for(i in 1:nyears){
      
      t.year <- c(t.year, rep(years[i], 360))
      
      t.day <- c(t.day, 1:360)
    }
    
    obj.clim$year <- t.year # vector with the years of each data set
    
    obj.clim$day <- t.day   # vector with the days of each data set
  }
  else{
    
    obj.clim$year <- rep(years, each = 12) # vector with the years of each data set
    
    obj.clim$month <- months               # vector with the months of each data set
  }
  
  obj.clim
}

# End of Function 1 ----

# **********************************************************************
# FUNCTIONS 2: find relevant filenames for ChessScape climate data ----
# **********************************************************************
# Written Sept 2023 

## Filenames for monthly data:
future.ceda.filename.mon <- function(path.root, dataset, vn){ paste0(path.root, dataset, "_",vn,"_uk_1km_monthly_19801201-20801130.nc") }

## Filenames for daily data:
future.ceda.filename.day <- function(path.root, dataset, vn, yr, mn, dyimv){ 
  
  mnl <- paste0((c("", "0"))[1 + (substr(mn,2,2) == "")], mn)
  
  dyim <- (dyimv)[mn]
  
  #if(is.leap(yr) & (mn == 2)){ dyim <- 29 }
  
  paste0(path.root, dataset, "_" ,vn, "_uk_1km_daily_", yr, mnl, "01-", yr, mnl, dyim, ".nc") 
}

# End of functions 2 ----


# ******************************************************************************
# FUNCTIONS 3: Read monthly climate data from ChessScape netCDF file ----
# ******************************************************************************
#' @title Read monthly climate data from ChessScape netCDF file
#' @description Read in a ChessScape netCDF file, and extract a variable 'varname'
#' @param filename Name of netCDF file to read; a character string
#' @param varname Name of variable to read; a character string. Must match to 
#' one of the variable names within the netCDF file.
#' @param bigval Threshold above which values of the variable should be converted 
#' to be missing (NA).
#' @return An R object containing the variable 'varname', as extracted from the netCDF file.
#' @export

chessscape.mon.read <- function(filename, varname, bigval, yr.st, yr.end){
  # The ChessScape data include daily and monthly data from 1st December 1980
  # to 30th November 2080, i.e. a total of 1200 months and 36000 days.
  # Here we extract only the monthly data for the selected years
  # Index:    Month:    Year:
  # 1         12        1980
  # 2         01        1981
  # 3         02        1981
  # ...
  # 1200      11        2080
  
  tmp.each <- nc_open(filename)  # open NetCDF file
  v1       <- tmp.each$var[[varname]]
  varsize  <- v1$varsize
  ns       <- 1+(yr.st-1981)*12+1  # index for first month
  ne       <- 1+(yr.end-1981)*12+12  # index for last month
  nm       <- length(ns:ne)
  start    <- c(1,1,ns)
  count    <- varsize
  count[3] <- nm
  tmp.each <- ncvar_get(tmp.each, varid = varname, start=start, count=count) #  get variable
  tmp.each[tmp.each > bigval] <- NA
  tmp.each
}


chessscape.read.specific.mon <- function(filename, varname, regfn, bigval, yr, mm){
  # function for extracting chessScape for a given month and year
  # The ChessScape data include daily and monthly data from 1st December 1980
  # to 30th November 2080, i.e. a total of 1200 months and 36000 days.

  # Index:    Month:    Year:
  # 1         12        1980
  # 2         01        1981
  # 3         02        1981
  # ...
  # 1200      11        2080
  
  tmp.each <- nc_open(filename)  # open NetCDF file
  v1       <- tmp.each$var[[varname]]
  varsize  <- v1$varsize
  # find index for given month
  mi       <- 1+(yr-1981)*12+mm  # index for first month
  start    <- c(1,1,mi)
  count    <- varsize
  count[3] <- 1
  tmp.each <- ncvar_get(tmp.each, varid = varname, start=start, count=count) #  get variable
  tmp.each[tmp.each > bigval] <- NA
  tmp.each <- regfn(tmp.each, mt = FALSE)
  tmp.each
}

# End of Functions 3 ----


# ******************************************************************************
# FUNCTIONS 4: Calculate future PET using the Penman-Monteith equation ----
# ******************************************************************************
#' @title Calculate future PET using the Penman-Monteith equation
#' @description Calculate future PET using the Penman-Monteith equation. 
#' Modified Penman Monteith equation to calculate future PET based 
#' on climate data projection from CHESS-SCAPE. The modification was needed 
#' because CHESS-SCAPE does not include exactly the same climatic variables 
#' as HAD-UK which was used for calculating PET for the historical NIRAMS runs.
#' Created September 2023.
#' @param obj.clim.mon A list containing monthly climate data.
#'    The list contains monthly spatio-temporal climate data for
#'    climate variables, as named in the CHESS-SCAPE gridded data. Each of these
#'    variables is a three-dimensional array: 'tasmax' (maximum temperature),
#'    'taxmin' (minimum temperature), 'tas' (mean temperature), 
#'    'sfcWind' (mean windspeed at 10m), 'hurs' (relative humidity), 'rsds' (shortwave radiation).
#'    Also contains 'lat' (spatial grid, with latitude values), 'evel' (spatial grid with 
#'    elevation values) and 'month'.
#' @param pars.pet.pm A data frame, with a single row, containing parameters to use for P-M
#'    equation calculations; contains columns: 
#'    'a.s' (Angstrom values: Fraction of et_rad reaching the earth on overcast days), 
#'    'b.s' (Angstrom values: (a_s + b_s) is fraction of et_rad reaching the earth on a clear day) 
#'    and 'alb.coef' (Albedo coefficient)'
#' @return    
## ********************************************
#' @export

calc.pet.pm.future <- function(obj.clim.mon, pars.pet.pm){
  
  zoo <- function(z,k){ if(is.array(z)){ if(length(dim(z)) > 2){ z <- z[,,k] } } ; z }
  
  nm <- dim(obj.clim.mon$tas)[3] # total number of months
  
  pet.pm <- array(dim=dim(obj.clim.mon$tas))
  
  for(k in 1:nm){
    
    ## ***********************************************************
    ## Extract a single month of data from monthly ChessScape data
    
    datk <- lapply(obj.clim.mon, zoo, k = k) 
    
    datk$year <- datk$year[k] ; datk$month <- datk$month[k] 
    
    ## ***********************************************************
    ## PM calculation of PET
    
    pet.pm[,,k] <- calc.singlemonth.pet.pm.future(datk, pars.pet.pm = pars.pet.pm)
    
    ## ***********************************************************
  }
  
  pet.pm
}

calc.singlemonth.pet.pm.future <- function(dat, pars.pet.pm){
  
  
  ## *********************************************************************
  ## Extract data 
  
  a.s <- pars.pet.pm$a.s
  b.s <- pars.pet.pm$b.s
  alb.coef <- pars.pet.pm$alb.coef
  
  tmean.cur  <- dat$tas     # mean temp [degC]
  tmax.cur   <- dat$tasmax  # max temp [degC]
  tmin.cur   <- dat$tasmin  # min temp [degC]
  tmean.diff <- dat$tasdiff # temp difference [degC]
  u.10       <- dat$sfcWind # wind speed at 10m [m/s]
  rel.hum    <- dat$hurs    # relative humidity [%]
  shortwave.rad <- dat$rsds # shortwave radiation [W m-2]
  
  latitude <- dat$lat
  elev     <- dat$elev
  
  month <- dat$month
  
  ## ***********************************************************
  ## Step 0 - preliminaries
  
  ## Step 0.1: month length (as in Python version)
  #months.dict <- c(30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30)
  months.dict <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  month.dict <- months.dict[month]
  
  ## Step 0.2: calculate Julian day of year midpoint
  days.of.yr <- c(0, cumsum(months.dict)[-12]) + 15
  
  day.of.yr <- days.of.yr[month]
  
  # Step 0.3: Convert latitude decimal degrees to radians
  lat.rad <- pi*latitude / 180.0
  
  # Step 0.4: Atmos P at altitude (in kPa)
  atm.p <- 101.3*(((293-0.0065*elev)/293)^5.26)
  
  # Step 0.5: Psychrometric constant (kPa/degC)
  psy.c <- 6.65E-4 * atm.p
  
  ## *******************************************************************
  ## Step 1 - Convert units
  
  ## 1.2: vap.p to kPa
  #vap.p <- vap.p / 10
  
  ## 1.3: u.10 to m/s
  #u.10 <- 0.5144 * u.10   ## MT 2023: Step removed. Wind speed data is already in m/s, not in knots
  
  ## ******************************************************************
  ## Step 2 - Estimate wind speed at 2m above ground from u.10
  
  ## 2.1
  u.2 <- (4.87*u.10)/(log(678-5.42))
  
  ## ******************************************************************
  ## Step 3 - Estimate key atmospheric parameters
  
  # 3.3: Slope of saturation vapour pressure curve (kPa/degC)
  slp.sat.vp<-((4098*0.6108*exp(17.27*tmean.cur/(tmean.cur+237.3)))/
                 ((tmean.cur+237.3)^2))
  
  # 3.4: Saturation vapour pressure for min and max T
  e.min <- 0.6108*exp(17.27*tmin.cur/(tmin.cur+237.3))
  e.max <- 0.6108*exp(17.27*tmax.cur/(tmax.cur+237.3))
  
  # 3.5a: Modified 27/11/2014. Average saturation vapour pressure (kPa)
  e.s <- (e.min+e.max)/2.0
  
  # 3.5b: Added 27/11/2014. Calculate e.a (kPa)
  e.a <- e.s*rel.hum/100.
  
  # 3.6: Modified 27/11/2014. Vapour pressure deficit (kPa)
  vap.p.def <- e.s - e.a
  
  ## ******************************************************************
  # Step 4 - Estimate key radiative parameters
  
  # 4.2: Inverse relative Earth-Sun distance
  rel.dist <- 1+0.033*cos(2*pi*day.of.yr/365)
  
  # 4.3: Solar declination (radians)
  sol.dec <- 0.409*sin((2*pi*day.of.yr/365)-1.39)
  
  # 4.4: Sunset hour angle
  set.hr <- acos(-1.*tan(lat.rad)*tan(sol.dec))
  
  # 4.5: Extra-terrestrial radiation (MJ m-2 day-1)
  cetr <- 0.082 * 24 * 60 / pi
  et.rad <- cetr * rel.dist * (set.hr*sin(lat.rad)*sin(sol.dec) +
                                 cos(lat.rad)*cos(sol.dec)*sin(set.hr))
  
  # 4.6: Daylight hours
  day.hrs <- 24*set.hr/pi
  
  # 4.7: Solar radiation (MJ m-2 day-1)  
  #      MT 2023: This is an input from ChessScape. Units need to be coverted
  #               from W m-2 to MJ m-2 day-1
  # sol.rad <- (a.s+(b.s*sun.hrs/day.hrs))*et.rad
  # sol.rad <- (a.s+(b.s*(sun.hrs/month.dict)/day.hrs))*et.rad
  # 
  sol.rad <- shortwave.rad*3600*24/1000000
  
  
  # 4.8: Clear-sky solar radiation (MJ m-2 day-1) 
  #      MT 2023: Small change so fraction of Ra reaching earth on clear days 
  #      equals as+bs
  #clr.sky.rad <- (0.75+(elev*2E-5))*et.rad
  clr.sky.rad <- (a.s+b.s+(elev*2E-5))*et.rad
  
  # 4.9: Net shortwave radiation (MJ m-2 day-1) 
  sw.rad <- (1-alb.coef)*sol.rad
  
  # 4.10: Modified 27/11/2014. Net longwave radiation (MJ m-2 day-1) 
  # First convert Tmin and Tmax to kelvin
  tmin.K <- tmin.cur + 273.15
  tmax.K <- tmax.cur + 273.15
  lw.rad <- (4.903E-9 * ((tmin.K^4+tmax.K^4)/2) *
               (0.34 - 0.14*(e.a^0.5)) *
               ((1.35*sol.rad/clr.sky.rad)-0.35))
  
  # 4.11: Net radiation (MJ m-2 day-1) 
  net.rad <- sw.rad - lw.rad
  
  ## ******************************************************************
  # Step 5 - Estimate soil heat flux (MJ m-2 day-1)
  
  # 5.1: Soil heat flux. Assumes soil heat cap <- 2.1 MJ/m3/C
  #      and effective soil depth of 2m
  
  soil.heat <- 0.07*tmean.diff
  
  ## ******************************************************************
  # Step 6 - Estimate PET for grassland in mm/day
  
  # 6.1: The Penman-Monteith equation (av. PET over month in mm/day)
  # Currently includes a fudge to set any values less than zero to
  # zero. Need to investigate this!
  # MT 2023 notes on units: 
  #   0.408 (kg MJ-1). Inverse of the latent heat of vaporization (i.e. energy required to evaporate 1kg of water; 2.45 MJ/kg)
  #   slp.sat.vp (kPa degC-1)
  #   net.rad (MJ m-2 day-1)
  #   soil.heat (MJ m-2 day-1)
  #   psy.c (kPa degC-1)
  #   900 (s kg K day-1 kJ-1)
  #   tmean.cur (degC)
  #   u.2 (m s-1)
  #   vap.p.def (kPa)
  #   0.34 (s m-1). This is based on rs/ra = rs*u.2/208 and assuming rs=70 s m-1
  
  pet <- (((0.408*slp.sat.vp*(net.rad-soil.heat)) +
             (psy.c*(900/(tmean.cur+273.16))*u.2*vap.p.def)) /
            (slp.sat.vp + psy.c*(1+0.34*u.2)))
  
  pet[pet<0] <- 0
  
  ## ******************************************************************
  ## Step 7 - Convert to monthly totals i.e. mm/month
  
  pet.mon <- pet * month.dict
  
  ## ******************************************************************
  
  pet.mon
}


# End of Functions 4 ----











