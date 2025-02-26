## This function loads the 2km gridded AGCENSUS raster data for Scotland, and 
## then re-scales and re-samples the rasters to match the extent and resolution
## of the NIRAMS. 
## The AGCENSUS data are needed for estimating nitrogen inputs and outputs.
## Re-written from scratch. MT 2023.
## 
## The AGCENSUS data contain land use and livestock data summarized at 2k grid. 
## The input data have been extracted for Scotland only. Raw AGCENSUS data files  
## from 1969-2015 were available from IP previous work in .csv format. More 
## recent files (up until 2019) can be downloaded from Digimap.edina.ac.uk which 
## JHI has access to. However, these data are in a different format, i.e. land
## use and livestock are available as separate files. AGCENSUS data have been
## downloaded in raster format for the years 2003-2019 and processed/combined 
## into one raster stack for each year.
## Note that the extent of the raster files before 2016 is smaller than from 
## 2016+. Also, note NIRAMS has a larger extent in the x direction than agcensus
##  
## Before 2016: 
##   dimensions : 345, 208, 71760  (nrow, ncol, ncell)
##   resolution : 2000, 2000  (x, y)
##   extent     : 54000, 470000, 530000, 1220000  (xmin, xmax, ymin, ymax)
## After 2016:
##   dimensions : 360, 240, 86400  (nrow, ncol, ncell)
##   resolution : 2000, 2000  (x, y)
##   extent     : 0, 480000, 520000, 1240000  (xmin, xmax, ymin, ymax)
##
## NIRAMS:
##   dimensions : 715, 485, 346775  (nrow, ncol, ncell)
##   resolution : 1000, 1000  (x, y)
##   extent     : 0, 485000, 520000, 1235000  (xmin, xmax, ymin, ymax)


get.agcensus.data <- function(dir.path, st.yr=2010, end.yr=2019, nirams.grid) {
  # INPUT:
  # dir.path          path to folder with agcensus raster files
  # st.yr             the first year for which agcensus data should be loaded
  # end.yr            the last year for which agcensus data should be loaded
  # nirams.grid       SpatialPixels object with the NIRAMS grid cell center 
  #                   coordinates (or for whatever gridded area that the emep 
  #                   data need to be extracted for)
  # OUTPUT:
  # agcensus_raster_resampled
  
  
  # ------------------------------------------------------------------------- #
  # L O A D   &   S A V E   C E N S U S   R A S T E R   F I L E S
  # ------------------------------------------------------------------------- #
  # AGCensus data has been downloaded in raster format from 2003-2019. 
  # Here agcensus data are loaded for a given year and saved as a raster stack.
  # Each raster stack is saved in a list. 
  
  pat     = ".tif"
  
  # List the folders with the annual agcensus data 
  yrs <- seq(st.yr,end.yr, by=1)
  nyr <- length(yrs)
  all.folders <- list.dirs(path=dir.path, recursive=FALSE)
  all.yrs     <- as.numeric(str_extract(basename(all.folders), "[0-9]{4}"))
  my.folders  <- all.folders[all.yrs %in% yrs]
  
  agcensus_raster <- list()
  curr.folder     <- 0
  
  # loop over years
  for (folder in my.folders){
    curr.folder <- curr.folder + 1
    agcensus.yr <- str_extract(basename(folder), "[0-9]{4}")
    rs <- stack()
    
    # first list all files matching the pattern (pat) in the directory and its subfolders
    myfiles=list.files(folder, pattern = pat, ignore.case=TRUE, full.names = FALSE,recursive=TRUE)
    # ignore farmers-and-workers and hives files
    myfiles <- myfiles[ !grepl("farmers-and-workers/",myfiles)]
    myfiles <- myfiles[ !grepl("hives/",myfiles)]
  
    # loop over files
    for (file in myfiles){
    
      my.filepath <- paste(folder,file,sep="/")
      my.filename <- sub(".*-2km-(.*?)\\..*", "\\1", file)
      r <- raster(my.filepath)
      names(r) <- my.filename

      rs <- stack(rs,r)
    
    }
  
    agcensus_raster[[curr.folder]] <- rs  # add to raster stack to list
    names(agcensus_raster)[[curr.folder]] <- agcensus.yr 
  
  }
  
  #check_r <- as.data.frame(rs, xy=TRUE)
  # ------------------------------------------------------------------------- #
  #  M A T C H   E X T E N T   T O   N I R A M S   &   R E S A M P L E   
  # ------------------------------------------------------------------------- #
  # Note that NIRAMS has a larger extent in the x direction than agcensus.
  # A new extent (e.new) is therefore defined, which covers the extent of both 
  # agcensus and nirams. The agcensus rasters are first extended to the new 
  # extent, then resampled to 1km resolution, and finally the resampled agcensus
  # rasters are cropped to match the nirams extent
  e.new     <- extent(0, 486000, 520000, 1240000)  # 
  e.nirams  <- extent(nirams.grid)# extent(0, 485000, 520000, 1235000)
  n.rasters <- length(agcensus_raster)
  ag.crs    <- proj4string(agcensus_raster[[n.rasters]]) 
  raster1km <- raster(resolution=c(1000,1000), crs=ag.crs, ext=e.new)
  
  agcensus_raster_resampled <- list()
  
  year.index      <- 0
  curr.year       <- as.character(st.yr) # "2003"
  for (curr.year in names(agcensus_raster)){
    year.index <- year.index + 1
    # extend raster
    raster.ext <- extend(agcensus_raster[[curr.year]],e.new)
    print("extended")
    # resample
    raster.rs  <- resample(raster.ext, raster1km, method='ngb')
    print("resampled")
    raster.rs  <- raster.rs/4  # divide by 4 as each 2km cell is now divided into four 1km cells
    print("corrected")
    # crop to match nirams
    agcensus_raster_resampled[[year.index]] <- crop(raster.rs,e.nirams)
    print("cropped")
    names(agcensus_raster_resampled)[[year.index]] <- curr.year 
    rm(raster.ext,raster.rs)
    print(paste(curr.year,"completed"))
  }
  return(agcensus_raster_resampled)
}


# # ---------------------------------------------------------------------------
# # G E T   A L L   C E N S U S   C A T E G O R I E S 
# # ---------------------------------------------------------------------------
# # This only had to be done once!
# # Each agcensus land use and livestock category is associated with a NIRAMS 
# # category that determines the amount of N input (organic or inorganic) and
# # N uptake.
# # The number and the names of the census categories have changed over time. 
# # Here, all the unique agcensus category names are extracted and saved in a
# # .csv file. These have subsequently manually been assigned NIRAMS categories
# 
# # list all filenames in the directory and its subfolders:
# myfiles <- list.files(dir.path, pattern = pat, ignore.case=TRUE, full.names = FALSE,recursive=TRUE)
# myfiles <- myfiles[ !grepl("farmers-and-workers/",myfiles)]
# myfiles <- myfiles[ !grepl("hives/",myfiles)]
# myfiles <- basename(myfiles)
# census_yr   = str_extract(myfiles, "[0-9]{4}")
# census_code = str_extract(myfiles, "c[0-9]{1,3}")
# census_cat  = sub(".*-2km-(.*?)\\..*", "\\1", myfiles)
# tmp.df      = data.frame(census_yr=census_yr,census_code=census_code, census_cat=census_cat)
# 
# unique.names <- unique(tmp.df$census_cat)
# write.table(unique.names, file = "EDINA_NAMES_upd.txt")





