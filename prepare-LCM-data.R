# This function reads the land cover map (LCM) raster files (tiffs). 
# The script was modified from IP 2017.
#
# LCM raster files from 2007 to 2021 have been downloaded from CEH data portal.
#
# The LCM data include both aggregate land cover classes (10 in total) and more 
# detailed (target) land cover classes (21 in total). The classes are summarized 
# for each raster pixel, either as the dominant land cover or as percentages of 
# different land covers within the pixel.Therefore, for each available year,
# the following 4 land cover data files, each covering all of GB and at 1km 
# resolution, are available:
#  a) dominant target (tiff file with one band with the dominant target LC class)
#  b) dominant aggregate (tiff file with one band with the dominant aggregate LC class)
#  c) percentage target (tiff file with 21 bands with the percentage of LC classes)
#  d) percentage aggregate (tiff file with 10 bands with the percentage of aggregate LC classes)
#
# The classes are summarized below. However, note in 2007 there are 23 target 
# classes. The semi-natural grassland aggregate class also includes rough 
# grassland, and the mountain, heath & bog aggregate class also includes
# montane habitats.
#
# AGGREGATE CLASS         	    TARGET CLASS
# 1 Broadleaf woodland  	      1 Deciduous woodland	
# 2 Coniferous woodland	    	  2 Coniferous woodland	
# 3 Arable	                	  3 Arable	            
# 4 Improved grassland 	    	  4 Improved grassland	  
# 5 Semi-natural grassland		  5 Neutral grassland	    
#                               6 Calcareous grassland	
#                               7 Acid grassland      	
#                               8 Fen                  	
# 6 Mountain, heath & bog 		  9 Heather              	
#                               10 Heather grassland   	
#                               11 Bog	                
#                               12 Inland rock         	
# 7 Saltwater	                  13 Saltwater            
# 8 Freshwater             		  14 Freshwater          	
# 9 Coastal                		  15 Supralittoral rock   
#                               16 Supralittoral sediment
#                               17 Littoral  rock      	
#                               18 Littoral sediment 	  
#                               19 Saltmarsh 	          
# 10 Built-up areas & gardens	  20 Urban 	
#                               21 Suburban	
#
# NOTE: All the LCM raster files are saved using the following naming convention:
# gbYYYlcm1km_percentage_target
# gbYYYYlcm1km_percentage_aggregate
# gbYYYYlcm1km_dominant_target
# gbYYYYlcm1km_dominant_aggregate


get.lcm.data <- function(dir.path, st.yr=2010, end.yr=2021, nirams.grid) {
  # ------------------------------------------------------------------------- #
  # L O A D   &   S A V E   L C M   R A S T E R   F I L E S
  # ------------------------------------------------------------------------- #
  # LCM data has been downloaded in raster format for the years 2007, 2015, and
  # 2017-2021. Here the target percentage LCM data for a given year are loaded 
  # and saved in a list. The LCM data are cropped to fit the extent of NIRAMS.

  e.nirams <- extent(nirams.grid)   #extent(0, 485000, 520000, 1235000)
  #dir.path = "C:/Users/mt40375/Documents/NIRAMS/input-files/Raw-data-LCM/"
  pat     = ".tif|.img"
  
  yrs <- seq(st.yr,end.yr, by=1)
  
  # List the folders with the annual lcm data (only for the years later than 
  # or equal to the specified starting year)
  all.lcm.folders <- list.dirs(path=dir.path, recursive=FALSE)
  all.yrs         <- as.numeric(str_extract(basename(all.lcm.folders), "[0-9]{4}"))
  lcm.folders     <- all.lcm.folders[all.yrs>st.yr]
  
  lcm.raster <- list()
  curr.folder <- 0
  
  # loop over years
  for (folder in lcm.folders){
    curr.folder <- curr.folder + 1
    lcm.yr <- str_extract(basename(folder), "[0-9]{4}")
    
    # list all .tif files in this folder and its subfolders
    myfiles=list.files(folder, pattern = pat, ignore.case=TRUE, full.names = FALSE,recursive=TRUE)
    myfiles <- myfiles[ !grepl(".xml",myfiles)]
    myfiles.names <- sub(".*lcm1km_(.*?)\\..*", "\\1", myfiles)
    my.id <- which(myfiles.names == "percentage_target")
    
    # load raster
    my.filepath <- paste(folder,myfiles[my.id],sep="/")
    r <- brick(my.filepath) # or use stack
    
    # crop
    r <- crop(r,e.nirams)
    
    # add to list
    lcm.raster[[curr.folder]] <- r
    names(lcm.raster)[[curr.folder]] <- lcm.yr 
    }
 
  
  # change LCM category names to be consistent for all years
  for (i in names(lcm.raster)){
    lcm.oldnames <- names(lcm.raster[[i]])
    lcm.newnames <- paste0("LC",seq(1,length(lcm.oldnames),by=1))
    names(lcm.raster[[i]]) <- lcm.newnames
    }

return(lcm.raster)

}

