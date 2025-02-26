## This function loads the gridded CRAFTY raster data on future land use 
## for the UK, and then re-projects, re-scales and re-samples the rasters to 
## match the extent and resolution of the NIRAMS. 
## The CRAFTY data are needed for estimating future nitrogen inputs and outputs.
## Re-written from scratch. MT 2024.
## 
## CRAFTY-GB simulates land use change between 2020 and 2080 as the result of 
## climatic and socio-economic impacts derived from the Representative 
## Concentration Pathways (RCPs) and Shared Socioeconomic Pathways (SSPs) 
## global scenarios. These have been tailored to the British context through a 
## combination of computational modelling and stakeholder engagement. 
## For more details about CRAFTY-GB: https://doi.org/10.1029/2022EF002905 
##
## The CRAFTY files loaded here have been pre-processed by Zisis Gagkas. The 
## data have already been reprojected and resampled to BNG at 1k resolution.
## 
## Data have been downloaded for 2020, 2040, 2060 and 2080.
## The CRAFTY data contain spatial predictions of 17 different land use agents. 
## These agents have been matched to the land use classes in NIRAMS.
##
## CRAFTY_LU_agent                AGWENT	   LU_INDEX   LANDUSE (P)         LANDUSE (NIRAMS)
## Unmanaged                      Lazy FR	         -1	  Wildscape           Montane, bog, heathland         
## Agroforestry                   AF	              0   Woodland            Woodland and forest
## Bioenergy                      Bioenergy	        1	  Arable intensive    Arable (ware and seed potatoes)
## Extensive arable               EA	              2	  Arable extensive    Arable (cereals,fruit,veg, other)
## Extensive pastoral             EP	              3	  Rough grazing       Rough grazing
## Intensive arable (fodder)      IAfodder          4	  Arable intensive    Arable (ware and seed potatoes)
## Intensive arable (food)        IAfood	          5	  Arable intensive    Arable (ware and seed potatoes)
## Intensive pastoral             IP              	6	  Grassland improved  Improved grassland
## Multifunctional mixed woodland MW              	7	  Woodland            Woodland and forest
## Native woodland (conservation) NWCons          	8 	Woodland            Woodland and forest
## Productive native broadleaf    PNB	              9	  Woodland            Woodland and forest
## Productive native conifer      PNC             	10	Forestry            Woodland and forest
## Productive nonnative broadleaf PNNB	            11	Woodland            Woodland and forest
## Productive nonnative conifer   PNNC  	          12	Forestry            Woodland and forest
## Sustainable arable             SusAr 	          13	Arable extensive    Arable (cereals,fruit,veg, other)
## Very extensive pastoral        VEP         	    14	Wildscape           Fen, marsh and swamp, bog, heathland, montane
## Urban                          Urban	            15	Other               Bare and built-up
##
## INPUT
## dir.path       path to folder containing the CRAFTY raster files
## nirams.grid    the NIRAMS grid
##
## OUTPUT
## 

get.crafty.data <- function(dir.path=paste0("C:/Users/mt40375/Documents/NIRAMS/input-files/","Raw-data-future-landuse/Reproject/RCP6_0-SSP3/"), 
                                 nirams.grid) {

  pat       = ".tif"
  myfiles   <- list.files(dir.path, pattern = pat, ignore.case=TRUE, full.names = FALSE,recursive=TRUE)
  e.nirams  <- extent(nirams.grid)# extent(0, 485000, 520000, 1235000)
  #scenario <- basename(dir.path)  # name of scenario (i.e. name of folder)
  
  crafty_raster <- list()
  curr.file     <- 0
  
  # loop over files
  for (file in myfiles){
    curr.file <- curr.file + 1
    crafty.yr <- regmatches(file, regexpr("(\\d{4})(?=_([a-zA-Z]+))",file, perl=T)) # get year (assumes year is followed by _ and then text)
   
    my.filepath <- paste0(dir.path,file)
    r <- raster(my.filepath)
    names(r) <- crafty.yr
    
    #my.filename <- paste(scenario,crafty.yr,sep="_") 
    #names(r) <- my.filename
    
    # crop raster and add to list
    crafty_raster[[curr.file]] <- crop(r,e.nirams)  # 
    names(crafty_raster)[[curr.file]] <- crafty.yr 
  
    
  }
  # # Check some land use fractions
  # crafty_df   <- as.data.frame(crafty_raster$"2040", xy=TRUE)
  # names(crafty_df)[3] <- "crafty"
  # nlu_tot = length(crafty_df$crafty)-sum(is.na(crafty_df$crafty))
  # n_rough = sum(crafty_df$crafty==3, na.rm=TRUE)
  # n_impgr = sum(crafty_df$crafty==6, na.rm=TRUE)
  # n_arabl = sum(crafty_df$crafty==1, na.rm=TRUE) + sum(crafty_df$crafty==2, na.rm=TRUE) +
  #           sum(crafty_df$crafty==4, na.rm=TRUE) + sum(crafty_df$crafty==5, na.rm=TRUE) + 
  #           sum(crafty_df$crafty==13, na.rm=TRUE)
  # n_wood  = sum(crafty_df$crafty==0, na.rm=TRUE) + sum(crafty_df$crafty==7, na.rm=TRUE) +
  #           sum(crafty_df$crafty==8, na.rm=TRUE) + sum(crafty_df$crafty==9, na.rm=TRUE) + 
  #           sum(crafty_df$crafty==10, na.rm=TRUE) + sum(crafty_df$crafty==11, na.rm=TRUE) +
  #           sum(crafty_df$crafty==12, na.rm=TRUE)
  # n_built = sum(crafty_df$crafty==15, na.rm=TRUE)
  # n_wild  = sum(crafty_df$crafty==-1, na.rm=TRUE) + sum(crafty_df$crafty==14, na.rm=TRUE)
  # 
  # n_built/nlu_tot
  # n_wood/nlu_tot
  # n_arabl/nlu_tot
  # n_rough/nlu_tot
  # n_impgr/nlu_tot
  # n_wild/nlu_tot
  
  # # ====================================== #
  # # Use below if resampling is required. 
  # # ====================================== #
  # # Resampling done based on lcm grid
  # path2lcm <- paste0("C:/Users/mt40375/Documents/NIRAMS/input-files/Raw-data-future-landuse/Reproject/")
  # lcm_1km  <- raster(paste0(path2lcm,"lcm2015_gb_1km_dominant_aggregate_class.tif"))
  # 
  # dir.path=paste0("C:/Users/mt40375/Documents/NIRAMS/input-files/","Raw-data-future-landuse/Reproject/RCP2_6-SSP1/")
  # myfiles   <- list.files(dir.path, pattern = pat, ignore.case=TRUE, full.names = FALSE,recursive=TRUE)
  # 
  # crafty_raster <- list()
  # curr.file     <- 0
  # 
  # for (file in myfiles){
  #   curr.file <- curr.file + 1
  #   crafty.yr <- regmatches(file, regexpr("(\\d{4})(?=.([a-zA-Z]+))",file, perl=T)) # get year (assumes year is followed by . and then text)
  #   
  #   my.filepath <- paste0(dir.path,file)
  #   r <- raster(my.filepath)
  #   names(r) <- crafty.yr
  #   
  #   r.bng <-  projectRaster(r,crs = crs(lcm_1km))      # set projection
  #   r.rs  <-  resample(r.bng, lcm_1km, method="ngb")   # resample
  # 
  #   # crop raster and add to list
  #   crafty_raster[[curr.file]] <- crop(r.rs,e.nirams)  # 
  #   names(crafty_raster)[[curr.file]] <- crafty.yr 
  #   r.c <- crop(r.rs,e.nirams)
  #   writeRaster(r.c, paste0("crafty_",crafty.yr,"_bng.tif"))
  # }

  return(crafty_raster)
  
}



  # ## using a 1km lcm grid as a template for reprojection
  # lcm_1km <- raster("lcm2015_gb_1km_dominant_aggregate_class.tif")
  # 
  # ## using one of the crafty datasets as an example (for UK; projection is WGS84) 
  # r20 <- raster(paste0(dir.path,"RCP6_0-SSP3-0-99-UK-Cell-2020.csv_LandUseIndex.tif"))
  # r40 <- raster(paste0(dir.path,"RCP6_0-SSP3-0-99-UK-Cell-2040.csv_LandUseIndex.tif"))
  # r60 <- raster(paste0(dir.path,"RCP6_0-SSP3-0-99-UK-Cell-2060.csv_LandUseIndex.tif"))
  # r80 <- raster(paste0(dir.path,"RCP6_0-SSP3-0-99-UK-Cell-2080.csv_LandUseIndex.tif"))
  # 
  # r40_2 <- raster(paste0(dir.path,"RCP6_0_SSP3_LandUseIndex_2040_OS.tif"))
  # 
  # extent(r40)
  # extent(lcm_1km)
  # 
  # r40BNG <- projectRaster(r40,crs = crs(nirams.grid))
  # r40_OS <- resample(r40BNG, lcm_1km, method="ngb") 
  # 
  # extent(r40BNG)
  # extent(lcm_1km)
  # extent(nirams.grid)
  # 
  # ## resample crafty data onto lcm 1km BNG
  # r40_OS <- resample(r20, lcm_1km, method="ngb") 
  # r40_OS
  # 
  # ## output
  # writeRaster(r80_OS, "crafty_80_OS.tif")






