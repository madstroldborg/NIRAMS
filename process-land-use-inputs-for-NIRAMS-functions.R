# This script contains the functions to process the prepared AGCensus and LCM
# data to calculate the spatial N inputs and the PET factor required for NIRAMS.
# The script contains the following functions:
# (1) Assign NIRAMS land use codes to each grid cell in the NIRAMS grid based on
#     the distribution of land uses within each cell as determined from both 
#     AGCensus (for the different arable categories) and LCM. 
# (2) Calculate the total organic N input to each grid cell based on AGCensus 
#     livestock numbers.
# (3) Distribute annual organic N inputs to broad land use classes
# (4) Calculate the inorganic N application and uptake based on the distribution
#     of NIRAMS land uses within each cell.
# (3) calculate the weighted AET/PET factor for each grid cell based on the
#     distribution of NIRAMS land uses within each cell.
#
# MT 2023

# **************************************************************************** #
# FUNCTION 1: Assign NIRAMS land use categories based on AGCensus and LCM ----
# **************************************************************************** #

assign_nirams_landuse_categories <- function(my.agcensus, my.lcm, Categories_AGCENSUS_NIRAMS){
  # AGCensus contain a lot of detailed land use classes (as well as livestock
  # categories). Each of these land uses are here assigned a corresponding 
  # NIRAMS land use class. The Categories_AGCENSUS_NIRAMS file shows which 
  # NIRAMS classes are assigned to which AGCensus classes. The AGCensus is 
  # mainly used to give the proportion of arable subclasses (e.g. SB, WB etc) 
  # whereas LCM is used to give the distribution of the broader land use 
  # categories (e.g., Arable, Improved grassland etc.)
  # The NIRAMS land use codes are as follows:
  # CODE    DETAILED_CLASS	        BROAD_CLASS
  #   1     Spring barley	          Arable
  #   2     Winter barley	          Arable
  #   3     Spring wheat	          Arable
  #   4     Winter wheat	          Arable
  #   5     Spring oil seed rape	  Arable
  #   6     Winter oil seed rape	  Arable
  #   7     Spring oats	            Arable
  #   8   	Winter oats	            Arable
  #   9     Seed potatoes 	        Arable
  #   10    Ware potatoes	          Arable
  #   11  	Fruit	                  Arable
  #   12	  Vegetables	            Arable
  #   13	  Other arable	          Arable
  #   14	  Set aside	              Arable
  #   15	  Rough grassland	        Rough grassland
  #   16  	Improved grassland  	  Improved grassland
  #   17	  Woodland and forest	    Woodland and forest
  #   18  	Short rotation coppice	Woodland and forest
  #   19	  Bare and built up	      Bare and built up
  #   20  	Water	                  Water
  #   21  	Grazed woodland	        Woodland and forest
  #   22	  Fen, marsh and swamp  	Fen, marsh and swamp
  #   23  	Heathland	              Heathland
  #   24	  Bog	                    Bog
  #   25  	Montane	                Montane
  #   26	  Saltmarsh	              Saltmarsh
  #   27	  Other	                  Other
  #   28  	None	                  LCM2007 arable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # INPUT:
  # my.agcensus  Raster brick containing the agcensus data for a given year. 
  #              Extent and resolution must match that of NIRAMS
  # my.lcm       Raster stack containing the lcm data for a given year.
  #              Extent and resolution must match that of NIRAMS
  # OUTPUT
  # landusecats  Data.frame with the areas of each NIRAMS land use category 
  #              within each cell of the NIRAMS grid
  
  # COMBINE AGCENSUS AND LCM DATA IN ONE DATA.FRAME
  agcensus.df <- as.data.frame(my.agcensus, xy=TRUE)
  lcm.df      <- as.data.frame(my.lcm, xy=TRUE)
  matched     <- merge(agcensus.df, lcm.df, by=c("x","y"),all=TRUE)
  
  landusecats <- matched[,c(1:2)]  # data.frame for storing final values
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - -
  # ASSIGN NIRAMS LAND USE CODES 
  # - - - - - - - - - - - - - - - - - - - - - - - - - -
  # LUC 1-14: ARABLE
  # Land in ha is (specific class / total agcensus arable) * lcm arable

  # Extract census areas (ha) for all of the arable categories
  ARABLE       <- matched[,names(matched) %in% Categories_AGCENSUS_NIRAMS$AGCENSUS_cat[Categories_AGCENSUS_NIRAMS$Summary_Class == "Arable"]]
  # Calculate total arable area based on census data
  ARABLE$TOTAL <- rowSums(matched[,names(matched) %in% Categories_AGCENSUS_NIRAMS$AGCENSUS_cat[Categories_AGCENSUS_NIRAMS$Summary_Class == "Arable"]], na.rm = TRUE)

  # get arable areas as reported in LCM
  ARABLE$LCM  <- matched[,names(matched) == "LC3"]

  for (i in c(1:14)){
    # Extract agcensus area values associated NIRAMS land use code i
    agric <- matched[,names(matched) %in% Categories_AGCENSUS_NIRAMS$AGCENSUS_cat[Categories_AGCENSUS_NIRAMS$LU_CODE == i]]
    
    # normalize extracted agcensus land areas by total arable area from lcm
    if (is.null(dim(agric))){
      agric.landuse <- (agric / ARABLE$TOTAL) * ARABLE$LCM
      } else {
        agric.landuse <- (rowSums( agric , na.rm = TRUE ) / ARABLE$TOTAL) * ARABLE$LCM
      }
    # add normalized values to the landusecats df
    landusecats[,(i + 2)]     <- agric.landuse
    names(landusecats)[i + 2] <- paste("LUC",i,sep = "_")
    }
  #set NA to zero
  #landusecats[is.na(landusecats)] <- 0
  
  # If no agricultural land use according to agcensus, but according to lcm, 
  # then set to "other - 13"
  landusecats$LUC_13[rowSums(landusecats[,(c(3:16))], na.rm = TRUE) == 0] <- ARABLE$LCM[rowSums(landusecats[,(c(3:16))],na.rm=TRUE) == 0]
  
  #landusecats[is.na(landusecats)] <- 0
  rm(ARABLE)
  
  
  # LUC 15: ROUGH GRAZING (SEMINATURAL)
  # Here rough grazing from lcm is used as these appear flawed in agcensus
  # Note: IP 2017 used the LCM2007 which has 23 target classes, including an 
  # additional seminatural class (called rough grassland). This is not included
  # in LCM2015+
  LUC_15 <- rowSums(matched[,names(matched) %in% c("LC5","LC6","LC7")], na.rm = TRUE)
  landusecats$LUC_15 <- LUC_15
  rm(LUC_15)
  
  # LUC 16: IMPROVED GRASSLAND
  # Improved grassland from lcm 
  landusecats$LUC_16 <- matched[,names(matched) %in% c("LC4")]
  
  # LUC 17: WOORDLAND & FOREST
  # Values from lcm
  LUC_17 <- rowSums(matched[,names(matched) %in% c("LC1","LC2")], na.rm = TRUE)
  landusecats$LUC_17 <-  LUC_17
  rm(LUC_17)
  
  # LUC 18: SHORT ROTATION COPPICE
  # We do not have data on this one
  landusecats$LUC_18 <- 0
  
  # LUC 19: BARE & BUILT UP 
  # Values from lcm. Note that IP 2017 included inland rock as "bare & 
  # built up". However, here it is assumed inland rock belongs to "montane" (see 
  # later). This is done because LCM2015+ does no longer include a separate 
  # montane target class, unlike the LCM2007 (which IP used). 
  LUC_19 <- rowSums(matched[,names(matched) %in% c("LC20", "LC21")], na.rm = TRUE)
  #LUC_19 <- rowSums(matched[,names(matched) %in% c("LC14","LC22", "23")], na.rm = TRUE)
  landusecats$LUC_19 <-  LUC_19
  rm(LUC_19)
  
  # LUC 20: WATER 
  # Values from lcm
  LUC_20 <- rowSums(matched[,names(matched) %in% c("LC13","LC14")], na.rm = TRUE)
  landusecats$LUC_20 <-  LUC_20
  
  # LUC 21: GRAZED WOODLAND
  # Values from agcensus 
  # LUC_21 <- 0 * LUC_20  # ?????
  LUC_21 <- matched[,names(matched) %in% Categories_AGCENSUS_NIRAMS$AGCENSUS_cat[Categories_AGCENSUS_NIRAMS$LU_CODE == 21]]
  LUC_21[is.na(LUC_21)] <- 0 #set NA to zero
  landusecats$LUC_21 <- LUC_21
  rm(LUC_21)
  rm(LUC_20)
  
  # LUC 22: FEN MARSH & SWAMP 
  # Values from lcm
  landusecats$LUC_22 <- matched[,names(matched) %in% c("LC8")]
  
  # LUC 23: HEATHLAND
  # Values from lcm
  LUC_23 <- rowSums(matched[,names(matched) %in% c("LC9","LC10")], na.rm = TRUE)
  landusecats$LUC_23 <- LUC_23 
  rm(LUC_23)
  
  # LUC 24: BOG 
  # Values from lcm
  landusecats$LUC_24 <- matched[,names(matched) %in% c("LC11")]
  
  # LUC 25: MONTANE 
  # Values from lcm. However, note the montane target class does not exist in 
  # LCM after 2007.This is here assumed to be inland rock
  landusecats$LUC_25 <- matched[,names(matched) %in% c("LC12")]
  
  # LUC 26: SALTMARSH
  # Values from lcm. IP assumed saltmarsh to be based on supra-littoral rock and 
  # sediment (LC17 and 18 in LCM2007). However, there is a separate target class 
  # in LCM called saltmarsh and overall it seems like littoral sediments would 
  # be a better representation of saltmarsh. Here, it is therefore assumed that
  # saltmarsh is constituted by littoral sediment (LC18) and saltmarsh (LC19) 
  LUC_26 <- rowSums(matched[,names(matched) %in% c("LC18","LC19")], na.rm = TRUE)
  landusecats$LUC_26 <- LUC_26
  rm(LUC_26)
  
  # MT2023 NOTE: Other (LUC_27) and None (LUC_28) not included in landusecats. 
  # Also, LCM classes 15 (supra-littoral rock), 16 (supra-littoral sediment) and 
  # 17 (littoral rock) are not included. Should these just be included as saltmarsh?
  
  # Set NA to zero
  #landusecats[is.na(landusecats)] <- 0
  # Limit all category values to 100%
  landusecats[,-c(1:2)][landusecats[,-c(1:2)] > 100] <- 100
  
  # Calculate sum of land use categories
  landusecats$sums <- rowSums(  landusecats[,-c(1:2)], na.rm = TRUE)
  # Set all zeros to NA (i dont want to divide by zero) -
  # ONLY RELEVANT IF THE NEXT STEP IS INCLUDED
  landusecats$sums[landusecats$sums == 0] <- NA
  summary(landusecats$sums)
  
  ## For land use cover not 100% redistribute in relation to land use of all
  ## types in land use class.
  ##  MT 20223 comment: Why is this done? Is it not possible for some cells not 
  ##                    to 100% land cover, such as cells close to the coast or 
  ##                    just by the fact that not all LCM land use categories 
  ##                    are included in the data frame?
  ## STEP BELOW REMOVED!!!
  # landusecats[,-c(1,2,29)] <- landusecats[,-c(1,2,29)] * 100 / landusecats[,29] 
  # landusecats$sums <- rowSums(  landusecats[,-c(1,2,29)], na.rm = TRUE)
  
  # # EXTREMELY SLOW TO RUN
  # for (i in c(3:28)){
  #   for (y in c(1:length(landusecats[,i]))){
  #     if (is.na(landusecats[y,i]) & !is.na(landusecats$sums[y]))
  #       landusecats[y,i] <- 0
  #   }}
  for (i in c(3:28)){
    
    naid <- which(is.na(landusecats[,i]) & !is.na(landusecats$sums))
    landusecats[naid,i] <- 0
    }
  
  landusecats$sums[landusecats$sums == 0] <- NA
  
  return(landusecats)
}
# ---- End of Function 1 ----


# **************************************************************************** #
# FUNCTION 2: Calculate N output (kg/ha) from all livestock to grid cells ----
# **************************************************************************** #

Livestock_N <- function(my.agcensus,Categories_AGCENSUS_NIRAMS,
                        JAC_Livestock_N_Per_Animal_InputMatrix){
  # Calculates the total kg org N/ha from all livestock for each grid cell

  resampled.agcensus <- as.data.frame(my.agcensus, xy=TRUE)
  
  # -----------------------------------------
  # PIGS
  # -----------------------------------------
  # Extract number of pigs from census data
  Pigs <- resampled.agcensus[,names(resampled.agcensus) %in% 
                               Categories_AGCENSUS_NIRAMS$AGCENSUS_cat[Categories_AGCENSUS_NIRAMS$AGCENSUS_group == "pigs" & 
                                                                       Categories_AGCENSUS_NIRAMS$AGCENSUS_subclass == "individual"]]

  # Get annual N output (kg/stock unit/yr) for each individual pig class
  pigs_QID <- Categories_AGCENSUS_NIRAMS$QID[match(names(Pigs), 
                   Categories_AGCENSUS_NIRAMS$AGCENSUS_cat)]                                    
  N_output.pigs <- JAC_Livestock_N_Per_Animal_InputMatrix$Annual.N.output..kg.stock.unit.year.[match(pigs_QID,JAC_Livestock_N_Per_Animal_InputMatrix$QID)] * 
                   JAC_Livestock_N_Per_Animal_InputMatrix$Default[match(pigs_QID,JAC_Livestock_N_Per_Animal_InputMatrix$QID)]
  
  # Determine annual N output (kg) from pigs to each grid cell 
  # Pigs.N  <- as.data.frame(t(t(Pigs)*(N_output.pigs))) 
  # NPIG    <- rowSums(Pigs.N, na.rm = TRUE)
  Pigs[is.na(Pigs)] = 0 # set na to 0
  NPIG <- as.matrix(Pigs) %*% as.vector(N_output.pigs)
  
  # ## To do calculations in raster format
  # raster.id <- names(my.agcensus) %in% Categories_AGCENSUS_NIRAMS$AGCENSUS_cat[Categories_AGCENSUS_NIRAMS$AGCENSUS_group == "pigs" & 
  #                                                                               Categories_AGCENSUS_NIRAMS$AGCENSUS_subclass == "individual"]
  # Pigs.raster <- my.agcensus[[names(my.agcensus)[raster.id]]]
  # pigs_QID <- Categories_AGCENSUS_NIRAMS$QID[match(names(Pigs.raster), 
  #                                                  Categories_AGCENSUS_NIRAMS$AGCENSUS_cat)]                                    
  # N_output.pigs <- JAC_Livestock_N_Per_Animal_InputMatrix$Annual.N.output..kg.stock.unit.year.[match(pigs_QID,JAC_Livestock_N_Per_Animal_InputMatrix$QID)] * 
  #                  JAC_Livestock_N_Per_Animal_InputMatrix$Default[match(pigs_QID,JAC_Livestock_N_Per_Animal_InputMatrix$QID)]
  # 
  # fun1 <- function(x) { x[is.na(x)] <- 0; return(x)} 
  # rs1  <- calc(Pigs.raster, fun1)
  # rs2  <- calc(rs1, fun=function(x){x * N_output.pigs})
  # NPIG <- calc(rs2, sum)
  
  # - - - - - - - - - - - - - - - - - - - - - 
  # HORSES
  # - - - - - - - - - - - - - - - - - - - - - 
  # Extract number of horses from census data
  Horses <- resampled.agcensus[,names(resampled.agcensus) %in% 
                                 Categories_AGCENSUS_NIRAMS$AGCENSUS_cat[Categories_AGCENSUS_NIRAMS$AGCENSUS_group == "horses" & 
                                                                         Categories_AGCENSUS_NIRAMS$AGCENSUS_subclass == "individual"]]
  # Get annual N output (kg/stock unit/yr) for each individual horse class
  horses_QID <- Categories_AGCENSUS_NIRAMS$QID[match(names(Horses), 
                                                     Categories_AGCENSUS_NIRAMS$AGCENSUS_cat)]                                     
  N_output.horses <- JAC_Livestock_N_Per_Animal_InputMatrix$Annual.N.output..kg.stock.unit.year.[match(horses_QID,JAC_Livestock_N_Per_Animal_InputMatrix$QID)]* 
                     JAC_Livestock_N_Per_Animal_InputMatrix$Default[match(horses_QID,JAC_Livestock_N_Per_Animal_InputMatrix$QID)]
  
  # Determine annual N output (kg) from horses to each grid cell  
  #Horses.N  <- as.data.frame(t(t(Horses)*(N_output.horses))) 
  #NHORSE <- rowSums(Horses.N, na.rm = TRUE)
  Horses[is.na(Horses)] = 0 # set na to 0
  NHORSE <- as.matrix(Horses) %*% as.vector(N_output.horses)
  
  # - - - - - - - - - - - - - - - - - - - - - 
  # SHEEP
  # - - - - - - - - - - - - - - - - - - - - - 
  # Extract number of sheep from census data
  Sheep <- resampled.agcensus[,names(resampled.agcensus) %in% 
                                Categories_AGCENSUS_NIRAMS$AGCENSUS_cat[Categories_AGCENSUS_NIRAMS$AGCENSUS_group == "sheep" & 
                                                                        Categories_AGCENSUS_NIRAMS$AGCENSUS_subclass == "individual"]]
  # Get annual N output (kg/stock unit/yr) for each individual sheep class
  sheep_QID <- Categories_AGCENSUS_NIRAMS$QID[match(names(Sheep), 
                                                     Categories_AGCENSUS_NIRAMS$AGCENSUS_cat)]                                   
  N_output.sheep <- JAC_Livestock_N_Per_Animal_InputMatrix$Annual.N.output..kg.stock.unit.year.[match(sheep_QID,JAC_Livestock_N_Per_Animal_InputMatrix$QID)]* 
                    JAC_Livestock_N_Per_Animal_InputMatrix$Default[match(sheep_QID,JAC_Livestock_N_Per_Animal_InputMatrix$QID)]
  
  # Determine annual N output (kg) from sheep to each grid cell 
  Sheep[is.na(Sheep)] = 0 # set na to 0
  NSHEEP <- as.matrix(Sheep) %*% as.vector(N_output.sheep)
  #Sheep.N  <- as.data.frame(t(t(Sheep)*(N_output.sheep))) 
  #NSHEEP <- rowSums(Sheep.N, na.rm = TRUE)
  
  # - - - - - - - - - - - - - - - - - - - - - 
  # POULTRY
  # - - - - - - - - - - - - - - - - - - - - - 
  # Extract number of poultry from census data
  Poultry <- resampled.agcensus[,names(resampled.agcensus) %in% 
                                  Categories_AGCENSUS_NIRAMS$AGCENSUS_cat[Categories_AGCENSUS_NIRAMS$AGCENSUS_group == "poultry" & 
                                                                          Categories_AGCENSUS_NIRAMS$AGCENSUS_subclass == "individual"]]
  # Get annual N output (kg/stock unit/yr) for each individual poultry class
  poultry_QID <- Categories_AGCENSUS_NIRAMS$QID[match(names(Poultry), 
                                                    Categories_AGCENSUS_NIRAMS$AGCENSUS_cat)]                                   
  N_output.poultry <- JAC_Livestock_N_Per_Animal_InputMatrix$Annual.N.output..kg.stock.unit.year.[match(poultry_QID,JAC_Livestock_N_Per_Animal_InputMatrix$QID)]* 
                      JAC_Livestock_N_Per_Animal_InputMatrix$Default[match(poultry_QID,JAC_Livestock_N_Per_Animal_InputMatrix$QID)]
  
  # Determine annual N output (kg) from poultry to each grid cell 
  Poultry[is.na(Poultry)] = 0 # set na to 0
  NPOULTRY <- as.matrix(Poultry) %*% as.vector(N_output.poultry)
  #Poultry.N  <- as.data.frame(t(t(Poultry)*(N_output.poultry))) 
  #NPOULTRY   <- rowSums(Poultry.N, na.rm = TRUE)
  
  # - - - - - - - - - - - - - - - - - - - - - 
  # Goat - include totals only
  # - - - - - - - - - - - - - - - - - - - - - 
  # Extract number of goats from census data
  Goats <- resampled.agcensus[,names(resampled.agcensus) %in% 
                                Categories_AGCENSUS_NIRAMS$AGCENSUS_cat[Categories_AGCENSUS_NIRAMS$AGCENSUS_group == "goats" & 
                                                                        Categories_AGCENSUS_NIRAMS$AGCENSUS_subclass == "total"]]
  # For total goats, there is only one match and the extracted values 
  if(is.null(dim(Goats))) { 
    goatNames <- Categories_AGCENSUS_NIRAMS$AGCENSUS_cat[Categories_AGCENSUS_NIRAMS$AGCENSUS_group == "goats" & 
                                                           Categories_AGCENSUS_NIRAMS$AGCENSUS_subclass == "total"];
    Goats <- as.data.frame(Goats);
    names(Goats) <- goatNames[1];
    
  } 
  # Get annual N output (kg/stock unit/yr) for each individual goat class
  goat_QID <- Categories_AGCENSUS_NIRAMS$QID[match(names(Goats), 
                                                      Categories_AGCENSUS_NIRAMS$AGCENSUS_cat)]                                   
  N_output.goats <- JAC_Livestock_N_Per_Animal_InputMatrix$Annual.N.output..kg.stock.unit.year.[match(goat_QID,JAC_Livestock_N_Per_Animal_InputMatrix$QID)]* 
                    JAC_Livestock_N_Per_Animal_InputMatrix$Default[match(goat_QID,JAC_Livestock_N_Per_Animal_InputMatrix$QID)]
  
  # Determine annual N output (kg) from goats for each grid cell 
  Goats[is.na(Goats)] = 0 # set na to 0
  NGOATS <- as.matrix(Goats) %*% as.vector(N_output.goats)
  #Goats.N  <- as.data.frame(t(t(Goats)*(N_output.goats))) 
  #NGOATS <- rowSums(Goats.N, na.rm = TRUE)
  
  # - - - - - - - - - - - - - - - - - - - - - 
  # CATTLE
  # - - - - - - - - - - - - - - - - - - - - - 
  # Extract number of cattle from census data
  Cattle <- resampled.agcensus[,names(resampled.agcensus) %in% 
                                 Categories_AGCENSUS_NIRAMS$AGCENSUS_cat[Categories_AGCENSUS_NIRAMS$AGCENSUS_group == "cattle" & 
                                                                         Categories_AGCENSUS_NIRAMS$AGCENSUS_subclass == "individual"]]
  # Get annual N output (kg/stock unit/yr) for each individual sheep class
  cattle_QID <- Categories_AGCENSUS_NIRAMS$QID[match(names(Cattle), 
                                                      Categories_AGCENSUS_NIRAMS$AGCENSUS_cat)]                                   
  N_output.cattle <- JAC_Livestock_N_Per_Animal_InputMatrix$Annual.N.output..kg.stock.unit.year.[match(cattle_QID,JAC_Livestock_N_Per_Animal_InputMatrix$QID)] 
                     #*JAC_Livestock_N_Per_Animal_InputMatrix$Default[match(cattle_QID,JAC_Livestock_N_Per_Animal_InputMatrix$QID)]
  # Do not divide the cow numbers by 3!!! This is to mimic cows being outside all year
  
  # Determine annual N output (kg) from cattle for each grid cell 
  Cattle[is.na(Cattle)] = 0 # set na to 0
  NCATTLE <- as.matrix(Cattle) %*% as.vector(N_output.cattle)
  #Cattle.N  <- as.data.frame(t(t(Cattle)*(N_output.cattle))) 
  #NCATTLE <- rowSums(Cattle.N, na.rm = TRUE)
  
  # Divide by 100 as we want kg N/ha
  NANIMALS <- (NPIG + NCATTLE + NPOULTRY + NSHEEP + NGOATS + NHORSE)/100
  
  NLIVESTOCK  <- cbind(resampled.agcensus[,c(1:2)], NANIMALS) 
  
  return(NLIVESTOCK)
  
}
# End of function 2 ----

# **************************************************************************** #
# FUNCTION 3: Distribute organic N to broad land use classes ----
# **************************************************************************** #
# This function distributes the total annual org N per km2 (as calculated using 
# the Livestock_N function) to 4 broad land use classes (grassland, spring
# crops, winter crops and other). This done because for each of these broad
# land use classes, NIRAMS defines time-series functions that allow to  
# distribute the annual N values to daily input values.
#
# The distribution of org N is firstly assigned to rough grazing (up to a maximum 
# threshold of 10 kg/ha), then to improved grassland (up to a maximum improved 
# grazing threshold of 250 kg/ha), then to arable (up to threshold of 170 kg/ha).
# The arable N is subsequently distributed to winter crops, spring crops and 
# other in line with the NIRAMS broad land use classes. NOTE: 
# In some cases, the annual org N for a given cell exceeds the total threshold. 
# In these cases, the excess N (above the total threshold) is redistributed 
# evenly to other cells that still have "N holding capacity" (i.e. those cells 
# that have rough grazing, improved grassland and/or arable, but where the 
# livestock numbers are low enough so that the threshold values have not been 
# exceeded). 
# Script initially developed by IP 2017. Modified by MT 2023
#
# INPUT:
# livestock_n     data.frame with the total kg organic N/ha for each grid cell
# landcover       data.frame with the areas (ha) of each NIRAMS land use class  
#                 for each cell in the NIRAMS grid
# roughgraz_max   max org N/ha threshold for rough grazing
# impgraz_max     max org N/ha threshold for improved grassland
# arable_max      max org N/ha threshold for arable
#
# OUTPUT:
# organicn_kgha   data.frame with gridded annual org N input by four broad 
#                 land use classes (grassland, spring crops, winter crops, other)

org_n_distrib <- function(livestock_n, landcover, roughgraz_max = 10,
                          impgraz_max = 250, arable_max = 170){
  
  lc_and_n <- merge(landcover, livestock_n, by = c("x","y"), all = TRUE)
  
  # get total kg of N per grid cell (1 km2 = 100 ha)
  lc_and_n$kg <-  lc_and_n$NANIMALS * 100
  
  # FIRST DISTRIBUTE KG/HA VALUES TO ROUGH GRAZING
  # identify cells with rough grazing area values > 0
  # it is assumed that grazed woodland (LUC21) is included in the rough grazing category
  #roughgrazing <- lc_and_n[,names(lc_and_n) == "LUC_15"] > 0 
  roughgrazing <- rowSums(lc_and_n[,names(lc_and_n) == "LUC_15" | 
                                     names(lc_and_n) == "LUC_21"], na.rm = TRUE) > 0
  roughgrazing[is.na(roughgrazing)] <- FALSE
  # determine area of rough grazing
  roughgrazing.area <- rowSums(lc_and_n[,names(lc_and_n) == "LUC_15" | 
                                         names(lc_and_n) == "LUC_21"], na.rm = TRUE)
  
  # determine kg/ha by dividing total org N with rough grazing area
  lc_and_n$kgha_roughgrazing[roughgrazing] <- lc_and_n$kg[roughgrazing]/ roughgrazing.area[roughgrazing]
  #lc_and_n$kgha_roughgrazing[roughgrazing] <- lc_and_n$kg[roughgrazing]/ lc_and_n[,names(lc_and_n) == "LUC_15"][roughgrazing]
  
  # adjust for max threshold for rough grazing
  lc_and_n$kgha_roughgrazing[lc_and_n$kgha_roughgrazing > roughgraz_max] <- roughgraz_max
  lc_and_n$kgha_roughgrazing[is.na(lc_and_n$kgha_roughgrazing)] <- 0
  
  # kg that are assigned to rough grazing
  lc_and_n$kg_roughgrazing[roughgrazing] <- lc_and_n$kgha_roughgrazing[roughgrazing]*roughgrazing.area[roughgrazing]
  #lc_and_n$kg_roughgrazing[roughgrazing] <- lc_and_n$kgha_roughgrazing[roughgrazing]*lc_and_n[,names(lc_and_n) == "LUC_15"][roughgrazing]
  lc_and_n$kg_roughgrazing[is.na(lc_and_n$kg_roughgrazing)] <- 0
  
  # NEXT DISTRIBUTE KG/HA VALUES TO IMPROVED GRASSLAND FROM ORG N THAT IS LEFT
  # identify cells with improved grassland areas > 0 
  improvedgrazing <- lc_and_n[,names(lc_and_n) == "LUC_16"] > 0
  improvedgrazing[is.na(improvedgrazing)] <- FALSE
  
  # kg N remaining after allocation to rough grazing. 
  # Note: rounding necessary as some results are minus 0.0000000000000000000001 which would ruin the balance
  lc_and_n$kg_after_roughgrazing <- round(lc_and_n$kg - lc_and_n$kg_roughgrazing,10)
  
  # determine kg/ha by dividing remaining org N with improved grazing area
  lc_and_n$kgha_improvedgrazing[improvedgrazing] <- lc_and_n$kg_after_roughgrazing[improvedgrazing] / lc_and_n[,names(lc_and_n) == "LUC_16"][improvedgrazing]
  #lc_and_n$kgha_improvedgrazing[(lc_and_n$kgha_improvedgrazing) < 0] <- 0
  lc_and_n$kgha_improvedgrazing[is.na(lc_and_n$kgha_improvedgrazing)] <- 0
  
  # adjust for max threshold for improved grazing
  lc_and_n$kgha_improvedgrazing[lc_and_n$kgha_improvedgrazing > impgraz_max] <- impgraz_max
  
  # calculate kg N assigned to improved grazing 
  lc_and_n$kg_improvedgrazing[improvedgrazing] <- lc_and_n$kgha_improvedgrazing[improvedgrazing]*lc_and_n[,names(lc_and_n) == "LUC_16"][improvedgrazing]
  lc_and_n$kg_improvedgrazing[is.na(lc_and_n$kg_improvedgrazing)] <- 0
  
  # determine kg N remaining after allocation to rough grazing and improved grazing. 
  lc_and_n$kg_after_improvedgrazing <- round(lc_and_n$kg - (lc_and_n$kg_roughgrazing + lc_and_n$kg_improvedgrazing) ,9)
  
  lc_and_n$kg_sum_grazing    <- lc_and_n$kg_roughgrazing + lc_and_n$kg_improvedgrazing
  lc_and_n$kg_ha_sum_grazing <- lc_and_n$kgha_roughgrazing + lc_and_n$kgha_improvedgrazing
  
  # NEXT DISTRIBUTE KG/HA VALUES TO ARABLE FROM ORG N THAT IS LEFT
  #  Note: set aside (LUC_14) was previously included in this and then assigned
  #  to the other crops broad land use category (see later). However, 'other crops'  
  #  is not considered when distributing annual orgN to daily values, i.e. there 
  #  is no idealized time-series for orgN and 'other crops'. Therefore set aside 
  #  has been removed here. Alternatively set aside could be included in the grass
  #  category if set aside is actually being grazed; something I'm not sure about.
  #
  # identify cells with arable area > 0
  arable <- rowSums(lc_and_n[,c(3:15)], na.rm = TRUE) > 0 # rowSums(lc_and_n[,c(3:16)], na.rm = TRUE) > 0
  arable[is.na(arable)] <- FALSE
  # calculate arable area for each cell
  arable.land <- rowSums(lc_and_n[,c(3:15)], na.rm = TRUE)
  
  # determine kg/ha by dividing remaining org N with arable land area
  lc_and_n$kgha_arable[arable] <- round(lc_and_n$kg_after_improvedgrazing[arable]/ arable.land[arable],10)
  lc_and_n$kgha_arable[is.na(lc_and_n$kgha_arable)] <- 0
  #lc_and_n$kgha_arable[lc_and_n$kgha_arable < 0] <- 0
  
  # adjust for max threshold for arable
  lc_and_n$kgha_arable[lc_and_n$kgha_arable > arable_max] <- arable_max
  lc_and_n$kg_arable[arable] <- lc_and_n$kgha_arable[arable]*arable.land[arable]
  lc_and_n$kg_arable[is.na(lc_and_n$kg_arable)] <- 0
  
  # Distribute arable values to spring crops, winter crops and other
  # NOTE: 'other' not considered for orgN in NIRAMS, hence this step has now been removed.
  # Spring crops
  lc_and_n$springcrops.area <- rowSums(lc_and_n[,c(1,3,5,7,9,10,11,12,13)+2], na.rm = TRUE)
  springcrops <- lc_and_n$springcrops.area > 0
  springcrops[is.na(springcrops)] <- FALSE
  # Winter crops
  lc_and_n$wintercrops.area <- rowSums(lc_and_n[,c(2,4,6,8)+2], na.rm = TRUE)
  wintercrops <- lc_and_n$wintercrops.area > 0
  wintercrops[is.na(wintercrops)] <- FALSE
  # Other 
  lc_and_n$othercrops.area <- lc_and_n[,16] # corresponds to LUC_14 (set aside)
  othercrops <- lc_and_n$othercrops.area > 0
  othercrops[is.na(othercrops)] <- FALSE
  
  # Distribute arable kg/ha to the different broad land classes 
  lc_and_n$kgha_springcrops[springcrops] <- lc_and_n$kgha_arable[springcrops]
  lc_and_n$kgha_wintercrops[wintercrops] <- lc_and_n$kgha_arable[wintercrops]
  lc_and_n$kgha_othercrops[othercrops]   <- lc_and_n$kgha_arable[othercrops]
  
  lc_and_n$kg_springcrops[springcrops]  <- lc_and_n$kgha_arable[springcrops] * lc_and_n$springcrops.area[springcrops]
  lc_and_n$kg_wintercrops[wintercrops]  <- lc_and_n$kgha_arable[wintercrops] * lc_and_n$wintercrops.area[wintercrops]
  lc_and_n$kg_othercrops <- 0  # ?
  #lc_and_n$kg_othercrops[othercrops]    <- lc_and_n$kgha_arable[othercrops] * lc_and_n$othercrops.area[othercrops]
  
  # Determine excess org N after distributing to rough grazing, improved grassland and arable
  lc_and_n$kg_excess <- round(lc_and_n$kg_after_improvedgrazing - lc_and_n$kg_arable,0)
  
  # Sum to get total excess 
  #   Note: Ina originally only summed the excess for cells that had some kind 
  #         of land cover in them (i.e., lc_and_n$sums > 0). For some cells, the 
  #         sum of land covers is 0 (NA) but the N application rate is >0. This 
  #         mainly occurs because the census data on N application is at 2km 
  #         resolution whereas the land cover data is at 1 km resolution. It 
  #         therefore seems fair to determine the total excess by summing the 
  #         excess for all cells (incl. those for which lc_and_n$sums=NA)
  #       
  #excess_total_kg <-  sum(lc_and_n$kg_excess[lc_and_n$sums > 0], na.rm = TRUE)
  excess_total_kg <-  sum(lc_and_n$kg_excess, na.rm = TRUE)
  
  # # find cells for which land cover is NA but N application is >0.
  # cell_id <- lc_and_n[which(is.na(lc_and_n$sums) & lc_and_n$kg>0),c(1,2,29:31)]
  # nid <- length(cell_id[,1])
  # cell_vec <- vector(length=nid)
  # 
  # for (i in 1:nid) {
  #   cell_ext <-   lc_and_n[which(lc_and_n$x>=cell_id$x[i]-1000 & lc_and_n$x<=cell_id$x[i]+1000 & 
  #                                  lc_and_n$y>=cell_id$y[i]-1000 & lc_and_n$y<=cell_id$y[i]+1000 &
  #                                  lc_and_n$kg==cell_id$kg[i]),
  #                          c(1,2,29:31)]
  #   cell_vec[i] <- sum(cell_ext$sums>0,na.rm=TRUE)#cell_ext$sums
  # }

  
  # RESDISTRIBUTE EXCESS N
  # Redistribute excess N according to similar routine as for original 
  # distribution within the 1 km2 cells.  
  # In the first round of calculations, the total amount of N was divided by 
  # the area of rough grazing; if the resulting kg/ha was bigger than the 
  # threshold, then the remaining N was passed on to improved grassland. Hence,
  # any cell that had a kg/ha above the threshold (and thus was set equal to the 
  # threshold), will not have any further capacity to receive more N. Only cells
  # with rough razing that have been assigned kg/ha values below the threshold  
  # have capacity to received excess N. This means that any excess N will be 
  # allocated to cells that still have capacity to receive N. So effectively, 
  # this means reducing livestock numbers in cells with excess N and increase 
  # them in cells with additional capacity. 
  #
  # Case 1: Area of rough grazing is 20 ha. Total N applied to cell is 400 kg N
  #   Application rate = 400 kg/20 ha = 20 kg/ha 
  #    i.e. above the threshold so application rate set to threshold of 10 kg/ha 
  #    and the remaining 200 kg N is distributed to improved grass. There is no
  #    further capacity to increase rough grazing N. 
  #
  # Case 2: Area of rough grazing is 20 ha. Total N applied to cell is 100 kg N
  #   Application rate = 100 kg/20 ha = 5 kg/ha 
  #    i.e. below the threshold so all N is assigned to rough grazing; there is 
  #    no excess N. There is further capacity to increase rough grazing N by  
  #    100 kg before reaching capacity
  # 
  
  # UPDATE ROUGH GRAZING 
  # First add excess N to cells with rough grazing and with capacity
  lc_and_n$kgha_roughgrazing_available <- roughgraz_max - lc_and_n$kgha_roughgrazing # the remaining "capacity" for rough grazing N
  lc_and_n$kg_roughgrazing_available   <- lc_and_n$kgha_roughgrazing_available * roughgrazing.area #* lc_and_n[,names(lc_and_n) == "LUC_15"]
  
  # Evenly distribute and update
  total_kg_roughgrazing_available       <- sum(lc_and_n$kg_roughgrazing_available, na.rm = TRUE) # N amount that rough grazing land can still receive
  total_to_redistribute_to_roughgrazing <- min(total_kg_roughgrazing_available,excess_total_kg)
  
  # fraction of kg N to redistribute to rough grazing vs the capacity (total_to_redistribute_to_roughgrazing_vs_available)
  f_rgm2rgc  <- total_to_redistribute_to_roughgrazing/total_kg_roughgrazing_available
  
  # Adjust kg N for rough grazing 
  lc_and_n$kg_roughgrazing <- lc_and_n$kg_roughgrazing + lc_and_n$kg_roughgrazing_available * f_rgm2rgc
  
  # Update kg/ha for rough grazing 
  lc_and_n$kgha_roughgrazing <- lc_and_n$kg_roughgrazing/roughgrazing.area #/lc_and_n[,names(lc_and_n) == "LUC_15"]
  
  # Update excess kg 
  excess_total_kg <- excess_total_kg - total_to_redistribute_to_roughgrazing
  
  # UPDATE IMPROVED GRAZING
  # Next add any further excess N to cells with improved grazing and with capacity
  lc_and_n$kgha_improvedgrazing_available <- impgraz_max - lc_and_n$kgha_improvedgrazing # determine remaining capacity for improved grazing
  lc_and_n$kg_improvedgrazing_available   <- lc_and_n$kgha_improvedgrazing_available * lc_and_n[,names(lc_and_n) == "LUC_16"]
  
  # Determine N amount that improved grazing land can still receive
  total_kg_improvedgrazing_available         <- sum(lc_and_n$kg_improvedgrazing_available, na.rm = TRUE)
  total_to_redistribute_to_improvedgrassland <- min(total_kg_improvedgrazing_available,excess_total_kg)   
  
  # fraction of kg N to redistribute to improved grazing vs the capacity
  f_igm2igc <- total_to_redistribute_to_improvedgrassland/total_kg_improvedgrazing_available
  
  # update kg N for improved grazing  
  lc_and_n$kg_improvedgrazing <- lc_and_n$kg_improvedgrazing + lc_and_n$kg_improvedgrazing_available * f_igm2igc
  
  # update kg/ha for improved grazing
  lc_and_n$kgha_improvedgrazing <- lc_and_n$kg_improvedgrazing/lc_and_n[,names(lc_and_n) == "LUC_16"]  
  
  lc_and_n[lc_and_n == "NaN"] <- NA
  
  # kg for total grassland
  lc_and_n$kg_grass <- lc_and_n$kg_roughgrazing + lc_and_n$kg_improvedgrazing
  
  # Arable not updated as no excess N

  
  # save annual organic kg N grid-by-grid values for the 3 broad land use classes
  # note: other crops removed so only 3 broad land use classes
  organicn <- lc_and_n[names(lc_and_n) %in% c("x","y","kg_grass","kg_springcrops","kg_wintercrops","kg_othercrops")]
  #organicn[is.na(organicn)] <- 0
  #organicn <- lc_and_n[names(lc_and_n) %in% c("x","y","kg_grass","kg_springcrops","kg_wintercrops")]
  organicn[is.na(lc_and_n$sums),-c(1,2)] <- NA
  
  for (i in c(3:6)){ 
    na_id <- which(is.na(organicn[,i]) & !is.na(lc_and_n$sums))
    if(length(na_id)>0) organicn[na_id,i] <- 0
  }
  
  # determine org N application rate in kg/ha 
  organicn_kgha <- organicn
  
  # overall kg per ha for the raster cell not just for the individual land use
  #organicn_kgha[,c(3:6)] <- organicn_kgha[,c(3:6)]/100
  organicn_kgha[,c("kg_grass","kg_springcrops","kg_wintercrops","kg_othercrops")] <- organicn_kgha[,c("kg_grass","kg_springcrops","kg_wintercrops", "kg_othercrops")]/100
  names(organicn_kgha) <- gsub("kg_","kg_ha_",names(organicn))# c("kg_ha_wintercrops","kg_ha_othercrops","kg_ha_grass","kg_ha_springcrops")
  
  
  return(organicn_kgha)
}
# End of function 3 ----




# **************************************************************************** #
# FUNCTION 3A: Distribute organic N to broad land use classes CRAFTY ----
# **************************************************************************** #
# The distribution of organic N is modified here to be based on CRAFTY land use.
# The total annual org N per km2 is distributed to 4 broad land use classes: 
# - Grassland:    CRAFTY class 3 (Extensive pasture = rough grazing) and class 6 
#                 (intensive pasture = improved grassland), 
# - Spring crops: CRAFTY classes 1, 2, 4, 5 and 13 (Bioenergy, Extensive arable, 
#                 Intensive arable (fodder), Intensive arable (food), Sustainable arable
# - Winter crops: CRAFTY classes 1, 2, 4, 5 and 13 (Bioenergy, Extensive arable, 
#                 Intensive arable (fodder), Intensive arable (food), Sustainable arable
# - Other:        NA.
#
# The distribution of org N is firstly assigned to rough grazing (CRAFTY class 3)
# up to a maximum threshold of 10 kg/ha, then to improved grassland (CRAFTY 
# class 6) up to a maximum threshold of 250 kg/ha, then to arable up to threshold 
# of 170 kg/ha. The arable N is subsequently distributed to winter crops, spring 
# crops and other. Because CRAFTY does not estimate actual crop distributions, 
# it is assumed that 70% of arable is spring crops and 30% is winter crops 
# (based on AGCensus data 2015-2019). 
# NOTE: It is assumed that the future livestock numbers remain constant and are 
# the same as in 2019. This is quie a big assumption and possible not very 
# consistent with the land use scenario/SSP - NEED TO CHECK BETTER APPROACH
# In some cases, the annual org N for a given cell exceeds the total threshold. 
# In these cases, the excess N (above the total threshold) is redistributed 
# evenly to other cells that still have "N holding capacity" (i.e. those cells 
# that have rough grazing, improved grassland and/or arable, but where the 
# livestock numbers are low enough so that the threshold values have not been 
# exceeded). 
# MT 2024
#
# INPUT:
# livestock_n     data.frame with the total kg organic N/ha for each grid cell for a given year
# my.crafty       raster with dominant future land use classes from CRAFTY for a given year
# roughgraz_max   max org N/ha threshold for rough grazing
# impgraz_max     max org N/ha threshold for improved grassland
# arable_max      max org N/ha threshold for arable
#
# OUTPUT:
# organicn_kgha   data.frame with gridded annual org N input by four broad 
#                 land use classes (grassland, spring crops, winter crops, other)


future_org_n_distrib <- function(livestock_n, my.crafty, roughgraz_max = 10,
                          impgraz_max = 250, arable_max = 170){
  
  # Combine CRAFTY and org N data in one data.frame
  crafty.df   <- as.data.frame(my.crafty, xy=TRUE)
  names(crafty.df)[3] <- "crafty"
  lc_and_n <- merge(crafty.df, livestock_n, by = c("x","y"), all = TRUE)
  
  # Get total kg of N per grid cell (1 km2 = 100 ha). 
  # Not sure this is really needed here because all cells only have one
  # dominant land use (not a distribution/fraction of land uses). 
  lc_and_n$kg <-  lc_and_n$NANIMALS * 100   
  
  # ----------------------------------------------------- #
  # 1: FIRST DISTRIBUTE KG/HA VALUES TO ROUGH GRAZING
  # ----------------------------------------------------- #
  # Identify cells with rough grazing area values > 0
  # It is assumed that extensive pastoral (class 3) in CRAFTY is rough grazing
  # CONSIDER WHETHER CRAFTY CLASS 14 SHOULD BE ROUGH GRAZING
  roughgrazing <- (lc_and_n$crafty == 3)  # identify cells with rough grazing
  roughgrazing[is.na(roughgrazing)] <- FALSE
  
  # Get org N kg/ha for rough grazing
  lc_and_n$kgha_roughgrazing[roughgrazing] <- lc_and_n$NANIMALS[roughgrazing]
  
  # Determine area of rough grazing 
  roughgrazing.area <- roughgrazing*100 # cells are either 0 ha or 100 ha rough grazing
 
  # Adjust for max threshold for rough grazing
  lc_and_n$kgha_roughgrazing[lc_and_n$kgha_roughgrazing > roughgraz_max] <- roughgraz_max
  lc_and_n$kgha_roughgrazing[is.na(lc_and_n$kgha_roughgrazing)] <- 0
  
  # kg that are assigned to rough grazing - IS THIS ONLY NEEDED FOR MASS BALANCE?
  lc_and_n$kg_roughgrazing[roughgrazing] <- lc_and_n$kgha_roughgrazing[roughgrazing]*100 # roughgrazing.area[roughgrazing]
  lc_and_n$kg_roughgrazing[is.na(lc_and_n$kg_roughgrazing)] <- 0
  
  # kg N remaining after allocation to rough grazing. 
  # Note: rounding necessary as some results are minus 0.0000000000000000000001 which would ruin the balance
  lc_and_n$kg_after_roughgrazing <- round(lc_and_n$kg - lc_and_n$kg_roughgrazing,10)
  
  # ------------------------------------------------------ #
  # 2: NEXT DISTRIBUTE KG/HA VALUES TO IMPROVED GRASSLAND 
  # ------------------------------------------------------ #
  # Identify cells with improved grassland areas > 0 
  improvedgrazing <- (lc_and_n$crafty == 6)  
  improvedgrazing[is.na(improvedgrazing)] <- FALSE
  
  # Determine org N kg/ha for improved grazing
  lc_and_n$kgha_improvedgrazing[improvedgrazing] <- lc_and_n$NANIMALS[improvedgrazing]
  lc_and_n$kgha_improvedgrazing[is.na(lc_and_n$kgha_improvedgrazing)] <- 0
  
  # Area of improved grazing
  improvedgrazing.area <- improvedgrazing*100 
  
  # Adjust for max threshold for improved grazing
  lc_and_n$kgha_improvedgrazing[lc_and_n$kgha_improvedgrazing > impgraz_max] <- impgraz_max
  
  # Calculate kg N assigned to improved grazing 
  lc_and_n$kg_improvedgrazing[improvedgrazing] <- lc_and_n$kgha_improvedgrazing[improvedgrazing]*100 
  lc_and_n$kg_improvedgrazing[is.na(lc_and_n$kg_improvedgrazing)] <- 0
  
  # Determine kg N remaining after allocation to rough grazing and improved grazing. 
  lc_and_n$kg_after_improvedgrazing <- round(lc_and_n$kg - (lc_and_n$kg_roughgrazing + lc_and_n$kg_improvedgrazing) ,9)
  
  # ------------------------------------------------------ #
  # 3: NEXT DISTRIBUTE KG/HA VALUES TO ARABLE 
  # ------------------------------------------------------ #
  # Identify cells with arable area (CRAFTY class 1, 2, 4, 5 and 13)
  arable <- (lc_and_n$crafty == 1 | lc_and_n$crafty == 2 | lc_and_n$crafty == 4 | lc_and_n$crafty == 5 | lc_and_n$crafty == 13)   
  arable[is.na(arable)] <- FALSE
  
  # Determine kg/ha for arable land 
  lc_and_n$kgha_arable[arable] <- lc_and_n$NANIMALS[arable]
  
  arable.area <- arable*100 
  
  # Adjust for max threshold for arable
  lc_and_n$kgha_arable[lc_and_n$kgha_arable > arable_max] <- arable_max
  lc_and_n$kg_arable[arable] <- lc_and_n$kgha_arable[arable]*100 # arable.land[arable]
  lc_and_n$kg_arable[is.na(lc_and_n$kg_arable)] <- 0
  

  # Determine excess org N after distributing to rough grazing, improved grassland and arable
  lc_and_n$kg_excess <- round(lc_and_n$kg_after_improvedgrazing - lc_and_n$kg_arable,0)
  
  # Sum to get total excess 
  #excess_total_kg <-  sum(lc_and_n$kg_excess[!is.na(lc_and_n$crafty)], na.rm = TRUE) # use this if N application on cells with no crafty land use should be ignored
  excess_total_kg <-  sum(lc_and_n$kg_excess, na.rm = TRUE)
  
  # ========================================================================= #
  # RESDISTRIBUTE EXCESS N
  # ========================================================================= #
  # Redistribute excess N according to similar routine as for original 
  # distribution within the 1 km2 cells.  
  # In the first round of calculations, the total amount of N was divided by 
  # the area of rough grazing; if the resulting kg/ha was bigger than the 
  # threshold, then the remaining N was passed on to improved grassland. Hence,
  # any cell that had a kg/ha above the threshold (and thus was set equal to the 
  # threshold), will not have any further capacity to receive more N. Only cells
  # with rough razing that have been assigned kg/ha values below the threshold  
  # have capacity to received excess N. This means that any excess N will be 
  # allocated to cells that still have capacity to receive N. So effectively, 
  # this means reducing livestock numbers in cells with excess N and increase 
  # them in cells with additional capacity. 
  
  # ------------------------------------------------------ #
  # 5: UPDATE ROUGH GRAZING 
  # ------------------------------------------------------ #
  # First add excess N to cells with rough grazing and with capacity
  lc_and_n$kgha_roughgrazing_available <- roughgraz_max - lc_and_n$kgha_roughgrazing # the remaining "capacity" for rough grazing N
  lc_and_n$kg_roughgrazing_available   <- lc_and_n$kgha_roughgrazing_available * roughgrazing.area 
  
  # Evenly distribute and update
  total_kg_roughgrazing_available       <- sum(lc_and_n$kg_roughgrazing_available, na.rm = TRUE) # N amount that rough grazing land can still receive
  total_to_redistribute_to_roughgrazing <- min(total_kg_roughgrazing_available,excess_total_kg)
  
  # fraction of kg N to redistribute to rough grazing vs the capacity (total_to_redistribute_to_roughgrazing_vs_available)
  f_rgm2rgc  <- total_to_redistribute_to_roughgrazing/total_kg_roughgrazing_available
  
  # Adjust kg N for rough grazing - if total_kg_roughgrazing_available is less
  # than the total excess, all rough grazing is set to maximum application rate
  lc_and_n$kg_roughgrazing <- lc_and_n$kg_roughgrazing + lc_and_n$kg_roughgrazing_available * f_rgm2rgc
  
  # Update kg/ha for rough grazing 
  lc_and_n$kgha_roughgrazing <- lc_and_n$kg_roughgrazing/roughgrazing.area #/lc_and_n[,names(lc_and_n) == "LUC_15"]
  
  # Update excess kg 
  excess_total_kg <- excess_total_kg - total_to_redistribute_to_roughgrazing
  
  # ------------------------------------------------------ #
  # 6: UPDATE IMPROVED GRAZING
  # ------------------------------------------------------ #
  # Next add any further excess N to cells with improved grazing and with capacity
  lc_and_n$kgha_improvedgrazing_available <- impgraz_max - lc_and_n$kgha_improvedgrazing # determine remaining capacity for improved grazing
  lc_and_n$kg_improvedgrazing_available   <- lc_and_n$kgha_improvedgrazing_available * improvedgrazing.area
  
  # Determine N amount that improved grazing land can still receive
  total_kg_improvedgrazing_available         <- sum(lc_and_n$kg_improvedgrazing_available, na.rm = TRUE)
  total_to_redistribute_to_improvedgrassland <- min(total_kg_improvedgrazing_available,excess_total_kg)   
  
  # fraction of kg N to redistribute to improved grazing vs the capacity
  f_igm2igc <- total_to_redistribute_to_improvedgrassland/total_kg_improvedgrazing_available
  
  # update kg N for improved grazing  
  lc_and_n$kg_improvedgrazing <- lc_and_n$kg_improvedgrazing + lc_and_n$kg_improvedgrazing_available * f_igm2igc
  
  # update kg/ha for improved grazing
  lc_and_n$kgha_improvedgrazing <- lc_and_n$kg_improvedgrazing/improvedgrazing.area  

  # # ------------------------------------------------------ #
  # # 7: UPDATE ARABLE 
  # # ------------------------------------------------------ #
  # # Update excess kg 
  # excess_total_kg <- excess_total_kg - total_to_redistribute_to_improvedgrassland
  # 
  # # Next add any further excess N to cells with improved grazing and with capacity
  # lc_and_n$kgha_arable_available <- arable_max - lc_and_n$kgha_arable # determine remaining capacity for arable
  # lc_and_n$kg_arable_available   <- lc_and_n$kgha_arable_available * arable.area
  # 
  # # Determine N amount that improved grazing land can still receive
  # total_kg_arable_available         <- sum(lc_and_n$kg_arable_available, na.rm = TRUE)
  # total_to_redistribute_to_arableland <- min(total_kg_arable_available,excess_total_kg)   
  # 
  # # fraction of kg N to redistribute to improved grazing vs the capacity
  # f_am2ac <- total_to_redistribute_to_arableland/total_kg_arable_available
  # 
  # # update kg N for improved grazing  
  # lc_and_n$kg_arable <- lc_and_n$kg_arable + lc_and_n$kg_arable_available * f_am2ac
  # 
  # # update kg/ha for improved grazing
  # lc_and_n$kgha_arable <- lc_and_n$kg_arable/arable.area  
  
  # ========================================================================= #
  # ASSIGN N VALUES TO BROAD LAND CLASSES
  # ========================================================================= #

  lc_and_n[lc_and_n == "NaN"] <- NA
  
  # # Not needed? 
  # lc_and_n$kg_sum_grazing    <- lc_and_n$kg_roughgrazing + lc_and_n$kg_improvedgrazing
  # lc_and_n$kg_ha_sum_grazing <- lc_and_n$kgha_roughgrazing + lc_and_n$kgha_improvedgrazing
  
  # GRASSLAND
  # kg for total grassland
  lc_and_n$kg_grass <- lc_and_n$kg_roughgrazing + lc_and_n$kg_improvedgrazing
  
  # SPRING CROPS, WINTER CROPS, OTHER CROPS
  # Distribute arable kg/ha values to spring crops, winter crops and other crops
  # Application per hectare is assumed the same for all land classes.
  # Assumed that area of spring crop is 70% and winter crop 30%
  lc_and_n$kgha_springcrops[arable]  <- lc_and_n$kgha_arable[arable]
  lc_and_n$kgha_wintercrops[arable]  <- lc_and_n$kgha_arable[arable]
  lc_and_n$kgha_othercrops[arable]   <- lc_and_n$kgha_arable[arable]
  
  lc_and_n$kg_springcrops[arable]  <- lc_and_n$kgha_arable[arable] * 70
  lc_and_n$kg_wintercrops[arable]  <- lc_and_n$kgha_arable[arable] * 30
  lc_and_n$kg_othercrops <- 0  # 
  
  
  # Save annual organic kg N grid-by-grid values for the broad land use classes
  organicn <- lc_and_n[names(lc_and_n) %in% c("x","y","kg_grass","kg_springcrops","kg_wintercrops","kg_othercrops")]
  organicn[is.na(lc_and_n$crafty),-c(1,2)] <- NA
  
  for (i in c(3:6)){ 
    #na_id <- which(is.na(organicn[,i]) & !is.na(lc_and_n$sums))
    na_id <- which(is.na(organicn[,i]) & !is.na(lc_and_n$crafty))
    if(length(na_id)>0) organicn[na_id,i] <- 0
  }
  
  # determine org N application rate in kg/ha 
  organicn_kgha <- organicn
  
  # Overall kg per ha for the raster cell not just for the individual land use
  organicn_kgha[,c("kg_grass","kg_springcrops","kg_wintercrops","kg_othercrops")] <- organicn_kgha[,c("kg_grass","kg_springcrops","kg_wintercrops", "kg_othercrops")]/100
  names(organicn_kgha) <- gsub("kg_","kg_ha_",names(organicn)) 
  
  
  return(organicn_kgha)
}
# End of function 3A ----




# **************************************************************************** #
# FUNCTION 4: Distribute inorganic N and uptake to broad land use classes ----
# **************************************************************************** #
# This function calculates the gridded total annual inorganic N and N uptake 
# (kg/ha) and distributes these to 4 broad land use classes (grassland, spring
# crops, winter crops and other). For each of these broad land uses, NIRAMS 
# defines time-series functions that allow to distribute the annual N values to
# daily input values.
#
# INPUT:
# landcover       data.frame with the areas (ha) of each NIRAMS land use class  
#                 for each cell in the NIRAMS grid
# Inorg_uptake    data.frame with annual inorganic N application rates and 
#                 N uptake rates by NIRAMS land use classes
#
# OUTPUT:
# res             list with gridded annual inorganic N input and N uptake by 4  
#                 broad land use classes (grassland, spring crops, winter crops, other)
#
# Script by IP 2017

Inorg_N_Application_Uptake <- function (landcover, Inorg_uptake){
  
  uptake <- landcover[,-c(1,2,29)]
  uptake <- as.data.frame(t(t(uptake)*Inorg_uptake$N_Uptake_kgpyr[1:26]))
  
  application <- landcover[,-c(1,2,29)]
  application <- as.data.frame(t(t(application)*Inorg_uptake$Inorg_N_Applied_kgpyr[1:26]))
  
  springcrops.uptake      <- uptake[,c(1,3,5,7,9,10,11,12,13)]  # removed short rotation coppice (LUC18) from spring crops
  springcrops.uptake$kgha <- rowSums(springcrops.uptake, na.rm = TRUE) / 100
  springcrops.application      <- application[,c(1,3,5,7,9,10,11,12,13)]
  springcrops.application$kgha <- rowSums(springcrops.application, na.rm = TRUE) / 100
  
  wintercrops.uptake      <- uptake[,c(2,4,6,8)]
  wintercrops.uptake$kgha <- rowSums(wintercrops.uptake, na.rm = TRUE) / 100
  wintercrops.application      <- application[,c(2,4,6,8)]
  wintercrops.application$kgha <- rowSums(wintercrops.application, na.rm = TRUE) / 100
  
  other.uptake      <- uptake[,c(14,17:26)] # added short rotation coppice (LUC18) and grazed woodland (LUC21) to other
  other.uptake$kgha <- rowSums(other.uptake, na.rm = TRUE) / 100
  other.application <- application[,c(14,17,19:26)]   # inorg N application not considered for other land use
  other.application$kgha <- rowSums(other.application, na.rm = TRUE) / 100
  
  grass.uptake      <- uptake[,c(15,16)]
  grass.uptake$kgha <- rowSums(grass.uptake, na.rm = TRUE) / 100
  grass.application      <- application[,c(15,16,18)]  # added short rotation coppice (LUC18) to grass application
  grass.application$kgha <- rowSums(grass.application, na.rm = TRUE) / 100
  
  in_sp <- cbind(landcover[,c(1:2)],springcrops.application$kgha)
  in_wi <- cbind(landcover[,c(1:2)],wintercrops.application$kgha)
  in_ot <- cbind(landcover[,c(1:2)],other.application$kgha)
  in_gr <- cbind(landcover[,c(1:2)],grass.application$kgha)
  
  up_sp <- cbind(landcover[,c(1:2)],springcrops.uptake$kgha)
  up_wi <- cbind(landcover[,c(1:2)],wintercrops.uptake$kgha)
  up_ot <- cbind(landcover[,c(1:2)],other.uptake$kgha)
  up_gr <- cbind(landcover[,c(1:2)],grass.uptake$kgha)
  
  res <- list(in_sp,in_wi,in_ot,in_gr,up_sp,up_wi,up_ot,up_gr)
  names(res) <- c("in_sp","in_wi","in_ot","in_gr","up_sp","up_wi","up_ot","up_gr")
  
  for (i in c(1:length(res))){
    res[[i]][,3][is.na(landcover$sums)] <- NA
  }
  return(res)
  
  
}
# end of function 4 ----


# **************************************************************************** #
# FUNCTION 4A: Add inorganic N and uptake to broad land use classes CRAFTY ----
# **************************************************************************** #
# This function calculates the gridded total annual inorganic N and N uptake 
# (kg/ha) and distributes these to 4 broad land use classes (grassland, spring
# crops, winter crops and other). For each of these broad land uses, NIRAMS 
# defines time-series functions that allow to distribute the annual N values to
# daily input values.
#
# INPUT:
# my.crafty       raster with dominant future land use classes from CRAFTY for 
#                 a given year
# Inorg_uptake    data.frame with annual inorganic N application rates and 
#                 N uptake rates by NIRAMS land use classes
#
# OUTPUT:
# res             list with gridded annual inorganic N input and N uptake by 4  
#                 broad land use classes (grassland, spring crops, winter crops, other)
#
# Script by MT 2024


Future_Inorg_N_Application_Uptake <- function (my.crafty, Inorg_uptake){

  # CRAFTY raster to data.frame
  crafty.df   <- as.data.frame(my.crafty, xy=TRUE)
  names(crafty.df)[3] <- "crafty"
  
  # Because the CRAFTY data only gives the dominant land use class, this is just 
  # a task of assigning uptake and fertiliser application rates to the right
  # CRAFTY classes as defined in Inorg_uptake.
  
  # ------------------------------------------------------- #
  # 1: IDENTIFY THE CELLS WITH GRASS, ARABLE OR OTHER 
  # ------------------------------------------------------- #
  # Identify cells with rough grazing. It is assumed that extensive pastoral
  #  (class 3) in CRAFTY is rough grazing. SHOULD CRAFTY CLASS 14 ALSO BE ROUGH GRAZING?
  roughgrazing <- (crafty.df$crafty == 3)  # identify cells with rough grazing
  roughgrazing[is.na(roughgrazing)] <- FALSE

  # Identify cells with improved grassland areas > 0 
  improvedgrazing <- (crafty.df$crafty == 6)  
  improvedgrazing[is.na(improvedgrazing)] <- FALSE
  
  # Identify cells with arable area (CRAFTY class 1, 2, 4, 5 and 13)
  arable <- (crafty.df$crafty == 1 | crafty.df$crafty == 2 | crafty.df$crafty == 4 | crafty.df$crafty == 5 | crafty.df$crafty == 13)   
  arable[is.na(arable)] <- FALSE
  
  # Identify cells with woodland area (CRAFTY class 0, 7-12)
  woodland <- (crafty.df$crafty == 0 | crafty.df$crafty == 7 | crafty.df$crafty == 8 | crafty.df$crafty == 9 | 
                 crafty.df$crafty == 10 | crafty.df$crafty == 11 | crafty.df$crafty == 12 )   
  woodland[is.na(woodland)] <- FALSE
  
  # Identify cells with "wildscape" area (CRAFTY class -1 and 14)
  wild <- (crafty.df$crafty == -1 | crafty.df$crafty == 14)   
  wild[is.na(wild)] <- FALSE
  
  # Identify cells with built area (CRAFTY class 15)
  built <- (crafty.df$crafty == 15 )   
  built[is.na(built)] <- FALSE
  
  # -------------------------------------------------------- #
  # 2: EXTRACT UPTAKE AND FERTILSER RATES FOR THE LAND USES
  # -------------------------------------------------------- #
  # Note, all CRAFTY arable classes are currently given the same values 
  roughgrass_uptake <- Inorg_uptake$N_Uptake_kgpyr[Inorg_uptake$CRAFTY_INDEX == 3]
  roughgrass_inorg  <- Inorg_uptake$Inorg_N_Applied_kgpyr[Inorg_uptake$CRAFTY_INDEX == 3]
  
  imprgrass_uptake <- Inorg_uptake$N_Uptake_kgpyr[Inorg_uptake$CRAFTY_INDEX == 6]
  imprgrass_inorg  <- Inorg_uptake$Inorg_N_Applied_kgpyr[Inorg_uptake$CRAFTY_INDEX == 6]
    
  arable_uptake <- Inorg_uptake$N_Uptake_kgpyr[Inorg_uptake$CRAFTY_INDEX == 2]          
  arable_inorg  <- Inorg_uptake$Inorg_N_Applied_kgpyr[Inorg_uptake$CRAFTY_INDEX == 2]
  
  wood_uptake <- Inorg_uptake$N_Uptake_kgpyr[Inorg_uptake$CRAFTY_INDEX == 7]          
  wood_inorg  <- Inorg_uptake$Inorg_N_Applied_kgpyr[Inorg_uptake$CRAFTY_INDEX == 7]
  
  wild_uptake <- Inorg_uptake$N_Uptake_kgpyr[Inorg_uptake$CRAFTY_INDEX == 14]          
  wild_inorg  <- Inorg_uptake$Inorg_N_Applied_kgpyr[Inorg_uptake$CRAFTY_INDEX == 14]
  
  built_uptake <- Inorg_uptake$N_Uptake_kgpyr[Inorg_uptake$CRAFTY_INDEX == 15]          
  built_inorg  <- Inorg_uptake$Inorg_N_Applied_kgpyr[Inorg_uptake$CRAFTY_INDEX == 15]
  
  # ---------------------------------------------------------------- #
  # 3: ASSIGN UPTAKE AND FERTILISER RATES TO BROAD LAND USE CLASSES
  # ---------------------------------------------------------------- #
  # Note: Assumed 70% of arable is spring crops and 30% is winter crops
  # GRASS
  grass.uptake      <- crafty.df[,c(1,2)]
  grass.uptake$kgha <- 0
  grass.uptake$kgha[roughgrazing]    <-  roughgrass_uptake
  grass.uptake$kgha[improvedgrazing] <-  imprgrass_uptake
    
  grass.application      <- crafty.df[,c(1,2)]
  grass.application$kgha <- 0
  grass.application$kgha[roughgrazing]    <-  roughgrass_inorg
  grass.application$kgha[improvedgrazing] <-  imprgrass_inorg
  
  # SPRING CROPS
  springcrops.uptake      <- crafty.df[,c(1,2)]
  springcrops.uptake$kgha <- 0
  springcrops.uptake$kgha[arable]    <-  arable_uptake*0.7

  springcrops.application      <- crafty.df[,c(1,2)]
  springcrops.application$kgha <- 0
  springcrops.application$kgha[arable]  <-  arable_inorg*0.7
  
  # WINTER CROPS
  wintercrops.uptake      <- crafty.df[,c(1,2)]
  wintercrops.uptake$kgha <- 0
  wintercrops.uptake$kgha[arable]    <-  arable_uptake*0.3
  
  wintercrops.application      <- crafty.df[,c(1,2)]
  wintercrops.application$kgha <- 0
  wintercrops.application$kgha[arable]  <-  arable_inorg*0.3
  
  # OTHER
  other.uptake      <- crafty.df[,c(1,2)]
  other.uptake$kgha <- 0
  other.uptake$kgha[woodland]    <-  wood_uptake
  other.uptake$kgha[wild]        <-  wild_uptake
  other.uptake$kgha[built]       <-  built_uptake
  
  other.application      <- crafty.df[,c(1,2)]
  other.application$kgha <- 0
  other.application$kgha[woodland]  <-  wood_inorg
  other.application$kgha[wild]      <-  wild_inorg
  other.application$kgha[built]     <-  built_inorg
  
  # FINAL PROCESSING
  in_sp <- springcrops.application# cbind(landcover[,c(1:2)],springcrops.application$kgha)
  in_wi <- wintercrops.application# cbind(landcover[,c(1:2)],wintercrops.application$kgha)
  in_ot <- other.application# cbind(landcover[,c(1:2)],other.application$kgha)
  in_gr <- grass.application# cbind(landcover[,c(1:2)],grass.application$kgha)
  
  up_sp <- springcrops.uptake# cbind(landcover[,c(1:2)],springcrops.uptake$kgha)
  up_wi <- wintercrops.uptake# cbind(landcover[,c(1:2)],wintercrops.uptake$kgha)
  up_ot <- other.uptake# cbind(landcover[,c(1:2)],other.uptake$kgha)
  up_gr <- grass.uptake# cbind(landcover[,c(1:2)],grass.uptake$kgha)
  
  res <- list(in_sp,in_wi,in_ot,in_gr,up_sp,up_wi,up_ot,up_gr)
  names(res) <- c("in_sp","in_wi","in_ot","in_gr","up_sp","up_wi","up_ot","up_gr")
  
  for (i in c(1:length(res))){
    #res[[i]][,3][is.na(landcover$sums)] <- NA
    res[[i]][,3][is.na(crafty.df$crafty)] <- NA
  }
  return(res)
  
  
}
# end of function 4A ----



# **************************************************************************** #
# FUNCTION 5: Calculate weighted averaged PET facts ----
# **************************************************************************** #
# Calculation of PET facts by weighted average of NIRAMS land covers

pet_facts_calc <- function(landcover, PET_facts){
  # Calculates weighted average PET/AET factor for each grid cell. 
  # This is done by multiplying the land cover areas within each grid cell with
  # their corresponding land-cover-specific PET factor, sum and then divide by 100 
  #
  # INPUT:
  # landcover     data.frame with the areas (ha) of each NIRAMS land use class  
  #               for each cell in the NIRAMS grid
  # PET_facts     data.frame with PET factors for each of the 28 NIRAMS land 
  #               use classes
  #
  # OUTPUT: 
  #

  PF_weighted <- as.data.frame(t(t(landcover[,c(3:28)])*(PET_facts[1:26])))
  PF_weighted[PF_weighted == 0] <- NA
  PF_weighted_mean <- cbind(landcover[,c(1:2)],round(rowSums(PF_weighted,na.rm = TRUE),2)/100)
  
  names(PF_weighted_mean)[3] <- "petfac"
  PF_weighted_mean$petfac[PF_weighted_mean$petfac == 0] <- NA
  
  return(PF_weighted_mean)
}

# end of function 5 ----


# **************************************************************************** #
# FUNCTION 5A: Calculate weighted averaged PET facts CRAFTY ----
# **************************************************************************** #
# Calculation of PET facts by weighted average of NIRAMS land covers

future_pet_facts_calc <- function(my.crafty, CRAFTY_LU){
  # Calculates weighted average PET/AET factor for each grid cell. 
  # This is done by multiplying the land cover areas within each grid cell with
  # their corresponding land-cover-specific PET factor, sum and then divide by 100 
  #
  # INPUT:
  # my.crafty     raster with future land use from CRAFTY
  # CRAFTY_LU     data.frame with PET factors for each of the 27 crafty land 
  #               use classes
  #
  # OUTPUT: 
  #
  
  # CRAFTY raster to data.frame
  crafty.df           <- as.data.frame(my.crafty, xy=TRUE)
  names(crafty.df)[3] <- "CRAFTY_INDEX"
  
  crafty_pet <- CRAFTY_LU[,c(1,5)]  # get PET_facts for CRAFTY land uses
  
  PF_weighted = merge(crafty.df, crafty_pet, by = "CRAFTY_INDEX", all.x=TRUE)
  names(PF_weighted)[4] <- "petfac"
  PF_weighted <- PF_weighted[,c("x","y", "petfac")] 
  return(PF_weighted)
  
  # # ANOTHER OPTION:
  # # Because the CRAFTY data only gives the dominant land use class, this is just 
  # # a task of assigning uptake and fertiliser application rates to the right
  # # CRAFTY classes as defined in Inorg_uptake.
  # 
  # # ------------------------------------------------------- #
  # # 1: IDENTIFY THE CELLS WITH GRASS, ARABLE OR OTHER 
  # # ------------------------------------------------------- #
  # names(crafty.df)[3] <- "crafty"
  # # Identify cells with rough grazing. It is assumed that extensive pastoral
  # #  (class 3) in CRAFTY is rough grazing. SHOULD CRAFTY CLASS 14 ALSO BE ROUGH GRAZING?
  # roughgrazing <- (crafty.df$crafty == 3)  # identify cells with rough grazing
  # roughgrazing[is.na(roughgrazing)] <- FALSE
  # 
  # # Identify cells with improved grassland areas > 0 
  # improvedgrazing <- (crafty.df$crafty == 6)  
  # improvedgrazing[is.na(improvedgrazing)] <- FALSE
  # 
  # # Identify cells with arable area (CRAFTY class 1, 2, 4, 5 and 13)
  # arable <- (crafty.df$crafty == 1 | crafty.df$crafty == 2 | crafty.df$crafty == 4 | crafty.df$crafty == 5 | crafty.df$crafty == 13)   
  # arable[is.na(arable)] <- FALSE
  # 
  # # Identify cells with woodland area (CRAFTY class 0, 7-12)
  # woodland <- (crafty.df$crafty == 0 | crafty.df$crafty == 7 | crafty.df$crafty == 8 | crafty.df$crafty == 9 | 
  #                crafty.df$crafty == 10 | crafty.df$crafty == 11 | crafty.df$crafty == 12 )   
  # woodland[is.na(woodland)] <- FALSE
  # 
  # # Identify cells with "wildscape" area (CRAFTY class -1 and 14)
  # wild <- (crafty.df$crafty == -1 | crafty.df$crafty == 14)   
  # wild[is.na(wild)] <- FALSE
  # 
  # # Identify cells with built area (CRAFTY class 15)
  # built <- (crafty.df$crafty == 15 )   
  # built[is.na(built)] <- FALSE
  # 
  # # -------------------------------------------------------- #
  # # 2: EXTRACT PETFACT FOR THE LAND USES
  # # -------------------------------------------------------- #
  # # Note, all CRAFTY arable classes are currently given the same values 
  # roughgrass_pet <- CRAFTY_LU_CODE$PET_Fact[CRAFTY_LU_CODE$CRAFTY_INDEX == 3]
  # imprgrass_pet  <- CRAFTY_LU_CODE$PET_Fact[CRAFTY_LU_CODE$CRAFTY_INDEX == 6]
  # arable_pet     <- CRAFTY_LU_CODE$PET_Fact[CRAFTY_LU_CODE$CRAFTY_INDEX == 2]          
  # wood_pet       <- CRAFTY_LU_CODE$PET_Fact[CRAFTY_LU_CODE$CRAFTY_INDEX == 7]          
  # wild_pet       <- CRAFTY_LU_CODE$PET_Fact[CRAFTY_LU_CODE$CRAFTY_INDEX == 14]          
  # built_pet      <- CRAFTY_LU_CODE$PET_Fact[CRAFTY_LU_CODE$CRAFTY_INDEX == 15]          
  # 
  # # ---------------------------------------------------------------- #
  # # 3: ASSIGN UPTAKE AND FERTILISER RATES TO  LAND USE CLASSES
  # # ---------------------------------------------------------------- #
  # # Note: Assumed 70% of arable is spring crops and 30% is winter crops
  # # GRASS
  # PF.LU <-  crafty.df[,c(1,2)]
  # PF.LU$pet_fact <- 0
  # PF.LU$pet_fact[roughgrazing]    <-  roughgrass_pet
  # PF.LU$pet_fact[improvedgrazing] <-  imprgrass_pet
  # PF.LU$pet_fact[arable]          <-  arable_pet
  # PF.LU$pet_fact[woodland]        <-  wood_pet
  # PF.LU$pet_fact[wild]            <-  wild_pet
  # PF.LU$pet_fact[built]           <-  built_pet
  # PF.LU$pet_fact[PF.LU$pet_fact == 0] <- NA
  # PF_weighted <- PF.LU
  # return(PF_weighted)
  

}

# end of function 5A ----

# library(reshape2)
# 


