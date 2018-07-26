library(cffdrs)

### PARAMETERS

## FUEL TYPES PARAMTERS FROM LANDIS
maxcol <- 21 #max(count.fields(file.path(getPaths()$inputPath, "dynamic-biomass-fuels.txt"), sep = ""))
dynamicBiomassFuels <- prepInputs(targetFile = "dynamic-biomass-fuels.txt", 
                                   url = "https://raw.githubusercontent.com/CeresBarros/Extension-Dynamic-Biomass-Fuels/master/testings/version-tests/v6.0-2.0/dynamic-biomass-fuels.txt", ## change to extractURL later 
                                   destinationPath = getPaths()$inputPath,   ## change to module data/ later
                                   fun = "utils::read.table", 
                                   fill = TRUE, row.names = NULL,
                                   sep = "",
                                   header = FALSE,
                                   blank.lines.skip = TRUE,
                                   col.names = c(paste("col",1:maxcol, sep = "")),
                                   stringsAsFactors = FALSE)
dynamicBiomassFuels <- data.table(dynamicBiomassFuels[, 1:15])
dynamicBiomassFuels <- dynamicBiomassFuels[!(col1 %in% c("LandisData", "Timestep", "MapFileNames",
                                                     "PctConiferFileName", "PctDeadFirFileName"))]
dynamicBiomassFuels <- dynamicBiomassFuels[!(col2 %in% c("The", "Users", "Optional"))]
dynamicBiomassFuels <- dynamicBiomassFuels[!(col1 == ">>" & grepl("---", col2))]
dynamicBiomassFuels[col1 == ">>"] <- data.table(dynamicBiomassFuels[col1 == ">>", col2:col15], col15 = "NA")
dynamicBiomassFuels[,col15:=NULL]
dynamicBiomassFuels[col1 == "Fuel" & col2 == "Type"] <- data.table(dynamicBiomassFuels[col1 == "Fuel" & col2 == "Type", col2:col14], 
                                                                               col14 = "NA")
dynamicBiomassFuels[col1 == "Type", col1 := "FuelType"] 



## SPECIES NAMES CORREPONDENCES ------------------------------------
## get species codes from dynamic fuels example table
# speciesNames <- data.table(LANDISNames = dynamicBiomassFuels[(which(col1=="Fuel") + 2) : (which(col1 == "HardwoodMaximum") - 1),
# speciesNames[, LANDISNames1:=as.character(substring(LANDISNames, 1, 4))]
# speciesNames[, LANDISNames2:=as.character(substring(LANDISNames, 5, 7))]
# 
# speciesNames[, ':='(LandRNames = paste0(toupper(substring(LANDISNames1, 1, 1)),
#                                         tolower(as.character(substring(LANDISNames1, 2, 4))),
#                                         "_", tolower(substring(LANDISNames2, 1, 1)),
#                                         tolower(as.character(substring(LANDISNames2, 2, 3)))),
#                     LANDISNames1 = NULL,
#                     LANDISNames2 = NULL)]
                                                             # col1])
# speciesNames[grepl("spp", LANDISNames2), LANDISNames2:="sp"]

## OR get species codes speciestable (only after running Land-R)
speciesNames <- data.table(LANDISNames = unique(LBMR_testSim@.envir$speciesTable$LandisCode))
speciesNames[, LANDISNames := tolower(gsub("\\.", "", LANDISNames))]
speciesNames[, LANDISNames1:=as.character(substring(LANDISNames, 1, 4))]
speciesNames[, LANDISNames2:=as.character(substring(LANDISNames, 5, 7))]
speciesNames[LANDISNames2 == "spp" | LANDISNames2 == "all", LANDISNames2:="sp"]


speciesNames[, ':='(LandRNames = paste0(toupper(substring(LANDISNames1, 1, 1)),
                                        tolower(as.character(substring(LANDISNames1, 2, 4))),
                                        "_", tolower(substring(LANDISNames2, 1, 1)),
                                        tolower(as.character(substring(LANDISNames2, 2, 3)))),
                    LANDISNames1 = NULL,
                    LANDISNames2 = NULL)]

if(any(!tolower(speciesList[,1]) %in% tolower(speciesNames$LandRNames)))
  warning(cat("\nFollowing selected species not found in the LANDIS species list.\nCheck if this is correct:\n",
                paste0(speciesList[!tolower(speciesList[,1]) %in%
                                     tolower(speciesNames$LandRNames), 1],
                collapse = ", ")))

## Convert species names for selected species, creating a 'species' column that matches other tables
## append species codes only after init of LBMR
tempList <- speciesList
rownames(tempList) <- tolower(tempList[,1])
commonSpp <- tolower(speciesNames$LandRNames[tolower(LandRNames) %in% rownames(tempList)])

speciesNames[tolower(LandRNames) %in% rownames(tempList), species := tempList[commonSpp,2]]

speciesNames <- merge(speciesNames, LBMR_testSimout@.envir$species[, c("species", "speciesCode"), with = FALSE],
                      by = "species", all = TRUE)
#sim$speciesNames

## SPECIES COEFFICIENTS PARAM TABLE
sppMultipliers <- dynamicBiomassFuels[(which(col1=="Fuel") + 1) : (which(col1 == "HardwoodMaximum") - 1),
                                      col1:col2]

names(sppMultipliers) <- as.character(sppMultipliers[1,])
sppMultipliers <- sppMultipliers[-1]

## remove last character to match other tables
sppMultipliers[, Species := substring(Species, 1, 7)]

## merge species names and codes
## keep all species, even thoe that do not match between fuels inputs and LandR tables
sppMultipliers <- merge(speciesNames, sppMultipliers, by.x = "LANDISNames", by.y="Species", all = TRUE)
sppMultipliers[, LANDISNames := NULL]
sppMultipliers$Coefficient <- as.numeric(sppMultipliers$Coefficient)

## LANDIS advises keeping one for all species or values close to 1
sppMultipliers[is.na(Coefficient), Coefficient := 1.0]

## exclude lines that have no spp codes
sppMultipliers <- sppMultipliers[, sum(is.na(.SD)) < 3, by = 1:nrow(sppMultipliers), .SDcols = 1:3] %>%
  .$V1 %>%
  sppMultipliers[.]

# sim$sppMultipliers <- sppMultipliers

## FUEL TYPES PARAM TABLE
## if a fuel types paramter table exists in inputs, don't use LANDIS example
if(!file.exists(file.path(getPaths()$inputPath, "FuelTypeCorresp.csv"))) {
  FuelTypes <- dynamicBiomassFuels[(which(col1=="FuelTypes") + 1) : (which(col1 == ">>EcoregionsTable") - 1),
                                   col1:col14]
  ## rename columns
  FuelTypes[1, `:=`(col3 = "minAge",
                    col4 = "NA",
                    col5 = "maxAge",
                    col6 = "Species")]
  FuelTypes[, col4 := NULL]
  names(FuelTypes) <- c(as.character(FuelTypes[1, 1:5, with = FALSE]), paste0("Species", 2:9))
  FuelTypes <- FuelTypes[-1]
  
  ## melt table to get all species in one column
  FuelTypes <- melt(FuelTypes, id.vars = 1:4, variable.name = "Species") %>%
    .[value != ""] %>%
    .[, Species := value] %>%
    .[, value := NULL]
  
  ## create a column with negative switch
  ## species with a negative switch negatively contribute to a fuel type
  FuelTypes[grepl("-", Species), negSwitch := -1L]
  FuelTypes[is.na(negSwitch), negSwitch := 1]
  FuelTypes[, Species := sub("-", "", Species)]
  
  ## remove last character to match other tables
  FuelTypes[, Species := substring(Species, 1, 7)]
  
} else {
  FuelTypes <- prepInputs(targetFile = "FuelTypeCorresp.csv",
                          destinationPath = getPaths()$inputPath,
                          fun = "utils::read.csv",
                          header = TRUE) %>%
    data.table
  
  ## remove last character to match other tables
  FuelTypes[, Species := substring(Species, 1, 7)]
}

## merge species names and codes, and only keep those that match other tables
## note that this WILL not drop unselected spcies  (which are important when detecting fuel type)
FuelTypes <- merge(FuelTypes, speciesNames, by.x = "Species", by.y = "LANDISNames", all = TRUE)
FuelTypes[, Species := NULL]

## convert some columns to numeric
numCols <- c("minAge", "maxAge", "negSwitch")
FuelTypes[, (numCols) := lapply(.SD, function(x) as.numeric(x)), .SDcols = numCols]

#sim$FuelTypes <- FuelTypes

## HARDWOOD CUTOFF
# hardwoodMax <- as.integer(dynamicBiomassFuels[col1 == "HardwoodMaximum", col2])
hardwoodMax <- 0    ## manuel recommendation
# sim$hardwoodMax <- hardwoodMax


## CALCULATE SPP VALUES FOR EACH FUEL TYPE IN EACH PIXEL ------------------------------
## (spp influence towards fuel type classif; the fuel type with the largest value "wins") 
## subset cohortData according to species age range for fuel types
## for each sp in each pixel group/ecoregion, sum biomass across ages

## subset cohort data and non-na fuel types
sppValues <- copy(LBMR_testSim@.envir$cohortData[, pixelGroup:B])
tempFT <- na.omit(copy(FuelTypes[, -c("FuelTypeDesc", "species", "LandRNames"), with = FALSE]))  ## keep only complete lines with spp codes

## merge the two tables and add sppMultipliers
sppValues <- sppValues[tempFT, on = .(speciesCode), allow.cartesian = TRUE, nomatch = 0] %>%
  .[order(pixelGroup)]
sppValues <- sppMultipliers[!is.na(speciesCode), speciesCode:Coefficient] %>%
  .[!duplicated(.)] %>%
  sppValues[., on = .(speciesCode), nomatch =0] %>%
  .[order(pixelGroup)]

## calculate sppValues per pixelGroup, ecoregion, and fuel type
## species biomass is weighted by the coeff and becomes 
## negative if the species has a negative contribution (negSwitch) to that fuel type
sppValues <- sppValues %>%
  group_by(pixelGroup, ecoregionGroup, FuelType, BaseFuel) %>%
  filter(age>= minAge & age <= maxAge) %>%
  summarise(fuelTypeVal = sum(B*Coefficient*negSwitch)) %>%
  data.table()

## ASSIGN FINAL FUEL TYPES
## pixelGroups get the fuel type that corresponds to maxValue
## two fuel types may share maxValue, assignment priority is:
## Conifer > Deciduous > ConiferPlantation > Open > Slash

## get max spp value in each pixelGroup
cols <- c("pixelGroup", "ecoregionGroup")
sppValues[, maxValue := max(fuelTypeVal),
          by = cols]

## TODO move this function to sourced script
calcFinalFuels <- function(BaseFuel, FuelType, fuelTypeVal, maxValue) {
  # browser()
  finalFuelType  <- FuelType[which(fuelTypeVal == unique(maxValue))]
  finalBaseFuel <- BaseFuel[which(fuelTypeVal == unique(maxValue))]
  
  ## in case more that one fuel type has maxValue:
  if(length(finalFuelType) > 1) {
    if(any(grepl("Conifer$", finalBaseFuel))){
      finalFuelType <- finalFuelType[finalBaseFuel == "Conifer"]
      finalBaseFuel <- "Conifer"
    } else if(any(grepl("Deciduous", finalBaseFuel))) {
      finalFuelType <- finalFuelType[finalBaseFuel == "Deciduous"]
      finalBaseFuel <- "Deciduous"
    } else if(any(grepl("Plantation", finalBaseFuel))) {
      finalFuelType <- finalFuelType[finalBaseFuel == "ConiferPlantation"]
      finalBaseFuel <- "ConiferPlantation" 
    } else if(any(grepl("Open", finalBaseFuel))) {
      finalFuelType <- finalFuelType[finalBaseFuel == "Open"]
      finalBaseFuel <- "Open"  
    } else if(any(grepl("Slash", finalBaseFuel))) {
      finalFuelType <- finalFuelType[finalBaseFuel == "Slash"]
      finalBaseFuel <- "Slash"
    }
  }
  
  list(finalFuelType = as.integer(finalFuelType), finalBaseFuel = as.character(finalBaseFuel))
}

sppValues[, finalFuelType := as.integer()]
sppValues[, finalBaseFuel := as.character()]
sppValues[, c("finalFuelType", "finalBaseFuel") := calcFinalFuels(BaseFuel, 
                                                                    FuelType,
                                                                    fuelTypeVal,
                                                                    maxValue),
          by = cols]

## CALCULATE CONIFEROUS/DECIDUOUS DOMINANCE ------------
## sum sppValues across conifer/deciduous BaseFuels 
## for each pixelGroup

coniferDom <- sppValues[grepl("Conifer", BaseFuel)] %>%
  group_by(pixelGroup, ecoregionGroup) %>%
  summarise(sumCon = sum(fuelTypeVal)) %>%
  data.table

deciduousDom <- sppValues[grepl("Deciduous", BaseFuel)] %>%
  group_by(pixelGroup, ecoregionGroup) %>%
  summarise(sumDec = sum(fuelTypeVal)) %>%
  data.table

## merge and clean WS
cols <- c("pixelGroup", "ecoregionGroup")
sppValues <- merge(coniferDom, deciduousDom, by = cols, all = TRUE) %>%
  merge(sppValues, ., by = cols)

rm(coniferDom, deciduousDom)

## ConiferPlantation, Open and Slash have their own rules 
## for conifer/deciduos dominance, so sumCon and sumDec need to be overriden 
cols <- c("sumCon", "sumDec")

sppValues[finalBaseFuel == "ConiferPlantation", 
          (cols) := list(100, 0)]
sppValues[finalBaseFuel == "Slash", 
          (cols) := list(0, 0)]
sppValues[finalBaseFuel == "Open", 
          (cols) := list(0, 0)]

## For Conifer and Deciduous, the conifer vs hardwood dominance
## are calculated for each as their percent dominance + 0.5
## NOTE: this is not = to LANDIS source code, but approaches the idea of the manual

calcDominance <- function(sumCon, sumDec, finalBaseFuel){
  ## get the initial fuel type attribute dominances
  finalBaseFuel <- finalBaseFuel
  if(finalBaseFuel == "ConiferPlantation") {
    coniferDom <- 100
    hardwoodDom <- 0
  } else {
    ## initial values will be kept for Open and Slash
    coniferDom <- 0
    hardwoodDom <- 0
  }
  
  
  if(sumCon > 0 | sumDec > 0) {
    coniferDom <- ceiling(sumCon/(sumCon + sumDec) * 100)
    hardwoodDom <- ceiling(sumDec/(sumCon + sumDec) * 100)
    
    if(coniferDom > hardwoodMax & hardwoodDom > hardwoodMax)
      finalBaseFuel <- "Mixed"
    
    if(coniferDom <= hardwoodMax & hardwoodDom <= hardwoodMax) {
      finalBaseFuel <- if(coniferDom > 0 & hardwoodDom > 0)
        "Mixed" else if(coniferDom > 0) 
          "Conifer" else
            "Deciduous"
    }
    
    if(coniferDom <= hardwoodMax & hardwoodDom > hardwoodMax) {
      finalBaseFuel <- "Deciduous"
      coniferDom <- 0
      hardwoodDom <-  100
    }
    
    if(coniferDom > hardwoodMax & hardwoodDom <= hardwoodMax) {
      finalBaseFuel <- "Conifer"
      coniferDom <- 100
      hardwoodDom <-  0
    }
  }
  
  list(coniferDom = coniferDom, hardwoodDom = hardwoodDom, finalBaseFuel = finalBaseFuel)
}

## replace NAs by zeros
sppValues[is.na(sumCon), sumCon := 0]
sppValues[is.na(sumDec), sumDec := 0]

sppValues[, c("coniferDom", "hardwoodDom", "finalBaseFuel") := calcDominance(sumCon, sumDec, finalBaseFuel),
     by = 1:nrow(sppValues)]

## HERE ##
## TODO test integration with FBP

## dataframe of minimum FBP inputs set to their defaults. See ?fbp for info
FBPinputs <- data.frame(id = NA,
                        FuelType = "C2",
                        LAT = 55,
                        LONG = -120, 
                        FFMC = 90,
                        BUI = 60,
                        WS = 10,
                        GS = 0,
                        Dj = 180,
                        Aspect = 0)
