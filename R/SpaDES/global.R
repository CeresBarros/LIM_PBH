## ------------------------------------------------------
## FIRE MODELLING WITH SpaDES
## TESTS
##
## Ceres: Nov 2017
## ------------------------------------------------------

## clean workspace
rm(list=ls()); gc(reset = TRUE)

## requires
# devtools::install_github("PredictiveEcology/SpaDES.core@development")    ## Feb2018 
# devtools::install_github("PredictiveEcology/reproducible@development")   ## Feb2018 
library(SpaDES)

## define paths
setPaths(cachePath = file.path("D:/Cherry/GitHub/LandscapesInMotion/R/SpaDES/cache"),
         modulePath = file.path("D:/Cherry/GitHub/LandscapesInMotion/R/SpaDES/m"),
         inputPath = file.path("D:/Cherry/GitHub/LandscapesInMotion/R/SpaDES/inputs"),
         outputPath = file.path("D:/Cherry/GitHub/LandscapesInMotion/R/SpaDES/outputs"))

## NOTE REVISE INPUTE OBJECTS TO ENSURE MODULARITY
## HERE
## check if external modules exist if not download them
if(!checkModuleLocal("cropReprojectLccAge", path = getPaths()$modulePath)){
  downloadModule("cropReprojectLccAge", repo = "PredictiveEcology/SpaDES-modules")
}
if(!checkModuleLocal("LccToBeaconsReclassify", path = getPaths()$modulePath)){
  downloadModule("LccToBeaconsReclassify", repo = "PredictiveEcology/SpaDES-modules")
}

## Get study area, and make a smaller region
foothills <- raster::shapefile("D:/Cherry/GitHub/LandscapesInMotion/R/SpaDES/inputs/Alberta_study_area/Alberta_study_area")
foothillsSMALL <- rgeos::gBuffer(foothills, width = -0.3)

## Get vegetation data - if checksums do not exist, then download manually and do checksums
if(!file.exists(file.path(moduleDir, "fire_spreadSTSM", "data", "CHECKSUMS.txt"))) {
  download.file(url = "http://www.cec.org/sites/default/files/Atlas/Files/Land_Cover_2010/Land_Cover_2010_TIFF.zip", 
                destfile = file.path(moduleDir, "fire_spreadSTSM", "data", "Land_Cover_2010_TIFF.zip"))
  checksums("fire_spreadSTSM", moduleDir, write = TRUE)
} else {
  downloadData(module = "fire_spreadSTSM", path = moduleDir)   
}

if(!file.exists(file.path(moduleDir, "fire_spreadSTSM", "data", "NA_LandCover_2010_25haMMU.tif"))){
  LCC2010 <- grep("NA_LandCover_2010_25haMMU.tif$", 
                  unzip(file.path(moduleDir, "fire_spreadSTSM", "data", "Land_Cover_2010_TIFF.zip"), list = TRUE)$Name,
                  value = TRUE)
  unzip(zipfile = file.path(moduleDir, "fire_spreadSTSM", "data", "Land_Cover_2010_TIFF.zip"),
        files = LCC2010, exdir = file.path(moduleDir, "fire_spreadSTSM", "data"), junkpaths = TRUE)

  LCC2010 <- raster::raster(file.path(moduleDir, "fire_spreadSTSM", "data", basename(LCC2010)))
} else {
  LCC2010 <- raster::raster(file.path(moduleDir, "fire_spreadSTSM", "data", "NA_LandCover_2010_25haMMU.tif"))
}


## simulation parameters

modules <- list("simplifyLCCVeg", "fire_spreadSTSM", "fireSeverity")

times <- list(start = 1.0, end = 5, timeunit = "year")

parameters <- list(
  .globals = list(.useCache = TRUE),
  fire_spreadSTSM = list(fireSize = 1000, noStartPix = 100),
  fireSeverity = list(.plotMaps = TRUE),
  fireStats = list(.plotStats = TRUE)
)

objects <- list("studyArea" = foothillsSMALL, "vegetationRas" = LCC2010)

mySim <- simInit(times = times, params = parameters, modules = modules,
                 objects = objects, 
                 loadOrder = c("simplifyLCCVeg", "fire_spreadSTSM", "fireSeverity"))

dev()
clearPlot()
mySim <- spades(mySim, debug = TRUE)




