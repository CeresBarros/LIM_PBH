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
# devtools::install_github("PredictiveEcology/SpaDES.tools@development")    ## Feb2018

library(SpaDES)
source("R/R_tools/Useful_functions.R")

## define paths
setPaths(cachePath = file.path("R/SpaDES/cache"),
         modulePath = file.path("R/SpaDES/m"),
         inputPath = file.path("R/SpaDES/inputs"),
         outputPath = file.path("R/SpaDES/outputs"))


## STUDY AREA(S) ---------------------------------------

## SUB REGION
## Foothills and a smaller region for testing
# foothills <- raster::shapefile(file.path(getPaths()$inputPath, "Alberta_study_area/Alberta_study_area"))
# foothillsSMALL <- rgeos::gBuffer(foothills, width = -0.3)

## Dave's AB and SK fire datasets - buffer around the polygons will be the study area
## this should be cached
firesABSK <- Cache(bindSpatialObjs, 
                   files = c("albertafires1_postfire", "albertafires2_postfire", "saskatchewanfires_postfire"),
                   folder = "data/fires_Dave/Projected_renamed", 
                   cacheRepo = getPaths()$cachePath)

firesABSK.buf <- Cache(outerBuffer, 
                       x = firesABSK, 
                       cacheRepo = getPaths()$cachePath)

crs(firesABSK.buf) = crs(firesABSK)
plot(firesABSK) ; plot(firesABSK.buf, add =TRUE)

## FULL REGION - From LandWeb 
source("R/R_tools/inputMaps.R")   ## functions for input maps
studyArea = firesABSK.buf

shpStudyRegions <- Cache(loadStudyRegion,
                         shpPath = asPath(file.path("../LandWeb/inputs", "shpLandWEB.shp")),
                         studyArea = studyArea,
                         crsKNNMaps = crs(studyArea), cacheRepo = getPaths()$cachePath)
list2env(shpStudyRegions, envir = environment())


## simulation parameters
# modules <- list("simplifyLCCVeg", "fireSpread", "fireSeverity", "simpleLCCSuccession")
modules <- list("Boreal_LBMRDataPrep", "LBMR")

times <- list(start = 1.0, end = 50, timeunit = "year")

parameters <- list(
  # .globals = list(.useCache = TRUE),
  # fireSpread = list(fireSize = 1000, noStartPix = 100),
  # fireSeverity = list(.plotMaps = TRUE),
  # fireStats = list(.plotStats = TRUE)
  .plotInitialTime = times$start,
  .saveInitialTime = times$start
)

objects <- list("shpStudyRegionFull" = shpStudyRegionFull,
                "shpStudySubRegion" = shpStudyRegion,
                "successionTimestep" = 10,
                "useParallel" = FALSE)

showCache(getPaths()$cachePath)
# clearCache(getPaths$cachePath, userTags = "LBMR")
mySim <- simInit(times = times, modules = modules,
                 objects = objects)
# moduleDiagram(mySim)
# objectDiagram(mySim)
# events(mySim)

dev()
clearPlot()
mySim2 <- spades(mySim, debug = TRUE)





## NEW LBMR VERSION TO INTEGRATE WITH CODE ABOVE:

## as of April 10th
# devtools::install_github("PredictiveEcology/reproducible@development")
# devtools::install_github("PredictiveEcology/quickPlot@development")
# devtools::install_github("PredictiveEcology/SpaDES.core@development")
# devtools::install_github("PredictiveEcology/SpaDES.tools@development")
# devtools::install_github("PredictiveEcology/SpaDES.tools@prepInputs")
rm(list = ls()); gc(reset = TRUE)

library(SpaDES)

setwd("C:/Ceres/GitHub/LandscapesInMotion/")
setPaths(cachePath = file.path("R/SpaDES/cache"),
         modulePath = file.path("R/SpaDES/m"),
         inputPath = file.path("R/SpaDES/inputs"),
         outputPath = file.path("R/SpaDES/outputs"))
paths <- getPaths()
subStudyRegionName <- "RIA"

## FULL REGION - From LandWeb 
source("R/R_tools/inputMaps.R")   ## functions for input maps

# These are used in inputTables.R for filling the tables of parameters in
landisInputs <- readRDS("../LandWeb/inputs/landisInputs.rds")
spEcoReg <- readRDS("../LandWeb/inputs/SpEcoReg.rds")

# The CRS for the Study -- spTransform converts this first one to the second one, they are identical geographically
# crsStudyRegion <- CRS(paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0",
#                         "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
crsStudyRegion <- sp::CRS(paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0",
                                "+ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

studyRegionFilePath <- {
  studyRegionFilename <- if ("RIA" %in% subStudyRegionName) {
    "RIA_SE_ResourceDistricts_Clip.shp"
  } else {
    "studyarea-correct.shp"
  }
  file.path(paths$inputPath, studyRegionFilename)
}

studyRegionsShps <- Cache(loadStudyRegions, shpStudyRegionCreateFn = shpStudyRegionCreate,
                          asPath(studyRegionFilePath),
                          fireReturnIntervalMap = asPath(file.path(paths$inputPath, "ltfcmap correct.shp")),
                          subStudyRegionName = subStudyRegionName,
                          crsStudyRegion = crsStudyRegion, cacheRepo = paths$cachePath)
list2env(studyRegionsShps, envir = environment()) # shpStudyRegion & shpStudyRegion

# Time steps
fireTimestep <- 1
successionTimestep <- 10 # was 2

## Simulation setup
eventCaching <- c(".inputObjects", "init")
timesSim <- list(start = 0, end = 100)
modulesSim <- list("Boreal_LBMRDataPrep", "LBMR")
objectsSim <- list("shpStudyRegionFull" = shpStudyRegion,
                   "shpStudySubRegion" = shpSubStudyRegion,
                   "useParallel" = 2)
paramsSim <- list(Boreal_LBMRDataPrep = list(.useCache = FALSE),
                  LBMR = list(successionTimestep = successionTimestep,
                              .plotInitialTime = timesSim$start,
                              .saveInitialTime = NA,
                              .useCache = FALSE))

LBMR_testSim <- simInit(times = timesSim, params = paramsSim, modules = modulesSim,
                        objects = objectsSim, paths = paths)

moduleDiagram(LBMR_testSim)
objectDiagram(LBMR_testSim)

LBMR_testSim <- spades(LBMR_testSim, cache = FALSE, debug = FALSE)







