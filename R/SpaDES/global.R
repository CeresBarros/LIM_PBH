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






