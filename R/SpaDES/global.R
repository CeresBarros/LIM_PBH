## ------------------------------------------------------
## FIRE MODELLING WITH SpaDES
## TESTS
##
## Ceres: Nov 2017
## ------------------------------------------------------

## clean workspace
rm(list=ls())

## requires
# devtools::install_github("PredictiveEcology/SpaDES.core@development")    ## Feb2018
# devtools::install_github("PredictiveEcology/reproducible@development")   ## Feb2018
# devtools::install_github("PredictiveEcology/SpaDES.tools@development")    ## Feb2018

library(SpaDES)

## define paths
setPaths(cachePath = file.path("R/SpaDES/cache"),
         modulePath = file.path("R/SpaDES/m"),
         inputPath = file.path("R/SpaDES/inputs"),
         outputPath = file.path("R/SpaDES/outputs"))


## Get study area, and make a smaller region
foothills <- raster::shapefile(file.path(getPaths()$inputPath, "Alberta_study_area/Alberta_study_area"))
foothillsSMALL <- rgeos::gBuffer(foothills, width = -0.3)

## simulation parameters
# modules <- list("simplifyLCCVeg", "fireSpread", "fireSeverity", "simpleLCCSuccession")
modules <- list("LBMRDataPrep", "LBMR")

times <- list(start = 1.0, end = 30, timeunit = "year")

# parameters <- list(
#   .globals = list(.useCache = TRUE),
#   fireSpread = list(fireSize = 1000, noStartPix = 100),
#   fireSeverity = list(.plotMaps = TRUE),
#   fireStats = list(.plotStats = TRUE)
# )

objects <- list("studyArea" = foothillsSMALL)

mySim <- simInit(times = times, modules = modules,
                 objects = objects)
# moduleDiagram(mySim)
# objectDiagram(mySim)
# events(mySim)

dev()
clearPlot()
mySim <- spades(mySim, debug = TRUE)




