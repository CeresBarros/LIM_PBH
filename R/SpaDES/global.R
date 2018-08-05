## ------------------------------------------------------
## FIRE MODELLING WITH SpaDES
## TESTS
##
## Ceres: Nov 2017
## ------------------------------------------------------

## clean workspace
rm(list=ls()); gc(reset = TRUE)

## requires as of july 23rd 2018
# loading reproducible     0.2.2
# loading quickPlot        0.1.4
# loading SpaDES.core      0.2.1
# loading SpaDES.tools     0.3.0
# loading SpaDES.addins    0.1.1
# devtools::install_github("CeresBarros/reproducible@development") ## July 29th
# devtools::install_github("PredictiveEcology/quickPlot@development")
# devtools::install_github("PredictiveEcology/SpaDES.core@development")
# devtools::install_github("Predictive/SpaDES.tools@development")
library(SpaDES)

## testing packages
# try(detach("package:SpaDES.core", unload = TRUE)); try(detach("package:SpaDES.tools", unload = TRUE)); try(detach("package:reproducible", unload = TRUE)); devtools::load_all("E:/GitHub/reproducible"); devtools::load_all("E:/GitHub/SpaDES.tools"); devtools::load_all("E:/GitHub/SpaDES.core")

source("R/R_tools/Useful_functions.R")

## define paths
setPaths(cachePath = file.path("R/SpaDES/cache"),
         modulePath = file.path("R/SpaDES/m"),
         inputPath = file.path("R/SpaDES/inputs"),
         outputPath = file.path("R/SpaDES/outputs"))


## STUDY AREA(S) ---------------------------------------

## Foothills and a smaller region for testing
foothills <- raster::shapefile("data/maps/Foothills_study_area.shp")
foothillsSMALL <- raster::buffer(foothills, width = -0.3)

## Betu_pap was causing spades to hang on writing CASFRI raster to disk
speciesList <- as.matrix(read.csv(file.path(getPaths()$inputPath, "speciesList.csv"), header = TRUE))

## SIMULATION SETUP ------------------------------

## simulation parameters
pathsSim <- getPaths()
timesSim <- list(start = 1, end = 50)

# modulesSimtoy <- list("simplifyLCCVeg", "simpleLCCSuccession", "fireSpread", "fireSeverity")
# objectsSimtoy <-  list(studyArea = foothillsSMALL)
# paramsSimtoy <- list(
#   .globals = list(.useCache = TRUE),
#   fireSpread = list(fireSize = 1000, noStartPix = 100, fireFreq = 3),
#   fireSeverity = list(.plotMaps = TRUE, .plotInterval = 1),
#   fireStats = list(.plotStats = TRUE)
# )


# eventCaching <- c(".inputObjects")
modulesSimLBMR <- list("BiomassSpeciesData", "Boreal_LBMRDataPrep",
                       "LandR_BiomassFuels", "fireSpread", "LBMR", "fireSeverity")
                       
objectsSimLBMR <- list("shpStudyRegionFull" = foothillsSMALL,
                       "shpStudySubRegion" = foothillsSMALL,
                       "speciesList" = speciesList)
outputs <- data.frame(expand.grid(objectName = c("cohortData","rstCurrentBurn"),
                                  saveTime = seq(2, 50, by = 5),
                                  stringsAsFactors = FALSE))


paramsSimLBMR <- list(
  # Boreal_LBMRDataPrep = list(crsUsed = projCRS), 
  # BiomassSpeciesData = list(crsUsed = projCRS),
  LBMR = list(successionTimestep = 1,
              fireTimestep = 5L,
              seedingAlgorithm = "wardDispersal",
              .saveInitialTime = 1,
              .plotInitialTime = timesSim$start#, 
              # .useCache = eventCaching
              ),
  fireSpread = list(fireTimestep = 5L,
                    vegFeedback = FALSE),
  fireSeverity = list(fireTimestep = 5L,
                      .plotMaps = FALSE,
                      .saveInitialTime = 1)
)

pathsSim$outputPath <- "R/SpaDES/outputs/vegFB_0"
# pathsSim$outputPath <- file.path(pathsSim$outputPath, "vegFB_1")

showCache(getPaths()$cachePath)
# clearCache(getPaths()$cachePath#, userTags = "Boreal_LBMRDataPrep"
# )

LBMR_testSim <- simInit(times = timesSim, params = paramsSimLBMR, modules = modulesSimLBMR,
                        objects = objectsSimLBMR, outputs = outputs, paths = pathsSim)

graphics.off()
dev()
clearPlot()
LBMR_testSimout <- spades(LBMR_testSim, cache = TRUE, debug = TRUE)   ## debug = TRUE activates automatic browsing when errors occur
# events(LBMR_testSim)


paramsSimLBMR <- list(
  # Boreal_LBMRDataPrep = list(crsUsed = projCRS), 
  # BiomassSpeciesData = list(crsUsed = projCRS),
  LBMR = list(successionTimestep = 1,
              fireTimestep = 5L,
              seedingAlgorithm = "wardDispersal",
              .saveInitialTime = 1,
              .plotInitialTime = timesSim$start#, 
              # .useCache = eventCaching
  ),
  fireSpread = list(fireTimestep = 5L,
                    vegFeedback = TRUE),
  fireSeverity = list(fireTimestep = 5L,
                      .plotMaps = FALSE,
                      .saveInitialTime = 1)
)

pathsSim$outputPath <- "R/SpaDES/outputs/vegFB_1"
# pathsSim$outputPath <- file.path(pathsSim$outputPath, "vegFB_1")

showCache(getPaths()$cachePath)
# clearCache(getPaths()$cachePath#, userTags = "Boreal_LBMRDataPrep"
# )

LBMR_testSim2 <- simInit(times = timesSim, params = paramsSimLBMR, modules = modulesSimLBMR,
                        objects = objectsSimLBMR, paths = pathsSim)

graphics.off()
dev()
clearPlot()
LBMR_testSimout2 <- spades(LBMR_testSim2, cache = TRUE, debug = TRUE)   ## debug = TRUE activates automatic browsing when errors occur
# events(LBMR_testSim)







# expParams <- list(
#   # Boreal_LBMRDataPrep = list(crsUsed = projCRS),
#   # BiomassSpeciesData = list(crsUsed = projCRS),
#   LBMR = list(successionTimestep = 1,
#               fireTimestep = 5L,
#               seedingAlgorithm = "wardDispersal",
#               .saveInitialTime = 1,
#               .plotInitialTime = timesSim$start#,
#               # .useCache = eventCaching
#   ),
#   fireSpread = list(fireTimestep = 5L,
#                     vegFeedback = c(FALSE,TRUE)),
#   fireSeverity = list(fireTimestep = 5L,
#                       .plotMaps = FALSE,
#                       .saveInitialTime = 1)
# )
# 
# graphics.off()
# dev()
# clearPlot()
# LBMR_experiment <- experiment(LBMR_testSim, replicates = 5, params = expParams,
#                               experimentFile = TRUE)

## error when veg feeback are true
# endCluster()

# moduleDiagram(LBMR_testSim)
# objectDiagram(LBMR_testSim)
# events(LBMR_testSim)


