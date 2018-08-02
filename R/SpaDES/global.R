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
timesSim <- list(start = 1, end = 100)

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
## projection for FBP system compatibility -  not tolerated for dispersal leave default parameter for now
# projCRS = "+proj=longlat +datum=WGS84"

paramsSimLBMR <- list(
  # Boreal_LBMRDataPrep = list(crsUsed = projCRS), 
  # BiomassSpeciesData = list(crsUsed = projCRS),
  LBMR = list(successionTimestep = 1,
              seedingAlgorithm = "wardDispersal",
              .saveInitialTime = NA, 
              .plotInitialTime = timesSim$start
              # seedingAlgorithm = "universalDispersal",
              # .useCache = eventCaching
              ),
  fireSpread = list(fireFreq = 2L,
                    vegFeedback = FALSE),
  fireSeverity = list(.plotMaps = FALSE)
)


showCache(getPaths()$cachePath)
# clearCache(getPaths()$cachePath#, userTags = "Boreal_LBMRDataPrep"
           # )

LBMR_testSim <- simInit(times = timesSim, params = paramsSimLBMR, modules = modulesSimLBMR,
                        objects = objectsSimLBMR, paths = pathsSim)

expParams <- list(
  # Boreal_LBMRDataPrep = list(crsUsed = projCRS), 
  # BiomassSpeciesData = list(crsUsed = projCRS),
  LBMR = list(successionTimestep = 1,
              seedingAlgorithm = "wardDispersal",
              .saveInitialTime = NA, 
              .plotInitialTime = NA
              # seedingAlgorithm = "universalDispersal",
              # .useCache = eventCaching
  ),
  fireSpread = list(fireFreq = 2L,
                    vegFeedback = c(FALSE,TRUE)),
  fireSeverity = list(.plotMaps = FALSE)
)

graphics.off()
dev()
clearPlot()
LBMR_experiment <- experiment(LBMR_testSim, replicates = 10, params = expParams,
                              experimentFile = FALSE)
endCluster()

# moduleDiagram(LBMR_testSim)
# objectDiagram(LBMR_testSim)
# events(LBMR_testSim)

graphics.off()
dev()
clearPlot()
LBMR_testSimout <- spades(LBMR_testSim, cache = TRUE, debug = TRUE)   ## debug = TRUE activates automatic browsing when errors occur
# events(LBMR_testSim)


