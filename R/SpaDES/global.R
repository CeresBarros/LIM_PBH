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
# devtools::install_github("CeresBarros/reproducible@development", force = TRUE)
# devtools::install_github("PredictiveEcology/quickPlot@development")
# devtools::install_github("PredictiveEcology/SpaDES.core@development")  ## August 15th
# devtools::install_github("PredictiveEcology/SpaDES.tools@devel-amc") ## August 27th
library(SpaDES)

## testing packages
# try(detach("package:SpaDES.core", unload = TRUE))
# try(detach("package:SpaDES.tools", unload = TRUE))
# try(detach("package:reproducible", unload = TRUE))
# devtools::load_all("E:/GitHub/reproducible")
# devtools::load_all("E:/GitHub/SpaDES.tools")
# devtools::load_all("E:/GitHub/SpaDES.core")

source("R/R_tools/Useful_functions.R")

## define paths
setPaths(modulePath = file.path("R/SpaDES/m"),
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
timesSim <- list(start = 1, end = 3)

# eventCaching <- c(".inputObjects")
modulesSim <- list("BiomassSpeciesData", "Boreal_LBMRDataPrep",   ## biomassSpeciesData needs a data prep -can't cope with LBMR defaults
                   "LBMR", "LandR_BiomassGMOrig"#,
                   # "LandR_BiomassRegen", "LandR_BiomassFuels", "fireSpread", 
                   #"fireSeverity"
                   )

objectsSim <- list("shpStudyRegionFull" = foothillsSMALL,
                   "shpStudySubRegion" = foothillsSMALL#,
                   # "speciesList" = speciesList
                   )

outputs <- data.frame(expand.grid(objectName = c("cohortData"),
                                  saveTime = seq(2, 50, by = 5),
                                  stringsAsFactors = FALSE))
outputs <- rbind(outputs, data.frame(objectName = "rstCurrentBurn", 
                                     saveTime = tail(seq(2, 50, by = 5), 1)))

paramsSim <- list(
  LBMR = list(successionTimestep = 1,
              seedingAlgorithm = "wardDispersal",
              .saveInitialTime = 1,
              .plotInitialTime = timesSim$start#, 
              # .purge = 7
              # .useCache = eventCaching
              )#,
  # LandR_BiomassRegen = list(fireTimestep = 2L)#,
  # fireSpread = list(fireTimestep = 5L,
  #                   vegFeedback = TRUE)#,
  # fireSeverity = list(fireTimestep = 5L,
  #                     .plotMaps = FALSE,
  #                     .saveInitialTime = 1)
)

# pathsSim$outputPath <- "R/SpaDES/outputs/vegFB_0"
pathsSim$outputPath <- file.path(pathsSim$outputPath, "vegFB_1/tests")
pathsSim$cachePath <- file.path("R/SpaDES/cache/LIM_tests")

showCache(pathsSim$cachePath, userTags = "simList")
# reproducible::clearCache(pathsSim$cachePath, userTags = c("simList"))
LBMR_testSim <- simInit(times = timesSim, params = paramsSim, modules = modulesSim,
                        objects = objectsSim, outputs = outputs, paths = pathsSim)

graphics.off()
dev()
clearPlot()
LBMR_testSimout <- spades(LBMR_testSim, cache = TRUE, debug = TRUE)   ## debug = TRUE activates automatic browsing when errors occur

# events(LBMR_testSimout)
# year <- 5L
# LBMR_testSimout <- scheduleEvent(LBMR_testSimout, eventTime = year, eventType = "save", moduleName = "LBMR")
# events(LBMR_testSimout)
# 
# year <- 5.0
# LBMR_testSimout <- scheduleEvent(LBMR_testSimout, eventTime = year, eventType = "save", moduleName = "LBMR")
# events(LBMR_testSimout)
# 
# year <- time(LBMR_testSimout) + 2
# LBMR_testSimout <- scheduleEvent(LBMR_testSimout, eventTime = year, eventType = "save", moduleName = "LBMR")
# events(LBMR_testSimout)


## TEST WITH FAKE FIRE MAP
## make fake fire map
rstCurrentBurn <- LBMR_testSimout@.envir$pixelGroupMap
rstCurrentBurn[rstCurrentBurn[]>0] <- 1
rstCurrentBurn[rstCurrentBurn[] <= 0] <- NA
IDs <- which(rstCurrentBurn[] == 1)
rstCurrentBurn[IDs[1:round(length(IDs)/2)]] <- NA

objectsSim["rstCurrentBurn"] <- rstCurrentBurn

modulesSim <- list("BiomassSpeciesData", "Boreal_LBMRDataPrep",   ## biomassSpeciesData needs a data prep -can't cope with LBMR defaults
                   "LBMR", "LandR_BiomassGMOrig",
                   "LandR_BiomassFuels",
                   "LandR_BiomassRegen", "fireSpread",
                   "fireSeverity")

LBMR_testSim <- simInit(times = timesSim, params = paramsSim, modules = modulesSim,
                        objects = objectsSim, outputs = outputs, paths = pathsSim)

graphics.off()
dev()
clearPlot()
LBMR_testSimout <- spades(LBMR_testSim, cache = TRUE, debug = TRUE)   ## debug = TRUE activates automatic browsing when errors occur
completed(LBMR_testSimout)


## test species layers in landweb
# LandWebSA <- raster::shapefile("../LandWeb/inputs/studyarea-correct")
# 
# modulesSim <- list("BiomassSpeciesData")
# pathsSim$cachePath <- "R/SpaDES/cache/tempLandWeb"
# paramsSim <- list(BiomassSpeciesData = list(.crsUsed = as.character(raster::crs(LandWebSA))))
# objectsSim <- list("shpStudyRegionFull" = LandWebSA,
#                    "shpStudySubRegion" = LandWebSA)
# 
# # clearCache(pathsSim$cachePath, userTags = "Pickell")
# LandWebSpp <- simInit(times = timesSim, params = paramsSim, modules = modulesSim,
#                         objects = objectsSim, paths = pathsSim)
# 
# LandWebSppout <- spades(LandWebSpp, cache = TRUE, debug = TRUE) 
# 
