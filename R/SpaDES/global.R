## ------------------------------------------------------
## FIRE MODELLING WITH SpaDES
## TESTS
##
## Ceres: Nov 2017
## ------------------------------------------------------

## clean workspace
rm(list=ls()); amc::.gc()

## Get packages ----------------------
## requires as of Jan 2nd 2020
# loading reproducible     0.2.11.9000
# loading quickPlot        0.1.6.9000
# loading SpaDES.core      0.2.8
# loading SpaDES.tools     0.3.4.9000
# loading SpaDES.addins    0.1.2
# loading LandR            0.0.3.9000

# devtools::install_github("PredictiveEcology/reproducible@development", upgrade = "always", type = "binary")
# devtools::install_github("achubaty/amc@development", upgrade = "always", type = "binary")
# devtools::install_github("PredictiveEcology/pemisc@development", upgrade = "always", type = "binary")
# devtools::install_github("PredictiveEcology/map@development", upgrade = "always", type = "binary")
# devtools::install_github("PredictiveEcology/LandR@development", upgrade = "always", type = "binary")
# devtools::install_github("PredictiveEcology/quickPlot@development", upgrade = "always", type = "binary")
# devtools::install_github("PredictiveEcology/SpaDES.tools@development", upgrade = "always", type = "binary")
# devtools::install_github("PredictiveEcology/SpaDES.core@development", upgrade = "always", type = "binary")
library(SpaDES)
library(LandR)

## testing packages
# try(detach("package:LandR", unload = TRUE))
# try(detach("package:pemisc", unload = TRUE))
# try(detach("package:amc", unload = TRUE))
# try(detach("package:map", unload = TRUE))
# try(detach("package:SpaDES", unload = TRUE))
# try(detach("package:SpaDES.addins", unload = TRUE))
# try(detach("package:pemisc", unload = TRUE))
# try(detach("package:amc", unload = TRUE))
# try(detach("package:map", unload = TRUE))
# try(detach("package:SpaDES.core", unload = TRUE))
# try(detach("package:SpaDES.tools", unload = TRUE))
# try(detach("package:pemisc", unload = TRUE))
# try(detach("package:amc", unload = TRUE))
# try(detach("package:map", unload = TRUE))
# try(detach("package:reproducible", unload = TRUE))
# devtools::load_all("../reproducible")
# devtools::load_all("../SpaDES.tools")
# library(SpaDES.addins)
# library(logging)
# devtools::load_all("../SpaDES.core")
# library(map)
# library(amc)
# library(pemisc)
# devtools::load_all("../LandR")

options("reproducible.useNewDigestAlgorithm" = TRUE)
options("spades.moduleCodeChecks" = FALSE)
options("reproducible.useCache" = TRUE)
options("reproducible.inputPaths" = normPath("R/SpaDES/inputs"))  ## store everything in data/ so that there are no duplicated files across modules
options("reproducible.destinationPath" = normPath("R/SpaDES/inputs"))

source("R/R_tools/Useful_functions.R")

## -----------------------------------------------
## SIMULATION SETUP
## -----------------------------------------------

## Set up simulation name  ---------------------------
# runName <- "blogSep2019_PM"
# runName <- "blogSep2019_noPM"

# runName <- "blogSep2019_PM_oneFire"
# runName <- "blogSep2019_noPM_oneFire"

runName <- "PM_oneFire_newSppParams"
# runName <- "noPM_oneFire_newSppParams"
eventCaching <- c(".inputObjects", "init")
useParallel <- FALSE

## paths define simulation paths
simPaths <- list(cachePath = file.path("R/SpaDES/cache/LIM_tests", runName),
                 modulePath = file.path("R/SpaDES/m"),
                 inputPath = file.path("R/SpaDES/inputs"),
                 outputPath = file.path("R/SpaDES/outputs", runName))

## Get necessary objects -----------------------
source("R/SpaDES/1_simObjects.R")

## Run Biomass_speciesData to get species layers
source("R/SpaDES/2_speciesLayers.R")

source("R/SpaDES/3_fireWeather.R")

## Define simulation params --------------------
simTimes <- list(start = 0, end = 65)
vegLeadingProportion <- 0 # indicates what proportion the stand must be in one species group for it to be leading.
# If all are below this, then it is a "mixed" stand
fireInitialTime <- 5L
fireTimestep <- if (grepl("oneFire", runName)) 100L else 2L
successionTimestep <- 1L

## Make simulation module list, parameters objects, objects and outputs accoding to run
## name and the parameters above
source("R/SpaDES/4_simulationSetup.R")


fireWeatherSim <- simInitAndSpades(times = simTimes
                                   , params = simParams
                                   # , modules = simModules[c(1:3, 5:6, 8)] ## for blog post
                                   , modules = simModules[c(1:5, 7)]
                                   , objects = simObjects
                                   , paths = simPaths
                                   , outputs = outputs
                                   , debug = TRUE
                                   # , .plotInitialTime = NA
)

## -----------------------------------------------
## SIMULATION RUN
## -----------------------------------------------

# showCache(simPaths$cachePath, after = "2018-09-26 00:00:00")
# reproducible::clearCache(simPaths$cachePath, userTags = c("prepInputsLCC2005_rtm", "Biomass_borealDataPrep"))

## TODO CHANGE FIRE MODULES TO USE COHORT DATA RATHER THAN SUMMARY BMG OUTPUTS, LIKE BIOMASSMAP
graphics.off()

# reproducible::clearCache(simPaths$cachePath, userTags = c("^Biomass_core$", "init"), ask = FALSE)
## TODO RUN SIMUALTIONS W/ AND W/O PM for blog
# set.seed(524326)

## TODO: implement LANDIS pixel fire severity calculation:
## Each fire event has an associated mean fire severity which is the average of the severities at all of the event’s sites. (LANDIS-II DNFS v3)
# reproducible::clearCache(simPaths$cachePath, userTags = c("statsModel"))
# Biomass_core_testSim <- simInitAndSpades(times = simTimes
#                                  , params = simParams
#                                  , modules = simModules[1:6]
#                                  , objects = simObjects
#                                  , paths = simPaths
#                                  , outputs = outputs
#                                  , debug = TRUE
#                                  , .plotInitialTime = NA
# )
# saveRDS(Biomass_core_testSim, file.path(simPaths$outputPath, paste0("simList_", runName, ".rds")))
# simTimes$end <- 21

Biomass_core_testSim <- simInitAndSpades(times = simTimes
                                 , params = simParams
                                 # , modules = simModules[c(1:3, 5:6, 8)] ## for blog post
                                 , modules = simModules[c(1:4, 5:6, 8)]
                                 , objects = simObjects
                                 , paths = simPaths
                                 , outputs = outputs
                                 , debug = TRUE
                                 # , .plotInitialTime = NA
)

saveSimList(Biomass_core_testSim,
            file.path(simPaths$outputPath, paste0("simList_fakeRstCurrentBurn", runName, ".RData")))
dev.print(tiff, file.path(simPaths$outputPath, paste0("simPlots_", runName, ".tiff")),
          res = 300, units = "in")

