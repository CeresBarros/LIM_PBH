## ------------------------------------------------------
## FIRE MODELLING WITH SpaDES
## TESTS
##
## Ceres: Nov 2017
## ------------------------------------------------------

## clean workspace
rm(list=ls()); amc::.gc()

## Get packages ----------------------
## requires as of Dec 18th 2020
# loading reproducible     1.0.0.9010
# loading quickPlot        0.1.6.9000
# loading SpaDES.core      1.0.0.9
# loading SpaDES.tools     0.3.4.9000
# loading SpaDES.addins    0.1.2
# loading LandR            0.0.3.9006
# loading LandR.CS         0.0.1.0001

# devtools::install_github("PredictiveEcology/reproducible@development", dependencies = FALSE)
# devtools::install_github("achubaty/amc@development", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/pemisc@development", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/map@development", dependencies = FALSE)
# devtools::install_github("CeresBarros/LandR@dispersalRcpp", dependencies = FALSE)
# devtools::install_github("ianmseddy/LandR.CS", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/quickPlot@development", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/SpaDES.tools@development", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/SpaDES.core@development", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/SpaDES.experiment@development", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/LandWebUtils@development", dependencies = FALSE)  # maybe not needed

## test packages
# devtools::install_local("../LandR", dependencies = FALSE, force = TRUE)
# devtools::install_local("../reproducible", dependencies = FALSE, force = TRUE)
library(SpaDES)
library(SpaDES.experiment)
library(LandR)

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

# runName <- "PM_oneFire_newSppParams"
# runName <- "noPM_oneFire_newSppParams"
# runName <- "PM_newSppParams"
# runName <- "noPM_newSppParams"
runName <- c("noPM_newSppParams_fullSA", "PM_newSppParams_fullSA")

eventCaching <- c(".inputObjects", "init")
useParallel <- FALSE

## paths define simulation paths
# simDirName <- "AI_report"
simDirName <- "dec2020Tests"
simPaths <- list(cachePath = file.path("R/SpaDES/cache", simDirName)
                 , modulePath = file.path("R/SpaDES/m")
                 , inputPath = file.path("R/SpaDES/inputs")
                 , outputPath = file.path("R/SpaDES/outputs", simDirName))

## Get necessary objects -----------------------
source("R/SpaDES/1_simObjects.R")

## Run Biomass_speciesData to get species layers
source("R/SpaDES/2_speciesLayers.R")

## maybe drop some species - Black spruce, and Ponderosa pine have v. few occurrences
plot(simOutSpeciesLayers$speciesLayers)
keepSpp <- sapply(unstack(simOutSpeciesLayers$speciesLayers), FUN = function(ras) {
  propPres <- sum(ras[] > 0, na.rm = TRUE)/sum(!is.na(simOutSpeciesLayers$speciesLayers[[2]][]))
  propPres > 0.05  ## species need to be in at least 5% of the landscaoe
})

keepSpp <- names(simOutSpeciesLayers$speciesLayers)[keepSpp]
simOutSpeciesLayers$speciesLayers <- subset(simOutSpeciesLayers$speciesLayers, keepSpp)
sppEquivalencies_CA <- sppEquivalencies_CA[LIM %in% keepSpp]
sppColorVect <- sppColorVect[keepSpp]
plot(simOutSpeciesLayers$speciesLayers)

## Prepare fire weather tables --------------------
source("R/SpaDES/3_fireWeather.R")

## Run more data prep -----------------------------
# Biomass_borealDataPrep, LandR_speciesParameters, Biomass_core (just init and year 0) and Biomass_fuelsPFG
## to prepare objects for simulation and FireSense ignition/fire frquency fits
## Define simulation params
simTimes <- list(start = 1, end = 100)
vegLeadingProportion <- 0 # indicates what proportion the stand must be in one species group for it to be leading.
# If all are below this, then it is a "mixed" stand
fireInitialTime <- 5L
fireTimestep <- if (sum(grepl("oneFire", runName))) 100000L else 1L
successionTimestep <- 1L
source("R/SpaDES/4_preSimulation.R")

## -----------------------------------------------
## SIMULATION RUN
## -----------------------------------------------
# plotInitialTime <- simTimes$start
plotInitialTime <- NA

## Make actuaL simulation module list, parameters objects, objects and outputs accoding to run
## name and the parameters above
source("R/SpaDES/5_simulationSetup.R")

# simPaths$cachePath <- file.path(simPaths$cachePath, runName[1])
graphics.off()
# simInitOut <- simInitAndSpades(times = simTimes
#                                , params = simParams
#                                , modules = simModulesNoPM
#                                # , modules = simModulesPM
#                                , objects = simObjects
#                                , paths = simPaths
#                                , outputs = outputs
#                                , debug = TRUE
#                                , .plotInitialTime = plotInitialTime
# )

# saveSimList(simOut, file.path(simPaths$outputPath, paste0("simList_", runName, ".qs")))
# if (!is.na(plotInitialTime))
#   dev.print(tiff, file.path(simPaths$outputPath, paste0("simPlots_", runName, ".tiff")),
#             res = 300, units = "in")

## using experiment:
simInitList <- mapply(function(modList, pathSim) {
  simPaths$cachePath <- file.path(simPaths$cachePath, pathSim)
  simInit(times = simTimes
          , params = simParams
          , modules = modList
          , objects = simObjects
          , paths = simPaths
          , outputs = outputs)
},
modList = list(simModulesNoPM, simModulesPM),
pathSim = runName, SIMPLIFY = FALSE)

amc::.gc()
library(future)
plan("multiprocess", workers = 2)
simExperimentOut <- experiment2(noPM = simInitList[[1]], PM = simInitList[[2]],
                                clearSimEnv = TRUE,
                                replicates = 1)
lapply(names(simExperimentOut), FUN = function(simName) {
  saveSimList(simExperimentOut[[simName]],
              filename = file.path(outputPath(simExperimentOut[[simName]]), paste0("simList_", simName, ".qs")))
})

future:::ClusterRegistry("stop")

q("no")
