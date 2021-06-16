## ------------------------------------------------------
## FIRE MODELLING WITH SpaDES
## TESTS
##
## Ceres: Nov 2017
## ------------------------------------------------------

## clean workspace
rm(list = ls()); amc::.gc()

## Get packages ----------------------
## requires as of June 16th 2021
# loading reproducible      1.2.7.9006
# loading quickPlot         0.1.7.9002
# loading SpaDES.core       1.0.8.9000
# loading SpaDES.tools      0.3.7.9007
# loading SpaDES.addins     0.1.2
# loading SpaDES.experiment 0.0.2.9002
# loading LandR             1.0.4.0001
# loading LandR.CS          0.0.2.0002
## loading fireSenseUtils   0.0.4.9080

# devtools::install_github("PredictiveEcology/reproducible@DotsBugFix", dependencies = FALSE)
# devtools::install_github("achubaty/amc@development", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/pemisc@development", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/map@development", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/LandR@modelBiomass", dependencies = FALSE)
# devtools::install_github("ianmseddy/LandR.CS", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/quickPlot@development", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/SpaDES.tools@development", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/SpaDES.core@development")
# devtools::install_github("PredictiveEcology/SpaDES.experiment@development", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/fireSenseUtils@development", dependencies = FALSE)

if (!require("Require")) {
  devtools::install_github("PredictiveEcology/Require@development")
  library(Require)
}

SpaDES.install::makeSureAllPackagesInstalled("R/SpaDES/m")

Require(c("SpaDES",
          "raster",
          "data.table",
          "CeresBarros/ToolsCB",
          "PredictiveEcology/SpaDES.experiment",
          "PredictiveEcology/LandR@modelBiomass",
          "PredictiveEcology/reproducible@DotsBugFix"),
        upgrade = FALSE)

options("reproducible.useNewDigestAlgorithm" = 2,
        "spades.moduleCodeChecks" = FALSE,
        "reproducible.useCache" = TRUE,
        "reproducible.inputPaths" = normPath("R/SpaDES/inputs"),  ## store everything in data/ so that there are no duplicated files across modules
        "reproducible.destinationPath" = normPath("R/SpaDES/inputs"),
        "reproducible.useGDAL" = FALSE,
        "reproducible.cacheSaveFormat" = "qs",
        "reproducible.useMemoise" = TRUE)

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
runName <- c("noPM", "PM")

eventCaching <- c(".inputObjects", "init")
useParallel <- FALSE

## paths define simulation paths
# simDirName <- "AI_report"
# simDirName <- "mar2021Runs"
simDirName <- "jun2021Runs"
simPaths <- list(cachePath = file.path("R/SpaDES/cache", simDirName)
                 , modulePath = file.path("R/SpaDES/m")
                 , inputPath = file.path("R/SpaDES/inputs")
                 , outputPath = file.path("R/SpaDES/outputs", simDirName))

## Get necessary objects -----------------------
source("R/SpaDES/1_simObjects.R")

## Run Biomass_speciesData to get species layers
source("R/SpaDES/2_speciesLayers.R")

## maybe drop some species - Black spruce, and Ponderosa pine have v. few occurrences
# plot(simOutSpeciesLayers$speciesLayers)
keepSpp <- sapply(unstack(simOutSpeciesLayers$speciesLayers), FUN = function(ras) {
  propPres <- sum(ras[] > 0, na.rm = TRUE)/sum(!is.na(ras[]))
  propPres > 0.05  ## species need to be in at least 5% of the landscape
})

keepSpp <- names(simOutSpeciesLayers$speciesLayers)[keepSpp]
simOutSpeciesLayers$speciesLayers <- subset(simOutSpeciesLayers$speciesLayers, keepSpp)
sppEquivalencies_CA <- sppEquivalencies_CA[LIM %in% keepSpp]
sppColorVect <- sppColorVect[keepSpp]
# raster::plot(simOutSpeciesLayers$speciesLayers)

## Prepare fire weather tables --------------------
source("R/SpaDES/3_fireWeather.R")

## Run more data prep -----------------------------
# Biomass_borealDataPrep, LandR_speciesParameters, Biomass_core (just init and year 0) and Biomass_fuelsPFG
## to prepare objects for simulation and FireSense ignition/fire frquency fits
## Define simulation params
simTimes <- list(start = 2011L, end = 2111L)
vegLeadingProportion <- 0 # indicates what proportion the stand must be in one species group for it to be leading.
# If all are below this, then it is a "mixed" stand
fireInitialTime <- simTimes$start + 5L
fireTimestep <- if (sum(grepl("oneFire", runName))) 100000L else 1L
successionTimestep <- 10L
plotInitialTime <- NA

# reproducible::clearCache(file.path(simPaths$cachePath, "noPM"), userTags = "simInitAndSpades", ask = FALSE)
# reproducible::clearCache(file.path(simPaths$cachePath, "PM"), userTags = "simInitAndSpades", ask = FALSE)
source("R/SpaDES/4_preSimulation.R")

## tests
# LIM_simInitList <- lapply(list.files(simPaths$outputPath, pattern = "LIM_simInit_", full.names = TRUE, recursive = TRUE),
#                           loadSimList)
# end(LIM_simInitList[[1]]) <- 2020
# end(LIM_simInitList[[2]]) <- 2020
# simOut1 <- spades(LIM_simInitList[[1]])
# simOut2 <- spades(LIM_simInitList[[2]])


## -----------------------------------------------
## SIMULATION RUN
## -----------------------------------------------
# graphics.off()

## using experiment:
library(future)
if (Sys.info()["sysname"] == "Windows") {
  plan("multisession", workers = 2)   ## each worker consuming roughly 16Gb
} else {
  plan("multicore", workers = 2)
}
simExperimentOut <- experiment2(noPM = LIM_simInitList[["noPM"]],
                                PM = LIM_simInitList[["PM"]],
                                clearSimEnv = TRUE,
                                replicates = 10)
future:::ClusterRegistry("stop")

## save simLists object.
qs::qsave(simExperimentOut, file.path(simPaths$outputPath, paste0("LIM_simLists_noPM_PM")))

q("no")
