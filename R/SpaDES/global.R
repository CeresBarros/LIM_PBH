## ------------------------------------------------------
## FIRE MODELLING WITH SpaDES - global simulation scrip
##
## Ceres: Nov 2017, last updated: July 2025
## ------------------------------------------------------

source("R/SpaDES/0_packages.R")

library(Require)
Require(c("raster",
          "SpaDES",
          "data.table",
          "ToolsCB",
          "SpaDES.experiment",
          "LandR",
          "reproducible",
          "future"
), upgrade = FALSE, install = FALSE)

## to prevent overflow to threads that aren't actually available
data.table::setDTthreads(threads = 1)

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
# eventCaching <- NULL ## there's a bug
useParallel <- FALSE

## paths define simulation paths
# simDirName <- "AI_report"
# simDirName <- "jun2021Runs"
simDirName <- "mar2022Runs"

if (Sys.info()["nodename"] == "W-VIC-A127584") {
  simPaths <- list(cachePath = file.path("F:", basename(getwd()), "R/SpaDES/cache", simDirName)
                   , modulePath = file.path("R/SpaDES/m")
                   , inputPath = file.path("F:", basename(getwd()), "R/SpaDES/inputs", simDirName)
                   , outputPath = file.path("F:", basename(getwd()), "R/SpaDES/outputs", simDirName))
} else if (grepl("for-cast", Sys.info()["nodename"])) {
  simPaths <- list(cachePath = file.path("R/SpaDES/cache", simDirName)
                   , modulePath = file.path("R/SpaDES/m")
                   , inputPath = file.path("R/SpaDES/inputs")
                   , outputPath = file.path("R/SpaDES/outputs", simDirName)
                   , rasterPath = file.path("/mnt/scratch/cbarros", basename(getwd()), "R/SpaDES/scratch/raster")
                   , scratchPath = file.path("/mnt/scratch/cbarros", basename(getwd()), "R/SpaDES/scratch"))
} else {
  simPaths <- list(cachePath = file.path("R/SpaDES/cache", simDirName)
                   , modulePath = file.path("R/SpaDES/m")
                   , inputPath = file.path("R/SpaDES/inputs")
                   , outputPath = file.path("R/SpaDES/outputs", simDirName))
}

if (grepl("for-cast", Sys.info()["nodename"]) ||
    grepl("4458e1a42ddc", Sys.info()["nodename"])) {
  options(bitmapType="cairo")
}


options("reproducible.useNewDigestAlgorithm" = 2,
        "spades.moduleCodeChecks" = FALSE,
        "reproducible.useCache" = TRUE,
        "reproducible.inputPaths" = simPaths$inputPath,
        "reproducible.destinationPath" = simPaths$inputPath,
        "reproducible.useGDAL" = FALSE,
        "reproducible.cacheSaveFormat" = "qs",
        "reproducible.useMemoise" = TRUE,
        "mc.cores" = 1,
        "spades.useRequire" = FALSE,
        "reproducible.useTerra" = FALSE,
        "reproducible.rasterRead" = "raster::raster")

# ## Get necessary objects -----------------------
# source("R/SpaDES/1_simObjects.R")

# ## Run Biomass_speciesData to get species layers
# source("R/SpaDES/2_speciesLayers.R")

## maybe drop some species - Black spruce, and Ponderosa pine have v. few occurrences
# plot(simOutSpeciesLayers$speciesLayers)
# keepSpp <- sapply(unstack(simOutSpeciesLayers$speciesLayers), FUN = function(ras) {
#   propPres <- sum(ras[] > 0, na.rm = TRUE)/sum(!is.na(ras[]))
#   propPres > 0.05  ## species need to be in at least 5% of the landscape
# })

# keepSpp <- names(simOutSpeciesLayers$speciesLayers)[keepSpp]
# speciesLayers <- subset(simOutSpeciesLayers$speciesLayers, keepSpp)
# sppEquivalencies_CA <- sppEquivalencies_CA[get(sppEquivCol) %in% keepSpp]
# sppColorVect <- simOutSpeciesLayers$sppColorVect[keepSpp]
# sppNameVector <- intersect(simOutSpeciesLayers$sppNameVector, names(sppColorVect))
# raster::plot(speciesLayers)

## Prepare fire weather tables --------------------
# source("R/SpaDES/3_fireWeather.R")

## Run more data prep -----------------------------
## Biomass_borealDataPrep, LandR_speciesParameters, Biomass_core (just init and year 0) and Biomass_fuelsPFG
## to prepare objects for simulation and FireSense ignition/fire frequency fits
## Define simulation params
simTimes <- list(start = 2011L, end = 4011L)   ## capture several cycles of Doug-fir dynamics
vegLeadingProportion <- 0 # indicates what proportion the stand must be in one species group for it to be leading.
## If all are below this, then it is a "mixed" stand
fireInitialTime <- simTimes$start + 5L
fireTimestep <- if (sum(grepl("oneFire", runName))) 100000L else 1L
successionTimestep <- 10L

# source("R/SpaDES/4_preSimulation.R")
# LIM_preSimulation <- loadSimList(file.path(simPaths$outputPath, "LIM_preSimulation.qs"))

## re-initialize simulation modules (overcoming seed problem)
# reproducible::clearCache(file.path(simPaths$cachePath, "noPM"), ask = FALSE)
# reproducible::clearCache(file.path(simPaths$cachePath, "PM"), ask = FALSE)
# source("R/SpaDES/4_preSimulationPART2.R")

simListFiles <- list.files(simPaths$outputPath, pattern = "LIM_simInit_", full.names = TRUE, recursive = TRUE)
simListFiles <- grep("PART2.qs$", simListFiles, value = TRUE)
names(simListFiles) <- sub(paste0(simPaths$outputPath, "/(.*)/.*"), "\\1", simListFiles)
LIM_simInitList <- lapply(simListFiles, loadSimList)

## tests
# end(LIM_simInitList[[1]]) <- 2020
# end(LIM_simInitList[[2]]) <- 2020

## try increasing max spread to simulate a greater spread of fire sizes
# P(LIM_simInitList$noPM, "spreadProbRange", "FavierFireSpread") <- c(0.23, 0.26)
# P(LIM_simInitList$PM, "spreadProbRange", "FavierFireSpread") <- c(0.23, 0.25)

## schedule final plots to new end
# LIM_simInitList <- lapply(LIM_simInitList, function(sim) {
#   sim <- scheduleEvent(sim, end(sim),
#                        "Biomass_core", "plotSummaryBySpecies", eventPriority = 9.00)
#   sim <- scheduleEvent(sim, end(sim),
#                        "Biomass_core", "plotAvgs", eventPriority = 9.00 + 0.5)
# })
#
# simOut1 <- spades(LIM_simInitList[["noPM"]])
# simOut2 <- spades(LIM_simInitList[["PM"]])


# saveSimList(simOut1, filename = file.path(simPaths$outputPath, "noPM", "simOut1.qs"))
# saveSimList(simOut2, filename = file.path(simPaths$outputPath, "PM", "simOut2.qs"))

# plotnoPM <- qs::qread(file.path(simPaths$outputPath, "noPM", "figures", "biomass_by_species_gg.qs"))

# LIM_simInitList$PM <- loadSimList(simListFiles[2])
# LIM_simInitList$PM$cohortData <- readRDS("R/SpaDES/outputs/mar2022Runs/PM/PM_rep1/cohortData_year2371.rds")
# LIM_simInitList$PM$pixelGroupMap <- readRDS("R/SpaDES/outputs/mar2022Runs/PM/PM_rep1/pixelGroupMap_year2371.rds")
# LIM_simInitList$PM$vegTypeMap <- readRDS("R/SpaDES/outputs/mar2022Runs/PM/PM_rep1/vegTypeMap_year2371.rds")
# LIM_simInitList$PM$fireYear <- 2011
#
# params(LIM_simInitList$PM)$Biomass_regenerationPM$fireInitialTime <- 2011
# params(LIM_simInitList$PM)$fireProperties$fireInitialTime <- 2011
# params(LIM_simInitList$PM)$FavierFireSpread$fireInitialTime <- 2011
# params(LIM_simInitList$PM)$Biomass_fuelsPFG$fireInitialTime <- 2011
#
# reproducible::clearCache(cachePath(LIM_simInitList$PM), userTags = c("init", "FavierFireSpread"))
# spades(LIM_simInitList[["PM"]], debug = "fireDisturbance")
## -----------------------------------------------
## SIMULATION RUN
## -----------------------------------------------
# graphics.off()

## clean unnecessary simLists
# rm(simOutSpeciesLayers, simOutFireWeather);
# gc()

LIM_simInitList$noPM$._randomSeed <- NULL
LIM_simInitList$PM$._randomSeed <- NULL

## using experiment:
## multicore no longer available from RStudio
plan("multicore", workers = 2)   ##
# plan("multisession", workers = 2)   ## each worker consuming more than 70GB
# plan("multisession", workers = 5)   ## add interactive check for this one/
# plan("sequential")

# out <- future.apply::future_replicate(1, spades(LIM_simInitList[["PM"]], .saveInitialTime = NA))  ## no errors
# future:::ClusterRegistry("stop")

clearSimEnv <- TRUE
simExperimentOut <- experiment2(#noPM = LIM_simInitList[["noPM"]],
  PM = LIM_simInitList[["PM"]],
  clearSimEnv = clearSimEnv,
  replicates = 5,
  useCache = FALSE,
  meanStaggerIntervalInSecs = 60)
future:::ClusterRegistry("stop")

## save simLists object.
if (isFALSE(clearSimEnv)) {  ## we have a caching bug so need to clear the env before saving
  for (i in seq_along(simExperimentOut)) {
    rm(list = ls(simExperimentOut[[i]], all.names = TRUE), envir = envir(s))
  }
}
qs::qsave(simExperimentOut, file.path(outputPath(LIM_simInitList[["PM"]]), paste0("LIM_simLists_PM", ".qs")))

rm(simExperimentOut); gc(reset = TRUE)

simExperimentOut <- experiment2(noPM = LIM_simInitList[["noPM"]],
                                #PM = LIM_simInitList[["PM"]],
                                clearSimEnv = clearSimEnv,
                                replicates = 5,
                                useCache = FALSE,
                                meanStaggerIntervalInSecs = 60)
future:::ClusterRegistry("stop")

## save simLists object.
if (isFALSE(clearSimEnv)) {  ## we have a caching bug so need to clear the env before saving
  for (i in seq_along(simExperimentOut)) {
    rm(list = ls(simExperimentOut[[i]], all.names = TRUE), envir = envir(s))
  }
}
qs::qsave(simExperimentOut, file.path(outputPath(LIM_simInitList[["noPM"]]), paste0("LIM_simLists_noPM", ".qs")))

q("no")
