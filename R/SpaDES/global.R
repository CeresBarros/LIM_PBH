## ------------------------------------------------------
## FIRE MODELLING WITH SpaDES
## TESTS
##
## Ceres: Nov 2017
## ------------------------------------------------------

## clean workspace
rm(list = ls()); amc::.gc()

if (!exists("pkgDir")) {
  pkgDir <- file.path("packages", version$platform,
                      paste0(version$major, ".", strsplit(version$minor, "[.]")[[1]][1]))

  if (!dir.exists(pkgDir)) {
    dir.create(pkgDir, recursive = TRUE)
  }
  .libPaths(pkgDir)
}

# devtools::install_github("PredictiveEcology/reproducible@cb41d78c2cdcaa06d5a98412302c1f4d01850e78", dependencies = FALSE)
# devtools::install_github("achubaty/amc@15c5229951700f9a638fd186f176f0e793d76c10", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/pemisc@dd2be4a9a15981d0d6f3740d8b3d4de07f255b95", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/map@9b401b88ac4d2ceef6de821d718f66a525599d74", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/LandR@466836778f4050d752bc8348d021655226cd504e", dependencies = FALSE)
# devtools::install_github("ianmseddy/LandR.CS@2b056a5d9efea150f3145c8497b33b7fbb726488", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/quickPlot@878dcb2a421c239adfc6d65de37edde58689492b", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/SpaDES.tools@e4add9495d8e9c38d31e5325f37282140d38d8af", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/SpaDES.core@8a7886a6afd7f3b90df10ea6b87caae8661f8709")
# devtools::install_github("PredictiveEcology/SpaDES.experiment@5a23c40f8aa9a9efc6dc16e040f8771561059152", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/fireSenseUtils@4a23a2071599b0c8322d07997871eb53c77db5ff", dependencies = FALSE)

if (!require("Require")) {
  if (!require("devtools")) {
    install.packages("devtools")
  }
  devtools::install_github("PredictiveEcology/Require@development")
}

if (FALSE) {
  Require::pkgSnapshot("packages/pkgSnapshot.txt", standAlone = TRUE, exact = TRUE)
  # Much later on a different or same machine
  # Require::Require(packageVersionFile = "packages/pkgSnapshot.txt", standAlone = TRUE)
}

# devtools::install_local("../LandR/", force = TRUE)
# devtools::install_github("ianmseddy/PSPclean@development")

Require::Require("PredictiveEcology/SpaDES.install@development (>= 0.0.7.9000)")

SpaDES.install::makeSureAllPackagesInstalled("R/SpaDES/m")

Require::Require(c("SpaDES",
                   "raster",
                   "data.table",
                   "CeresBarros/ToolsCB",
                   "PredictiveEcology/SpaDES.experiment",
                   "PredictiveEcology/LandR@LANDISinitialB (>= 1.0.7.9013)",
                   "PredictiveEcology/reproducible"
                   ),
                 upgrade = FALSE)

options("reproducible.useNewDigestAlgorithm" = 2,
        "spades.moduleCodeChecks" = FALSE,
        "reproducible.useCache" = TRUE,
        "reproducible.inputPaths" = normPath("R/SpaDES/inputs"),  ## store everything in data/ so that there are no duplicated files across modules
        "reproducible.destinationPath" = normPath("R/SpaDES/inputs"),
        "reproducible.useGDAL" = FALSE,
        "reproducible.cacheSaveFormat" = "qs",
        "reproducible.useMemoise" = TRUE,
        "spades.useRequire" = TRUE)

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

# eventCaching <- c(".inputObjects", "init")
eventCaching <- NULL ## there's a bug
useParallel <- FALSE

## paths define simulation paths
# simDirName <- "AI_report"
# simDirName <- "jun2021Runs"
simDirName <- "mar2022Runs"

if (Sys.info()["nodename"] == "W-VIC-A127584") {
simPaths <- list(cachePath = file.path("F:", basename(getwd()), "R/SpaDES/cache", simDirName)
                 , modulePath = file.path("R/SpaDES/m")
                 , inputPath = file.path("R/SpaDES/inputs")
                 , outputPath = file.path("F:", basename(getwd()), "R/SpaDES/outputs", simDirName))
} else if (grepl("for-cast", Sys.info()["nodename"])) {
  simPaths <- list(cachePath = file.path("/mnt/scratch/cbarros", basename(getwd()), "R/SpaDES/cache", simDirName)
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
simTimes <- list(start = 2011L, end = 2611L)
vegLeadingProportion <- 0 # indicates what proportion the stand must be in one species group for it to be leading.
# If all are below this, then it is a "mixed" stand
fireInitialTime <- simTimes$start + 5L
fireTimestep <- if (sum(grepl("oneFire", runName))) 100000L else 1L
successionTimestep <- 10L

# reproducible::clearCache(file.path(simPaths$cachePath, "noPM"), userTags = "simInitAndSpades", ask = FALSE)
# reproducible::clearCache(file.path(simPaths$cachePath, "PM"), userTags = "simInitAndSpades", ask = FALSE)
# source("R/SpaDES/4_preSimulation.R")
simListFiles <- list.files(simPaths$outputPath, pattern = "LIM_simInit_", full.names = TRUE, recursive = TRUE)
names(simListFiles) <- sub(paste0(simPaths$outputPath, "/(.*)/.*"), "\\1", simListFiles)
LIM_simInitList <- lapply(simListFiles, loadSimList)

## tests
# end(LIM_simInitList[[1]]) <- 2025
# end(LIM_simInitList[[2]]) <- 2020
# end(simOut1) <- 2112
# simOut1 <- spades(LIM_simInitList[["noPM"]])
# simOut2 <- spades(LIM_simInitList[["PM"]])

# saveSimList(simOut1, filename = file.path(simPaths$outputPath, "noPM", "simOut1.qs"))
# saveSimList(simOut2, filename = file.path(simPaths$outputPath, "PM", "simOut2.qs"))

# plotnoPM <- qs::qread(file.path(simPaths$outputPath, "noPM", "figures", "biomass_by_species_gg.qs"))

## -----------------------------------------------
## SIMULATION RUN
## -----------------------------------------------
# graphics.off()

## using experiment:
library(future)
if (Sys.info()["sysname"] == "Windows") {
  # plan("multisession", workers = 2)   ## each worker consuming roughly 16Gb
  plan("sequential")
} else {
  plan("multicore", workers = 2)
}
simExperimentOut <- experiment2(noPM = LIM_simInitList[["noPM"]],
                                PM = LIM_simInitList[["PM"]],
                                clearSimEnv = TRUE,
                                replicates = 10,
                                useCache = TRUE)
future:::ClusterRegistry("stop")

## save simLists object.
qs::qsave(simExperimentOut, file.path(simPaths$outputPath, paste0("LIM_simLists_noPM_PM", ".qs")))

q("no")
