## ------------------------------------------------------
## FIRE MODELLING WITH SpaDES
## TESTS
##
## Ceres: Nov 2017
## ------------------------------------------------------

## clean workspace
rm(list=ls()); amc::.gc()

## Get packages ----------------------
## requires as of June 18th 2019
# loading reproducible     0.2.10.9000
# loading quickPlot        0.1.6.9000
# loading SpaDES.core      0.2.6.9000
# loading SpaDES.tools     0.3.2.9002
# loading SpaDES.addins    0.1.2

# devtools::install_github("PredictiveEcology/reproducible@development", upgrade = "always")
# devtools::install_github("achubaty/amc@development", upgrade = "always")
# devtools::install_github("PredictiveEcology/pemisc@development", upgrade = "always")
# devtools::install_github("PredictiveEcology/map@development", upgrade = "always")
# devtools::install_github("PredictiveEcology/LandR@development", upgrade = "always")
# devtools::install_github("PredictiveEcology/quickPlot@development", upgrade = "always")
# devtools::install_github("PredictiveEcology/SpaDES.tools@development", upgrade = "always")
# devtools::install_github("PredictiveEcology/SpaDES.core@development", upgrade = "always")
library(SpaDES)
library(LandR)

## testing packages
# try(detach("package:LandR", unload = TRUE))
# try(detach("package:SpaDES.core", unload = TRUE))
# try(detach("package:SpaDES.tools", unload = TRUE))
# try(detach("package:reproducible", unload = TRUE))
# devtools::load_all("../reproducible")
# devtools::load_all("../SpaDES.tools")
# devtools::load_all("../SpaDES.core")
# devtools::load_all("../LandR")

options('reproducible.useNewDigestAlgorithm' = TRUE)

source("R/R_tools/Useful_functions.R")

## -----------------------------------------------
## SIMULATION SETUP
## -----------------------------------------------

## Get necessary objects -----------------------
source("R/SpaDES/1_simObjects.R")

## Set up modelling parameters  ---------------------------
# runName <- "blogSep2019_PM"
# runName <- "blogSep2019_noPM"

# runName <- "blogSep2019_PM_oneFire"
runName <- "blogSep2019_noPM_oneFire"
eventCaching <- c(".inputObjects", "init")
useParallel <- FALSE

## Run Biomass_speciesData to get species layers
source("R/SpaDES/2_speciesLayers.R")

## paths define simulation paths
simPaths <- list(cachePath = file.path("R/SpaDES/cache/LIM_tests", runName),
                 modulePath = file.path("R/SpaDES/m"),
                 inputPath = file.path("R/SpaDES/inputs"),
                 outputPath = file.path("R/SpaDES/outputs", runName))

## simulation params
simTimes <- list(start = 0, end = 65)
vegLeadingProportion <- 0 # indicates what proportion the stand must be in one species group for it to be leading.
# If all are below this, then it is a "mixed" stand
fireInitialTime <- 5L
fireTimestep <- if (grepl("oneFire", runName)) 100L else 2L
successionTimestep <- 1L

if (grepl("blogSep2019_noPM", runName)) {
  simModules <- list("Biomass_borealDataPrep"
                     , "Biomass_fireProperties"
                     , "Biomass_fireWeather"
                     , "Biomass_core"
                     , "Biomass_fuels"
                     , "fireSpread"
                     , "Biomass_regeneration"
                     , "fireSeverity"
  )

  simParams <- list(
    Biomass_borealDataPrep = list(
      "sppEquivCol" = sppEquivCol
      , "forestedLCCClasses" = c(1:15, 34:36)
      # next two are used when assigning pixelGroup membership; what resolution for
      #   age and biomass
      , "pixelGroupAgeClass" = successionTimestep * 10L
      , "pixelGroupBiomassClass" = 100
      , "useCloudCacheForStats" = FALSE
      , "cloudFolderID" = NA
      , ".useCache" = eventCaching
    )
    , Biomass_core = list(
      "calcSummaryBGM" = c("start")
      , "initialBiomassSource" = "cohortData" # can be 'biomassMap' or "spinup" too
      , ".plotInitialTime" = simTimes$start
      , "plotOverstory" = TRUE
      , "seedingAlgorithm" = "wardDispersal"
      , "sppEquivCol" = sppEquivCol
      , "successionTimestep" = successionTimestep * 10L
      , "vegLeadingProportion" = vegLeadingProportion
      , ".plotInterval" = 1
      , ".plotMaps" = FALSE
      , ".saveInitialTime" = NA
      , ".useCache" = eventCaching[eventCaching] # seems slower to use Cache for both
      , ".useParallel" = useParallel
    )
    , Biomass_fuels = list(
      "fireInitialTime" = fireInitialTime
      , "fireTimestep" = fireTimestep
      , "nonForestFire" = TRUE
      , "sppEquivCol" = sppEquivCol
      , ".useCache" = eventCaching
    )
    , Biomass_regeneration = list(
      "fireInitialTime" = fireInitialTime
      , "fireTimestep" = fireTimestep
      , "successionTimestep" = successionTimestep
    )
    , Biomass_fireWeather = list(
      ".useCache" = eventCaching
    )
    , Biomass_fireProperties = list(
      "fireInitialTime" = fireInitialTime
      , "fireTimestep" = fireTimestep
      , "vegFeedback" = TRUE
      , ".useCache" = eventCaching
    )
    , fireSpread = list(
      "fireInitialTime" = fireInitialTime
      , "fireTimestep" = fireTimestep
      , "fireSize" = ncell(simOutSpeciesLayers$rasterToMatchLarge)   ## try allowing fires to spread beyond SA
      , "noStartPix" = 10
      , ".useCache" = eventCaching
    )
    , fireSeverity = list(
      "fireTimestep" = fireTimestep
      , ".plotMaps" = TRUE
      , ".saveInitialTime" = 1
      , ".useCache" = eventCaching
    )
  )
}

if (grepl("blogSep2019_PM", runName)) {
  simModules <- list("Biomass_borealDataPrep"
                     , "Biomass_fireProperties"
                     , "Biomass_fireWeather"
                     , "Biomass_core"
                     , "Biomass_fuels"
                     , "fireSpread"
                     , "Biomass_regenerationPM"
                     , "fireSeverity"
  )

  simParams <- list(
    Biomass_borealDataPrep = list(
      "sppEquivCol" = sppEquivCol
      , "forestedLCCClasses" = c(1:15, 34:36)
      # next two are used when assigning pixelGroup membership; what resolution for
      #   age and biomass
      , "pixelGroupAgeClass" = successionTimestep * 10L
      , "pixelGroupBiomassClass" = 100
      , "useCloudCacheForStats" = FALSE
      , "cloudFolderID" = NA
      , ".useCache" = eventCaching
    )
    , Biomass_core = list(
      "calcSummaryBGM" = c("start")
      , "initialBiomassSource" = "cohortData" # can be 'biomassMap' or "spinup" too
      , ".plotInitialTime" = simTimes$start
      , "plotOverstory" = TRUE
      , "seedingAlgorithm" = "wardDispersal"
      , "sppEquivCol" = sppEquivCol
      , "successionTimestep" = successionTimestep * 10L
      , "vegLeadingProportion" = vegLeadingProportion
      , ".plotInterval" = 1
      , ".plotMaps" = FALSE
      , ".saveInitialTime" = NA
      , ".useCache" = eventCaching[eventCaching] # seems slower to use Cache for both
      , ".useParallel" = useParallel
    )
    , Biomass_fuels = list(
      "fireInitialTime" = fireInitialTime
      , "fireTimestep" = fireTimestep
      , "nonForestFire" = TRUE
      , "sppEquivCol" = sppEquivCol
      , ".useCache" = eventCaching
    )
    , Biomass_regenerationPM = list(
      "fireInitialTime" = fireInitialTime
      , "fireTimestep" = fireTimestep
      , "successionTimestep" = successionTimestep
    )
    , Biomass_fireWeather = list(
      ".useCache" = eventCaching
    )
    , Biomass_fireProperties = list(
      "fireInitialTime" = fireInitialTime
      , "fireTimestep" = fireTimestep
      , "vegFeedback" = TRUE
      , ".useCache" = eventCaching
    )
    , fireSpread = list(
      "fireInitialTime" = fireInitialTime
      , "fireTimestep" = fireTimestep
      , "fireSize" = ncell(simOutSpeciesLayers$rasterToMatchLarge)  ## try allowing fires to spread to the whole SA
      , "noStartPix" = 10
      , ".useCache" = eventCaching
    )
    , fireSeverity = list(
      "fireTimestep" = fireInitialTime
      , ".plotMaps" = TRUE
      , ".saveInitialTime" = 1
      , ".useCache" = eventCaching
    )
  )
}

simObjects <- list("studyArea" = foothillsSMALL
                   , "studyAreaLarge" = foothillsMED
                   , "sppEquiv" = sppEquivalencies_CA
                   , "sppColorVect" = sppColorVect
                   , "speciesLayers" = simOutSpeciesLayers$speciesLayers
                   , "treed" =  simOutSpeciesLayers$treed
                   , "numTreed" =  simOutSpeciesLayers$numTreed
                   , "nonZeroCover" =  simOutSpeciesLayers$nonZeroCover
)

if (grepl("oneFire", runName)) {
  outputs <- data.frame(expand.grid(objectName = c("cohortData"),
                                    saveTime = sort(c(1, seq(simTimes$start, simTimes$end,
                                                   by = 5))),
                                    eventPriority = 10,
                                    stringsAsFactors = FALSE))
  outputs[1, "eventPriority"] <- 5.5  ## after init events, before mortalityAndGrowth
  outputs <- rbind(outputs, data.frame(objectName = "rstCurrentBurn",
                                       saveTime = seq(simParams$fireSpread$fireInitialTime,
                                                      simTimes$end, by = simParams$fireSpread$fireTimestep),
                                       eventPriority = 10))
  outputs <- rbind(outputs, data.frame(objectName = "fireCFBRas",
                                       saveTime = seq(simParams$fireSpread$fireInitialTime,
                                                      simTimes$end, by = simParams$fireSpread$fireTimestep),
                                       eventPriority = 10))
  outputs <- rbind(outputs, data.frame(objectName = "vegTypeMap",
                                       saveTime = sort(c(1, seq(simTimes$start, simTimes$end,
                                                                by = 5))),
                                       eventPriority = 10))
  outputs <- rbind(outputs, data.frame(objectName = "pixelGroupMap",
                                       saveTime = sort(c(1, seq(simTimes$start, simTimes$end,
                                                                by = 5))),
                                       eventPriority = 10))
} else {
  outputs <- data.frame(expand.grid(objectName = c("cohortData"),
                                    saveTime = seq(simTimes$start, simTimes$end,
                                                   by = simParams$fireSpread$fireTimestep),
                                    eventPriority = 10,
                                    stringsAsFactors = FALSE))
  outputs[1, "eventPriority"] <- 5.5  ## after init events, before mortalityAndGrowth
  outputs <- rbind(outputs, data.frame(objectName = "rstCurrentBurn",
                                       saveTime = seq(simParams$fireSpread$fireInitialTime,
                                                      simTimes$end, by = simParams$fireSpread$fireTimestep),
                                       eventPriority = 10))
  outputs <- rbind(outputs, data.frame(objectName = "fireCFBRas",
                                       saveTime = seq(simParams$fireSpread$fireInitialTime,
                                                      simTimes$end, by = simParams$fireSpread$fireTimestep),
                                       eventPriority = 10))
  outputs <- rbind(outputs, data.frame(objectName = "vegTypeMap",
                                       saveTime = seq(simTimes$start, simTimes$end,
                                                      by = simParams$fireSpread$fireTimestep),
                                       eventPriority = 10))
  outputs <- rbind(outputs, data.frame(objectName = "pixelGroupMap",
                                       saveTime = seq(simTimes$start, simTimes$end,
                                                      by = 5),
                                       eventPriority = 10))
}

# showCache(simPaths$cachePath, after = "2018-09-26 00:00:00")
# reproducible::clearCache(simPaths$cachePath, userTags = c("prepInputsLCC2005_rtm", "Biomass_borealDataPrep"))

## TODO CHANGE FIRE MODULES TO USE COHORT DATA RATHER THAN SUMMARY BMG OUTPUTS, LIKE BIOMASSMAP
options(spades.moduleCodeChecks = TRUE)
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

## TEST WITH FAKE FIRE MAP
## make fake fire map
rstCurrentBurn <- raster(file.path(sub("_oneFire", "", simPaths$cache), "rasterToMatch.tif"))
IDs <- which(!is.na(rstCurrentBurn[]))
rstCurrentBurn[IDs] <- 1
rstCurrentBurn[IDs[1:round(length(IDs)/2)]] <- NA

simObjects$rstCurrentBurn <- rstCurrentBurn

# simTimes$end <- 21

Biomass_core_testSim <- simInitAndSpades(times = simTimes
                                 , params = simParams
                                 , modules = simModules[c(1:5, 7)]
                                 , objects = simObjects
                                 , paths = simPaths
                                 , outputs = outputs
                                 , debug = TRUE
                                 # , .plotInitialTime = NA
)
saveRDS(Biomass_core_testSim, file.path(simPaths$outputPath, paste0("simList_fakeRstCurrentBurn", runName, ".rds")))
dev.print(tiff, file.path(simPaths$outputPath, paste0("simPlots_", runName, ".tiff")),
          res = 300, units = "in")

