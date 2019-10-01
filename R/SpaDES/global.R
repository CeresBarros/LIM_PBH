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
runName <- "blogSep2019_noPM"
eventCaching <- c(".inputObjects", "init")
useParallel <- FALSE

## paths define simulation paths
simPaths <- list(cachePath = file.path("R/SpaDES/cache/LIM_tests", runName),
                 modulePath = file.path("R/SpaDES/m"),
                 inputPath = file.path("R/SpaDES/inputs"),
                 outputPath = file.path("R/SpaDES/outputs", runName))

## simulation params
simTimes <- list(start = 0, end = 100)
vegLeadingProportion <- 0 # indicates what proportion the stand must be in one species group for it to be leading.
# If all are below this, then it is a "mixed" stand
fireTimestep <- 2L
successionTimestep <- 1L

if (runName == "blogSep2019_noPM") {
  simModules <- list("Boreal_LBMRDataPrep"
                     , "Biomass_fireProperties"
                     , "LBMR"
                     , "Biomass_fuels"
                     , "fireSpread"
                     , "Biomass_regeneration"
                     , "fireSeverity"
  )

  simParams <- list(
    Boreal_LBMRDataPrep = list(
      "sppEquivCol" = sppEquivCol
      , "forestedLCCClasses" = c(1:15, 34:36)
      # next two are used when assigning pixelGroup membership; what resolution for
      #   age and biomass
      , "pixelGroupAgeClass" = successionTimestep * 10L
      , "pixelGroupBiomassClass" = 100
      , "runName" = runName
      , "useCloudCacheForStats" = FALSE
      , "cloudFolderID" = NA
      , ".useCache" = eventCaching
    )
    , LBMR = list(
      "calcSummaryBGM" = c("start")
      , "initialBiomassSource" = "cohortData" # can be 'biomassMap' or "spinup" too
      , ".plotInitialTime" = simTimes$start
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
      "fireInitialTime" = fireTimestep
      , "fireTimestep" = fireTimestep
      , "sppEquivCol" = sppEquivCol
      , ".useCache" = eventCaching
    )
    , Biomass_regeneration = list(
      "fireInitialTime" = fireTimestep
      , "fireTimestep" = fireTimestep
      , "successionTimestep" = successionTimestep
    )
    , Biomass_fireProperties = list(
      "fireInitialTime" = fireTimestep
      , "fireTimestep" = fireTimestep
      , "vegFeedback" = TRUE
      , ".useCache" = eventCaching
    )
    , fireSpread = list(
      "fireInitialTime" = fireTimestep
      , "fireTimestep" = fireTimestep
      , "fireSize" = 1000L
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

if (runName == "blogSep2019_PM") {
  simModules <- list("Boreal_LBMRDataPrep"
                     , "Biomass_fireProperties"
                     , "LBMR"
                     , "Biomass_fuels"
                     , "fireSpread"
                     , "Biomass_regenerationPM"
                     , "fireSeverity"
  )

  simParams <- list(
    Boreal_LBMRDataPrep = list(
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
    , LBMR = list(
      "calcSummaryBGM" = c("start")
      , "initialBiomassSource" = "cohortData" # can be 'biomassMap' or "spinup" too
      , ".plotInitialTime" = simTimes$start
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
      "fireInitialTime" = fireTimestep
      , "fireTimestep" = fireTimestep
      , "sppEquivCol" = sppEquivCol
      , ".useCache" = eventCaching
    )
    , Biomass_regenerationPM = list(
      "fireInitialTime" = fireTimestep
      , "fireTimestep" = fireTimestep
      , "successionTimestep" = successionTimestep
    )
    , Biomass_fireProperties = list(
      "fireInitialTime" = fireTimestep
      , "fireTimestep" = fireTimestep
      , "vegFeedback" = TRUE
      , ".useCache" = eventCaching
    )
    , fireSpread = list(
      "fireInitialTime" = fireTimestep
      , "fireTimestep" = fireTimestep
      , "fireSize" = 1000L
      , "noStartPix" = 10
      , "vegFeedback" = TRUE
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

## Run Biomass_speciesData to get species layers
source("R/SpaDES/2_speciesLayers.R")

simObjects <- list("studyArea" = foothillsSMALL
                   , "studyAreaLarge" = foothillsMED
                   , "sppEquiv" = sppEquivalencies_CA
                   , "sppColorVect" = sppColorVect
                   , "speciesLayers" = simOutSpeciesLayers$speciesLayers
                   , "treed" =  simOutSpeciesLayers$treed
                   , "numTreed" =  simOutSpeciesLayers$numTreed
                   , "nonZeroCover" =  simOutSpeciesLayers$nonZeroCover
)

outputs <- data.frame(expand.grid(objectName = c("cohortData"),
                                  saveTime = seq(2, 50, by = 5),
                                  stringsAsFactors = FALSE))
outputs <- rbind(outputs, data.frame(objectName = "rstCurrentBurn",
                                     saveTime = seq(2, 50, by = 5)))
outputs <- rbind(outputs, data.frame(objectName = "fireCFBRas",
                                     saveTime = seq(2, 50, by = 5)))

# showCache(simPaths$cachePath, after = "2018-09-26 00:00:00")
# reproducible::clearCache(simPaths$cachePath, userTags = c("prepInputsLCC2005_rtm", "Boreal_LBMRDataPrep"))

## TODO CHANGE FIRE MODULES TO USE COHORT DATA RATHER THAN SUMMARY BMG OUTPUTS, LIKE BIOMASSMAP
options(spades.moduleCodeChecks = TRUE)
graphics.off()

# reproducible::clearCache(simPaths$cachePath, userTags = c("^LBMR$", "init"), ask = FALSE)
## TODO TEST PARAMETRISING WITH SA LARGE
# set.seed(524326)

## TODO RUN SIMUALTIONS W/ AND W/O PM for blog
# set.seed(524326)

## TODO: implement LANDIS pixel fire severity calculation:
## Each fire event has an associated mean fire severity which is the average of the severities at all of the event’s sites. (LANDIS-II DNFS v3)
# reproducible::clearCache(simPaths$cachePath, userTags = c("Boreal", ".inputObjects", "rasterToMatch"))
LBMR_testSim <- simInitAndSpades(times = simTimes
                                 , params = simParams
                                 , modules = simModules[1:6]
                                 , objects = simObjects
                                 , paths = simPaths
                                 , debug = TRUE
                                 , .plotInitialTime = NA
)


## TEST WITH FAKE FIRE MAP
## make fake fire map
rstCurrentBurn <- SpaDES.core:::.pkgEnv$.sim$pixelGroupMap
rstCurrentBurn[rstCurrentBurn[]>0] <- 1
rstCurrentBurn[rstCurrentBurn[] <= 0] <- NA
IDs <- which(rstCurrentBurn[] == 1)
rstCurrentBurn[IDs[1:round(length(IDs)/2)]] <- NA

simObjects$rstCurrentBurn <- rstCurrentBurn

LBMR_testSim <- simInitAndSpades(times = simTimes
                                 , params = simParams
                                 , modules = simModules[1:6]
                                 , objects = simObjects
                                 , paths = simPaths
                                 , debug = TRUE
                                 , .plotInitialTime = NA
)

