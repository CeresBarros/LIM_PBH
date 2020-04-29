## --------------------------------------------------
## DEFINE SIMULATION PARAMETERS AND INPUT OBJECTS
## --------------------------------------------------

## this script defines the modules, parameters, objects and outputs of
## SpaDES simulations in the Foothills of SW Alberta, according to
## the run name defined in global.R. Run names can differ in their
## parametrisation, modules used and necessary objects

## SIM PARAMS ------------------------------------------------

if (grepl("noPM", runName)) {
  simModules <- list("Biomass_fireProperties"
                     , "Biomass_core"
                     , "Biomass_fuelsPFG"
                     , "fireSpread"
                     , "Biomass_regeneration"
                     , "fireSeverity"
  )

  simParams <- list(
    Biomass_core = list(
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
      , ".useCache" = eventCaching[1] # seems slower to use Cache for both
      , ".useParallel" = useParallel
    )
    , Biomass_fuelsPFG = list(
      "fireInitialTime" = fireInitialTime
      , "fireTimestep" = fireTimestep
      , "nonForestFire" = TRUE
      , "sppEquivCol" = sppEquivCol
      , ".plotMaps" = FALSE
      , ".useCache" = eventCaching
    )
    , Biomass_regeneration = list(
      "fireInitialTime" = fireInitialTime
      , "fireTimestep" = fireTimestep
      , "successionTimestep" = successionTimestep
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
      , "noStartPix" = NA  ## NA to make sure this isn't used to randomly draw fires.
      , ".useCache" = eventCaching
    )
    , fireSeverity = list(
      "fireTimestep" = fireTimestep
      , ".plotMaps" = TRUE
      , ".saveInitialTime" = 1
      , ".useCache" = eventCaching
    )
  )
} else {
  simModules <- list("Biomass_fireProperties"
                     , "Biomass_core"
                     , "Biomass_fuelsPFG"
                     , "fireSpread"
                     , "Biomass_regenerationPM"
                     , "fireSeverity"
  )

  simParams <- list(
    Biomass_core = list(
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
      , ".useCache" = eventCaching # seems slower to use Cache for both
      , ".useParallel" = useParallel
    )
    , Biomass_fuelsPFG = list(
      "fireInitialTime" = fireInitialTime
      , "fireTimestep" = fireTimestep
      , "nonForestFire" = TRUE
      , "sppEquivCol" = sppEquivCol
      , ".plotMaps" = FALSE
      , ".useCache" = eventCaching
    )
    , Biomass_regenerationPM = list(
      "fireInitialTime" = fireInitialTime
      , "fireTimestep" = fireTimestep
      , "successionTimestep" = successionTimestep
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
      , "noStartPix" = NA  ## NA to make sure this isn't used to randomly draw fires.
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

## SIM OBJECTS ------------------------------------------------
## make base object list
simObjects <- lapply(ls(simOutPreSim@.xData), FUN = function(x) {
  get(x, envir = simOutPreSim@.xData)
})

names(simObjects) <- ls(simOutPreSim@.xData)

simObjects <- c(simObjects,
                list("weatherData" = simOutFireWeather$weatherData
                     , "weatherDataCRS" = simOutFireWeather$weatherDataCRS
                     , "fireIgnitionProb" = simOutFireFreq$fireSense_IgnitionPredicted
                )
)

## change the site shade level relative biomass.
simObjects$minRelativeB[, X1 := 0.1]
simObjects$minRelativeB[, X2 := 0.2]
simObjects$minRelativeB[, X3 := 0.3]  ## this and next were added in second round
simObjects$minRelativeB[, X4 := 0.5]
simObjects$minRelativeB[, X5 := 0.7]


## add fake fire map if need be
if (grepl("oneFire", runName)) {
  rstCurrentBurn <- raster(list.files("R/SpaDES/cache/LIM_tests", recursive = TRUE,
                                      pattern = "rasterToMatch.tif", full.names = TRUE)[1])
  IDs <- which(!is.na(rstCurrentBurn[]))
  rstCurrentBurn[IDs] <- 1
  rstCurrentBurn[IDs[1:round(length(IDs)/2)]] <- NA

  simObjects$rstCurrentBurn <- rstCurrentBurn
}


## SIM OUTPUTS ------------------------------------------------
outputs <- data.frame(expand.grid(objectName = c("cohortData"),
                                  saveTime = unique(sort(c(1, simTimes$end,
                                                           seq(simTimes$start, simTimes$end, by = 5)))),
                                  eventPriority = 10,
                                  stringsAsFactors = FALSE))
outputs <- rbind(outputs, data.frame(objectName = "rstCurrentBurn",
                                     saveTime = seq(simParams$fireSpread$fireInitialTime,
                                                    simTimes$end, by = simParams$fireSpread$fireTimestep),
                                     eventPriority = 10))
# outputs <- rbind(outputs, data.frame(objectName = "fireCFBRas",
#                                      saveTime = seq(simParams$fireSpread$fireInitialTime,
#                                                     simTimes$end, by = simParams$fireSpread$fireTimestep),
#                                      eventPriority = 10))
outputs <- rbind(outputs, data.frame(objectName = "vegTypeMap",
                                     saveTime = unique(sort(c(1, simTimes$end,
                                                              seq(simTimes$start, simTimes$end, by = 5)))),
                                     eventPriority = 10))
outputs <- rbind(outputs, data.frame(objectName = "pixelGroupMap",
                                     saveTime = unique(sort(c(1, simTimes$end,
                                                              seq(simTimes$start, simTimes$end, by = 5)))),
                                     eventPriority = 10))
## on the first year save after init events, but before mortalityAndGrowth
outputs[outputs$saveTime == 0, "eventPriority"] <- 5.5

