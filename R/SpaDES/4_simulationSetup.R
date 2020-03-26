## --------------------------------------------------
## DEFINE SIMULATION PARAMETERS AND INPUT OBJECTS
## --------------------------------------------------

## this script defines the modules, parameters, objects and outputs of
## SpaDES simulations in the Foothills of SW Alberta, according to
## the run name defined in global.R. Run names can differ in their
## parametrisation, modules used and necessary objects

## SIM PARAMS ------------------------------------------------

if (grepl("noPM", runName)) {
  simModules <- list("Biomass_borealDataPrep"
                     , "Biomass_fireProperties"
                     , "Biomass_core"
                     , "Biomass_fuelsPFG"
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
      , ".useCache" = eventCaching # seems slower to use Cache for both
      , ".useParallel" = useParallel
    )
    , Biomass_fuelsPFG = list(
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
} else {
  simModules <- list("Biomass_borealDataPrep"
                     , "Biomass_fireProperties"
                     , "Biomass_core"
                     , "Biomass_fuelsPFG"
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
      , ".useCache" = eventCaching # seems slower to use Cache for both
      , ".useParallel" = useParallel
    )
    , Biomass_fuelsPFG = list(
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

## add LandR_speciesParameters module if need be
if (grepl("newSppParam", runName)) {
  simModules <- c("LandR_speciesParameters", simModules)

  simParams[["LandR_speciesParameters"]] <- list(
    "sppEquivCol" = sppEquivCol
    , ".useCache" = eventCaching
    )
}

## SIM OBJECTS ------------------------------------------------
## make base object list
simObjects <- list("studyArea" = foothillsSMALL
                   , "studyAreaLarge" = foothillsMED
                   , "sppEquiv" = sppEquivalencies_CA
                   , "sppColorVect" = sppColorVect
                   , "speciesLayers" = simOutSpeciesLayers$speciesLayers
                   , "treed" =  simOutSpeciesLayers$treed
                   , "numTreed" =  simOutSpeciesLayers$numTreed
                   , "nonZeroCover" =  simOutSpeciesLayers$nonZeroCover
                   , "weatherData" = simOutFireWeather$weatherData
                   , "weatherDataCRS" = simOutFireWeather$weatherDataCRS
)

## add PSP data if need be
if (grepl("newSppParams", runName)) {
  simObjects$PSPgis <- PSPgis
  simObjects$PSPmeasure <- PSPmeasure
  simObjects$PSPplot <- PSPplot
}

## add fake fire map if need be
if (grepl("oneFire", runName)) {
  rstCurrentBurn <- raster("R/SpaDES/cache/LIM_tests/blogSep2019_noPM_oneFire/rasterToMatch.tif")
  IDs <- which(!is.na(rstCurrentBurn[]))
  rstCurrentBurn[IDs] <- 1
  rstCurrentBurn[IDs[1:round(length(IDs)/2)]] <- NA

  simObjects$rstCurrentBurn <- rstCurrentBurn
}


## SIM OUTPUTS ------------------------------------------------
if (grepl("oneFire", runName)) {
  outputs <- data.frame(expand.grid(objectName = c("cohortData"),
                                    saveTime = sort(c(1, seq(simTimes$start, simTimes$end,
                                                             by = 5))),
                                    eventPriority = 10,
                                    stringsAsFactors = FALSE))
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
  ## on the first year save after init events, but before mortalityAndGrowth
  outputs[outputs$saveTime == 0, "eventPriority"] <- 5.5
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
