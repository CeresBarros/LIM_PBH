## --------------------------------------------------
## DEFINE PRE-SIMULATION PARAMETERS AND INPUT OBJECTS
## --------------------------------------------------

## this script defines the modules, parameters, objects and outputs of
## SpaDES PRE- simulations in the Foothills of SW Alberta, according to
## the run name defined in global.R. Run names can differ in their
## parametrisation, modules used and necessary objects

## pre-simulations are used to prepare objects necessary to estimate
## fire frequency using FireSense

## SIM PARAMS ------------------------------------------------
preSimPaths <- list(cachePath = file.path("R/SpaDES/cache/LIM_tests/preSim"),
                    modulePath = file.path("R/SpaDES/m"),
                    inputPath = file.path("R/SpaDES/inputs"),
                    outputPath = file.path("R/SpaDES/outputs/preSim"))

preSimModules <- list("Biomass_borealDataPrep"
                      , "Biomass_core"
                      , "Biomass_fuelsPFG"
                      , "fireSense_dataPrep"
                      , "fireSense_IgnitionFit"
)

preSimParams <- list(
  Biomass_borealDataPrep = list(
    "sppEquivCol" = sppEquivCol
    , "forestedLCCClasses" = c(1:15, 34:36)
    , "LCCClassesToReplaceNN" = c(34:36)
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
    "fireInitialTime" = 0
    , "fireTimestep" = 1
    , "nonForestFire" = TRUE
    , "sppEquivCol" = sppEquivCol
    , ".plotMaps" = FALSE
    , ".useCache" = eventCaching
  )
  , fireSense_dataPrep = list(
    "fireInitialTime" = 0
    , "fireTimestep" = 1
    , ".useCache" = eventCaching
  )
  , fireSense_IgnitionFit = list(
    formula = formula(n_fires ~ coniferous:julMDC + D2:julMDC +
                        M2:julMDC + O1b:julMDC + NF:julMDC - 1)
    , family = poisson(link = "identity")
    , ub = list(coef = 1)
    , data = "dataFireSense_IgnitionFit"
    , trace = 1
    , iterDEoptim = 60
    , iterNlminb = 100
    , cores = 1
    , ".useCache" = c(eventCaching, "run")
  )
)

if (grepl("newSppParam", runName)) {
  preSimModules <- c("Biomass_speciesParameters", preSimModules)

  preSimParams[["Biomass_speciesParameters"]] <- list(
    "sppEquivCol" = sppEquivCol
    , ".useCache" = eventCaching
  )
}

## SIM OBJECTS ------------------------------------------------
## make base object list
preSimObjects <- list("studyArea" = foothillsSMALL
                      , "studyAreaLarge" = foothillsMED
                      , "sppEquiv" = sppEquivalencies_CA
                      , "sppColorVect" = sppColorVect
                      , "speciesLayers" = simOutSpeciesLayers$speciesLayers
                      , "treed" =  simOutSpeciesLayers$treed
                      , "numTreed" =  simOutSpeciesLayers$numTreed
                      , "nonZeroCover" =  simOutSpeciesLayers$nonZeroCover
                      , "weatherData" = simOutFireWeather$weatherData
                      , "weatherDataMDC" = simOutFireWeather$weatherDataMDC
                      , "weatherDataMDCCRS" = simOutFireWeather$weatherDataCRS
)

## add PSP data if need be
if (grepl("newSppParams", runName)) {
  preSimObjects$PSPgis <- PSPgis
  preSimObjects$PSPmeasure <- PSPmeasure
  preSimObjects$PSPplot <- PSPplot
}


## SIM OUTPUTS ------------------------------------------------
## on the first year save after init events, but before mortalityAndGrowth
outputs <- data.frame(expand.grid(objectName = c("cohortData"),
                                  saveTime = 0,
                                  eventPriority = 5.5,
                                  stringsAsFactors = FALSE))
outputs <- rbind(outputs, data.frame(objectName = "fuelTypesMaps",
                                     saveTime = 0,
                                     eventPriority = 5.5))
outputs <- rbind(outputs, data.frame(objectName = "pixelGroupMap",
                                     saveTime = 0,
                                     eventPriority = 5.5))

simOutPreSim <- Cache(simInitAndSpades
                      , times = list(start = 0, end = 0)
                      , params = preSimParams
                      , modules = preSimModules
                      , paths = preSimPaths
                      , objects = preSimObjects
                      , outputs = outputs
                      , debug = TRUE
                      , .plotInitialTime = NA
                      # , useCache = "overwrite"
                      , cacheRepo = preSimPaths$cachePath
                      , userTags = "preSim"
                      , omitArgs = c("userTags"))

## this is not loading properly:
saveSimList(simOutPreSim,
            file.path(preSimPaths$outputPath, "preSimList.qs"))
# saveRDS(simOutPreSim,
#             file.path(simPaths$outputPath, paste0("preSimList_fakeRstCurrentBurn", runName, ".rds")))

## ESTIMATE FIRE FREQUENCY------------------------------------------------
## fireSenseIgnitionPredict needs to be run separately as the objects can't
## be directly passed between modules

## Define fireSense_IgnitionPredict module inputs
## we will use an average MDC across the years and output a raster that will not change
## in time
objects <- list(
  "fireSense_IgnitionFitted" = simOutPreSim$fireSense_IgnitionFitted
  , "D2" = simOutPreSim$fuelTypesCoverStk$D2
  , "M2" = simOutPreSim$fuelTypesCoverStk$M2
  , "coniferous" = simOutPreSim$fuelTypesCoverStk$coniferous
  , "O1b" = simOutPreSim$fuelTypesCoverStk$O1b
  , "NF" = simOutPreSim$fuelTypesCoverStk$NF
  , "julMDC" = mean(simOutPreSim$weatherDataMDCStk)
)

# Define fireSense_IgnitionPredict module outputs
outputs <- rbind(
  data.frame(
    file = paste0("fireSense_IgnitionPredicted.tif"),
    fun = "writeRaster",
    objectName = "fireSense_IgnitionPredicted",
    package = "raster",
    saveTime = 1
  )
)

# Define fireSense_IgnitionPredict module parameters
parameters <- list(
  fireSense_IgnitionPredict = list(
    "data" = c("D2", "M2", "coniferous", "O1b", "NF", "julMDC"),
    "modelObjName" = "fireSense_IgnitionFitted" # This is the default
  )
)

# Run the simulation
simOutFireFreq <- Cache(simInitAndSpades
                        , times = list(start = 0, end = 0)
                        , modules = "fireSense_IgnitionPredict"
                        , paths = preSimPaths
                        , objects = objects
                        , outputs = outputs
                        , params = parameters
                        , cacheRepo = preSimPaths$cachePath
                        , userTags = "preSim"
                        , omitArgs = c("userTags")
)
saveSimList(simOutFireFreq,
            file.path(preSimPaths$outputPath, "preSimList_FireSense_IgnitionPredict.qs"))



## DIAGNOSE FIRE IGNITION MODEL ----------------------------------
## to get predicted values we re-run the IgnitionPredict with the original data
objects <- list(
  "fireSense_IgnitionFitted" = simOutPreSim$fireSense_IgnitionFitted
  , "dataFireSense_IgnitionFit" = simOutPreSim$dataFireSense_IgnitionFit
)

parameters <- list(
  fireSense_IgnitionPredict = list(
    "data" = "dataFireSense_IgnitionFit",
    "modelObjName" = "fireSense_IgnitionFitted" # This is the default
  )
)
simOutFireFreqPredVals <- Cache(simInitAndSpades
                        , times = list(start = 0, end = 0)
                        , modules = "fireSense_IgnitionPredict"
                        , paths = preSimPaths
                        , objects = objects
                        , params = parameters
                        , cacheRepo = preSimPaths$cachePath
                        , userTags = "preSim"
                        , omitArgs = c("userTags")
)

## plot predicted and fitted values
ignitionsData <- na.omit(simOutPreSim$dataFireSense_IgnitionFit)
ignitionsData[, fittedVals := simOutFireFreqPredVals$fireSense_IgnitionPredicted]

ggplot(data = ignitionsData, aes(y = fittedVals, x = n_fires, col = year)) +
  geom_point()

