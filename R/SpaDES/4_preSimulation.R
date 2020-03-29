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
preSimModules <- list("Biomass_borealDataPrep"
                      , "Biomass_core"
                      , "Biomass_fuelsPFG"
                      , "fireSense_DataPrep"
                      , "fireSense_IgnitionFit"
                      , "fireSense_IgnitionPredict"
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
  , fireSense_IgnitionFit = list(
    formula = formula(n_fires ~ coniferous:julMDC + D2:julMDC +
                        M2:julMDC + O1b:julMDC + NF:julMDC - 1),
    family = poisson(link = "identity"),
    ub = list(coef = 1),
    data = "dataFireSense_IgnitionFit",
    trace = 1,
    iterDEoptim = 60,
    iterNlminb = 100,
    cores = 50
  )
)

if (grepl("newSppParam", runName)) {
  preSimModules <- c("LandR_speciesParameters", preSimModules)

  preSimParams[["LandR_speciesParameters"]] <- list(
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
                                       saveTime = 5.5,
                                       eventPriority = 1))
  outputs <- rbind(outputs, data.frame(objectName = "pixelGroupMap",
                                       saveTime = 5.5,
                                       eventPriority = 1))

simOutPreSim <- Cache(simInitAndSpades
                      , times = list(start = 0, end = 0)
                      , params = preSimParams
                      , modules = preSimModules
                      , paths = simPaths
                      , objects = preSimObjects
                      , outputs = outputs
                      , debug = TRUE
                      , .plotInitialTime = NA
                      # , useCache = "overwrite"
                      , cacheRepo = simPaths$cachePath
                      , userTags = "preSim"
                      , omitArgs = c("userTags"))

## this is not loading properly:
# saveSimList(simOutPreSim,
#             file.path(simPaths$outputPath, paste0("preSimList_fakeRstCurrentBurn", runName, ".RData")))

saveRDS(simOutPreSim,
            file.path(simPaths$outputPath, paste0("preSimList_fakeRstCurrentBurn", runName, ".rds")))

