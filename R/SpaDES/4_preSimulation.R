## --------------------------------------------------
## DEFINE PRE-SIMULATION PARAMETERS AND INPUT OBJECTS
## --------------------------------------------------

## this script defines the modules, parameters, objects and outputs of
## SpaDES PRE- simulations in the Foothills of SW Alberta, according to
## the run name defined in global.R. Run names can differ in their
## parametrisation, modules used and necessary objects

## pre-simulations are used to prepare objects necessary to estimate
## fire frequency using FireSense

## CUSTOMIZE SPECIES TRAIT VALUES -------
## pass a list of species traits parameter values to LandR::updateSpeciesTable
speciesParams <- list(
  "shadetolerance" = list(
    Abie_sp = 2.3,
    Pice_eng = 2.1,
    Pseu_men = 2.0,
    Pice_gla = 1.6,
    Pinu_sp = 1, Popu_sp = 1
  )
)

## SIM MODULES AND PARAMETERS -----------------------------------------
## make two lists of sim modules, one noPM and one PM
simModules <- list("noPM" = list(
  "Biomass_borealDataPrep"
  , "Biomass_speciesParameters"
  , "Biomass_core"
  , "Biomass_fuelsPFG"
  , "fireSense_dataPrep"
  , "fireSense_IgnitionFit"
  , "fireSense_IgnitionPredict"
  , "Biomass_fireProperties"
  , "fireSpread"
  , "Biomass_regeneration"
)
, "PM" = list(
  "Biomass_borealDataPrep"
  , "Biomass_speciesParameters"
  , "Biomass_core"
  , "Biomass_fuelsPFG"
  , "fireSense_dataPrep"
  , "fireSense_IgnitionFit"
  , "fireSense_IgnitionPredict"
  , "Biomass_fireProperties"
  , "fireSpread"
  , "Biomass_regenerationPM"
)
)

## noPM and PM moduels can both be in the same parameter list.
simParams <- list(
  Biomass_borealDataPrep = list(
    "sppEquivCol" = sppEquivCol
    , "forestedLCCClasses" = c(1:15, 34:36)
    , "LCCClassesToReplaceNN" = c(34:36)
    , "ecoregionLayerField" = "ecozoneCode"
    , "fitDeciduousCoverDiscount" = FALSE
    , "exportModels" = "all"
    , "fixModelBiomass" = TRUE
    , "speciesUpdateFunction" = list(
      quote(LandR::speciesTableUpdate(sim$species, sim$speciesTable, sim$sppEquiv, P(sim)$sppEquivCol)),
      quote(LandR::updateSpeciesTable(speciesTable = sim$species, params = sim$speciesParams)))
    # next two are used when assigning pixelGroup membership; what resolution for
    #   age and biomass
    , "pixelGroupAgeClass" = successionTimestep
    , "pixelGroupBiomassClass" = 100
    , "useCloudCacheForStats" = FALSE
    , "cloudFolderID" = NA
    , ".useCache" = eventCaching
  )
  , Biomass_speciesParameters = list(
    "sppEquivCol" = sppEquivCol
    , ".useCache" = eventCaching
  )
  , Biomass_core = list(
    "calcSummaryBGM" = c("start")
    , "initialBiomassSource" = "cohortData" # can be 'biomassMap' or "spinup" too
    , ".plotInitialTime" = plotInitialTime
    , "plotOverstory" = TRUE
    , "seedingAlgorithm" = "wardDispersal"
    , "sppEquivCol" = sppEquivCol
    , "successionTimestep" = successionTimestep
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
  , Biomass_fireProperties = list(
    "fireInitialTime" = fireInitialTime
    , "fireTimestep" = fireTimestep
    , "vegFeedback" = TRUE
    , ".useCache" = eventCaching
  )
  , Biomass_regeneration = list(
    "fireInitialTime" = fireInitialTime
    , "fireTimestep" = fireTimestep
    , "successionTimestep" = fireTimestep
  )
  , Biomass_regenerationPM = list(
    "fireInitialTime" = fireInitialTime
    , "fireTimestep" = fireTimestep
    , "successionTimestep" = fireTimestep
  )
  , fireSense_dataPrep = list(
    "averageWeather4Pred" = TRUE
    , "prepPredictionObjs" = TRUE
    , "propAbsences" = NA
    , "rescalePredictionObjs" = TRUE
    , ".useCache" = eventCaching
  )
  , fireSense_IgnitionFit = list(
    "fireSense_ignitionFormula" = paste0("n_fires ~ coniferous:julMDC + D2:julMDC +
                          M2:julMDC + O1b:julMDC + NF:julMDC +
                          coniferous:pw(julMDC, k_conif) + D2:pw(julMDC, k_D2) +
                          M2:pw(julMDC, k_M2) + O1b:pw(julMDC, k_O1b) + NF:pw(julMDC, k_NF) - 1")
    , "lb" = list(coef = 0,
                  knots = list("julMDC" = 19))   ## the rounded 5% quantile, pre scaling
    , "ub" = list(coef = 1,
                  knots = list("julMDC" = 21))   ## the rounded 80% quantile, pre scaling
    , "iterDEoptim" = 60
    , "iterNlminb" = 100
    , "cores" = 4
    , "rescaleVars" = TRUE
    , "rescalers" = NULL
    , ".plots" = "png"
    , ".useCache" = eventCaching
  )
  , fireSense_IgnitionPredict = list(
    ".runInterval" = NA    ## only run once at the start
  )
  , fireSpread = list(
    "fireInitialTime" = fireInitialTime
    , "fireTimestep" = fireTimestep
    , "noStartPix" = NA  ## NA to make sure this isn't used to randomly draw fires.
    , "spreadProbRange" = c(0.19, 0.24)
    , ".useCache" = eventCaching
  )
)


## SIM OBJECTS ------------------------------------------------
## make base object list
simObjects <- list(
  # "studyArea" = foothillsSMALL
  # , "studyAreaLarge" = foothillsMED
  "studyArea" = foothills
  , "studyAreaLarge" = foothills
  , "sppEquiv" = sppEquivalencies_CA
  , "sppColorVect" = sppColorVect
  , "speciesLayers" = simOutSpeciesLayers$speciesLayers
  , "speciesParams" = speciesParams
  , "ecoregionLayer" = ecoregionLayer
  , "treed" =  simOutSpeciesLayers$treed
  , "numTreed" =  simOutSpeciesLayers$numTreed
  , "nonZeroCover" =  simOutSpeciesLayers$nonZeroCover
  , "PSPgis" = PSPgis
  , "PSPmeasure" = PSPmeasure
  , "PSPplot" = PSPplot
  , "weatherData" = simOutFireWeather$weatherData
  , "weatherDataMDC" = simOutFireWeather$weatherDataMDC
  , "weatherDataMDCCRS" = simOutFireWeather$weatherDataCRS
  , "weatherDataCRS" = simOutFireWeather$weatherDataCRS
)

## SIM OUTPUTS ------------------------------------------------
## save objects at the end of each year
outputs <- data.frame(expand.grid(objectName = c("cohortData"),
                                  saveTime = unique(sort(c(1, simTimes$end,
                                                           seq(simTimes$start, simTimes$end, by = 5)))),
                                  eventPriority = 10,
                                  stringsAsFactors = FALSE))
outputs <- rbind(outputs, data.frame(objectName = "rstCurrentBurn",
                                     saveTime = seq(simParams$fireSpread$fireInitialTime,
                                                    simTimes$end, by = simParams$fireSpread$fireTimestep),
                                     eventPriority = 10))
outputs <- rbind(outputs, data.frame(objectName = "vegTypeMap",
                                     saveTime = unique(sort(c(1, simTimes$end,
                                                              seq(simTimes$start, simTimes$end, by = 5)))),
                                     eventPriority = 10))
outputs <- rbind(outputs, data.frame(objectName = "pixelGroupMap",
                                     saveTime = unique(sort(c(1, simTimes$end,
                                                              seq(simTimes$start, simTimes$end, by = 5)))),
                                     eventPriority = 10))
## on the first year save at the start.
outputs[outputs$saveTime == simTimes$start, "eventPriority"] <- 1

## -----------------------------------------------
## SIMULATION INITIALISATION
## --------------------------------------------

## Run two simInit calls, plus init events
names(runName) <- runName
LIM_simInitList <- lapply(runName, FUN = function(scenario, simPaths, simModules, simParams) {
  simPaths2 <- simPaths
  simPaths2$cachePath <- file.path(simPaths2$cachePath, scenario)
  simPaths2$outputPath <- file.path(simPaths2$outputPath, scenario)
  simModules2 <- simModules[[scenario]]
  Cache(simInitAndSpades
        , times = simTimes
        , params = simParams
        , modules = simModules2
        , loadOrder = unlist(simModules2)
        , paths = simPaths2
        , objects = simObjects
        , outputs = outputs
        , events = "init"
        , debug = TRUE
        , .plotInitialTime = NA
        , cacheRepo = simPaths2$cachePath
        , userTags = c("simInitAndInits", scenario)
        , omitArgs = c("userTags", ".plotInitialTime", "debug"))
}, simPaths = simPaths, simModules = simModules, simParams = simParams)

## save
lapply(names(LIM_simInitList), FUN = function(scenario) {
  saveSimList(LIM_simInitList[[scenario]],
              file.path(outputPath(LIM_simInitList[[scenario]]), paste0("LIM_simInit_", scenario)))
})

