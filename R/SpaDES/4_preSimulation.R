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
    Pinu_sp = 1,
    Popu_sp = 1
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
  , "fireProperties"
  , "FavierFireSpread"
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
  , "fireProperties"
  , "FavierFireSpread"
  , "Biomass_regenerationPM"
)
)

## noPM and PM moduels can both be in the same parameter list.
simParams <- list(
  .globals = list("dataYear" = 2011L    ## will not be used as the layers have been pre-preped, but just in case...
                  , "initialB" = NA     ## use LANDIS approach to estimate initial cohort B
                  , "sppEquivCol" = sppEquivCol
                  , "vegLeadingProportion" = vegLeadingProportion
                  , ".plots" = c("png", "object")
                  , ".plotInitialTime" = 1
                  , ".useCache" = eventCaching)
  , Biomass_borealDataPrep = list(
    # "forestedLCCClasses" = c(1:15, 34:36)
    # , "LCCClassesToReplaceNN" = c(34:36)
    "forestedLCCClasses" = c(1:6)   ## LCC 2010
    , "LCCClassesToReplaceNN" = numeric(0) ## no replacement - urban/cropland could be grassland, barren may or not support veg.
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
    , ".plots" = NA
    , ".seed" = list("init" = 123)
  )
  , Biomass_speciesParameters = list(
    "constrainGrowthCurve" = c(0, 1)
    , "constrainMaxANPP" = c(3.0, 3.5)
    , "constrainMortalityShape" = c(10, 25)
    , "quantileAgeSubset" = list(Abie_sp = 95,
                                 Pice_eng = 95,
                                 Pseu_men = 95,
                                 Pice_gla = 95,
                                 Pinu_sp = 95,
                                 Popu_sp = 99)
    , "speciesFittingApproach" = "focal"
    , ".plots" = NA
  )
  , Biomass_core = list(
    "calcSummaryBGM" = c("start")
    , "initialBiomassSource" = "cohortData" # can be 'biomassMap' or "spinup" too
    , "plotOverstory" = TRUE
    , "seedingAlgorithm" = "wardDispersal"
    , "successionTimestep" = successionTimestep
    , ".plotInterval" = 1
    , ".plotMaps" = FALSE
    , ".saveInitialTime" = NA
    , ".useParallel" = useParallel
  )
  , Biomass_fuelsPFG = list(
    "fireInitialTime" = fireInitialTime
    , "fireTimestep" = fireTimestep
    , "nonForestFire" = TRUE
    , ".plotMaps" = FALSE
    , ".plots" = NA
  )
  , fireProperties = list(
    "fireInitialTime" = fireInitialTime
    , "fireTimestep" = fireTimestep
    , "vegFeedback" = TRUE
    , ".plots" = NA
  )
  , Biomass_regeneration = list(
    "fireInitialTime" = fireInitialTime
    , "fireTimestep" = fireTimestep
    , "successionTimestep" = fireTimestep
    , ".plots" = NA
  )
  , Biomass_regenerationPM = list(
    "fireInitialTime" = fireInitialTime
    , "fireTimestep" = fireTimestep
    , "successionTimestep" = fireTimestep
    , ".plots" = NA
  )
  , fireSense_dataPrep = list(
    "averageWeather4Pred" = TRUE
    , "fitRes" = 1000
    , "parallel" = TRUE
    , "parallelCores" = 8
    , "prepPredictionObjs" = TRUE
    , "propAbsences" = 9
    , ".plots" = NA
  )
  , fireSense_IgnitionFit = list(
    "fireSense_ignitionFormula" = paste("n_fires ~ coniferous:meanMDC + D2:meanMDC +",
                                        "M2:meanMDC + O1b:meanMDC + NF:meanMDC +",
                                        "coniferous:pw(meanMDC, k_conif) + D2:pw(meanMDC, k_D2) +",
                                        "M2:pw(meanMDC, k_M2) + O1b:pw(meanMDC, k_O1b) + NF:pw(meanMDC, k_NF) - 1")
    , "lb" = list(coef = 0,
                  knots = list("meanMDC" = 19))   ## the rounded 5% quantile, pre scaling
    , "ub" = list(coef = 20,
                  knots = list("meanMDC" = 21))   ## the rounded 80% quantile, pre scaling
    , "iterDEoptim" = 60
    , "iterNlminb" = 500
    , "cores" = 1
    , "rescaleVars" = TRUE
    , "rescalers" = NULL
    , ".plots" = "png"
  )
  , fireSense_IgnitionPredict = list(
    ".runInterval" = NA    ## only run once at the start
  )
  , FavierFireSpread = list(
    "fireInitialTime" = fireInitialTime
    , "fireTimestep" = fireTimestep
    , "noStartPix" = NA  ## NA to make sure this isn't used to randomly draw fires.
    # , "spreadProbRange" = c(0.20, 0.25)   ## over-estimated small fires and under estimated medium/large
    , "spreadProbRange" = c(0.23, 0.25)     ## better alignment with fire size distributions -- see playing with spreadProb.pptx
    , ".plots" = NA
  )
)

## SIM OBJECTS ------------------------------------------------
## make base object list
simObjects <- list(
  # "studyArea" = foothillsSMALL
  # , "studyAreaLarge" = foothillsMED
  "ecoregionLayer" = ecoregionLayer
  , "nonForestFuelsTable" = nonForestFuelsTable
  , "nonZeroCover" =  simOutSpeciesLayers$nonZeroCover
  , "numTreed" =  simOutSpeciesLayers$numTreed
  , "PSPgis_sppParams" = PSPgis_sppParams
  , "PSPmeasure_sppParams" = PSPmeasure_sppParams
  , "PSPplot_sppParams" = PSPplot_sppParams
  , "rasterToMatch" = rasterToMatchLarge
  , "rasterToMatchLarge" = rasterToMatchLarge
  , "rawBiomassMap" = rawBiomassMap
  , "standAgeMap" = standAgeMap
  , "studyArea" = foothills
  , "studyAreaLarge" = foothills
  , "sppEquiv" = sppEquivalencies_CA
  , "sppColorVect" = sppColorVect
  , "sppNameVector" = sppNameVector
  , "speciesLayers" = speciesLayers
  , "speciesParams" = speciesParams
  , "treed" =  simOutSpeciesLayers$treed
  , "weatherData" = simOutFireWeather$weatherData
  , "weatherDataMDC" = simOutFireWeather$weatherDataMDC
  , "weatherDataMDCCRS" = simOutFireWeather$weatherDataCRS
  , "weatherDataCRS" = simOutFireWeather$weatherDataCRS
)

## SIM OUTPUTS ------------------------------------------------
## save objects at the end of each year
outputs <- data.frame(expand.grid(objectName = c("cohortData"),
                                  saveTime = unique(sort(c(simTimes$end,
                                                           seq(simTimes$start, simTimes$end, by = 5)))),
                                  eventPriority = 10,
                                  stringsAsFactors = FALSE))
outputs <- rbind(outputs, data.frame(objectName = "rstCurrentFires",
                                     saveTime = seq(simParams$FavierFireSpread$fireInitialTime,
                                                    simTimes$end, by = simParams$FavierFireSpread$fireTimestep),
                                     eventPriority = 10))
outputs <- rbind(outputs, data.frame(objectName = "severityData",
                                     saveTime = seq(simParams$FavierFireSpread$fireInitialTime,
                                                    simTimes$end, by = simParams$FavierFireSpread$fireTimestep),
                                     eventPriority = 10))
outputs <- rbind(outputs, data.frame(objectName = "vegTypeMap",
                                     saveTime = unique(sort(c(simTimes$end,
                                                              seq(simTimes$start, simTimes$end, by = 5)))),
                                     eventPriority = 10))
outputs <- rbind(outputs, data.frame(objectName = "pixelGroupMap",
                                     saveTime = unique(sort(c(simTimes$end,
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
        , cacheRepo = simPaths2$cachePath
        , userTags = c("simInitAndInits", scenario)
        , omitArgs = c("userTags", ".plotInitialTime", "debug"))
}, simPaths = simPaths, simModules = simModules, simParams = simParams)

## save
lapply(names(LIM_simInitList), FUN = function(scenario) {
  saveSimList(LIM_simInitList[[scenario]],
              file.path(outputPath(LIM_simInitList[[scenario]]), paste0("LIM_simInit_", scenario, ".qs")))
})

