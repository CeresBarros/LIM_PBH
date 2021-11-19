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

## Get land-cover raster 2001 now that we have a rasterToMatchLarge
# if (is.null(P(simOutSpeciesLayers)$.studyAreaName)) {
#   SAname <- reproducible::studyAreaName(simOutSpeciesLayers$studyAreaLarge)
# }
# rstLCC2005 <- LandR::prepInputsLCC(
#   year = 2005L,
#   destinationPath = simPaths$inputPath,
#   studyArea = simOutSpeciesLayers$studyAreaLarge,   ## Ceres: makePixel table needs same no. pixels for this, RTM rawBiomassMap, LCC.. etc
#   rasterToMatch = simOutSpeciesLayers$rasterToMatchLarge,
#   filename2 = .suffix("rstLCC.tif", paste0("_", SAname)),
#   overwrite = TRUE,
#   cacheRepo = simPaths$cachePath,
#   userTags = c("rstLCC", SAname),
#   omitArgs = c("userTags"))


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
  Biomass_borealDataPrep = list(
    "sppEquivCol" = sppEquivCol
    # , "forestedLCCClasses" = c(1:15, 34:36)
    # , "LCCClassesToReplaceNN" = c(34:36)
    , "forestedLCCClasses" = c(1:6)   ## LCC 2010
    , "LCCClassesToReplaceNN" = numeric(0) ## no replacement - urban/cropland could be grassland, barren may or not suppport veg.
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
    , ".plotInitialTime" = plotInitialTime
    , ".seed" = list("init" = 123)
    , ".useCache" = eventCaching
  )
  , Biomass_speciesParameters = list(
    "sppEquivCol" = sppEquivCol
    , ".plotInitialTime" = plotInitialTime
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
    , ".plotInitialTime" = plotInitialTime
    , ".plotMaps" = FALSE
    , ".useCache" = eventCaching
  )
  , fireProperties = list(
    "fireInitialTime" = fireInitialTime
    , "fireTimestep" = fireTimestep
    , "vegFeedback" = TRUE
    , ".plotInitialTime" = plotInitialTime
    , ".useCache" = eventCaching
  )
  , Biomass_regeneration = list(
    "fireInitialTime" = fireInitialTime
    , "fireTimestep" = fireTimestep
    , "successionTimestep" = fireTimestep
    , ".plotInitialTime" = plotInitialTime
  )
  , Biomass_regenerationPM = list(
    "fireInitialTime" = fireInitialTime
    , "fireTimestep" = fireTimestep
    , "successionTimestep" = fireTimestep
    , ".plotInitialTime" = plotInitialTime
  )
  , fireSense_dataPrep = list(
    "averageWeather4Pred" = TRUE
    , "fitRes" = 1000
    , "prepPredictionObjs" = TRUE
    , "propAbsences" = 9
    , ".plotInitialTime" = plotInitialTime
    , ".useCache" = eventCaching
  )
  , fireSense_IgnitionFit = list(
    "fireSense_ignitionFormula" = paste0("n_fires ~ coniferous:meanMDC + D2:meanMDC +
                          M2:meanMDC + O1b:meanMDC + NF:meanMDC +
                          coniferous:pw(meanMDC, k_conif) + D2:pw(meanMDC, k_D2) +
                          M2:pw(meanMDC, k_M2) + O1b:pw(meanMDC, k_O1b) + NF:pw(meanMDC, k_NF) - 1")
    , "lb" = list(coef = 0,
                  knots = list("meanMDC" = 19))   ## the rounded 5% quantile, pre scaling
    , "ub" = list(coef = 20,
                  knots = list("meanMDC" = 21))   ## the rounded 80% quantile, pre scaling
    , "iterDEoptim" = 60
    , "iterNlminb" = 500
    , "cores" = 4
    , "rescaleVars" = TRUE
    , "rescalers" = NULL
    , ".plots" = "png"
    , ".useCache" = eventCaching
  )
  , fireSense_IgnitionPredict = list(
    ".runInterval" = NA    ## only run once at the start
  )
  , FavierFireSpread = list(
    "fireInitialTime" = fireInitialTime
    , "fireTimestep" = fireTimestep
    , "noStartPix" = NA  ## NA to make sure this isn't used to randomly draw fires.
    , "spreadProbRange" = c(0.20, 0.25)
    , ".plotInitialTime" = plotInitialTime
    , ".useCache" = eventCaching
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
  , "PSPgis" = PSPgis
  , "PSPmeasure" = PSPmeasure
  , "PSPplot" = PSPplot
  , "rasterToMatch" = rasterToMatchLarge
  , "rasterToMatchLarge" = rasterToMatchLarge
  , "rawBiomassMap" = rawBiomassMap
  , "standAgeMap" = standAgeMap
  , "studyArea" = foothills
  , "studyAreaLarge" = foothills
  , "sppEquiv" = sppEquivalencies_CA
  , "sppColorVect" = sppColorVect
  , "speciesLayers" = simOutSpeciesLayers$speciesLayers
  , "speciesParams" = speciesParams
  , "treed" =  simOutSpeciesLayers$treed
  # , "rstLCC" = rstLCC2005
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
              file.path(outputPath(LIM_simInitList[[scenario]]), paste0("LIM_simInit_", scenario)))
})

