## --------------------------------------------------
## RE-INITIALIZE SIMULATION MODULES ONLY
## --------------------------------------------------

## this script re-runs a `simInit` call for the simulation
## modules to avoid inheriting a fixed seed across replicates
## Run names can differ in their
## parametrisation, modules used and necessary objects


## SIM MODULES AND PARAMETERS -----------------------------------------
## make two lists of sim modules, one noPM and one PM
simModules <- list("noPM" = list(
  "Biomass_core"
  , "Biomass_fuelsPFG"
  , "fireProperties"
  , "FavierFireSpread"
  , "Biomass_regeneration"
)
, "PM" = list(
  "Biomass_core"
  , "Biomass_fuelsPFG"
  , "fireProperties"
  , "FavierFireSpread"
  , "Biomass_regenerationPM"
)
)

## noPM and PM moduels can both be in the same parameter list.
simParams <- list(
  .globals = list("dataYear" = 2011L    ## will not be used as the layers have been pre-preped, but just in case...
                  , "initialB" = NA     ## use LANDIS approach to estimate initial cohort B
                  , "sppEquivCol" = params(LIM_preSimulation)$.globals$sppEquivCol
                  , "vegLeadingProportion" = params(LIM_preSimulation)$.globals$vegLeadingProportion
                  , ".plots" = c("png", "object")
                  , ".plotInitialTime" = 1
                  , ".useCache" = eventCaching)
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
    , ".useCache" = ".inputObjects"   ## don't cache init, because mod$ objects aren't being cached
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
  , FavierFireSpread = list(
    "fireInitialTime" = fireInitialTime
    , "fireTimestep" = fireTimestep
    , "noStartPix" = NA  ## NA to make sure this isn't used to randomly draw fires.
    # , "spreadProbRange" = c(0.20, 0.25)   ## over-estimated small fires and under estimated medium/large
    , "spreadProbRange" = c(0.23, 0.25)     ## better alignment with fire size distributions -- see playing with spreadProb.pptx
    , ".plots" = NA
  )
)

## SIM OUTPUTS ------------------------------------------------
## save objects at the end of each year
simOutputs <- data.frame(expand.grid(objectName = c("cohortData"),
                                  saveTime = unique(sort(c(simTimes$end,
                                                           seq(simTimes$start, simTimes$end, by = 5)))),
                                  eventPriority = 10,
                                  stringsAsFactors = FALSE))
simOutputs <- rbind(simOutputs, data.frame(objectName = "rstCurrentFires",
                                     saveTime = seq(simParams$FavierFireSpread$fireInitialTime,
                                                    simTimes$end, by = simParams$FavierFireSpread$fireTimestep),
                                     eventPriority = 10))
simOutputs <- rbind(simOutputs, data.frame(objectName = "severityData",
                                     saveTime = seq(simParams$FavierFireSpread$fireInitialTime,
                                                    simTimes$end, by = simParams$FavierFireSpread$fireTimestep),
                                     eventPriority = 10))
simOutputs <- rbind(simOutputs, data.frame(objectName = "vegTypeMap",
                                     saveTime = unique(sort(c(simTimes$end,
                                                              seq(simTimes$start, simTimes$end, by = 5)))),
                                     eventPriority = 10))
simOutputs <- rbind(simOutputs, data.frame(objectName = "pixelGroupMap",
                                     saveTime = unique(sort(c(simTimes$end,
                                                              seq(simTimes$start, simTimes$end, by = 5)))),
                                     eventPriority = 10))

## on the first year save at the start.
simOutputs[simOutputs$saveTime == simTimes$start, "eventPriority"] <- 1

## -----------------------------------------------
## SIMULATION INITIALISATION
## ----------------------------------------------

## Run two simInit calls, plus init events
names(runName) <- runName
LIM_simInitList <- Map(scenario = runName,
                        MoreArgs = list(simInitOut = LIM_preSimulation,
                                        simModules = simModules,
                                        simTimes = simTimes,
                                        simParams = simParams,
                                        simOutputs = simOutputs),
                        f = function(scenario, simInitOut, simModules, simTimes, simParams, simOutputs) {

                          simModules <- simModules[[scenario]]
                          simPaths2 <- paths(simInitOut)
                          simPaths2$cachePath <- file.path(simPaths2$cachePath, scenario)
                          simPaths2$outputPath <- file.path(simPaths2$outputPath, scenario)
                          simObjects <- mget(ls(simInitOut@.xData), as.environment(simInitOut))
                          Cache(simInit
                                , times = simTimes
                                , params = simParams
                                , modules = simModules
                                , loadOrder = unlist(simModules)
                                , paths = simPaths2
                                , objects = simObjects
                                , outputs = simOutputs
                                , cacheRepo = simPaths2$cachePath
                                , userTags = c("simInit", scenario)
                                , omitArgs = c("userTags", ".plotInitialTime", "debug"))
                        })


## zip with all objects in memory -- more foolproof for recovery
lapply(names(LIM_simInitList), FUN = function(scenario) {
  zipSimList(LIM_simInitList[[scenario]],
             zipfile = file.path(outputPath(LIM_simInitList[[scenario]]), paste0("LIM_simInit_", scenario, "PART2.zip")),
             filename = file.path(outputPath(LIM_simInitList[[scenario]]), paste0("LIM_simInit_", scenario, "PART2.qs")),
             fileBackend = 2)
})

## save without pulling file backed rasters to memory - these are the simlists to be used for simulation
lapply(names(LIM_simInitList), FUN = function(scenario) {
  saveSimList(LIM_simInitList[[scenario]],
              filename = file.path(outputPath(LIM_simInitList[[scenario]]), paste0("LIM_simInit_", scenario, "PART2.qs")))
})

