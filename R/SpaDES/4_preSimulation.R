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
    , "fitDeciduousCoverDiscount" = FALSE
    , "exportModels" = "all"
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
    , ".useCache" = eventCaching
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
preSimObjects <- list(# "studyArea" = foothillsSMALL
  # , "studyAreaLarge" = foothillsMED
  "studyArea" = foothills
  , "studyAreaLarge" = foothills
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

## IgnitionPredict can use finer scale rasters, as long as predictions are rescaled (P(sim$rescaleFactor)
fuelTypesCoverStkHR <- projectRaster(simOutPreSim$fuelTypesCoverStk, simOutPreSim$rasterToMatch,
                                     method = "bilinear")
weatherDataMDCStkHR <- projectRaster(simOutPreSim$weatherDataMDCStk, simOutPreSim$rasterToMatch,
                                     method = "bilinear")

## Define fireSense_IgnitionPredict module inputs
## we will use an average MDC across the years and output a raster that will not change
## in time
objects <- list(
  "fireSense_IgnitionFitted" = simOutPreSim$fireSense_IgnitionFitted
  , "D2" = fuelTypesCoverStkHR$D2
  , "M2" = fuelTypesCoverStkHR$M2
  , "coniferous" = fuelTypesCoverStkHR$coniferous
  , "O1b" = fuelTypesCoverStkHR$O1b
  , "NF" = fuelTypesCoverStkHR$NF
  , "julMDC" = mean(weatherDataMDCStkHR)
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
    "data" = c("D2", "M2", "coniferous", "O1b", "NF", "julMDC")
    , "modelObjName" = "fireSense_IgnitionFitted" # This is the default
    , "rescaleFactor" = res(fuelTypesCoverStkHR)[1]/res(simOutPreSim$fuelTypesCoverStk)[1] ^ 2
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
if (FALSE) {
  library(data.table)
  library(ggplot2)
  objects <- list(
    "fireSense_IgnitionFitted" = simOutPreSim$fireSense_IgnitionFitted
    , "dataFireSense_IgnitionFit" = simOutPreSim$dataFireSense_IgnitionFit
  )

  parameters <- list(
    fireSense_IgnitionPredict = list(
      "data" = "dataFireSense_IgnitionFit"
      , "modelObjName" = "fireSense_IgnitionFitted" # This is the default
      , "rescaleFactor" = 1
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
  ignitionsData <- simOutPreSim$dataFireSense_IgnitionFit
  ignitionsData[,  rows := 1:nrow(ignitionsData)]
  predVals <- data.table(rows = as.integer(names(simOutFireFreqPredVals$fireSense_IgnitionPredicted)),
                         fittedVals = simOutFireFreqPredVals$fireSense_IgnitionPredicted)
  ignitionsData <- predVals[ignitionsData, on = .(rows)]
  ignitionsData[, n_firesPred := rpois(.N, fittedVals)]

  plot1 <- ggplot(data = ignitionsData, aes(y = fittedVals, x = n_fires)) +
    geom_point() +
    labs(y = "lambda", x = "observed no. fires")

  plot2 <- ggplot(data = ignitionsData, aes(y = fittedVals, x = n_firesPred)) +
    geom_point() +
    labs(y = "lambda", x = "predicted (rpois) no. fires")

  plot3 <- ggplot(data = ignitionsData, aes(y = n_firesPred, x = n_fires)) +
    geom_point() +
    labs(y = "predicted (rpois) no. fires", x = "observed no. fires")

  plotData <- ignitionsData[, list(obsFires = sum(n_fires), predFires = sum(n_firesPred)),
                            by = year]
  plotData <- melt(plotData, id.var = "year")
  plot4 <- ggplot(data = plotData, aes(y = value, x = year, colour = variable)) +
    geom_line(size = 1) +
    scale_color_discrete(labels = c("obsFires" = "observed no. fires",
                                    "predFires" = "predicted (rpois) no. fires")) +
    theme(legend.position = "bottom") +
    labs(y = "no. fires", x = "year", colour = "")

  gridExtra::grid.arrange(plot1, plot2, plot3, plot4)

  plot(simOutFireFreq$fireSense_IgnitionPredicted)
  plot(simOutPreSim$fireLocations, add = TRUE)

}

## DIAGNOSE MODELBIOMASS ---------------------------------
## allFit is taking too long (>12h and didn't even with first optimizer)
if (FALSE) {
  library(lme4)
  modBiomass <- simOutPreSim$modelBiomass$mod

  pars <- unlist(getME(modBiomass, c("theta")))
  updateModBiomass <- update(modBiomass, devFunOnly = TRUE, data = modBiomass@frame)
  grad <- numDeriv::grad(updateModBiomass, pars)
  hess <- numDeriv::hessian(updateModBiomass, pars)
  sc_grad <- solve(hess, grad)

  if (length(modBiomass@optinfo$conv$lme4$messages) &
      max(pmin(abs(sc_grad), abs(grad))) > 0.001) {
    ## Nelder and L-BFGS-B methods didn't converge
    # , REML = FALSE,
    # control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    # ss <- getME(modBiomass, c("theta","fixef"))
    # modBiomass <- update(modBiomass, start = ss, data = modBiomass@frame,
    #                            control = lmerControl(calc.derivs = FALSE,  ## it's faster
    #                                                  optimizer = "optimx",
    #                                                  optCtrl = list(method='nlminb', maxit = 1e5)))  ## keeps hitting maxit...
    ncores <- detectCores()
    library(dfoptim)
    .specialData <<- modBiomass@frame
    diff_optims <- allFit(modBiomass$mod, parallel = 'multicore', ncpus = ncores)
    is.OK <- sapply(diff_optims, is, "merMod")  ## nlopt NELDERMEAD failed, others succeeded
    aa.OK <- diff_optims[is.OK]
    lapply(aa.OK,function(x) x@optinfo$conv$lme4$messages)

    pars <- unlist(getME(modBiomass, c("theta")))
    grad <- numDeriv::grad(update(modBiomass, devFunOnly = TRUE), pars)
    hess <- numDeriv::hessian(update(modBiomass, devFunOnly = TRUE), pars)
    sc_grad <- solve(hess, grad)
    max(pmin(abs(sc_grad),abs(grad)))
  }

}
