## ------------------------------------------------------
## FIRE MODELLING WITH SpaDES
## TESTS
##
## Ceres: Nov 2017
## ------------------------------------------------------

## clean workspace
rm(list=ls()); amc::.gc()

## Get packages ----------------------
## requires as of Jan 2nd 2020
# loading reproducible     1.0.0.9010
# loading quickPlot        0.1.6.9000
# loading SpaDES.core      1.0.0.9
# loading SpaDES.tools     0.3.4.9000
# loading SpaDES.addins    0.1.2
# loading LandR            0.0.3.9006
# loading LandR.CS         0.0.1.0001

# devtools::install_github("PredictiveEcology/reproducible@messagingOverhaul", dependencies = FALSE)
# devtools::install_github("achubaty/amc@development", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/pemisc@development", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/map@development", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/LandR@reworkCohorts", dependencies = FALSE)
# devtools::install_github("ianmseddy/LandR.CS", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/quickPlot@development", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/SpaDES.tools@development", dependencies = FALSE)
# devtools::install_github("PredictiveEcology/SpaDES.core@lowMemory", dependencies = FALSE)

## test packages
# devtools::install_local("../LandR", dependencies = FALSE, force = TRUE)
# devtools::install_local("../reproducible", dependencies = FALSE, force = TRUE)
library(SpaDES)
library(LandR)

options("reproducible.useNewDigestAlgorithm" = TRUE)
options("spades.moduleCodeChecks" = FALSE)
options("reproducible.useCache" = TRUE)
options("reproducible.inputPaths" = normPath("R/SpaDES/inputs"))  ## store everything in data/ so that there are no duplicated files across modules
options("reproducible.destinationPath" = normPath("R/SpaDES/inputs"))

source("R/R_tools/Useful_functions.R")

## -----------------------------------------------
## SIMULATION SETUP
## -----------------------------------------------

## Set up simulation name  ---------------------------
# runName <- "blogSep2019_PM"
# runName <- "blogSep2019_noPM"

# runName <- "blogSep2019_PM_oneFire"
# runName <- "blogSep2019_noPM_oneFire"

# runName <- "PM_oneFire_newSppParams"
runName <- "noPM_oneFire_newSppParams"
eventCaching <- c(".inputObjects", "init")
useParallel <- FALSE

## paths define simulation paths
simPaths <- list(cachePath = file.path("R/SpaDES/cache/LIM_tests", runName),
                 modulePath = file.path("R/SpaDES/m"),
                 inputPath = file.path("R/SpaDES/inputs"),
                 outputPath = file.path("R/SpaDES/outputs", runName))

## Get necessary objects -----------------------
source("R/SpaDES/1_simObjects.R")

## Run Biomass_speciesData to get species layers
source("R/SpaDES/2_speciesLayers.R")

## maybe drop some species - Black spruce, and Ponderosa pine have v. few occurrences
plot(simOutSpeciesLayers$speciesLayers)
keepSpp <- setdiff(names(simOutSpeciesLayers$speciesLayers), c("Pice_mar", "Pinu_pon"))
simOutSpeciesLayers$speciesLayers <- subset(simOutSpeciesLayers$speciesLayers, keepSpp)
sppEquivalencies_CA <- sppEquivalencies_CA[LIM %in% keepSpp]
sppColorVect <- sppColorVect[keepSpp]

source("R/SpaDES/3_fireWeather.R")

## Define simulation params --------------------
simTimes <- list(start = 0, end = 65)
vegLeadingProportion <- 0 # indicates what proportion the stand must be in one species group for it to be leading.
# If all are below this, then it is a "mixed" stand
fireInitialTime <- 5L
fireTimestep <- if (grepl("oneFire", runName)) 100L else 2L
successionTimestep <- 1L

## Make simulation module list, parameters objects, objects and outputs accoding to run
## name and the parameters above
source("R/SpaDES/4_simulationSetup.R")

## -----------------------------------------------
## SIMULATION RUN
## -----------------------------------------------

# showCache(simPaths$cachePath, after = "2018-09-26 00:00:00")
# reproducible::clearCache(simPaths$cachePath, userTags = c("prepInputsLCC2005_rtm", "Biomass_borealDataPrep"))

## TODO CHANGE FIRE MODULES TO USE COHORT DATA RATHER THAN SUMMARY BMG OUTPUTS, LIKE BIOMASSMAP
graphics.off()

# reproducible::clearCache(simPaths$cachePath, userTags = c("^Biomass_core$", "init"), ask = FALSE)
## TODO RUN SIMUALTIONS W/ AND W/O PM for blog
# set.seed(524326)

## TODO: implement LANDIS pixel fire severity calculation:
## Each fire event has an associated mean fire severity which is the average of the severities at all of the event’s sites. (LANDIS-II DNFS v3)
# reproducible::clearCache(simPaths$cachePath, userTags = c("statsModel"))
# Biomass_core_testSim <- simInitAndSpades(times = simTimes
#                                  , params = simParams
#                                  , modules = simModules[1:6]
#                                  , objects = simObjects
#                                  , paths = simPaths
#                                  , outputs = outputs
#                                  , debug = TRUE
#                                  , .plotInitialTime = NA
# )
# saveRDS(Biomass_core_testSim, file.path(simPaths$outputPath, paste0("simList_", runName, ".rds")))
# simTimes$end <- 21

Biomass_core_testSimInit <- simInitAndSpades(times = simTimes
                                             , params = simParams
                                             , modules = simModules[c(1:5, 7)]
                                             , objects = simObjects
                                             , paths = simPaths
                                             , outputs = outputs
                                             , debug = TRUE
                                             # , .plotInitialTime = NA
)

## CHECK CONVERGENCE OF MODELBIOMASS
## allFit is taking too long (>12h and didn't even with first optimizer)
modBiomass <- Biomass_core_testSimInit$modelBiomass$mod

pars <- unlist(getME(modBiomass, c("theta")))
grad <- numDeriv::grad(update(modBiomass, devFunOnly = TRUE), pars)
hess <- numDeriv::hessian(update(modBiomass, devFunOnly = TRUE), pars)
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



saveSimList(Biomass_core_testSimOut,
            file.path(simPaths$outputPath, paste0("simList_fakeRstCurrentBurn", runName, ".RData")))
if (!is.na(plotInitialTime))
  dev.print(tiff, file.path(simPaths$outputPath, paste0("simPlots_", runName, ".tiff")),
            res = 300, units = "in")

