## DIAGNOSE MODELBIOMASS ---------------------------------
## allFit is taking too long (>12h and didn't even with first optimizer)
library(lme4)
modBiomass <- LIM_simInitList$noPM$modelBiomass$mod

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


## DIAGNOSE FIRE IGNITION MODEL ----------------------------------
## to get predicted values we re-run the IgnitionPredict with the original data

library(data.table)
library(ggplot2)
library(ggspatial)
library(ggpubr)
library(SpaDES)

simDirName <- "mar2021Runs"
if (Sys.info()["nodename"] == "W-VIC-A127584") {
  simPaths <- list(cachePath = file.path("F:", basename(getwd()), "R/SpaDES/cache", simDirName, "noPM")
                   , modulePath = file.path("R/SpaDES/m")
                   , inputPath = file.path("R/SpaDES/inputs")
                   , outputPath = file.path("F:", basename(getwd()), "R/SpaDES/outputs", simDirName, "noPM"))
} else if (grepl("for-cast", Sys.info()["nodename"])) {
  simPaths <- list(cachePath = file.path("/mnt/scratch/cbarros", basename(getwd()), "R/SpaDES/cache", simDirName, "noPM")
                   , modulePath = file.path("R/SpaDES/m")
                   , inputPath = file.path("R/SpaDES/inputs")
                   , outputPath = file.path("R/SpaDES/outputs", simDirName, "noPM")
                   , rasterPath = file.path("/mnt/scratch/cbarros", basename(getwd()), "R/SpaDES/scratch/raster")
                   , scratchPath = file.path("/mnt/scratch/cbarros", basename(getwd()), "R/SpaDES/scratch"))
} else {
  simPaths <- list(cachePath = file.path("R/SpaDES/cache", simDirName, "noPM")
                   , modulePath = file.path("R/SpaDES/m")
                   , inputPath = file.path("R/SpaDES/inputs")
                   , outputPath = file.path("R/SpaDES/outputs", simDirName, "noPM"))
}
eventCaching <- c(".inputObjects", "init")


## old data
# oldSim <- qs::qread("C:/Users/cbarros/Desktop/oldDisp_noPM/preSimList.qs")

## get simList from Init - doesn't matter which scenario
LIM_simInitList <- loadSimList("R/SpaDES/outputs/mar2021Runs/noPM/LIM_simInit_noPM.qs")

parameters <- list(
  fireSense_IgnitionFit = list(
    "fireSense_ignitionFormula" = paste0("n_fires ~ coniferous:meanMDC + D2:meanMDC +",
                                         "M2:meanMDC + O1b:meanMDC + NF:meanMDC +",
                                         "coniferous:pw(meanMDC, k_conif) + D2:pw(meanMDC, k_D2) +",
                                         "M2:pw(meanMDC, k_M2) + O1b:pw(meanMDC, k_O1b) + NF:pw(meanMDC, k_NF) - 1")
    , "lb" = list(coef = 0,
                  knots = list("meanMDC" = 19))   ## the rounded 5% quantile, pre scaling
    , "ub" = list(coef = 20,
                  knots = list("meanMDC" = 21))   ## the rounded 80% quantile, pre scaling
    , "iterDEoptim" = 60
    , "iterNlminb" = 500
    , "family" = quote(MASS::negative.binomial(theta = 1, link = 'identity'))
    , "cores" = 4
    , "rescaleVars" = TRUE
    , "rescalers" = NULL
    , ".plots" = "png"
    , ".useCache" = eventCaching
  ),
  fireSense_IgnitionPredict = list(
    ".runInterval" = NA    ## only run once at the start
  ))

objects <- list(
  "ignitionFitRTM" = LIM_simInitList$ignitionFitRTM
  , "fireSense_ignitionCovariates" = LIM_simInitList$fireSense_ignitionCovariates
  , "fireSense_IgnitionAndEscapeCovariates" = LIM_simInitList$fireSense_IgnitionAndEscapeCovariates  ## FOR PRED AT 250M
)

simOutFireFreqPredVals <- simInitAndSpades(
  times = list(start = 0, end = 0)
  , modules = list("fireSense_IgnitionFit", "fireSense_IgnitionPredict")
  , paths = simPaths
  , objects = objects
  , params = parameters
  , cacheRepo = simPaths$cachePath
  , userTags = "fireSense_IgnitionFit&Predict"
  , omitArgs = c("userTags")
)

library(data.table)
library(ggplot2)
library(ggspatial)
library(ggpubr)

plot1 <- ggplot() +
  layer_spatial(data = simOutFireFreqPredVals$fireSense_IgnitionPredicted) +
  layer_spatial(data = LIM_simInitList$fireLocations, colour = "darkred") +
  annotation_north_arrow(style = north_arrow_minimal,
                         location = "tr", which_north = "true") +
  scale_fill_distiller(palette = "Greys", na.value = "transparent", direction = 1) +
  theme_pubr(margin = FALSE, legend = "right", base_size = 14) +
  labs(x = "longitude", y = "latitude", fill = expression(lambda))
plot1

## total no. fires (fitted)
fitted <- predictIgnition(simOutFireFreqPredVals$fireSense_IgnitionFitted[["formula"]][-2],
                          simOutFireFreqPredVals$fireSense_IgnitionFitted$data,
                          simOutFireFreqPredVals$fireSense_IgnitionFitted$coef,
                          1,
                          1,
                          simOutFireFreqPredVals$fireSense_IgnitionFitted$family$linkinv)
## fitted and observed
sum(fitted, na.rm = TRUE)
sum(simOutFireFreqPredVals$fireSense_IgnitionFitted$data$n_fires)

## average no. fires per year
fittedData <- data.table(simOutFireFreqPredVals$fireSense_IgnitionFitted$data, fittedVals = fitted)
mean(fittedData[, sum(fittedVals), by = year]$V1)
mean(fittedData[, sum(n_fires), by = year]$V1)
