## --------------------------------------------------
##  HYPERVOLUMES OF PYRODIVERSITY AND BIODIVERSITY
##  Analyses of results
## --------------------------------------------------

if (Sys.which("make") == "") {
  stop("Please install and setup RTools 4.0")
}

if (!exists("pkgDir")) {
  pkgDir <- file.path(
    if (Sys.info()[["user"]] == "rstudio") "packages_docker" else "packages",
    version$platform,
    paste0(version$major, ".", strsplit(version$minor, "[.]")[[1]][1])
  )

  if (!dir.exists(pkgDir)) {
    dir.create(pkgDir, recursive = TRUE)
  }
  .libPaths(pkgDir)
}

library(Require)

Require(
  c(
    "cowplot (>= 1.2.0)",
    "data.table",
    "dplyr",
    "ggplot2 (>= 4.0.1.9000)",
    "ggpubr",
    # "LandR (==1.0.7.9026)", ## needed but don't load -- here for install.
    "marginaleffects",
    "mgcv",
    "multcomp",
    "nlme",
    "performance",
    "rvlenth/emmeans (>= 2.0.0)",
    "reproducible",
    "SpaDES",
    "ToolsCB"
  ),
  install = FALSE
)


source("R/R_tools/Useful_functions.R")
source("R/R_tools/hypervolumesHelpers.R")
source("R/R_tools/glhtMethods.R")
source("R/R_tools/utilsForResultsAnalyses.R")
options(
  "reproducible.useNewDigestAlgorithm" = 2,
  "spades.moduleCodeChecks" = FALSE,
  "reproducible.useCache" = TRUE,
  "reproducible.destinationPath" = normPath("R/SpaDES/inputs"),
  "reproducible.useGDAL" = FALSE,
  "LandR.assertions" = TRUE
)

## general paths
# simDirName <- "mar2021Runs"
simDirName <- "mar2022Runs"
if (
  Sys.info()["nodename"] == "W-VIC-A127584" |
    Sys.info()["nodename"] == "L-VIC-A155348"
) {
  simPaths <- list(
    cachePath = file.path(
      "F:",
      basename(getwd()),
      "R/SpaDES/cache",
      simDirName,
      "postSimAnalyses"
    ),
    modulePath = file.path("R/SpaDES/m"),
    inputPath = file.path(
      "F:",
      basename(getwd()),
      "R/SpaDES/inputs",
      simDirName
    ),
    outputPath = file.path(
      "F:",
      basename(getwd()),
      "R/SpaDES/outputs",
      simDirName
    )
  )
} else if (grepl("for-cast", Sys.info()["nodename"])) {
  simPaths <- list(
    cachePath = file.path("R/SpaDES/cache", simDirName, "postSimAnalyses"),
    modulePath = file.path("R/SpaDES/m"),
    inputPath = file.path("R/SpaDES/inputs"),
    outputPath = file.path("R/SpaDES/outputs", simDirName),
    rasterPath = file.path(
      "/mnt/scratch/cbarros",
      basename(getwd()),
      "R/SpaDES/scratch/raster"
    ),
    scratchPath = file.path(
      "/mnt/scratch/cbarros",
      basename(getwd()),
      "R/SpaDES/scratch"
    )
  )
} else {
  simPaths <- list(
    cachePath = file.path("R/SpaDES/cache", simDirName, "postSimAnalyses"),
    modulePath = file.path("R/SpaDES/m"),
    inputPath = file.path("R/SpaDES/inputs"),
    outputPath = file.path("R/SpaDES/outputs", simDirName)
  )
}

if (grepl("for-cast|dc709f2508f1|45eafed436c8", Sys.info()["nodename"])) {
  data.table::setDTthreads(5)
  options(bitmapType = "cairo")
}


## path to figure folder and cache folder
figOutputPath <- file.path(simPaths$outputPath, "figuresAnalysis/Montane")
statsOutputPath <- file.path(simPaths$outputPath, "statsAnalysis")
HVoutputPath <- file.path(simPaths$outputPath, "hypervolumes")

## are we using the merged douglas-fir/dry-conifer stands?
mergeDMCPSME <- FALSE ## merge DMCPSME PSME dryPSME
mergePSME <- FALSE ## merge PSME dryPSME  -- 2026 not merging anymore
if (mergeDMCPSME) {
  HVoutputPathMergedVegType <- file.path(
    simPaths$outputPath,
    "hypervolumes",
    "mergeDMCPSME"
  )
  figOutputPath <- file.path(
    simPaths$outputPath,
    "figuresAnalysis/Montane/mergeDMCPSME"
  )
  statsOutputPath <- file.path(
    simPaths$outputPath,
    "statsAnalysis/mergeDMCPSME"
  )
}

if (mergePSME) {
  HVoutputPathMergedVegType <- file.path(
    simPaths$outputPath,
    "hypervolumes",
    "mergePSME"
  )
  figOutputPath <- file.path(
    simPaths$outputPath,
    "figuresAnalysis/Montane/mergePSME"
  )
  statsOutputPath <- file.path(simPaths$outputPath, "statsAnalysis/mergePSME")
}

dir.create(figOutputPath, recursive = TRUE)
dir.create(statsOutputPath, recursive = TRUE)

## PREP FIRE AND VEG DATA -------------------
yearSubset <- unique(as.integer(c(seq(3511, 4011, 5), 4011)))
runPrepResultsModule <- FALSE
source("R/R_tools/simResultsDataPrep.R")

opts <- options(reproducible.cachePath = simPaths$cachePath)
source("R/R_tools/prepFireData4HVs.R")

yearSamples <- sample5SimYears(allPixelCohortDataMnt[, .(year, rep)]) ## seed ensures same years are drawn
useFirstLastYear <- FALSE
source("R/R_tools/prepVegData4HVs.R")
options(opts)

## don't need these
rm(allPixelBurnData, allPixelCohortData, allPixelCohortDataMnt)
gc(reset = TRUE)

## LOAD HYPERVOLUMES RESULTS ----------------------
source("R/R_tools/prepHVData.R")

## LABELS AND COLOURS FOR PLOTTING ----------------
source("R/R_tools/plotLabels&Cols.R")
rm(preSimList) ## not needed anymore


## STATISTICS: EFFECT OF SCENARIO ON PYRODIVERSITY  -----------------------
## make separate datatable for stats so that we can change contrasts
modelData <- allHVData[
  year == max(yearSubset) & HVtype == "fireHV",
  .(
    HVtype,
    HV_noPM,
    HV_PM,
    Intersection,
    MinDist,
    CentroidDist,
    rep,
    repHV,
    vegType
  )
]

modelData <- melt.data.table(
  modelData,
  measure.vars = c("HV_noPM", "HV_PM"),
  variable.name = "scenario",
  value.name = "Volume"
)

modelData[, scenario := relevel(scenario, "HV_noPM")]
modelData[, vegType := relevel(vegType, "No veg.")]
modelData[, logVolume := log(Volume)]

if (
  isFALSE(
    nrow(modelData[, unique(repHV), by = .(vegType, rep)]) ==
      length(unique(modelData$vegType)) *
        length(unique(modelData$rep)) *
        length(unique(modelData$repHV))
  )
) {
  stop("Missing replicates")
}

### Landscape level ------
fireHVVolumeLandscape.lm <- lm(
  Volume ~ scenario,
  modelData[vegType == "landscape"]
)
## the data seems very dispersed - some reps are extreme outliers - logging helps -- Jul 14 2022 this is no longer true with 4000 years
hist(modelData[vegType == "landscape", Volume], breaks = 1000)
fireHVVolumeLandscape.lm2 <- lm(
  logVolume ~ scenario,
  modelData[vegType == "landscape"]
) ## a bit better even with 4000 years
summary(fireHVVolumeLandscape.lm)

check_model(fireHVVolumeLandscape.lm) ## looks better
check_model(fireHVVolumeLandscape.lm2)

caret::RMSE(fitted(fireHVVolumeLandscape.lm), modelData$Volume)
caret::RMSE(exp(fitted(fireHVVolumeLandscape.lm2)), exp(modelData$logVolume)) ## only marginally better, slightly worse resids


check_normality(fireHVVolumeLandscape.lm) ## it's ok
check_heteroscedasticity(fireHVVolumeLandscape.lm) ## it's ok

## March 2026: no indication that the model's below are necessary
fireHVVolumeLandscape.gls <- gls(
  Volume ~ scenario,
  weights = varIdent(form = ~ 1 | scenario),
  modelData[vegType == "landscape"]
)
fireHVVolumeLandscape.gls2 <- gls(
  logVolume ~ scenario,
  weights = varIdent(form = ~ 1 | scenario),
  modelData[vegType == "landscape"]
)
## no big improvements
check_model(fireHVVolumeLandscape.gls)
check_model(fireHVVolumeLandscape.gls2) # slightly worse
caret::RMSE(fitted(fireHVVolumeLandscape.lm), modelData$Volume)
caret::RMSE(exp(fitted(fireHVVolumeLandscape.lm2)), exp(modelData$logVolume))
caret::RMSE(fitted(fireHVVolumeLandscape.gls), modelData$Volume)
caret::RMSE(exp(fitted(fireHVVolumeLandscape.gls2)), exp(modelData$logVolume)) ## only marginally better, slightly worse resids

## can't directly compare transformed and non-transformed AICs
AIC(
  fireHVVolumeLandscape.lm, ## better
  fireHVVolumeLandscape.gls
)
AIC(
  fireHVVolumeLandscape.lm2, ## better
  fireHVVolumeLandscape.gls2
)


png(file.path(figOutputPath, "fireHVVolumeLandscapelmRESIDUALS.png"))
## a bit heteroskedastic, but not that much
check_model(fireHVVolumeLandscape.lm)
dev.off()

sink(file.path(statsOutputPath, "fireHVVolumeLandscapelmSUMMARY.txt"))
anova(fireHVVolumeLandscape.lm)
cat("\n*********************\n")
summary(fireHVVolumeLandscape.lm)
cat("\n*********************\n")
emmeans(fireHVVolumeLandscape.lm, specs = "scenario")
sink()

### vegType level --------
fireHVVolumeVegTypes.lm <- lm(
  Volume ~ scenario * vegType,
  modelData[vegType != "landscape"]
)
## same issue as before
fireHVVolumeVegTypes.lm2 <- lm(
  logVolume ~ scenario * vegType,
  modelData[vegType != "landscape"]
)

check_model(fireHVVolumeVegTypes.lm) ## (maybe?) better on homo. variance (Nov 2025)
check_model(fireHVVolumeVegTypes.lm2)

check_normality(fireHVVolumeVegTypes.lm) ## not ok
check_heteroscedasticity(fireHVVolumeVegTypes.lm) ## ok

check_normality(fireHVVolumeVegTypes.lm2) ## not ok
check_heteroscedasticity(fireHVVolumeVegTypes.lm2) ## not ok

caret::RMSE(
  fitted(fireHVVolumeVegTypes.lm),
  modelData[vegType != "landscape", Volume]
) ## marginally better
caret::RMSE(
  exp(fitted(fireHVVolumeVegTypes.lm2)),
  exp(modelData[vegType != "landscape", logVolume])
)

## glht can't see the missing levels, so we fit the model on a separate data table
modelData2 <- as.data.frame(modelData[vegType != "landscape"])
modelData2$scenario <- factor(modelData2$scenario)
modelData2$vegType <- factor(modelData2$vegType)
modelData2$scenario <- relevel(modelData2$scenario, "HV_noPM")
modelData2$vegType <- relevel(modelData2$vegType, "No veg.")

fireHVVolumeVegTypes.gls <- gls(
  Volume ~ scenario * vegType,
  weights = varIdent(form = ~ 1 | scenario * vegType),
  modelData2
)
check_model(fireHVVolumeVegTypes.gls) ## looks better
anova(fireHVVolumeVegTypes.gls) ## March 2026: scen not significant; after not merging PSME+dryPSME all effects significant
summary(fireHVVolumeVegTypes.gls)

## no longer relevant if not merging PSME+dryPSME
# fireHVVolumeVegTypes.gls2 <- gls(Volume ~ scenario + vegType,
#                                  weights = varIdent(form = ~ 1 | scenario * vegType),
#                                  data = modelData2)
# anova(fireHVVolumeVegTypes.gls, fireHVVolumeVegTypes.gls2)   ## interaction model is much better.

AICcmodavg::AICc(fireHVVolumeVegTypes.lm)
AICcmodavg::AICc(fireHVVolumeVegTypes.gls) ## best
# AICcmodavg::AICc(fireHVVolumeVegTypes.gls2)

caret::RMSE(
  fitted(fireHVVolumeVegTypes.lm),
  modelData[vegType != "landscape", Volume]
)
caret::RMSE(fitted(fireHVVolumeVegTypes.gls), modelData2$Volume) ## same as above
# caret::RMSE(fitted(fireHVVolumeVegTypes.gls2), modelData2$Volume)   ## worse

png(file.path(figOutputPath, "fireHVVolumeVegTypesglsRESIDUALS.png"))
check_model(fireHVVolumeVegTypes.gls)
dev.off()

sink(file.path(statsOutputPath, "fireHVVolumeVegTypesglsSUMMARY.txt"))
anova(fireHVVolumeVegTypes.gls)
cat("\n*********************\n")
summary(fireHVVolumeVegTypes.gls)
cat("\n*********************\n")
## tukey contrasts for scenarios within vegTypes
emmeans(
  fireHVVolumeVegTypes.gls,
  pairwise ~ scenario | vegType,
  adjust = "tukey"
)
sink()


## ------------------------------------------------------------------------
## STATISTICS: EFFECT OF SCENARIO ON BIODIVERSITY  ------------------------

### how did scenario affect overall diversity at the end of simulation? ----
## (i.e. HV size)
modelData <- allHVData[
  year == max(yearSubset) & HVtype == "vegHV",
  .(
    HVtype,
    HV_noPM,
    HV_PM,
    Intersection,
    MinDist,
    CentroidDist,
    rep,
    repHV,
    vegType
  )
]
modelData <- melt.data.table(
  modelData,
  measure.vars = c("HV_noPM", "HV_PM"),
  variable.name = "scenario",
  value.name = "Volume"
)

modelData[, scenario := relevel(scenario, "HV_noPM")]
modelData[, vegType := relevel(vegType, "No veg.")]
modelData[, logVolume := log(Volume)]

if (
  isFALSE(
    nrow(modelData[, unique(repHV), by = .(vegType, rep)]) ==
      length(unique(modelData$vegType)) *
        length(unique(modelData$rep)) *
        length(unique(modelData$repHV))
  )
) {
  stop("Missing replicates")
}

### Landscape level ------
modelData2 <- modelData[vegType == "landscape"]
## need to fit the model after dropping unused levels, otherwise `marginaleffects::predictions` will error
## when trying to get prediction intervals
cols <- names(modelData2)
modelData2[,
  (cols) := lapply(.SD, function(x) {
    if (is.factor(x)) droplevels(x) else x
  })
]

vegHVVolumeLandscape.lm <- lm(Volume ~ scenario, modelData2)
vegHVVolumeLandscape.lm2 <- lm(logVolume ~ scenario, modelData2)

caret::RMSE(fitted(vegHVVolumeLandscape.lm), modelData2$Volume)
caret::RMSE(exp(fitted(vegHVVolumeLandscape.lm2)), exp(modelData2$logVolume)) ## a little worse

check_model(vegHVVolumeLandscape.lm)
check_model(vegHVVolumeLandscape.lm2) ## better

check_normality(vegHVVolumeLandscape.lm) ##  not ok
check_heteroscedasticity(vegHVVolumeLandscape.lm) ##  not ok
check_normality(vegHVVolumeLandscape.lm2) ##  not ok
check_heteroscedasticity(vegHVVolumeLandscape.lm2) ## ok


vegHVVolumeLandscape.gls <- gls(
  Volume ~ scenario,
  weights = varIdent(form = ~ 1 | scenario),
  data = modelData2
)

vegHVVolumeLandscape.gls2 <- gls(
  logVolume ~ scenario,
  weights = varIdent(form = ~ 1 | scenario),
  data = modelData2
)

AICcmodavg::AICc(vegHVVolumeLandscape.lm)
AICcmodavg::AICc(vegHVVolumeLandscape.gls) ## better
check_model(vegHVVolumeLandscape.gls) ## March 2026: looks worse.
caret::RMSE(fitted(vegHVVolumeLandscape.lm), modelData2$Volume)
caret::RMSE(fitted(vegHVVolumeLandscape.gls), modelData2$Volume) ## identical

AICcmodavg::AICc(vegHVVolumeLandscape.lm2)
AICcmodavg::AICc(vegHVVolumeLandscape.gls2) ## worse
check_model(vegHVVolumeLandscape.gls2) ## equally bad as non-logged gls
caret::RMSE(fitted(vegHVVolumeLandscape.gls), modelData2$Volume)
caret::RMSE(exp(fitted(vegHVVolumeLandscape.gls2)), modelData2$Volume) ## slightly worse


png(file.path(figOutputPath, "vegHVVolumeLandscapelmRESIDUALS.png"))
check_model(vegHVVolumeLandscape.lm2)
dev.off()

sink(file.path(statsOutputPath, "vegHVVolumeLandscapelmSUMMARY.txt"))
anova(vegHVVolumeLandscape.lm2)
cat("\n*********************\n")
summary(vegHVVolumeLandscape.lm2)
cat("\n*********************\n")
emmeans(vegHVVolumeLandscape.lm2, specs = "scenario", data = modelData2)
sink()

### vegType level --------
modelData2 <- modelData[vegType != "landscape"]
## need to fit the model after dropping unused levels, otherwise `marginaleffects::predictions` will error
## when trying to get prediction intervals
cols <- names(modelData2)
modelData2[,
  (cols) := lapply(.SD, function(x) {
    if (is.factor(x)) droplevels(x) else x
  })
]

vegHVVolumeVegTypes.lm <- lm(Volume ~ scenario * vegType, modelData2)
vegHVVolumeVegTypes.lm2 <- lm(logVolume ~ scenario * vegType, modelData2)

anova(vegHVVolumeVegTypes.lm) ## interaction not sign.
anova(vegHVVolumeVegTypes.lm2) ## all sign.

check_normality(vegHVVolumeVegTypes.lm) ## not ok
check_heteroscedasticity(vegHVVolumeVegTypes.lm) ## not ok
check_normality(vegHVVolumeVegTypes.lm2) ## not ok
check_heteroscedasticity(vegHVVolumeVegTypes.lm2) ## ok

## some heteroscedasticity
check_model(vegHVVolumeVegTypes.lm) ## terrible
check_model(vegHVVolumeVegTypes.lm2) ## better
caret::RMSE(fitted(vegHVVolumeVegTypes.lm), modelData2$Volume)
caret::RMSE(exp(fitted(vegHVVolumeVegTypes.lm2)), modelData2$Volume) ## worse

modelData2 <- as.data.frame(modelData[vegType != "landscape"])
modelData2$scenario <- factor(modelData2$scenario)
modelData2$vegType <- factor(modelData2$vegType)
modelData2$scenario <- relevel(modelData2$scenario, "HV_noPM")
modelData2$vegType <- relevel(modelData2$vegType, "No veg.")

vegHVVolumeVegTypes.gls <- gls(
  Volume ~ scenario * vegType,
  weights = varIdent(form = ~ 1 | scenario * vegType),
  data = modelData2
)
vegHVVolumeVegTypes.gls2 <- gls(
  logVolume ~ scenario * vegType,
  weights = varIdent(form = ~ 1 | scenario * vegType),
  data = modelData2
)

check_model(vegHVVolumeVegTypes.gls)
check_model(vegHVVolumeVegTypes.gls2) ## better than lm2

anova(vegHVVolumeVegTypes.gls2) ## all sign.

png(file.path(figOutputPath, "vegHVVolumeVegTypesglsRESIDUALS.png"))
check_model(vegHVVolumeVegTypes.gls2) ## better
dev.off()

sink(file.path(statsOutputPath, "vegHVVolumeVegTypesglsSUMMARY.txt"))
anova(vegHVVolumeVegTypes.gls2)
cat("\n*********************\n")
summary(vegHVVolumeVegTypes.gls2)
cat("\n*********************\n")
emmeans(
  vegHVVolumeVegTypes.gls2,
  pairwise ~ scenario | vegType,
  adjust = "tukey",
  deata = modelData2
) ## tukey contrasts for scenarios within vegTypes
cat("\n*********************\n")
sink()

### ------------------------------------------------------------------------
### How much did communities change in time? ---------------------------------------
### as in Barros et al 2016, calculate the prop overlap
if (useFirstLastYear) {
  modelData <- allHVData[
    HVtype == "vegHV" &
      compare == paste(min(yearSubset), max(yearSubset), sep = "_"),
    .(overlap, MinDist, CentroidDist, rep, repHV, vegType, scenario)
  ]
  modelData[, scenario := relevel(factor(scenario), "noPM")]
  modelData[, vegType := relevel(vegType, "No veg.")]

  ## across the landscape
  vegHVOverlapLandscape.lm <- lm(
    overlap ~ scenario,
    modelData[vegType == "landscape"]
  )

  ## some heteroscedasticity, but not lots and gls is no big improvement and results are similar
  sets <- par(mfrow = c(2, 2))
  plot(vegHVOverlapLandscape.lm)
  par(sets)

  vegHVOverlapLandscape.gls <- gls(
    overlap ~ scenario,
    weights = varIdent(form = ~ 1 | scenario),
    data = modelData[vegType == "landscape"]
  )
  AIC(vegHVOverlapLandscape.lm, vegHVOverlapLandscape.gls) ## gls is worse and doens't improve residuals much

  png(file.path(figOutputPath, "vegHVOverlapLandscapelmRESIDUALS.png"))
  sets <- par(mfrow = c(2, 2))
  plot(vegHVOverlapLandscape.lm) ## still not great
  par(sets)
  dev.off()

  sink(file.path(statsOutputPath, "vegHVOverlapLandscapelmSUMMARY.txt"))
  anova(vegHVOverlapLandscape.lm)
  cat("\n*********************\n")
  summary(vegHVOverlapLandscape.lm)
  cat("\n*********************\n")
  emmeans(vegHVOverlapLandscape.lm, specs = "scenario")
  sink()

  ## by vegType
  vegHVOverlapVegTypes.lm <- lm(
    overlap ~ scenario * vegType,
    modelData[vegType != "landscape"]
  )
  ## heteroscedasticity
  sets <- par(mfrow = c(2, 2))
  plot(vegHVOverlapVegTypes.lm) ## not too bad for an anova
  par(sets)

  modelData2 <- as.data.frame(modelData[vegType != "landscape"])
  modelData2$scenario <- factor(modelData2$scenario)
  modelData2$vegType <- factor(modelData2$vegType)
  modelData2$scenario <- relevel(modelData2$scenario, "HV_noPM")
  modelData2$vegType <- relevel(modelData2$vegType, "No veg.")

  vegHVOverlapVegTypes.gls <- gls(
    overlap ~ scenario * vegType,
    weights = varIdent(form = ~ 1 | scenario * vegType),
    modelData2
  ) ## can't converge

  ## can't fix de variance on vegtype because some levels have 0 variance
  # vegHVOverlapVegTypes.gls <- gls(overlap ~ scenario * vegType,
  #                                 weights = varIdent(form = ~ 1 | scenario),
  #                                 modelData2)
  AIC(vegHVOverlapVegTypes.lm, vegHVOverlapVegTypes.gls)

  png(file.path(figOutputPath, "vegHVOvelapVegTypesglsRESIDUALS.png"))
  plot(vegHVOverlapVegTypes.gls) ## a bit better
  dev.off()

  sink(file.path(statsOutputPath, "vegHVOverlapVegTypesglsSUMMARY.txt"))
  anova(vegHVOverlapVegTypes.gls)
  cat("\n*********************\n")
  summary(vegHVOverlapVegTypes.gls)
  cat("\n*********************\n")
  emmeans(
    vegHVOverlapVegTypes.gls,
    pairwise ~ scenario | vegType,
    adjust = "tukey"
  ) ## tukey contrasts for scenarios within vegTypes
  cat("\n*********************\n")
  sink()
}


## ------------------------------------------------------------------------
## STATISTICS: EFFECT OF SCENARIO ON PYRO-BIODIVERSITY RELATIONSHIP  ------
modelData <- allHVData[
  year == max(yearSubset),
  .(HVtype, HV_noPM, HV_PM, rep, repHV, vegType)
]
modelData <- melt.data.table(
  modelData,
  measure.vars = c("HV_noPM", "HV_PM"),
  variable.name = "scenario",
  value.name = "Volume"
)

modelData <- dcast.data.table(modelData, ... ~ HVtype, value.var = "Volume")
modelData[, `:=`(logFireHV = log(fireHV), logVegHV = log(vegHV))]
modelData[, scenario := relevel(scenario, "HV_noPM")]
modelData[, vegType := relevel(vegType, "No veg.")]
modelData2 <- modelData[vegType == "landscape"]
## need to fit the model after dropping unused levels, otherwise `marginaleffects::predictions` will error
## when trying to get prediction intervals
cols <- names(modelData2)
modelData2[,
  (cols) := lapply(.SD, function(x) {
    if (is.factor(x)) droplevels(x) else x
  })
]

### Landscape level ------
pyroVSbiodiversityLandscape.lm <- lm(
  vegHV ~ fireHV * scenario,
  data = modelData2
)
## the data seems very dispersed for fire HVs - some reps are extreme outliers - logging helps -- no longer true Nov 2025
hist(modelData[vegType == "landscape", fireHV], breaks = 1000)
# hist(modelData[vegType == "landscape" & fireHV > 500, fireHV], breaks = 1000)   # no longer exists in 4000 years sims

pyroVSbiodiversityLandscape.lm2 <- lm(
  logVegHV ~ logFireHV * scenario,
  data = modelData2
)
check_model(pyroVSbiodiversityLandscape.lm)
check_model(pyroVSbiodiversityLandscape.lm2) ## a bit better

## no strong support for log
check_normality(pyroVSbiodiversityLandscape.lm) ## not ok
check_heteroscedasticity(pyroVSbiodiversityLandscape.lm) ## ok
check_normality(pyroVSbiodiversityLandscape.lm2) ## not ok
check_heteroscedasticity(pyroVSbiodiversityLandscape.lm2) ## ok


## in both: scenario main effect not significant, but interaction is (but close)
anova(pyroVSbiodiversityLandscape.lm)
anova(pyroVSbiodiversityLandscape.lm2)

## if not merging dryPSME + PSME, the following makes no sense bc interaction is signficant.
# pyroVSbiodiversityLandscape.lm3 <- lm(logVegHV ~ logFireHV + scenario, data = modelData2)
# anova(pyroVSbiodiversityLandscape.lm3)   ## scenario still not significant
# check_model(pyroVSbiodiversityLandscape.lm3)  ## probably as bad.
# anova(pyroVSbiodiversityLandscape.lm2, pyroVSbiodiversityLandscape.lm3) ##  sign., lm3 probably worse.

## March 2026: scenario effect/interaction is significant (using sevPropB as HV component now), so the model below
## no longer makes sense
# pyroVSbiodiversityLandscape.lm4 <- lm(logVegHV ~ logFireHV, data = modelData2)
# anova(pyroVSbiodiversityLandscape.lm4) ## scenarion effect not significant
# check_model(pyroVSbiodiversityLandscape.lm4)  ## probably as bad.
# anova(pyroVSbiodiversityLandscape.lm3, pyroVSbiodiversityLandscape.lm4) ## not sign. different (AIC is marginally better, see below)

## try quadratic (as per Steel et al 2021 and He et al 2019)
## note that by default, we're using orthogonal polinomials - AIC is the same, so are residuals
## coefficients differ
## due to a bug in nlme, gls doens't output the estimated polynomial (using poly) and therefore
## predictions with new data (necessary for emmeans) are wrong (see https://stackoverflow.com/questions/70746067/issues-predicting-with-nlmeglsquadratic-model-fitted-with-poly-2)
## instead of poly we now centre the data and use logFireHV + I(logFireHV^2) -- given the low degree polynomial this is okay.

## the bug was fixed

# modelData[, logFireHVcenter := scale(logFireHV, center = TRUE, scale = FALSE),
#           by = .(scenario, vegType)]
# pyroVSbiodiversityLandscape.lm3 <- lm(logVegHV ~ (logFireHVcenter + I(logFireHVcenter^2))*scenario,
#                                       data = modelData2)
# pyroVSbiodiversityLandscape.lm5 <- lm(logVegHV ~ poly(logFireHV, 2), data = modelData2)
# check_model(pyroVSbiodiversityLandscape.lm5)  ## as bad as previous
# anova(pyroVSbiodiversityLandscape.lm3, pyroVSbiodiversityLandscape.lm5)  ## not sign. different, AIC slightly worse

pyroVSbiodiversityLandscape.lm6 <- lm(
  vegHV ~ poly(fireHV, 2) * scenario,
  data = modelData2
)
anova(pyroVSbiodiversityLandscape.lm6)
anova(pyroVSbiodiversityLandscape.lm, pyroVSbiodiversityLandscape.lm6) # all significant

pyroVSbiodiversityLandscape.lm7 <- lm(
  logVegHV ~ poly(logFireHV, 2) * scenario,
  data = modelData2
)
anova(pyroVSbiodiversityLandscape.lm7) # only logFireHV is significant

check_model(pyroVSbiodiversityLandscape.lm6) ## better than lm, looks a bit better than lm7.
check_model(pyroVSbiodiversityLandscape.lm7) ## better than lm2

check_normality(pyroVSbiodiversityLandscape.lm6) ## not ok
check_heteroscedasticity(pyroVSbiodiversityLandscape.lm6) ## ok
check_normality(pyroVSbiodiversityLandscape.lm7) ## ok
check_heteroscedasticity(pyroVSbiodiversityLandscape.lm7) ## ok


## nov 2025, no support for the following model
# pyroVSbiodiversityLandscape.lm8 <- lm(logVegHV ~ poly(logFireHV, 2)+scenario,
#                                       data = modelData2)
# anova(pyroVSbiodiversityLandscape.lm8)  # only logFireHV is significant

AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm) ## can't compare to log models
AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm6) ## better than lm
AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm2)
AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm7) ## worse than lm2
# AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm3)   ## second best
# AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm4)
# AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm5)   ## best, but very marginally
# AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm8)

## and if we allow the variance to change?
pyroVSbiodiversityLandscape.gls <- nlme::gls(
  vegHV ~ fireHV * scenario,
  weights = nlme::varIdent(form = ~ 1 | scenario),
  data = modelData2
)
anova(pyroVSbiodiversityLandscape.gls) ## scenario main effect not significant

pyroVSbiodiversityLandscape.gls2 <- nlme::gls(
  logVegHV ~ logFireHV * scenario,
  weights = nlme::varIdent(form = ~ 1 | scenario),
  data = modelData2
)
anova(pyroVSbiodiversityLandscape.gls2) ## scenario main effect not significant

## if not merging dryPSME + PSME, the following makes no sense bc interaction is significant.
# pyroVSbiodiversityLandscape.gls3 <- nlme::gls(logVegHV ~ logFireHV + scenario,
#                                               weights = nlme::varIdent(form = ~ 1 | scenario),
#                                               data = modelData2)
# anova(pyroVSbiodiversityLandscape.gls3)
## March 2026: no support for following model
# pyroVSbiodiversityLandscape.gls3 <- nlme::gls(logVegHV ~ logFireHV,
#                                               weights = nlme::varIdent(form = ~ 1 | scenario),
#                                               data = modelData2)
# anova(pyroVSbiodiversityLandscape.gls3)

check_model(pyroVSbiodiversityLandscape.gls) ## still pretty bad
check_model(pyroVSbiodiversityLandscape.gls2) ## looks better on normality side

## try adding polinomial again
pyroVSbiodiversityLandscape.gls4 <- nlme::gls(
  logVegHV ~ poly(logFireHV, 2),
  weights = nlme::varIdent(form = ~ 1 | scenario),
  data = modelData2
)
pyroVSbiodiversityLandscape.gls5 <- nlme::gls(
  logVegHV ~ poly(logFireHV, 2) * scenario,
  weights = nlme::varIdent(form = ~ 1 | scenario),
  data = modelData2
)
pyroVSbiodiversityLandscape.gls6 <- nlme::gls(
  logVegHV ~ poly(logFireHV, 2) + scenario,
  weights = nlme::varIdent(form = ~ 1 | scenario),
  data = modelData2
)

## non log options
AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm)
AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm6)
AICcmodavg::AICc(pyroVSbiodiversityLandscape.gls) ## best

## log options
AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm2)
AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm7)
AICcmodavg::AICc(pyroVSbiodiversityLandscape.gls2)
AICcmodavg::AICc(pyroVSbiodiversityLandscape.gls5) ## best, and better residuals than non-log

# AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm3)   ## third best
# AICcmodavg::AICc(pyroVSbiodiversityLandscape.gls3)
# AICcmodavg::AICc(pyroVSbiodiversityLandscape.gls4)
# AICcmodavg::AICc(pyroVSbiodiversityLandscape.gls6)

anova(pyroVSbiodiversityLandscape.gls2, pyroVSbiodiversityLandscape.lm2) ## sign different, gls marginally worse
anova(pyroVSbiodiversityLandscape.gls2, pyroVSbiodiversityLandscape.gls5) ## sign different
anova(pyroVSbiodiversityLandscape.gls)
anova(pyroVSbiodiversityLandscape.gls5) ## interaction and scenario no longer significant

check_model(pyroVSbiodiversityLandscape.gls5) ## looks worse than gls2

## Nov 2025: in PM, replicate 4 seems to be an outlier, without it the relationship would be flatter
## also, scenario does seem to have an impact on mean volumes and the slopes.
## I think the GLS may be more appropriate given the outlier rep.
## I also don't see a strong support for polynomial relationships, especially removing the rep outlier in PM
## I think LM2/GLS2 are good options; gls2 didn't things improve much

plotData <- modelData2
plotData[, predLM2 := predict(pyroVSbiodiversityLandscape.lm2)]
plotData[, predLM7 := predict(pyroVSbiodiversityLandscape.lm7)]
plotData[, predGLS2 := predict(pyroVSbiodiversityLandscape.gls2)]
plotData[, predGLS5 := predict(pyroVSbiodiversityLandscape.gls5)]
ggplot(
  # plotData[scenario != "HV_noPM" & rep != 4]
  # plotData[scenario == "HV_noPM"]
  plotData,
  aes(x = logFireHV)
) +
  geom_line(aes(y = predLM2, linetype = scenario), col = "red") +
  geom_line(aes(y = predLM7, linetype = scenario), col = "green") +
  # geom_line(aes(y = predGLS2, linetype = scenario), col = "orange") + ## overlaps predLM2 completely.
  # geom_line(aes(y = predGLS5, linetype = scenario), col = "blue") +   ## overlaps predLM7 completely.
  # geom_smooth(aes(y = logVegHV, linetype = scenario), formula = y ~ poly(x, 2),
  #             colour = "red",
  #             method = "lm", se = FALSE) +
  # geom_smooth(aes(y = logVegHV, linetype = scenario), formula = y ~ x,
  #             colour = "blue",
  #             method = "lm", se = FALSE) +
  # geom_smooth(aes(y = logVegHV), formula = y ~ poly(x, 2),
  #             colour = "red",
  #             method = "lm", se = FALSE) +
  # geom_smooth(aes(y = logVegHV), formula = y ~ x,
  #             colour = "blue",
  #             method = "lm", se = FALSE) +
  geom_point(aes(y = logVegHV, colour = as.factor(rep), shape = scenario))


png(file.path(figOutputPath, "pyroVSbiodiversityLandscapelmRESIDUALS.png"))
plot(pyroVSbiodiversityLandscape.lm2)
dev.off()

sink(file.path(statsOutputPath, "pyroVSbiodiversityLandscapelmSUMMARY.txt"))
anova(pyroVSbiodiversityLandscape.lm2)
cat("\n*********************\n")
summary(pyroVSbiodiversityLandscape.lm2)
cat("\n*********************\n")
emtrends(
  pyroVSbiodiversityLandscape.lm2,
  specs = c("scenario"),
  var = "logFireHV",
  data = modelData2
)
sink()

### vegType level --------
modelData2 <- modelData[vegType != "landscape"]
## need to fit the model after dropping unused levels, otherwise `marginaleffects::predictions` will error
## when trying to get prediction intervals
cols <- names(modelData2)
modelData2[,
  (cols) := lapply(.SD, function(x) {
    if (is.factor(x)) droplevels(x) else x
  })
]

pyroVSbiodiversityVegTypes.lm <- lm(
  vegHV ~ fireHV * scenario * vegType,
  data = modelData2
)
pyroVSbiodiversityVegTypes.lm2 <- lm(
  logVegHV ~ logFireHV * scenario * vegType,
  data = modelData2
)
check_model(pyroVSbiodiversityVegTypes.lm)
check_model(pyroVSbiodiversityVegTypes.lm2) ## much better

check_normality(pyroVSbiodiversityVegTypes.lm) ## not ok
check_heteroscedasticity(pyroVSbiodiversityVegTypes.lm) ## not ok
check_normality(pyroVSbiodiversityVegTypes.lm2) ## not ok
check_heteroscedasticity(pyroVSbiodiversityVegTypes.lm2) ## not ok

pyroVSbiodiversityVegTypes.lm3 <- lm(
  vegHV ~ poly(fireHV, 2) * scenario * vegType,
  data = modelData2
)
check_model(pyroVSbiodiversityVegTypes.lm3) ## slightly better
anova(pyroVSbiodiversityVegTypes.lm, pyroVSbiodiversityVegTypes.lm3) # sign. different
anova(pyroVSbiodiversityVegTypes.lm3) ## interactions are significant, but not all

check_normality(pyroVSbiodiversityVegTypes.lm3) ## not ok
check_heteroscedasticity(pyroVSbiodiversityVegTypes.lm3) ## not ok

pyroVSbiodiversityVegTypes.lm4 <- lm(
  logVegHV ~ poly(logFireHV, 2) * scenario * vegType,
  data = modelData2
)
check_model(pyroVSbiodiversityVegTypes.lm4) ## slightly better
anova(pyroVSbiodiversityVegTypes.lm2, pyroVSbiodiversityVegTypes.lm4) # sign. different
anova(pyroVSbiodiversityVegTypes.lm4) ## all interactions are significant

check_normality(pyroVSbiodiversityVegTypes.lm4) ## not ok
check_heteroscedasticity(pyroVSbiodiversityVegTypes.lm4) ## not ok

## logging is best

## interactions are significant. no point in removing them.
# pyroVSbiodiversityVegTypes.lm4 <- lm(logVegHV ~ logFireHV*scenario + logFireHV*vegType + scenario*vegType,
#                                      data = modelData2)
# pyroVSbiodiversityVegTypes.lm5 <- lm(logVegHV ~ poly(logFireHV, 2)*scenario + poly(logFireHV, 2)*vegType + scenario*vegType,
#                                      data = modelData2)

## mnot logged
AICcmodavg::AICc(pyroVSbiodiversityVegTypes.lm) ##
AICcmodavg::AICc(pyroVSbiodiversityVegTypes.lm3) ## worse

## logged
AICcmodavg::AICc(pyroVSbiodiversityVegTypes.lm2)
AICcmodavg::AICc(pyroVSbiodiversityVegTypes.lm4) ## best but not by that much
# AICcmodavg::AICc(pyroVSbiodiversityVegTypes.lm5)

pyroVSbiodiversityVegTypes.gls <- nlme::gls(
  logVegHV ~ logFireHV * scenario * vegType,
  weights = nlme::varIdent(form = ~ 1 | scenario * vegType),
  data = modelData2
)
pyroVSbiodiversityVegTypes.gls2 <- nlme::gls(
  logVegHV ~ poly(logFireHV, 2) * scenario * vegType,
  weights = nlme::varIdent(form = ~ 1 | scenario * vegType),
  data = modelData2
) ## much better!
pyroVSbiodiversityVegTypes.gam <- mgcv::gam(
  logVegHV ~ s(logFireHV, k = 3, by = scenario) +
    s(logFireHV, k = 3, by = vegType),
  data = modelData2
)

check_model(pyroVSbiodiversityVegTypes.gls) ## pretty much the same as the LM equivalent

## not poly
AICcmodavg::AICc(pyroVSbiodiversityVegTypes.lm2)
AICcmodavg::AICc(pyroVSbiodiversityVegTypes.gls) ## better
## poly
AICcmodavg::AICc(pyroVSbiodiversityVegTypes.lm4)
AICcmodavg::AICc(pyroVSbiodiversityVegTypes.gls2) ## best
AICcmodavg::AICc(pyroVSbiodiversityVegTypes.gam)


## again, I don't see huge support for polynomial shapes in many cases. so option for linear effects
plotData <- modelData2
plotData[, predLM2 := predict(pyroVSbiodiversityVegTypes.lm2)]
plotData[, predLM4 := predict(pyroVSbiodiversityVegTypes.lm4)]
plotData[, predGLS := predict(pyroVSbiodiversityVegTypes.gls)]
plotData[, predGLS2 := predict(pyroVSbiodiversityVegTypes.gls2)]
ggplot(
  # plotData[scenario != "HV_noPM" & rep != 4]
  # plotData[scenario == "HV_noPM"]
  plotData,
  aes(x = logFireHV)
) +
  # geom_line(aes(y = predLM2, linetype = scenario), col = "red") +
  # geom_line(aes(y = predLM4, linetype = scenario), col = "green") +
  geom_line(aes(y = predGLS, linetype = scenario), col = "orange") + ## overlaps predLM2 completely.
  # geom_line(aes(y = predGLS2, linetype = scenario), col = "blue") +   ## overlaps predLM4 completely.
  # geom_smooth(aes(y = logVegHV, linetype = scenario), formula = y ~ poly(x, 2),
  #             colour = "red",
  #             method = "lm", se = FALSE) +
  # geom_smooth(aes(y = logVegHV, linetype = scenario), formula = y ~ x,
  #             colour = "blue",
  #             method = "lm", se = FALSE) +
  # geom_smooth(aes(y = logVegHV), formula = y ~ poly(x, 2),
  #             colour = "red",
  #             method = "lm", se = FALSE) +
  # geom_smooth(aes(y = logVegHV), formula = y ~ x,
  #             colour = "blue",
  #             method = "lm", se = FALSE) +
  geom_point(aes(y = logVegHV, colour = as.factor(rep), shape = scenario)) +
  facet_wrap(~vegType, scales = 'free')


png(file.path(figOutputPath, "pyroVSbiodiversityVegTypesglsRESIDUALS.png"))
check_model(pyroVSbiodiversityVegTypes.gls)
dev.off()

sink(file.path(statsOutputPath, "pyroVSbiodiversityVegTypesglsSUMMARY.txt"))
anova(pyroVSbiodiversityVegTypes.gls)
cat("\n*********************\n")
summary(pyroVSbiodiversityVegTypes.gls)
cat("\n*********************\n")
emtrends(
  pyroVSbiodiversityVegTypes.gls,
  pairwise ~ scenario | vegType,
  var = "logFireHV",
  max.degree = 2,
  mode = "df.error"
)
sink()

## PLOTS: EFFECT OF SCENARIO ON PYRODIVERSITY AND BIODIVERSITY ----------------------

### how does scenario affect hypervolume sizes at the end of the simulation?
### general plot

## reorder vegType for plotting
allHVData[, vegType := factor(vegType, levels = names(vegTypeCNLabels))]

## melt volume data
plotData <- allHVData[
  year == max(yearSubset),
  .(
    HVtype,
    HV_noPM,
    HV_PM,
    Intersection,
    MinDist,
    CentroidDist,
    rep,
    repHV,
    vegType
  )
]
plotData <- melt.data.table(
  plotData,
  measure.vars = c("HV_noPM", "HV_PM"),
  variable.name = "scenario",
  value.name = "Volume"
)
plotData[, logVolume := log(Volume)]

plotData[, scenario := sub("HV_", "", scenario)]

pyroHVvolumeVegTypesPlot <- HVBoxplots(
  plotData = plotData[vegType != "landscape" & HVtype == "fireHV"],
  yLab = "log-hypervolume size",
  xLab = "",
  fillLab = "",
  titleLab = "Pyrodiversity",
  xLabels = vegTypeCNLabels,
  fillLabels = scenLabels,
  fillVals = scenColours
)

bioHVvolumeVegTypesPlot <- HVBoxplots(
  plotData = plotData[vegType != "landscape" & HVtype == "vegHV"],
  yLab = "log-hypervolume size",
  xLab = "",
  fillLab = "",
  titleLab = "Forest diversity",
  xLabels = vegTypeCNLabels,
  fillLabels = scenLabels,
  fillVals = scenColours
)


pyroHVvolumeLandscapePlot <- HVBoxplots(
  plotData[vegType == "landscape" & HVtype == "fireHV"],
  yLab = "log-hypervolume size",
  xLab = "",
  titleLab = "Pyrodiversity",
  xLabels = vegTypeCNLabels,
  fillLabels = scenLabels,
  fillVals = scenColours
) +
  guides(fill = "none")


bioHVvolumeLandscapePlot <- HVBoxplots(
  plotData[vegType == "landscape" & HVtype == "vegHV"],
  yLab = "log-hypervolume size",
  xLab = "",
  titleLab = "Forest diversity",
  xLabels = vegTypeCNLabels,
  fillLabels = scenLabels,
  fillVals = scenColours
) +
  guides(fill = "none")


alignedVegPlots <- align_plots(
  pyroHVvolumeVegTypesPlot +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      legend.position = "none",
      axis.text.x = element_blank()
    ),
  bioHVvolumeVegTypesPlot +
    theme(plot.margin = margin(0, 0, 0, 0), legend.position = "none"),
  align = "v",
  axis = "l"
)
alignedLandPlots <- align_plots(
  pyroHVvolumeLandscapePlot +
    labs(y = "", title = "") +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      legend.position = "none",
      axis.text.x = element_blank()
    ),
  bioHVvolumeLandscapePlot +
    labs(y = "", title = "") +
    theme(plot.margin = margin(0, 0, 0, 0), legend.position = "none"),
  align = "v",
  axis = "l"
)

plotSave <- plot_grid(
  plot_grid(
    plot_grid(
      alignedVegPlots[[1]],
      alignedLandPlots[[1]],
      ncol = 2,
      nrow = 1,
      align = "h",
      axis = "b",
      rel_widths = c(1, 0.2),
      labels = c("a", "b")
    ),
    plot_grid(
      alignedVegPlots[[2]],
      alignedLandPlots[[2]],
      ncol = 2,
      nrow = 1,
      align = "v",
      axis = "b",
      rel_widths = c(1, 0.2),
      labels = c("c", "d")
    ),
    ncol = 1,
    nrow = 2,
    align = "v",
    axis = "l",
    rel_heights = c(0.76, 1)
  ),
  get_legend(bioHVvolumeVegTypesPlot + theme(legend.box = "horizontal")),
  ncol = 1,
  nrow = 2,
  rel_heights = c(1, 0.15)
)
ggsave(
  plot = plotSave,
  filename = file.path(figOutputPath, "HVVolumes.png"),
  width = 12,
  height = 8,
  bg = "white"
)

## were biodiversity volumes different at the start?
if (useFirstLastYear) {
  plotData <- allHVData[
    year == min(yearSubset),
    .(
      HVtype,
      HV_noPM,
      HV_PM,
      Intersection,
      MinDist,
      CentroidDist,
      rep,
      repHV,
      vegType
    )
  ]
  plotData <- melt.data.table(
    plotData,
    measure.vars = c("HV_noPM", "HV_PM"),
    variable.name = "scenario",
    value.name = "Volume"
  )
  HVvolumeVegTypesStartPlot <- ggplot(
    plotData[HVtype == "vegHV"],
    aes(x = vegType, y = Volume, alpha = scenario, fill = vegType)
  ) +
    geom_boxplot() +
    scale_x_discrete(labels = vegTypeCNLabels) +
    scale_fill_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
    scale_alpha_manual(
      values = c("HV_noPM" = 0.4, "HV_PM" = 1.0),
      labels = scenLabels
    ) +
    theme_pubr(
      base_size = 12,
      x.text.angle = 30,
      margin = FALSE,
      legend = "right"
    ) +
    theme(
      legend.box = "vertical",
      strip.background = element_blank(),
      panel.grid.major.y = element_line(
        colour = "grey",
        size = 11 / 22,
        linetype = "dotted"
      )
    ) +
    labs(x = "", y = "log-hypervolume size", fill = "", alpha = "") +
    guides(alpha = guide_legend(override.aes = list(fill = "grey50")))

  ggsave(
    plot = HVvolumeVegTypesStartPlot,
    filename = file.path(figOutputPath, "HVVolumesVegStart.png"),
    width = 7,
    height = 5,
    bg = "white"
  )

  ### How much did communities change in time? ---------------------------------------
  plotData <- allHVData[
    HVtype == "vegHV" &
      compare == paste(min(yearSubset), max(yearSubset), sep = "_"),
    .(overlap, MinDist, CentroidDist, rep, repHV, vegType, scenario)
  ]

  HVOverlapVegTypesPlot <- ggplot(
    plotData[vegType != "landscape"],
    aes(x = vegType, y = overlap, alpha = scenario, fill = vegType)
  ) +
    geom_boxplot() +
    scale_x_discrete(labels = vegTypeCNLabels) +
    scale_fill_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
    scale_alpha_manual(
      values = c("noPM" = 0.4, "PM" = 1.0),
      labels = scenLabels
    ) +
    theme_pubr(base_size = 10, x.text.angle = 30, margin = FALSE) +
    theme(
      legend.box = "vertical",
      strip.background = element_blank(),
      panel.grid.major.y = element_line(
        colour = "grey",
        size = 11 / 22,
        linetype = "dotted"
      )
    ) +
    labs(x = "", y = "Overlap", fill = "", alpha = "") +
    guides(alpha = guide_legend(override.aes = list(fill = "grey50")))

  HVOverlapLandscapePlot <- ggplot(
    plotData[vegType == "landscape"],
    aes(x = vegType, y = overlap, alpha = scenario, fill = vegType)
  ) +
    geom_boxplot() +
    scale_x_discrete(labels = vegTypeCNLabels) +
    scale_fill_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
    scale_alpha_manual(
      values = c("noPM" = 0.4, "PM" = 1.0),
      labels = scenLabels
    ) +
    theme_pubr(base_size = 10, margin = FALSE) +
    theme(
      strip.background = element_blank(),
      panel.grid.major.y = element_line(
        colour = "grey",
        size = 11 / 22,
        linetype = "dotted"
      )
    ) +
    labs(x = "", y = "Overlap", fill = "") +
    guides(alpha = "none", fill = "none")

  plotSave <- ggarrange(
    HVOverlapVegTypesPlot + theme(legend.box = "horizontal"),
    HVOverlapLandscapePlot + theme(legend.box = "horizontal") + labs(y = ""),
    ncol = 2,
    nrow = 1,
    align = "h",
    widths = c(1, 0.5),
    common.legend = TRUE,
    legend = "bottom",
    labels = "auto"
  )
  ggsave(
    plot = plotSave,
    filename = file.path(figOutputPath, "HVOverlap.png"),
    width = 8,
    height = 6,
    bg = "white"
  )
}

## PLOTS: BIODIVERSITY-PYRODIVERSITY RELATIONSHIP ----------------------
## melt volume data and dcast by volume type to relate the two
plotData <- allHVData[
  year == max(yearSubset),
  .(HVtype, HV_noPM, HV_PM, rep, repHV, vegType)
]
plotData <- melt.data.table(
  plotData,
  measure.vars = c("HV_noPM", "HV_PM"),
  variable.name = "scenario",
  value.name = "Volume"
)

plotData <- dcast(plotData, ... ~ HVtype, value.var = "Volume")
plotData[, `:=`(logFireHV = log(fireHV), logVegHV = log(vegHV))]
plotData[,
  logFireHVcenter := scale(logFireHV, center = TRUE, scale = FALSE),
  by = .(scenario, vegType)
]


## get predicted values and prediction intervals (CIs on predictions)
## predict for many more X values to smooth out intervals, and add the original values
## to facilitate plotting
## THIS DIDN'T WORK!!! Predictions came out flat.

# plotData2 <- plotData[, list(minLFH = min(logFireHV), maxLFH = max(logFireHV)),
#                       by = .(rep, scenario, vegType)]
#
# plotData2 <- plotData2[, list(logFireHV = seq(minLFH, maxLFH, length.out = 50)),
#                       by = .(rep, scenario, vegType)]
# plotData2 <- rbind(plotData2, plotData[, .(rep, scenario, vegType, logFireHV)])

lmPreds <- predict(
  pyroVSbiodiversityLandscape.lm2,
  # newdata = plotData2[vegType == "landscape"],
  se.fit = TRUE,
  interval = "confidence",
  type = "response"
)

if (paste(version$major, version$minor, sep = ".") < "4.4.0") {
  stop("marginaleffecs::predictions needs R version >= 4.4.0")
} else {
  ## make sure modelData is the same -- recommended in ?predictions
  modelData <- allHVData[
    year == max(yearSubset),
    .(HVtype, HV_noPM, HV_PM, rep, repHV, vegType)
  ]
  modelData <- melt.data.table(
    modelData,
    measure.vars = c("HV_noPM", "HV_PM"),
    variable.name = "scenario",
    value.name = "Volume"
  )

  modelData <- dcast.data.table(modelData, ... ~ HVtype, value.var = "Volume")
  modelData[, `:=`(logFireHV = log(fireHV), logVegHV = log(vegHV))]
  modelData2 <- modelData[vegType != "landscape"]
  ## need to fit the model after dropping unused levels, otherwise `marginaleffects::predictions` will error
  ## when trying to get prediction intervals
  cols <- names(modelData2)
  modelData2[,
    (cols) := lapply(.SD, function(x) {
      if (is.factor(x)) droplevels(x) else x
    })
  ]

  glsPreds <- predictions(
    pyroVSbiodiversityVegTypes.gls,
    # newdata = plotData3,
    vcov = TRUE,
    by = FALSE,
    type = "lp"
  )
}

plotData2 <- rbind(
  cbind(plotData[vegType == "landscape"], lmPreds$fit),
  cbind(
    plotData[vegType != "landscape"],
    data.table(
      fit = glsPreds$estimate,
      lwr = glsPreds$conf.low,
      upr = glsPreds$conf.high
    )
  )
)

setnames(plotData2, "fit", "pred")

## this works bcs original predictor values are in plotData2
# plotData3 <- plotData[plotData2, on = .(rep, scenario, vegType, logFireHV)]

## see influential replicates -- the same replicate is not influential across all veg types.
# ggplot(plotData3, aes(x = logFireHV, y = logVegHV, colour = as.factor(rep))) +
# geom_point() +
# facet_grid(scenario ~ vegType, scales = "free")

### VegTypes as colours; scenarios in separate plots

## the actual model is a gls, but because stat_smooth fits a separate model for each level,
## this is okay for visualisation.
## SEs are removed because they are not the same as the gls's
pyroVSbioDivVegTypesPlotnoPM <- plotBioPyroFun(
  plotData2[vegType != "landscape" & scenario == "HV_noPM"],
  # yPoints = "pred",  ## only to that preds match the smoother line
  ymin = "lwr",
  ymax = "upr",
  colourLabels = vegTypeCNLabels[names(vegTypeCNLabels) != "landscape"], ## needed for guides
  colourVals = vegTypeCNColours[names(vegTypeCNColours) != "landscape"], ## needed for guides
  titleLab = scenLabels["noPM"]
)
pyroVSbioDivLandscapePlotnoPM <- plotBioPyroFun(
  plotData2[vegType == "landscape" & scenario == "HV_noPM"],
  # yPoints = "pred",  ## only to that preds match the smoother line
  ymin = "lwr",
  ymax = "upr",
  colourLabels = vegTypeCNLabels[names(vegTypeCNLabels) == "landscape"], ## needed for guides
  colourVals = vegTypeCNColours[names(vegTypeCNColours) == "landscape"], ## needed for guides
  titleLab = scenLabels["noPM"]
)

pyroVSbioDivVegTypesPlotPM <- plotBioPyroFun(
  plotData2[vegType != "landscape" & scenario == "HV_PM"],
  # yPoints = "pred",  ## only to that preds match the smoother line
  ymin = "lwr",
  ymax = "upr",
  colourLabels = vegTypeCNLabels[names(vegTypeCNLabels) != "landscape"], ## needed for guides
  colourVals = vegTypeCNColours[names(vegTypeCNColours) != "landscape"], ## needed for guides
  titleLab = scenLabels["PM"]
)

pyroVSbioDivLandscapePlotPM <- plotBioPyroFun(
  plotData2[vegType == "landscape" & scenario == "HV_PM"],
  # yPoints = "pred",  ## only to that preds match the smoother line
  ymin = "lwr",
  ymax = "upr",
  colourLabels = vegTypeCNLabels[names(vegTypeCNLabels) == "landscape"], ## needed for guides
  colourVals = vegTypeCNColours[names(vegTypeCNColours) == "landscape"], ## needed for guides
  titleLab = scenLabels["PM"]
)

plotSave <- ggarrange(
  pyroVSbioDivVegTypesPlotnoPM,
  pyroVSbioDivLandscapePlotnoPM + labs(y = "", title = ""),
  ncol = 2,
  nrow = 1,
  widths = c(1, 0.4),
  labels = "auto",
  label.y = 0.95,
  common.legend = TRUE,
  legend = "bottom"
)
ggsave(
  plot = plotSave,
  filename = file.path(figOutputPath, "pyroVsbiodiversitynoPM.png"),
  width = 14,
  height = 7,
  bg = "white"
)

plotSave <- ggarrange(
  pyroVSbioDivVegTypesPlotPM,
  pyroVSbioDivLandscapePlotPM + labs(y = "", title = ""),
  ncol = 2,
  nrow = 1,
  widths = c(1, 0.4),
  labels = "auto",
  label.y = 0.95,
  common.legend = TRUE,
  legend = "bottom"
)
ggsave(
  plot = plotSave,
  filename = file.path(figOutputPath, "pyroVsbiodiversityPM.png"),
  width = 14,
  height = 7,
  bg = "white"
)


plotSave <- ggarrange(
  pyroVSbioDivVegTypesPlotnoPM +
    facet_wrap(
      ~vegType,
      scales = "free",
      labeller = as_labeller(vegTypeCNLabels)
    ) +
    scale_x_continuous(), ## free the scale,
  pyroVSbioDivLandscapePlotnoPM + labs(y = "", title = ""),
  ncol = 2,
  nrow = 1,
  widths = c(1, 0.4),
  labels = "auto",
  label.y = 0.95,
  common.legend = TRUE,
  legend = "bottom"
)
ggsave(
  plot = plotSave,
  filename = file.path(figOutputPath, "pyroVsbiodiversitynoPM2.png"),
  width = 14,
  height = 7,
  bg = "white"
)

plotSave <- ggarrange(
  pyroVSbioDivVegTypesPlotPM +
    facet_wrap(
      ~vegType,
      scales = "free",
      labeller = as_labeller(vegTypeCNLabels)
    ) +
    scale_x_continuous(), ## free the scale,
  pyroVSbioDivLandscapePlotPM + labs(y = "", title = ""),
  ncol = 2,
  nrow = 1,
  widths = c(1, 0.4),
  labels = "auto",
  label.y = 0.95,
  common.legend = TRUE,
  legend = "bottom"
)
ggsave(
  plot = plotSave,
  filename = file.path(figOutputPath, "pyroVsbiodiversityPM2.png"),
  width = 14,
  height = 7,
  bg = "white"
)


### scenarios as colours; vegTypes in facets
scenColours2 <- scenColours
names(scenColours2) <- paste0("HV_", names(scenColours2))
pyroVSbioDivVegTypesPlot <- plotBioPyroFun(
  plotData2[vegType != "landscape"],
  # yPoints = "pred",  ## only to that preds match the smoother line
  ymin = "lwr",
  ymax = "upr",
  colour = "scenario",
  colourLabels = scenLabels, ## needed for guides
  colourVals = scenColours2, ## needed for guides
  titleLab = ""
) +
  facet_wrap(
    ~vegType,
    scales = "free",
    labeller = as_labeller(vegTypeCNLabels)
  ) +
  scale_x_continuous() ## free the scale

pyroVSbioDivLandscapePlot <- plotBioPyroFun(
  plotData2[vegType == "landscape"],
  # yPoints = "pred",  ## only to that preds match the smoother line
  ymin = "lwr",
  ymax = "upr",
  colour = "scenario",
  colourLabels = scenLabels, ## needed for guides
  colourVals = scenColours2, ## needed for guides
  titleLab = ""
)


plotSave <- ggarrange(
  pyroVSbioDivVegTypesPlot,
  pyroVSbioDivLandscapePlot + labs(y = "", title = ""),
  ncol = 2,
  nrow = 1,
  widths = c(1, 0.4),
  labels = "auto",
  label.y = 0.95,
  common.legend = TRUE,
  legend = "bottom"
)
ggsave(
  plot = plotSave,
  filename = file.path(figOutputPath, "pyroVsbiodiversityAll.png"),
  width = 14,
  height = 7,
  bg = "white"
)
