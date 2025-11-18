## --------------------------------------------------
##  HYPERVOLUMES OF PYRODIVERSITY AND BIODIVERSITY
##  Analyses of results
## --------------------------------------------------
if (paste0(version$major, ".", strsplit(version$minor, "[.]")[[1]][1]) != "4.1") {
  stop("Please install and run R v4.1")
}

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

Require(c("SpaDES",
          "reproducible",
          "data.table",
          "ggplot2",
          "ggpubr",
          "nlme",
          "mgcv",
          "multcomp",
          "rvlenth/emmeans (>= 2.0.0)",
          "cowplot",
          "ToolsCB",
          "performance",
          "dplyr"), install = FALSE)


source("R/R_tools/Useful_functions.R")
source("R/R_tools/hypervolumesHelpers.R")
source("R/R_tools/glhtMethods.R")
options("reproducible.useNewDigestAlgorithm" = 2
        , "spades.moduleCodeChecks" = FALSE
        , "reproducible.useCache" = TRUE
        , "reproducible.destinationPath" = normPath("R/SpaDES/inputs")
        , "reproducible.useGDAL" = FALSE
        , "LandR.assertions" = TRUE)

## general paths
# simDirName <- "mar2021Runs"
simDirName <- "mar2022Runs"
if (Sys.info()["nodename"] == "W-VIC-A127584" |
    Sys.info()["nodename"] == "L-VIC-A155348") {
  simPaths <- list(cachePath = file.path("F:", basename(getwd()), "R/SpaDES/cache", simDirName, "postSimAnalyses")
                   , modulePath = file.path("R/SpaDES/m")
                   , inputPath = file.path("F:", basename(getwd()), "R/SpaDES/inputs", simDirName)
                   , outputPath = file.path("F:", basename(getwd()), "R/SpaDES/outputs", simDirName))
} else if (grepl("for-cast", Sys.info()["nodename"])) {
  simPaths <- list(cachePath = file.path("R/SpaDES/cache", simDirName, "postSimAnalyses")
                   , modulePath = file.path("R/SpaDES/m")
                   , inputPath = file.path("R/SpaDES/inputs")
                   , outputPath = file.path("R/SpaDES/outputs", simDirName)
                   , rasterPath = file.path("/mnt/scratch/cbarros", basename(getwd()), "R/SpaDES/scratch/raster")
                   , scratchPath = file.path("/mnt/scratch/cbarros", basename(getwd()), "R/SpaDES/scratch"))
} else {
  simPaths <- list(cachePath = file.path("R/SpaDES/cache", simDirName, "postSimAnalyses")
                   , modulePath = file.path("R/SpaDES/m")
                   , inputPath = file.path("R/SpaDES/inputs")
                   , outputPath = file.path("R/SpaDES/outputs", simDirName))
}

if (grepl("for-cast|dc709f2508f1|45eafed436c8",
          Sys.info()["nodename"])) {
  data.table::setDTthreads(5)
  options(bitmapType="cairo")
}


## path to figure folder and cache folder
figOutputPath <- file.path(simPaths$outputPath, "figuresAnalysis/Montane")
statsOutputPath <- file.path(simPaths$outputPath, "statsAnalysis")
HVoutputPath <- file.path(simPaths$outputPath, "hypervolumes")

## are we using the merged douglas-fir/dry-conifer stands?
mergeDMCPSME <- FALSE  ## merge DMCPSME PSME dryPSME
mergePSME <- TRUE ## merge PSME dryPSME
if (mergeDMCPSME) {
  HVoutputPathMergedVegType <- file.path(simPaths$outputPath, "hypervolumes", "mergeDMCPSME")
  figOutputPath <- file.path(simPaths$outputPath, "figuresAnalysis/Montane/mergeDMCPSME")
  statsOutputPath <- file.path(simPaths$outputPath, "statsAnalysis/mergeDMCPSME")
}

if (mergePSME) {
  HVoutputPathMergedVegType <- file.path(simPaths$outputPath, "hypervolumes", "mergePSME")
  figOutputPath <- file.path(simPaths$outputPath, "figuresAnalysis/Montane/mergePSME")
  statsOutputPath <- file.path(simPaths$outputPath, "statsAnalysis/mergePSME")
}

dir.create(figOutputPath, recursive = TRUE)
dir.create(statsOutputPath, recursive = TRUE)

## PREP FIRE AND VEG DATA -------------------
yearSubset <- unique(as.integer(c(seq(3511, 4011, 5), 4011)))
runPrepResultsModule <- FALSE
source("R/SpaDES/simResultsDataPrep.R")
source("R/R_tools/prepFireData4HVs.R")

yearSamples <- sample5SimYears(allPixelCohortDataMnt[, .(year, rep)])   ## seed ensures same years are drawn
useFirstLastYear <- FALSE
source("R/R_tools/prepVegData4HVs.R")

## don't need these
rm(allPixelBurnData, allPixelCohortData, allPixelCohortDataMnt)
gc(reset = TRUE)

## LOAD HYPERVOLUMES RESULTS ----------------------
source("R/R_tools/prepHVData.R")

## LABELS AND COLOURS FOR PLOTTING ----------------
source("R/R_tools/plotLabels&Cols.R")

## ------------------------------------------------------------------------
## STATISTICS: EFFECT OF SCENARIO ON PYRODIVERSITY  -----------------------
## make separate datatable for stats so that we can change contrasts
modelData <- allHVData[year == max(yearSubset) & HVtype == "fireHV",
                       .(HVtype, HV_noPM, HV_PM, Intersection, MinDist, CentroidDist,
                         rep, repHV, vegType)]

modelData <- melt.data.table(modelData, measure.vars = c("HV_noPM", "HV_PM"),
                             variable.name = "scenario", value.name = "Volume")

modelData[, scenario := relevel(scenario, "HV_noPM")]
modelData[, vegType := relevel(vegType, "No veg.")]
modelData[, logVolume := log(Volume)]

if (isFALSE(nrow(modelData[, unique(repHV), by = .(vegType, rep)]) == length(unique(modelData$vegType)) * length(unique(modelData$rep)) * length(unique(modelData$repHV)))) {
  stop("Missing replicates")
}

## model landscape separately
fireHVVolumeLandscape.lm <- lm(Volume ~ scenario, modelData[vegType == "landscape"])
## the data seems very dispersed - some reps are extreme outliers - logging helps -- Jul 14 2022 this is no longer true with 4000 years
hist(modelData[vegType == "landscape", Volume], breaks = 1000)
fireHVVolumeLandscape.lm2 <- lm(logVolume ~ scenario,
                                modelData[vegType == "landscape"])   ## a bit better even with 4000 years
summary(fireHVVolumeLandscape.lm)

check_model(fireHVVolumeLandscape.lm) ## looks better
check_model(fireHVVolumeLandscape.lm2)

caret::RMSE(fitted(fireHVVolumeLandscape.lm), modelData$Volume)
caret::RMSE(exp(fitted(fireHVVolumeLandscape.lm2)), exp(modelData$logVolume)) ## only marginally better, slightly worse resids

fireHVVolumeLandscape.gls <- gls(Volume ~ scenario, weights = varIdent(form = ~ 1 | scenario),
                                 modelData[vegType == "landscape"])
fireHVVolumeLandscape.gls2 <- gls(logVolume ~ scenario, weights = varIdent(form = ~ 1 | scenario),
                                  modelData[vegType == "landscape"])
## no big improvements
check_model(fireHVVolumeLandscape.gls)
check_model(fireHVVolumeLandscape.gls2)  # slightly worse
caret::RMSE(fitted(fireHVVolumeLandscape.lm), modelData$Volume)
caret::RMSE(exp(fitted(fireHVVolumeLandscape.lm2)), exp(modelData$logVolume))
caret::RMSE(fitted(fireHVVolumeLandscape.gls), modelData$Volume)
caret::RMSE(exp(fitted(fireHVVolumeLandscape.gls2)), exp(modelData$logVolume)) ## only marginally better, slightly worse resids

## can't directly compare transformed and non-transformed AICs
AIC(fireHVVolumeLandscape.lm,  ## better
    fireHVVolumeLandscape.gls)
AIC(fireHVVolumeLandscape.lm2,  ## better
    fireHVVolumeLandscape.gls2)


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

## model vegTypes
fireHVVolumeVegTypes.lm <- lm(Volume ~ scenario * vegType, modelData[vegType != "landscape"])
## same issue as before
fireHVVolumeVegTypes.lm2 <- lm(logVolume ~ scenario * vegType, modelData[vegType != "landscape"])

check_model(fireHVVolumeVegTypes.lm) ## better (Nov 2025)
check_model(fireHVVolumeVegTypes.lm2)

caret::RMSE(fitted(fireHVVolumeVegTypes.lm), modelData[vegType != "landscape", Volume])  ## marginally better
caret::RMSE(exp(fitted(fireHVVolumeVegTypes.lm2)), exp(modelData[vegType != "landscape", logVolume]))

## glht can't see the missing levels, so we fit the model on a separate data table
modelData2 <- as.data.frame(modelData[vegType != "landscape"])
modelData2$scenario <- factor(modelData2$scenario)
modelData2$vegType <- factor(modelData2$vegType)
modelData2$scenario <- relevel(modelData2$scenario, "HV_noPM")
modelData2$vegType <- relevel(modelData2$vegType, "No veg.")

fireHVVolumeVegTypes.gls <- gls(Volume ~ scenario * vegType,
                                weights = varIdent(form = ~ 1 | scenario * vegType),
                                modelData2)
check_model(fireHVVolumeVegTypes.gls)  ## looks better
anova(fireHVVolumeVegTypes.gls)    ## all significant
summary(fireHVVolumeVegTypes.gls)   ## all significant

fireHVVolumeVegTypes.gls2 <- gls(Volume ~ scenario + vegType,
                                 weights = varIdent(form = ~ 1 | scenario*vegType),
                                 data = modelData2)
anova(fireHVVolumeVegTypes.gls, fireHVVolumeVegTypes.gls2)   ## interaction model is much better.

AIC(fireHVVolumeVegTypes.lm,
    fireHVVolumeVegTypes.gls,  ## best
    fireHVVolumeVegTypes.gls2)

caret::RMSE(fitted(fireHVVolumeVegTypes.lm), modelData[vegType != "landscape", Volume])  ## marginally better
caret::RMSE(fitted(fireHVVolumeVegTypes.gls), modelData2$Volume)
caret::RMSE(fitted(fireHVVolumeVegTypes.gls2), modelData2$Volume)   ## marginally better, go with anova results above

png(file.path(figOutputPath, "fireHVVolumeVegTypesglsRESIDUALS.png"))
check_model(fireHVVolumeVegTypes.gls)
dev.off()

sink(file.path(statsOutputPath, "fireHVVolumeVegTypesglsSUMMARY.txt"))
anova(fireHVVolumeVegTypes.gls)
cat("\n*********************\n")
summary(fireHVVolumeVegTypes.gls)
cat("\n*********************\n")
## tukey contrasts for scenarios within vegTypes
emmeans(fireHVVolumeVegTypes.gls, pairwise ~ scenario | vegType,
        adjust = "tukey")
sink()


## ------------------------------------------------------------------------
## STATISTICS: EFFECT OF SCENARIO ON BIODIVERSITY  ------------------------
## how did scenario affect overall diversity at the end of simulation? ----
## (i.e. HV size)
modelData <- allHVData[year == max(yearSubset) & HVtype == "vegHV",
                       .(HVtype, HV_noPM, HV_PM, Intersection, MinDist, CentroidDist,
                         rep, repHV, vegType)]
modelData <- melt.data.table(modelData, measure.vars = c("HV_noPM", "HV_PM"),
                             variable.name = "scenario", value.name = "Volume")

modelData[, scenario := relevel(scenario, "HV_noPM")]
modelData[, vegType := relevel(vegType, "No veg.")]
modelData[, logVolume := log(Volume)]

if (isFALSE(nrow(modelData[, unique(repHV), by = .(vegType, rep)]) == length(unique(modelData$vegType)) * length(unique(modelData$rep)) * length(unique(modelData$repHV)))) {
  stop("Missing replicates")
}

## model landscape separately
modelData2 <- modelData[vegType == "landscape"]
vegHVVolumeLandscape.lm <- lm(Volume ~ scenario, modelData2)
vegHVVolumeLandscape.lm2 <- lm(logVolume ~ scenario, modelData2)

caret::RMSE(fitted(vegHVVolumeLandscape.lm), modelData2$Volume)
caret::RMSE(exp(fitted(vegHVVolumeLandscape.lm2)), exp(modelData2$logVolume)) ## a little worse

check_model(vegHVVolumeLandscape.lm)
check_model(vegHVVolumeLandscape.lm2) ## better

vegHVVolumeLandscape.gls <- gls(Volume ~ scenario,
                                weights = varIdent(form = ~ 1 | scenario),
                                data = modelData2)

vegHVVolumeLandscape.gls2 <- gls(logVolume ~ scenario,
                                 weights = varIdent(form = ~ 1 | scenario),
                                 data = modelData2)

AIC(vegHVVolumeLandscape.lm,
    vegHVVolumeLandscape.gls)  ## better
check_model(vegHVVolumeLandscape.gls) ## less deviation?
caret::RMSE(fitted(vegHVVolumeLandscape.lm), modelData2$Volume)
caret::RMSE(fitted(vegHVVolumeLandscape.gls), modelData2$Volume)  ## identical

AIC(vegHVVolumeLandscape.lm2,
    vegHVVolumeLandscape.gls2)  ## worse
check_model(vegHVVolumeLandscape.gls2) ## less deviation?
caret::RMSE(fitted(vegHVVolumeLandscape.gls), modelData2$Volume)
caret::RMSE(exp(fitted(vegHVVolumeLandscape.gls2)), modelData2$Volume)  ## slightly worse


png(file.path(figOutputPath, "vegHVVolumeLandscapeglsRESIDUALS.png"))
check_model(vegHVVolumeLandscape.gls2)
dev.off()

sink(file.path(statsOutputPath, "vegHVVolumeLandscapeglsSUMMARY.txt"))
anova(vegHVVolumeLandscape.gls2)
cat("\n*********************\n")
summary(vegHVVolumeLandscape.gls2)
cat("\n*********************\n")
emmeans(vegHVVolumeLandscape.gls2, specs = "scenario", data = modelData2)
sink()

## model vegTypes
modelData2 <- modelData[vegType != "landscape"]
vegHVVolumeVegTypes.lm <- lm(Volume ~ scenario * vegType, modelData2)
vegHVVolumeVegTypes.lm2 <- lm(logVolume ~ scenario * vegType, modelData2)

## some heteroscedasticity
check_model(vegHVVolumeVegTypes.lm)
check_model(vegHVVolumeVegTypes.lm2) ## better
caret::RMSE(fitted(vegHVVolumeVegTypes.lm), modelData2$Volume)
caret::RMSE(exp(fitted(vegHVVolumeVegTypes.lm2)), modelData2$Volume)  ## worse

modelData2 <- as.data.frame(modelData[vegType != "landscape"])
modelData2$scenario <- factor(modelData2$scenario)
modelData2$vegType <- factor(modelData2$vegType)
modelData2$scenario <- relevel(modelData2$scenario, "HV_noPM")
modelData2$vegType <- relevel(modelData2$vegType, "No veg.")

vegHVVolumeVegTypes.gls <- gls(Volume ~ scenario * vegType,
                               weights = varIdent(form = ~ 1 | scenario * vegType),
                               data = modelData2)
vegHVVolumeVegTypes.gls2 <- gls(logVolume ~ scenario * vegType,
                                weights = varIdent(form = ~ 1 | scenario * vegType),
                                data = modelData2)

check_model(vegHVVolumeVegTypes.gls)
check_model(vegHVVolumeVegTypes.gls2)   ## better than lm2

png(file.path(figOutputPath, "vegHVVolumeVegTypesglsRESIDUALS.png"))
check_model(vegHVVolumeVegTypes.gls2) ## better
dev.off()

sink(file.path(statsOutputPath, "vegHVVolumeVegTypesglsSUMMARY.txt"))
anova(vegHVVolumeVegTypes.gls2)
cat("\n*********************\n")
summary(vegHVVolumeVegTypes.gls2)
cat("\n*********************\n")
emmeans(vegHVVolumeVegTypes.gls2, pairwise ~ scenario | vegType, adjust = "tukey",
        deata = modelData2)    ## tukey contrasts for scenarios within vegTypes
cat("\n*********************\n")
sink()

## ------------------------------------------------------------------------
## How much did communities change in time? ---------------------------------------
## as in Barros et al 2016, calculate the prop overlap
if (useFirstLastYear) {
  modelData <- allHVData[HVtype == "vegHV" & compare == paste(min(yearSubset), max(yearSubset), sep = "_"),
                         .(overlap, MinDist, CentroidDist,
                           rep, repHV, vegType, scenario)]
  modelData[, scenario := relevel(factor(scenario), "noPM")]
  modelData[, vegType := relevel(vegType, "No veg.")]

  ## across the landscape
  vegHVOverlapLandscape.lm <- lm(overlap ~ scenario, modelData[vegType == "landscape"])

  ## some heteroscedasticity, but not lots and gls is no big improvement and results are similar
  sets <- par(mfrow = c(2,2))
  plot(vegHVOverlapLandscape.lm)
  par(sets)

  vegHVOverlapLandscape.gls <- gls(overlap ~ scenario, weights = varIdent(form = ~1 | scenario),
                                   data = modelData[vegType == "landscape"])
  AIC(vegHVOverlapLandscape.lm, vegHVOverlapLandscape.gls)  ## gls is worse and doens't improve residuals much

  png(file.path(figOutputPath, "vegHVOverlapLandscapelmRESIDUALS.png"))
  sets <- par(mfrow = c(2,2))
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
  vegHVOverlapVegTypes.lm <- lm(overlap ~ scenario * vegType,
                                modelData[vegType != "landscape"])
  ## heteroscedasticity
  sets <- par(mfrow = c(2,2))
  plot(vegHVOverlapVegTypes.lm) ## not too bad for an anova
  par(sets)

  modelData2 <- as.data.frame(modelData[vegType != "landscape"])
  modelData2$scenario <- factor(modelData2$scenario)
  modelData2$vegType <- factor(modelData2$vegType)
  modelData2$scenario <- relevel(modelData2$scenario, "HV_noPM")
  modelData2$vegType <- relevel(modelData2$vegType, "No veg.")

  vegHVOverlapVegTypes.gls <- gls(overlap ~ scenario * vegType,
                                  weights = varIdent(form = ~ 1 | scenario * vegType),
                                  modelData2)   ## can't converge

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
  emmeans(vegHVOverlapVegTypes.gls, pairwise ~ scenario | vegType, adjust = "tukey")    ## tukey contrasts for scenarios within vegTypes
  cat("\n*********************\n")
  sink()
}


## ------------------------------------------------------------------------
## STATISTICS: EFFECT OF SCENARIO ON PYRO-BIODIVERSITY RELATIONSHIP  ------
modelData <- allHVData[year == max(yearSubset),
                       .(HVtype, HV_noPM, HV_PM, rep, repHV, vegType)]
modelData <- melt.data.table(modelData, measure.vars = c("HV_noPM", "HV_PM"),
                             variable.name = "scenario", value.name = "Volume")

modelData <- dcast.data.table(modelData, ... ~ HVtype, value.var = "Volume")
modelData[, `:=`(logFireHV = log(fireHV),
                 logVegHV = log(vegHV))]
modelData2 <- modelData[vegType == "landscape"]

## model landscape separately
pyroVSbiodiversityLandscape.lm <- lm(vegHV ~ fireHV*scenario, data = modelData2)
## the data seems very dispersed for fire HVs - some reps are extreme outliers - logging helps -- no longer tru Nov 2025
hist(modelData[vegType == "landscape", fireHV], breaks = 1000)
# hist(modelData[vegType == "landscape" & fireHV > 500, fireHV], breaks = 1000)   # no longer exists in 4000 years sims

pyroVSbiodiversityLandscape.lm2 <- lm(logVegHV ~ logFireHV*scenario, data = modelData2)
check_model(pyroVSbiodiversityLandscape.lm)
check_model(pyroVSbiodiversityLandscape.lm2)  ## a bit better

anova(pyroVSbiodiversityLandscape.lm2) ## interaction not significant

pyroVSbiodiversityLandscape.lm3 <- lm(logVegHV ~ logFireHV + scenario, data = modelData2)
anova(pyroVSbiodiversityLandscape.lm3) ## scenarion effect not significant
check_model(pyroVSbiodiversityLandscape.lm3)  ## probably as bad.
anova(pyroVSbiodiversityLandscape.lm2, pyroVSbiodiversityLandscape.lm3) ## not sign. different (AIC is marginally better, see below)

pyroVSbiodiversityLandscape.lm4 <- lm(logVegHV ~ logFireHV, data = modelData2)
anova(pyroVSbiodiversityLandscape.lm4) ## scenarion effect not significant
check_model(pyroVSbiodiversityLandscape.lm4)  ## probably as bad.
anova(pyroVSbiodiversityLandscape.lm3, pyroVSbiodiversityLandscape.lm4) ## not sign. different (AIC is marginally better, see below)


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
pyroVSbiodiversityLandscape.lm5 <- lm(logVegHV ~ poly(logFireHV, 2), data = modelData2)
check_model(pyroVSbiodiversityLandscape.lm5)  ## as bad as previous
anova(pyroVSbiodiversityLandscape.lm4, pyroVSbiodiversityLandscape.lm5)  ## not sign. different, AIC slightly worse
pyroVSbiodiversityLandscape.lm6 <- lm(logVegHV ~ poly(logFireHV, 2)*scenario,
                                      data = modelData2)

## nov 2025, no support for the following models
# pyroVSbiodiversityLandscape.lm7 <- lm(logVegHV ~ logFireHV+scenario, data = modelData2)
# pyroVSbiodiversityLandscape.lm8 <- lm(logVegHV ~ poly(logFireHV, 2)+scenario,
#                                       data = modelData2)

# AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm)
AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm2)
AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm3)
AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm4)  ## best, but very marginally
AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm5)
AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm6)
# AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm7)
# AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm8)

## and if we allow the variance to change?
pyroVSbiodiversityLandscape.gls <- nlme::gls(logVegHV ~ logFireHV * scenario,
                                             weights = nlme::varIdent(form = ~ 1 | scenario),
                                             data = modelData2)
anova(pyroVSbiodiversityLandscape.gls)  ## interaction not significant
pyroVSbiodiversityLandscape.gls2 <- nlme::gls(logVegHV ~ logFireHV + scenario,
                                              weights = nlme::varIdent(form = ~ 1 | scenario),
                                              data = modelData2)
anova(pyroVSbiodiversityLandscape.gls2) ## senario effect not significant
pyroVSbiodiversityLandscape.gls3 <- nlme::gls(logVegHV ~ logFireHV,
                                              weights = nlme::varIdent(form = ~ 1 | scenario),
                                              data = modelData2)
anova(pyroVSbiodiversityLandscape.gls3)

check_model(pyroVSbiodiversityLandscape.gls)  ## still pretty bad
check_model(pyroVSbiodiversityLandscape.gls3) ## still pretty bad

## try adding polinomial again
pyroVSbiodiversityLandscape.gls4 <- nlme::gls(logVegHV ~ poly(logFireHV, 2),
                                              weights = nlme::varIdent(form = ~ 1 | scenario),
                                              data = modelData2)
pyroVSbiodiversityLandscape.gls5 <- nlme::gls(logVegHV ~ poly(logFireHV, 2)*scenario,
                                              weights = nlme::varIdent(form = ~ 1 | scenario),
                                              data = modelData2)
AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm4)   ## second best
AICcmodavg::AICc(pyroVSbiodiversityLandscape.gls)   ## third best
AICcmodavg::AICc(pyroVSbiodiversityLandscape.gls2)
AICcmodavg::AICc(pyroVSbiodiversityLandscape.gls3)
AICcmodavg::AICc(pyroVSbiodiversityLandscape.gls4)  ## worse
AICcmodavg::AICc(pyroVSbiodiversityLandscape.gls5)  ## best by less than 2 points
AICcmodavg::AICc(pyroVSbiodiversityLandscape.lm6)

anova(pyroVSbiodiversityLandscape.gls, pyroVSbiodiversityLandscape.gls5)  ## sign different
anova(pyroVSbiodiversityLandscape.gls5)  ## but interactions are not significant

## Nov 2025: in PM, replicate 4 seems to be an outlier, without it the relationship would be flatter
## also, scenario does seem to have an impact on mean volumes and the slopes.
## I think it seems reasonable to include the scenario interaction and with more data it'd prbably be a significant effect
plotData <- modelData2
plotData[, pred := predict(pyroVSbiodiversityLandscape.gls5)]
ggplot(
  # plotData[scenario != "HV_noPM" & rep != 4]
  # plotData[scenario == "HV_noPM"]
  plotData
  , aes(x = logFireHV, y = logVegHV)) +
  geom_smooth(aes(linetype = scenario), formula = y ~ poly(x, 2),
              colour = "red",
              method = "lm", se = FALSE) +
  geom_smooth(aes(linetype = scenario), formula = y ~ x,
              colour = "blue",
              method = "lm", se = FALSE) +
  geom_smooth(formula = y ~ poly(x, 2),
              colour = "red",
              method = "lm", se = FALSE) +
  geom_smooth(formula = y ~ x,
              colour = "blue",
              method = "lm", se = FALSE) +
  geom_point(aes(colour = as.factor(rep), shape = scenario))


png(file.path(figOutputPath, "pyroVSbiodiversityLandscapeglsRESIDUALS.png"))
plot(pyroVSbiodiversityLandscape.gls5)
dev.off()

sink(file.path(statsOutputPath, "pyroVSbiodiversityLandscapeglsSUMMARY.txt"))
anova(pyroVSbiodiversityLandscape.gls5)
cat("\n*********************\n")
summary(pyroVSbiodiversityLandscape.gls5)
cat("\n*********************\n")
emtrends(pyroVSbiodiversityLandscape.gls5, specs = c("scenario"),
         var = "logFireHV", data = modelData2)
sink()

## by vegType
modelData2 <- modelData[vegType != "landscape"]
## log variables
modelData2[, `:=`(logFireHV = log(fireHV),
                  logVegHV = log(vegHV))]
pyroVSbiodiversityVegTypes.lm <- lm(vegHV ~ fireHV*scenario*vegType, data = modelData2)
pyroVSbiodiversityVegTypes.lm2 <- lm(logVegHV ~ logFireHV*scenario*vegType,
                                     data = modelData2)
check_model(pyroVSbiodiversityVegTypes.lm)
check_model(pyroVSbiodiversityVegTypes.lm2)  ## much better

pyroVSbiodiversityVegTypes.lm3 <- lm(logVegHV ~ poly(logFireHV, 2)*scenario*vegType,
                                     data = modelData2)
check_model(pyroVSbiodiversityVegTypes.lm3)  ## maybe a bit better?
anova(pyroVSbiodiversityVegTypes.lm2, pyroVSbiodiversityVegTypes.lm3) # sign. different
anova(pyroVSbiodiversityVegTypes.lm3)   ## interactions are significant

## interactions are significant. no point in removing them.
# pyroVSbiodiversityVegTypes.lm4 <- lm(logVegHV ~ logFireHV*scenario + logFireHV*vegType + scenario*vegType,
#                                      data = modelData2)
# pyroVSbiodiversityVegTypes.lm5 <- lm(logVegHV ~ poly(logFireHV, 2)*scenario + poly(logFireHV, 2)*vegType + scenario*vegType,
#                                      data = modelData2)

# AICcmodavg::AICc(pyroVSbiodiversityVegTypes.lm)  ## don't compare non-logged
AICcmodavg::AICc(pyroVSbiodiversityVegTypes.lm2)
AICcmodavg::AICc(pyroVSbiodiversityVegTypes.lm3) ## best
# AICcmodavg::AICc(pyroVSbiodiversityVegTypes.lm4)
# AICcmodavg::AICc(pyroVSbiodiversityVegTypes.lm5)

pyroVSbiodiversityVegTypes.gls <- nlme::gls(logVegHV ~ logFireHV*scenario*vegType,
                                            weights = nlme::varIdent(form = ~ 1 | scenario * vegType),
                                            data = modelData2)
pyroVSbiodiversityVegTypes.gls2 <- nlme::gls(logVegHV ~ poly(logFireHV, 2)*scenario*vegType,
                                             weights = nlme::varIdent(form = ~ 1 | scenario * vegType),
                                             data = modelData2)   ## much better!
pyroVSbiodiversityVegTypes.gam <- mgcv::gam(logVegHV ~ s(logFireHV, k = 3, by = scenario) + s(logFireHV, k = 3, by = vegType),
                                            data =  modelData2)

AICcmodavg::AICc(pyroVSbiodiversityVegTypes.lm3)
AICcmodavg::AICc(pyroVSbiodiversityVegTypes.gls)
AICcmodavg::AICc(pyroVSbiodiversityVegTypes.gls2)  ## best
AICcmodavg::AICc(pyroVSbiodiversityVegTypes.gam)

png(file.path(figOutputPath, "pyroVSbiodiversityVegTypesglsRESIDUALS.png"))
check_model(pyroVSbiodiversityVegTypes.gls2)
dev.off()

sink(file.path(statsOutputPath, "pyroVSbiodiversityVegTypesglsSUMMARY.txt"))
anova(pyroVSbiodiversityVegTypes.gls2)
cat("\n*********************\n")
summary(pyroVSbiodiversityVegTypes.gls2)
cat("\n*********************\n")
emtrends(pyroVSbiodiversityVegTypes.gls2, pairwise ~ scenario | vegType,
         var = "logFireHV", max.degree = 2, mode = "df.error")
sink()

## PLOTS: EFFECT OF SCENARIO ON PYRODIVERSITY AND BIODIVERSITY ----------------------
## how does scenario affect hypervolume sizes at the end of the simulation?
## general plot

## reorder vegType for plotting
allHVData[, vegType := factor(vegType, levels = names(vegTypeCNLabels))]

## melt volume data
plotData <- allHVData[year == max(yearSubset),
                      .(HVtype, HV_noPM, HV_PM, Intersection, MinDist, CentroidDist,
                        rep, repHV, vegType)]
plotData <- melt.data.table(plotData, measure.vars = c("HV_noPM", "HV_PM"),
                            variable.name = "scenario", value.name = "Volume")
plotData[, logVolume := log(Volume)]


pyroHVvolumeVegTypesPlot <- ggplot(plotData[vegType != "landscape" & HVtype == "fireHV"],
                                   aes(x = vegType, y = logVolume
                                       , alpha = scenario
                                       , fill = vegType)) +
  geom_boxplot() +
  scale_x_discrete(labels = vegTypeCNLabels) +
  scale_fill_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
  scale_alpha_manual(values = c("HV_noPM" = 0.4, "HV_PM" = 1.0),
                     labels = scenLabels) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(legend.box = "vertical",
        strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", linewidth = 11/22, linetype = "dotted")) +
  labs(x = "", fill = "", y = "log-hypervolume size", alpha = "", title = "Pyrodiversity") +
  guides(alpha = guide_legend(override.aes = list(fill = "grey50")))
# facet_wrap(~ vegType == "landscape", nrow = 2, scales = "free_y",
#            labeller = labeller(HVtype = c("vegHV" = "Forest diversity",
#                                           "fireHV" = "Pyrodiversity")))

bioHVvolumeVegTypesPlot <- ggplot(plotData[vegType != "landscape" & HVtype == "vegHV"],
                                  aes(x = vegType, y = logVolume, alpha = scenario, fill = vegType)) +
  geom_boxplot() +
  scale_x_discrete(labels = vegTypeCNLabels) +
  scale_fill_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
  scale_alpha_manual(values = c("HV_noPM" = 0.4, "HV_PM" = 1.0),
                     labels = scenLabels) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(legend.box = "vertical",
        strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", linewidth = 11/22, linetype = "dotted")) +
  labs(x = "", y = "log-hypervolume size", fill = "", alpha = "", title = "Forest diversity") +
  guides(alpha = guide_legend(override.aes = list(fill = "grey50")))
# facet_wrap(~ HVtype, nrow = 2, scales = "free_y",
#            labeller = labeller(HVtype = c("vegHV" = "Forest diversity",
#                                           "fireHV" = "Pyrodiversity")))


pyroHVvolumeLandscapePlot <- ggplot(plotData[vegType == "landscape" & HVtype == "fireHV"],
                                    aes(x = vegType, y = logVolume, alpha = scenario, fill = vegType)) +
  geom_boxplot() +
  scale_x_discrete(labels = vegTypeCNLabels) +
  scale_fill_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
  scale_alpha_manual(values = c("HV_noPM" = 0.4, "HV_PM" = 1.0),
                     labels = scenLabels) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", linewidth = 11/22, linetype = "dotted")) +
  labs(x = "", y = "log-hypervolume size", fill = "", title = "Pyrodiversity") +
  guides(alpha = "none", fill = "none")
# facet_wrap(~ HVtype, nrow = 2, scales = "free_y",
#            labeller = labeller(HVtype = c("vegHV" = "Forest diversity",
#                                           "fireHV" = "Pyrodiversity")))

bioHVvolumeLandscapePlot <- ggplot(plotData[vegType == "landscape" & HVtype == "vegHV"],
                                   aes(x = vegType, y = logVolume, alpha = scenario, fill = vegType)) +
  geom_boxplot() +
  scale_x_discrete(labels = vegTypeCNLabels) +
  scale_fill_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
  scale_alpha_manual(values = c("HV_noPM" = 0.4, "HV_PM" = 1.0),
                     labels = scenLabels) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", linewidth = 11/22, linetype = "dotted")) +
  labs(x = "", y = "log-hypervolume size", fill = "", title = "Forest diversity") +
  guides(alpha = "none", fill = "none")
# facet_wrap(~ HVtype, nrow = 2, scales = "free_y",
#            labeller = labeller(HVtype = c("vegHV" = "Forest diversity",
#                                           "fireHV" = "Pyrodiversity")))

alignedVegPlots <- align_plots(pyroHVvolumeVegTypesPlot + theme(plot.margin = margin(0,0,0,0), legend.position = "none", axis.text.x = element_blank()),
                               bioHVvolumeVegTypesPlot + theme(plot.margin = margin(0,0,0,0), legend.position = "none"),
                               align = "v", axis = "l")
alignedLandPlots <- align_plots(pyroHVvolumeLandscapePlot + labs(y = "", title = "") +
                                  theme(plot.margin = margin(0,0,0,0), legend.position = "none", axis.text.x = element_blank()),
                                bioHVvolumeLandscapePlot + labs(y = "", title = "") + theme(plot.margin = margin(0,0,0,0), legend.position = "none"),
                                align = "v", axis = "l")

plotSave <- plot_grid(
  plot_grid(
    plot_grid(alignedVegPlots[[1]],
              alignedLandPlots[[1]],
              ncol = 2, nrow = 1, align = "h", axis = "b", rel_widths = c(1, 0.5),
              labels = c("a", "b")),
    plot_grid(alignedVegPlots[[2]],
              alignedLandPlots[[2]],
              ncol = 2, nrow = 1, align = "v", axis = "b", rel_widths = c(1, 0.5),
              labels = c("c", "d")),
    ncol = 1, nrow = 2, align = "v", axis = "l", rel_heights = c(0.76, 1)),
  get_legend(bioHVvolumeVegTypesPlot + theme(legend.box = "horizontal")),
  ncol = 1, nrow = 2, rel_heights = c(1, 0.15))
ggsave(plot = plotSave, filename = file.path(figOutputPath, "HVVolumes.png"),
       width = 10, height = 8, bg = "white")

## were biodiversity volumes different at the start?
if (useFirstLastYear) {
  plotData <- allHVData[year == min(yearSubset),
                        .(HVtype, HV_noPM, HV_PM, Intersection, MinDist, CentroidDist,
                          rep, repHV, vegType)]
  plotData <- melt.data.table(plotData, measure.vars = c("HV_noPM", "HV_PM"),
                              variable.name = "scenario", value.name = "Volume")
  HVvolumeVegTypesStartPlot <- ggplot(plotData[HVtype == "vegHV"],
                                      aes(x = vegType, y = Volume, alpha = scenario, fill = vegType)) +
    geom_boxplot() +
    scale_x_discrete(labels = vegTypeCNLabels) +
    scale_fill_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
    scale_alpha_manual(values = c("HV_noPM" = 0.4, "HV_PM" = 1.0),
                       labels = scenLabels) +
    theme_pubr(base_size = 12, x.text.angle = 30,
               margin = FALSE, legend = "right") +
    theme(legend.box = "vertical",
          strip.background = element_blank(),
          panel.grid.major.y = element_line(colour = "grey", size = 11/22, linetype = "dotted")) +
    labs(x = "", y = "log-hypervolume size", fill = "", alpha = "") +
    guides(alpha = guide_legend(override.aes = list(fill = "grey50")))

  ggsave(plot = HVvolumeVegTypesStartPlot, filename = file.path(figOutputPath, "HVVolumesVegStart.png"),
         width = 7, height = 5, bg = "white")

  ## How much did communities change in time? ---------------------------------------
  plotData <- allHVData[HVtype == "vegHV" & compare == paste(min(yearSubset), max(yearSubset), sep = "_"),
                        .(overlap, MinDist, CentroidDist, rep, repHV, vegType, scenario)]

  HVOverlapVegTypesPlot <- ggplot(plotData[vegType != "landscape"],
                                  aes(x = vegType, y = overlap,
                                      alpha = scenario, fill = vegType)) +
    geom_boxplot() +
    scale_x_discrete(labels = vegTypeCNLabels) +
    scale_fill_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
    scale_alpha_manual(values = c("noPM" = 0.4, "PM" = 1.0),
                       labels = scenLabels) +
    theme_pubr(base_size = 10, x.text.angle = 30, margin = FALSE) +
    theme(legend.box = "vertical",
          strip.background = element_blank(),
          panel.grid.major.y = element_line(colour = "grey", size = 11/22, linetype = "dotted")) +
    labs(x = "", y = "Overlap", fill = "", alpha = "") +
    guides(alpha = guide_legend(override.aes = list(fill = "grey50")))

  HVOverlapLandscapePlot <- ggplot(plotData[vegType == "landscape"],
                                   aes(x = vegType, y = overlap, alpha = scenario, fill = vegType)) +
    geom_boxplot() +
    scale_x_discrete(labels = vegTypeCNLabels) +
    scale_fill_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
    scale_alpha_manual(values = c("noPM" = 0.4, "PM" = 1.0),
                       labels = scenLabels) +
    theme_pubr(base_size = 10, margin = FALSE) +
    theme(strip.background = element_blank(),
          panel.grid.major.y = element_line(colour = "grey", size = 11/22, linetype = "dotted")) +
    labs(x = "", y = "Overlap", fill = "") +
    guides(alpha = "none", fill = "none")

  plotSave <- ggarrange(HVOverlapVegTypesPlot + theme(legend.box = "horizontal"),
                        HVOverlapLandscapePlot + theme(legend.box = "horizontal") + labs(y = ""),
                        ncol = 2, nrow = 1, align = "h", widths = c(1, 0.5),
                        common.legend = TRUE, legend = "bottom", labels = "auto")
  ggsave(plot = plotSave, filename = file.path(figOutputPath, "HVOverlap.png"),
         width = 8, height = 6, bg = "white")
}

## PLOTS: relationship between pyrodiversity and biodiversity ----------------------
## melt volume data and dcast by volume type to relate the two
plotData <- allHVData[year == max(yearSubset),
                      .(HVtype, HV_noPM, HV_PM, rep, repHV, vegType)]
plotData <- melt.data.table(plotData, measure.vars = c("HV_noPM", "HV_PM"),
                            variable.name = "scenario", value.name = "Volume")

plotData <- dcast(plotData, ... ~ HVtype, value.var = "Volume")
plotData[, `:=`(logFireHV = log(fireHV),
                logVegHV = log(vegHV))]
plotData[, logFireHVcenter := scale(logFireHV, center = TRUE, scale = FALSE),
         by = .(scenario, vegType)]

plotData[vegType == "landscape", pred := predict(pyroVSbiodiversityLandscape.gls5, type = "response")]   ## just to compare with smoother
plotData[vegType != "landscape", pred := predict(pyroVSbiodiversityVegTypes.gls2, type = "response")]   ## just to compare with smoother

## see influential replicates -- the same replicate is not influential across all veg types.
ggplot(plotData, aes(x = logFireHV, y = logVegHV, colour = as.factor(rep))) +
  geom_point() +
  facet_grid(scenario ~ vegType, scales = "free")

## the actual model is a gls, but because stat_smooth fits a separate model for each level,
## this is okay for visualisation.
## SEs are removed because they are not the same as the gls's
plotBioPyroFunSmooth <- function(plotData, title = "") {
  # if (all(plotData$vegType == "landscape")) {
  #   form <- quote(y ~ x)
  # } else {
  # form <- quote(y ~ x + I(x^2))   ## needs logFireHVcenter below
  form <- quote(y ~ poly(x, 2))
  # }

  ggplot(plotData,
         aes(
           # x = logFireHVcenter
           x = logFireHV
           , y = logVegHV
           # linetype = scenario,
           # shape = scenario,
           , colour = vegType
         )) +
    geom_point() +
    # geom_line(aes(y = pred)) +  ## just to check if it matches smoother
    stat_smooth(method = "lm", formula = form, se = FALSE) +
    scale_colour_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
    # scale_linetype_manual(labels = scenLabels,
    #                       values = scenLinetype) +
    # scale_shape_discrete(labels = scenLabels) +
    # scale_x_continuous(limits = range(plotData[, logFireHVcenter])) +
    scale_x_continuous(limits = range(plotData[, logFireHV])) +
    theme_pubr(base_size = 12, margin = TRUE) +
    theme(legend.box = "vertical",
          strip.background = element_blank(),
          panel.grid.major.y = element_line(colour = "grey", linewidth = 11/22, linetype = "dotted")) +
    labs(x = "Pyrodiversity", y = "Forest diversity", title = title, colour = ""
         # , linetype = "", shape = ""
    ) +
    # facet_wrap( ~ vegType, labeller = labeller(vegType = vegTypeCNLabels),
    #             scales = "free") +
    guides(linetype = guide_legend(override.aes = list(colour = "black")))
}

pyroVSbioDivVegTypesPlotnoPM <- plotBioPyroFunSmooth(plotData[vegType != "landscape" & scenario == "HV_noPM"],
                                                     title = scenLabels["noPM"])
pyroVSbioDivVegLandscapePlotnoPM <- plotBioPyroFunSmooth(plotData[vegType == "landscape" & scenario == "HV_noPM"],
                                                         title = scenLabels["noPM"])
pyroVSbioDivVegTypesPlotPM <- plotBioPyroFunSmooth(plotData[vegType != "landscape" & scenario == "HV_PM"],
                                                   title = scenLabels["PM"])
pyroVSbioDivVegLandscapePlotPM <- plotBioPyroFunSmooth(plotData[vegType == "landscape" & scenario == "HV_PM"],
                                                       title = scenLabels["PM"])

plotSave <- ggarrange(pyroVSbioDivVegTypesPlotnoPM, pyroVSbioDivVegLandscapePlotnoPM + labs(y = "", title = ""),
                      ncol = 2, nrow = 1, widths = c(1, 0.6), labels = "auto", label.y = 0.95,
                      common.legend = TRUE, legend = "bottom")
ggsave(plot = plotSave, filename = file.path(figOutputPath, "pyroVsbiodiversitySmoothnoPM.png"),
       width = 14, height = 7, bg = "white")

plotSave <- ggarrange(pyroVSbioDivVegTypesPlotPM, pyroVSbioDivVegLandscapePlotPM + labs(y = "", title = ""),
                      ncol = 2, nrow = 1, widths = c(1, 0.6), labels = "auto", label.y = 0.95,
                      common.legend = TRUE, legend = "bottom")
ggsave(plot = plotSave, filename = file.path(figOutputPath, "pyroVsbiodiversitySmoothPM.png"),
       width = 14, height = 7, bg = "white")


## without the smoother, with SEs, but using "smoothed" predictions
newData <- copy(plotData)
newData[, pred := NULL]
newData <- newData[, list(logFireHVcenter = seq(min(logFireHVcenter), max(logFireHVcenter), length.out = 1000),
                          logFireHV = seq(min(logFireHV), max(logFireHV), length.out = 1000)),
                   by = .(scenario, vegType)]
preds <- as.data.table(predict(pyroVSbiodiversityLandscape.lm2, newdata = newData[vegType == "landscape"], se.fit = TRUE)[c("fit", "se.fit")])
newData[vegType == "landscape", `:=`("pred" = preds$fit, "pred.se" = preds$se.fit)]

preds <- as.data.table(predict(pyroVSbiodiversityVegTypes.gls, newdata = newData[vegType != "landscape"], se.fit = TRUE)[c("fit", "se.fit")])
newData[vegType != "landscape", `:=`("pred" = preds$fit, "pred.se" = preds$se.fit)]

plotBioPyroFunPreds <- function(plotData, newData, title = "") {
  ggplot() +
    geom_ribbon(data = newData,
                aes(x = logFireHVcenter,
                    ymin = pred - pred.se, ymax = pred + pred.se),
                fill = "grey75", colour = "grey75") +
    geom_point(data = plotData,
               aes(x = logFireHVcenter, y = logVegHV,
                   colour = vegType)) +
    geom_line(data = newData,
              aes(x = logFireHVcenter, y = pred,
                  colour = vegType), size = 1) +  ## just to check if it matches smoother
    scale_colour_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
    scale_x_continuous(limits = range(plotData[, logFireHVcenter])) +
    theme_pubr(base_size = 12, margin = TRUE) +
    theme(legend.box = "vertical",
          strip.background = element_blank(),
          panel.grid.major.y = element_line(colour = "grey", size = 11/22, linetype = "dotted")) +
    labs(x = "Pyrodiversity", y = "Forest diversity", title = title, colour = "", linetype = "", shape = "") +
    facet_wrap( ~ vegType, labeller = labeller(vegType = vegTypeCNLabels),
                scales = "free") +
    guides(linetype = guide_legend(override.aes = list(colour = "black")))
}


pyroVSbioDivVegTypesPlotnoPM <- plotBioPyroFunPreds(plotData[scenario == "HV_noPM" & vegType != "landscape"],
                                                    newData[scenario == "HV_noPM" & vegType != "landscape"],
                                                    title = scenLabels["noPM"])
pyroVSbioDivVegLandscapePlotnoPM <- plotBioPyroFunPreds(plotData[scenario == "HV_noPM" & vegType == "landscape"],
                                                        newData[scenario == "HV_noPM" & vegType == "landscape"],
                                                        title = scenLabels["noPM"])
pyroVSbioDivVegTypesPlotPM <- plotBioPyroFunPreds(plotData[scenario == "HV_PM" & vegType != "landscape"],
                                                  newData[scenario == "HV_PM" & vegType != "landscape"],
                                                  title = scenLabels["PM"])
pyroVSbioDivVegLandscapePlotPM <- plotBioPyroFunPreds(plotData[scenario == "HV_PM" & vegType == "landscape"],
                                                      newData[scenario == "HV_PM" & vegType == "landscape"],
                                                      title = scenLabels["PM"])

plotSave <- ggarrange(pyroVSbioDivVegTypesPlotnoPM, pyroVSbioDivVegLandscapePlotnoPM + labs(y = "", title = ""),
                      ncol = 2, nrow = 1, widths = c(1, 0.6), labels = "auto", label.y = 0.95,
                      common.legend = TRUE, legend = "bottom")
ggsave(plot = plotSave, filename = file.path(figOutputPath, "pyroVsbiodiversityPrednoPM.png"),
       width = 14, height = 8, bg = "white")

plotSave <- ggarrange(pyroVSbioDivVegTypesPlotPM, pyroVSbioDivVegLandscapePlotPM + labs(y = "", title = ""),
                      ncol = 2, nrow = 1, widths = c(1, 0.6), labels = "auto", label.y = 0.95,
                      common.legend = TRUE, legend = "bottom")
ggsave(plot = plotSave, filename = file.path(figOutputPath, "pyroVsbiodiversityPredPM.png"),
       width = 14, height = 8, bg = "white")

