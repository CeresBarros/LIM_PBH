## --------------------------------------------------
##  HYPERVOLUMES OF PYRODIVERSITY AND BIODIVERSITY
##  Analyses of results
## --------------------------------------------------
library(SpaDES)
library(reproducible)
library(data.table)
library(ggplot2)
library(ggpubr)
library(nlme)
library(mgcv)
library(multcomp)
Require::Require("rvlenth/emmeans (>= 1.8.1-90002)")
library(cowplot)
library(ToolsCB)
library(performance)

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
if (Sys.info()["nodename"] == "W-VIC-A127584") {
  simPaths <- list(cachePath = file.path("F:", basename(getwd()), "R/SpaDES/cache", simDirName, "postSimAnalyses")
                   , modulePath = file.path("R/SpaDES/m")
                   , inputPath = file.path("R/SpaDES/inputs")
                   , outputPath = file.path("F:", basename(getwd()), "R/SpaDES/outputs", simDirName))
} else if (grepl("for-cast", Sys.info()["nodename"])) {
  simPaths <- list(cachePath = file.path("/mnt/scratch/cbarros", basename(getwd()), "R/SpaDES/cache", simDirName, "postSimAnalyses")
                   , modulePath = file.path("R/SpaDES/m")
                   , inputPath = file.path("R/SpaDES/inputs")
                   , outputPath = file.path("R/SpaDES/outputs", simDirName)
                   , rasterPath = file.path("/mnt/scratch/cbarros", basename(getwd()), "R/SpaDES/scratch/raster")
                   , scratchPath = file.path("/mnt/scratch/cbarros", basename(getwd()), "R/SpaDES/scratch"))
  data.table::setDTthreads(25)
  options(bitmapType="cairo")
} else {
  simPaths <- list(cachePath = file.path("R/SpaDES/cache", simDirName, "postSimAnalyses")
                   , modulePath = file.path("R/SpaDES/m")
                   , inputPath = file.path("R/SpaDES/inputs")
                   , outputPath = file.path("R/SpaDES/outputs", simDirName))
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

## model landscape separately
fireHVVolumeLandscape.lm <- lm(Volume ~ scenario, modelData[vegType == "landscape"])
## the data seems very dispersed - some reps are extreme outliers - logging helps -- Jul 14 2022 this is no longer true with 4000 years
hist(modelData[vegType == "landscape", Volume], breaks = 1000)
fireHVVolumeLandscape.lm2 <- lm(log(Volume) ~ scenario,
                                modelData[vegType == "landscape"])   ## a bit better even with 4000 years
summary(fireHVVolumeLandscape.lm)
AIC(fireHVVolumeLandscape.lm,
    fireHVVolumeLandscape.lm2) ## better

fireHVVolumeLandscape.gls <- gls(Volume ~ scenario, weights = varIdent(form = ~ 1 | scenario),
                                 modelData[vegType == "landscape"])
fireHVVolumeLandscape.gls2 <- gls(log(Volume) ~ scenario, weights = varIdent(form = ~ 1 | scenario),
                                  modelData[vegType == "landscape"])
## no big improvements, result is similar and AIC is slightly worse
plot(fireHVVolumeLandscape.gls2)
AIC(fireHVVolumeLandscape.lm,
    fireHVVolumeLandscape.lm2, #best
    fireHVVolumeLandscape.gls,
    fireHVVolumeLandscape.gls2)


png(file.path(figOutputPath, "fireHVVolumeLandscapelmRESIDUALS.png"))
## a bit heteroskedastic, but not that much
check_model(fireHVVolumeLandscape.lm2)
dev.off()

sink(file.path(statsOutputPath, "fireHVVolumeLandscapelmSUMMARY.txt"))
anova(fireHVVolumeLandscape.lm2)
cat("\n*********************\n")
summary(fireHVVolumeLandscape.lm2)
cat("\n*********************\n")
emmeans(fireHVVolumeLandscape.lm2, specs = "scenario")
sink()

## model vegTypes
fireHVVolumeVegTypes.lm <- lm(Volume ~ scenario * vegType, modelData[vegType != "landscape"])
## same issue as before
fireHVVolumeVegTypes.lm2 <- lm(log(Volume) ~ scenario * vegType, modelData[vegType != "landscape"])

sets <- par(mfrow = c(2,2))
plot(fireHVVolumeVegTypes.lm2) ## better
par(sets)

AIC(fireHVVolumeVegTypes.lm,
    fireHVVolumeVegTypes.lm2)

fireHVVolumeVegTypes.gls <- gls(log(Volume) ~ scenario * vegType,
                                weights = varIdent(form = ~ 1 | scenario * vegType),
                                modelData[vegType != "landscape"])
plot(fireHVVolumeVegTypes.gls)  ## looks better
summary(fireHVVolumeVegTypes.gls)   ## all significant (don't trust ANOVA... different results wrt summary)

## glht can't see the missing levels, so we fit the model on a separate data table
modelData2 <- as.data.frame(modelData[vegType != "landscape"])
modelData2$scenario <- factor(modelData2$scenario)
modelData2$vegType <- factor(modelData2$vegType)
modelData2$scenario <- relevel(modelData2$scenario, "HV_noPM")
modelData2$vegType <- relevel(modelData2$vegType, "No veg.")

fireHVVolumeVegTypes.gls2 <- gls(log(Volume) ~ scenario + vegType,
                                 weights = varIdent(form = ~ 1 | scenario*vegType),
                                 data = modelData2)
anova(fireHVVolumeVegTypes.gls, fireHVVolumeVegTypes.gls2)   ## interaction model is much better.

AIC(fireHVVolumeVegTypes.lm, fireHVVolumeVegTypes.lm2, fireHVVolumeVegTypes.gls)

png(file.path(figOutputPath, "fireHVVolumeVegTypesglsRESIDUALS.png"))
plot(fireHVVolumeVegTypes.gls)
dev.off()

sink(file.path(statsOutputPath, "fireHVVolumeVegTypesglsSUMMARY.txt"))
anova(fireHVVolumeVegTypes.gls)
cat("\n*********************\n")
summary(fireHVVolumeVegTypes.gls)
cat("\n*********************\n")
emmeans(fireHVVolumeVegTypes.gls, pairwise ~ scenario | vegType,
        adjust = "tukey", data = modelData2)    ## tukey contrasts for scenarios within vegTypes
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

## model landscape separately
modelData2 <- modelData[vegType == "landscape"]
vegHVVolumeLandscape.lm <- lm(log(Volume) ~ scenario, modelData2)

vegHVVolumeLandscape.gls <- gls(log(Volume) ~ scenario,
                                weights = varIdent(form = ~ 1 | scenario),
                                data = modelData2)

AIC(vegHVVolumeLandscape.lm,
    vegHVVolumeLandscape.gls)  ## slightly better

png(file.path(figOutputPath, "vegHVVolumeLandscapeglsRESIDUALS.png"))
plot(vegHVVolumeLandscape.gls) ## not too bad for an anova
dev.off()

sink(file.path(statsOutputPath, "vegHVVolumeLandscapeglsSUMMARY.txt"))
anova(vegHVVolumeLandscape.gls)
cat("\n*********************\n")
summary(vegHVVolumeLandscape.gls)
cat("\n*********************\n")
emmeans(vegHVVolumeLandscape.gls, specs = "scenario", data = modelData2)
sink()

## model vegTypes
modelData2 <- modelData[vegType != "landscape"]
vegHVVolumeVegTypes.lm <- lm(log(Volume) ~ scenario * vegType, modelData2)

## some heteroscedasticity
check_model(vegHVVolumeVegTypes.lm)

modelData2 <- as.data.frame(modelData[vegType != "landscape"])
modelData2$scenario <- factor(modelData2$scenario)
modelData2$vegType <- factor(modelData2$vegType)
modelData2$scenario <- relevel(modelData2$scenario, "HV_noPM")
modelData2$vegType <- relevel(modelData2$vegType, "No veg.")

vegHVVolumeVegTypes.gls <- gls(log(Volume) ~ scenario * vegType,
                               weights = varIdent(form = ~ 1 | scenario * vegType),
                               data = modelData2)

png(file.path(figOutputPath, "vegHVVolumeVegTypesglsRESIDUALS.png"))
plot(vegHVVolumeVegTypes.gls) ## better
dev.off()

sink(file.path(statsOutputPath, "vegHVVolumeVegTypesglsSUMMARY.txt"))
anova(vegHVVolumeVegTypes.gls)
cat("\n*********************\n")
summary(vegHVVolumeVegTypes.gls)
cat("\n*********************\n")
emmeans(vegHVVolumeVegTypes.gls, pairwise ~ scenario | vegType, adjust = "tukey",
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

## model landscape separately
pyroVSbiodiversityLandscape.lm <- lm(vegHV ~ fireHV*scenario, data = modelData[vegType == "landscape"])
## the data seems very dispersed for fire HVs - some reps are extreme outliers - logging helps
hist(modelData[vegType == "landscape", fireHV], breaks = 1000)
# hist(modelData[vegType == "landscape" & fireHV > 500, fireHV], breaks = 1000)   # no longer exists in 4000 years sims

## log variables
modelData[, `:=`(logFireHV = log(fireHV),
                 logVegHV = log(vegHV))]
pyroVSbiodiversityLandscape.lm2 <- lm(logVegHV ~ logFireHV*scenario, data = modelData[vegType == "landscape"])

## try quadratic (as per Steel et al 2021 and He et al 2019)
## note that by default, we're using orthogonal polinomials - AIC is the same, so are residuals
## coefficients differ
## due to a bug in nlme, gls doens't output the estimated polynomial (using poly) and therefore
## predictions with new data (necessary for emmeans) are wrong (see https://stackoverflow.com/questions/70746067/issues-predicting-with-nlmeglsquadratic-model-fitted-with-poly-2)
## instead of poly we now centre the data and use logFireHV + I(logFireHV^2) -- given the low degree polynomial this is okay.

## the bug was "seemingly" fixed but there's still an issue with estimates in gls.

# pyroVSbiodiversityLandscape.lm3 <- lm(logVegHV ~ poly(logFireHV, 2)*scenario,
#                                       data = modelData[vegType == "landscape"])
modelData[, logFireHVcenter := scale(logFireHV, center = TRUE, scale = FALSE),
          by = .(scenario, vegType)]
pyroVSbiodiversityLandscape.lm3 <- lm(logVegHV ~ (logFireHVcenter + I(logFireHVcenter^2))*scenario,
                                      data = modelData[vegType == "landscape"])

AIC(pyroVSbiodiversityLandscape.lm,
    pyroVSbiodiversityLandscape.lm2, ## best by 2.8, but the curve seems to match the data better
    pyroVSbiodiversityLandscape.lm3)

anova(pyroVSbiodiversityLandscape.lm2, pyroVSbiodiversityLandscape.lm3) ## not signficantly different

## for the noPM scenario the relationship seems like an exponential decay
pyroVSbiodiversityLandscape.lm2.2 <- lm(logVegHV ~ logFireHV,
                                        data = modelData[vegType == "landscape" & scenario == "HV_noPM"])
pyroVSbiodiversityLandscape.lm2.3 <- lm(logVegHV ~ fireHV,
                                        data = modelData[vegType == "landscape" & scenario == "HV_noPM"])
pyroVSbiodiversityLandscape.lm3.2 <- lm(logVegHV ~ (logFireHVcenter + I(logFireHVcenter^2)),
                                        data = modelData[vegType == "landscape" & scenario == "HV_noPM"])

AIC(pyroVSbiodiversityLandscape.lm2.2,
    pyroVSbiodiversityLandscape.lm2.3,
    pyroVSbiodiversityLandscape.lm3.2)   ## marginally better


png(file.path(figOutputPath, "pyroVSbiodiversityLandscapelmRESIDUALS.png"))
check_model(pyroVSbiodiversityLandscape.lm2)
dev.off()

sink(file.path(statsOutputPath, "pyroVSbiodiversityLandscapelmSUMMARY.txt"))
anova(pyroVSbiodiversityLandscape.lm2)
cat("\n*********************\n")
summary(pyroVSbiodiversityLandscape.lm2)
cat("\n*********************\n")
emtrends(pyroVSbiodiversityLandscape.lm2, specs = c("scenario"),
         var = "logFireHV")
sink()

## by vegType
pyroVSbiodiversityVegTypes.lm <- lm(vegHV ~ fireHV*scenario*vegType, data = modelData[vegType != "landscape"])
## the data seems very dispersed for fire HVs - some reps are extreme outliers - logging helps
hist(modelData[vegType != "landscape", fireHV], breaks = 1000)
hist(modelData[vegType != "landscape" & fireHV > 500, fireHV], breaks = 1000)
hist(modelData[vegType != "landscape" & fireHV < 500, fireHV], breaks = 1000)

pyroVSbiodiversityVegTypes.lm2 <- lm(logVegHV ~ logFireHV *scenario*vegType,
                                     data = modelData[vegType != "landscape"])
pyroVSbiodiversityVegTypes.lm3 <- lm(logVegHV ~ (logFireHVcenter + I(logFireHVcenter^2))*scenario*vegType,
                                     data = modelData[vegType != "landscape"])
AIC(pyroVSbiodiversityVegTypes.lm2,
    pyroVSbiodiversityVegTypes.lm3) ## better

check_model(pyroVSbiodiversityVegTypes.lm2)  ## a bit worse, but AIC is better

modelData2 <- as.data.frame(modelData[vegType != "landscape"])
modelData2$scenario <- factor(modelData2$scenario)
modelData2$vegType <- factor(modelData2$vegType)
modelData2$scenario <- relevel(modelData2$scenario, "HV_noPM")
modelData2$vegType <- relevel(modelData2$vegType, "No veg.")

pyroVSbiodiversityVegTypes.gls <- nlme::gls(logVegHV ~ (logFireHVcenter + I(logFireHVcenter^2))*scenario*vegType,
                                            weights = nlme::varIdent(form = ~ 1 | scenario * vegType),
                                            data = modelData2)   ## much better!
pyroVSbiodiversityVegTypes.gam <- mgcv::gam(logVegHV ~ s(logFireHV, k = 3, by = interaction(scenario, vegType)),
                                            data =  modelData2)

AIC(pyroVSbiodiversityVegTypes.lm,
    pyroVSbiodiversityVegTypes.lm2,
    pyroVSbiodiversityVegTypes.lm3,
    pyroVSbiodiversityVegTypes.gls,  ## best
    pyroVSbiodiversityVegTypes.gam)
plot(pyroVSbiodiversityVegTypes.gls)


png(file.path(figOutputPath, "pyroVSbiodiversityVegTypesglsRESIDUALS.png"))
plot(pyroVSbiodiversityVegTypes.gls)
dev.off()

sink(file.path(statsOutputPath, "pyroVSbiodiversityVegTypesglsSUMMARY.txt"))
anova(pyroVSbiodiversityVegTypes.gls)
cat("\n*********************\n")
summary(pyroVSbiodiversityVegTypes.gls)
cat("\n*********************\n")
emtrends(pyroVSbiodiversityVegTypes.gls, pairwise ~ scenario | vegType,
         var = "logFireHVcenter", max.degree = 2, mode = "df.error")
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

pyroHVvolumeVegTypesPlot <- ggplot(plotData[vegType != "landscape" & HVtype == "fireHV"],
                                   aes(x = vegType, y = log(Volume)
                                       , alpha = scenario
                                       , fill = vegType)) +
  geom_boxplot() +
  scale_x_discrete(labels = vegTypeCNLabels) +
  scale_fill_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
  scale_alpha_manual(values = c("HV_noPM" = 0.4, "HV_PM" = 1.0),
                     labels = scenLabels) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  theme_pubr(base_size = 12, x.text.angle = 30, margin = FALSE) +
  theme(legend.box = "vertical",
        strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size = 11/22, linetype = "dotted")) +
  labs(x = "", fill = "", y = "log-hypervolume size", alpha = "", title = "Pyrodiversity") +
  guides(alpha = guide_legend(override.aes = list(fill = "grey50")))
  # facet_wrap(~ vegType == "landscape", nrow = 2, scales = "free_y",
  #            labeller = labeller(HVtype = c("vegHV" = "Forest diversity",
  #                                           "fireHV" = "Pyrodiversity")))

bioHVvolumeVegTypesPlot <- ggplot(plotData[vegType != "landscape" & HVtype == "vegHV"],
                                  aes(x = vegType, y = log(Volume), alpha = scenario, fill = vegType)) +
  geom_boxplot() +
  scale_x_discrete(labels = vegTypeCNLabels) +
  scale_fill_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
  scale_alpha_manual(values = c("HV_noPM" = 0.4, "HV_PM" = 1.0),
                     labels = scenLabels) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  theme_pubr(base_size = 12, x.text.angle = 30, margin = FALSE) +
  theme(legend.box = "vertical",
        strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size = 11/22, linetype = "dotted")) +
  labs(x = "", y = "log-hypervolume size", fill = "", alpha = "", title = "Forest diversity") +
  guides(alpha = guide_legend(override.aes = list(fill = "grey50")))
# facet_wrap(~ HVtype, nrow = 2, scales = "free_y",
#            labeller = labeller(HVtype = c("vegHV" = "Forest diversity",
#                                           "fireHV" = "Pyrodiversity")))


pyroHVvolumeLandscapePlot <- ggplot(plotData[vegType == "landscape" & HVtype == "fireHV"],
                                    aes(x = vegType, y = log(Volume), alpha = scenario, fill = vegType)) +
  geom_boxplot() +
  scale_x_discrete(labels = vegTypeCNLabels) +
  scale_fill_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
  scale_alpha_manual(values = c("HV_noPM" = 0.4, "HV_PM" = 1.0),
                     labels = scenLabels) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size = 11/22, linetype = "dotted")) +
  labs(x = "", y = "log-hypervolume size", fill = "", title = "Pyrodiversity") +
  guides(alpha = "none", fill = "none")
# facet_wrap(~ HVtype, nrow = 2, scales = "free_y",
#            labeller = labeller(HVtype = c("vegHV" = "Forest diversity",
#                                           "fireHV" = "Pyrodiversity")))

bioHVvolumeLandscapePlot <- ggplot(plotData[vegType == "landscape" & HVtype == "vegHV"],
                                   aes(x = vegType, y = log(Volume), alpha = scenario, fill = vegType)) +
  geom_boxplot() +
  scale_x_discrete(labels = vegTypeCNLabels) +
  scale_fill_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
  scale_alpha_manual(values = c("HV_noPM" = 0.4, "HV_PM" = 1.0),
                     labels = scenLabels) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size = 11/22, linetype = "dotted")) +
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
              ncol = 2, nrow = 1, align = "h", axis = "b", rel_widths = c(1, 0.5),
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

## relationship between pyrodiversity and biodiversity --------------------
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

plotData[vegType == "landscape", pred := predict(pyroVSbiodiversityLandscape.lm2, type = "response")]   ## just to compare with smoother
plotData[vegType != "landscape", pred := predict(pyroVSbiodiversityVegTypes.gls, type = "response")]   ## just to compare with smoother


## the actual model is a gls, but because stat_smooth fits a separate model for each level,
## this is okay for visualisation.
## SEs are removed because they are not the same as the gls's
plotBioPyroFunSmooth <- function(plotData, title = "") {
  if (all(plotData$vegType == "landscape")) {
    form <- quote(y ~ x)
  } else {
    form <- quote(y ~ x + I(x^2))
  }

  ggplot(plotData,
         aes(x = logFireHVcenter, y = logVegHV,
             # linetype = scenario,
             # shape = scenario,
             colour = vegType)) +
    geom_point() +
    # geom_line(aes(y = pred)) +  ## just to check if it matches smoother
    stat_smooth(method = "lm", formula = form, se = FALSE) +
    scale_colour_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
    # scale_linetype_manual(labels = scenLabels,
    #                       values = scenLinetype) +
    # scale_shape_discrete(labels = scenLabels) +
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
       width = 14, height = 8, bg = "white")

plotSave <- ggarrange(pyroVSbioDivVegTypesPlotPM, pyroVSbioDivVegLandscapePlotPM + labs(y = "", title = ""),
                      ncol = 2, nrow = 1, widths = c(1, 0.6), labels = "auto", label.y = 0.95,
                      common.legend = TRUE, legend = "bottom")
ggsave(plot = plotSave, filename = file.path(figOutputPath, "pyroVsbiodiversitySmoothPM.png"),
       width = 14, height = 8, bg = "white")


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

