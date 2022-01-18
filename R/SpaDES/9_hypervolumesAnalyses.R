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
library(multcomp)
library(emmeans)

source("R/R_tools/Useful_functions.R")
source("R/R_tools/hypervolumesHelpers.R")
source("R/R_tools/glhtMethods.R")
options("reproducible.useNewDigestAlgorithm" = 2)
options("spades.moduleCodeChecks" = FALSE)
options("reproducible.useCache" = TRUE)
options("reproducible.destinationPath" = normPath("R/SpaDES/inputs"))
options("reproducible.useGDAL" = FALSE)

## general paths
simDirName <- "jun2021Runs"
simPaths <- list(cachePath = file.path("R/SpaDES/cache", simDirName, "postSimAnalyses")
                 , modulePath = file.path("R/SpaDES/m")
                 , inputPath = file.path("R/SpaDES/inputs")
                 , outputPath = file.path("R/SpaDES/outputs", simDirName))

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

## LOAD SIM LIST  ---------------------------------
preSimList <- loadSimList(file.path(simPaths$outputPath, "noPM", "LIM_simInit_noPM"))

## LOAD HYPERVOLUMES RESULTS ----------------------
source("R/R_tools/prepHVData.R")

## LABELS AND COLOURS FOR PLOTTING ----------------
source("R/R_tools/plotLabels&Cols.R")

## ------------------------------------------------------------------------
## STATISTICS: EFFECT OF SCENARIO ON PYRODIVERSITY  -----------------------
## make separate datatable for stats so that we can change contrasts
modelData <- allHVData[year == end(preSimList) & HVtype == "fireHV",
                       .(HVtype, HV_noPM, HV_PM, Intersection, MinDist, CentroidDist,
                         rep, repHV, vegType)]

modelData <- melt.data.table(modelData, measure.vars = c("HV_noPM", "HV_PM"),
                             variable.name = "scenario", value.name = "Volume")

modelData[, scenario := relevel(scenario, "HV_noPM")]
modelData[, vegType := relevel(vegType, "No veg.")]

## model landscape separately
fireHVVolumeLandscape.lm <- lm(Volume ~ scenario, modelData[vegType == "landscape"])
## the data seems very dispersed - some reps are extreme outliers - logging helps
hist(modelData[vegType == "landscape", Volume], breaks = 1000)
fireHVVolumeLandscape.lm2 <- lm(log(Volume) ~ scenario, modelData[vegType == "landscape"])
summary(fireHVVolumeLandscape.lm)

fireHVVolumeLandscape.gls <- gls(log(Volume) ~ scenario, weights = varIdent(form = ~ 1 | scenario),
                                 modelData[vegType == "landscape"])
## no big improvements, result is similar and AIC is slightly worse
plot(fireHVVolumeLandscape.gls)

tiff(file.path(figOutputPath, "fireHVVolumeLandscapelmRESIDUALS.tiff"))
## a bit heteroskedastic, but not that much
sets <- par(mfrow = c(2,2))
plot(fireHVVolumeLandscape.lm2)
par(sets)
dev.off()

sink(file.path(statsOutputPath, "fireHVVolumeLandscapelmSUMMARY.txt"))
summary(fireHVVolumeLandscape.lm2)   ## same p-value estimate.
sink()

## model vegTypes
fireHVVolumeVegTypes.lm <- lm(Volume ~ scenario * vegType, modelData[vegType != "landscape"])
## same issue as before
fireHVVolumeVegTypes.lm2 <- lm(log(Volume) ~ scenario * vegType, modelData[vegType != "landscape"])

sets <- par(mfrow = c(2,2))
plot(fireHVVolumeVegTypes.lm2) ## not much better
par(sets)


fireHVVolumeVegTypes.gls <- gls(log(Volume) ~ scenario * vegType, weights = varIdent(form = ~ 1 | scenario * vegType),
                                modelData[vegType != "landscape"])
plot(fireHVVolumeVegTypes.gls)  ## looks better
anova(fireHVVolumeVegTypes.gls)   ## no significant interaction between scenario and vegtype. let's drop it
# Denom. DF: 522
# numDF  F-value p-value
# (Intercept)          1 9599.827  <.0001
# scenario             1   96.784  <.0001
# vegType              8  153.635  <.0001
# scenario:vegType     8    1.500  0.1542


## glht can't see the missing levels, so we fit the model on a separate data table
modelData2 <- as.data.frame(modelData[vegType != "landscape"])
modelData2$scenario <- factor(modelData2$scenario)
modelData2$vegType <- factor(modelData2$vegType)
fireHVVolumeVegTypes.gls2 <- gls(log(Volume) ~ scenario + vegType,
                                 weights = varIdent(form = ~ 1 | scenario*vegType),
                                 data = modelData2)

tiff(file.path(figOutputPath, "fireHVVolumeVegTypesglsRESIDUALS.tiff"))
plot(fireHVVolumeVegTypes.gls2)
dev.off()

sink(file.path(statsOutputPath, "fireHVVolumeVegTypesglsSUMMARY.txt"))
anova(fireHVVolumeVegTypes.gls2)
cat("\n*********************\n")
summary(fireHVVolumeVegTypes.gls2)
cat("\n*********************\n")
emmeans(fireHVVolumeVegTypes.gls2, pairwise ~ vegType, adjust = "tukey")    ## tukey contrasts for vegTypes
sink()


## ------------------------------------------------------------------------
## STATISTICS: EFFECT OF SCENARIO ON BIODIVERSITY  ------------------------
## how did scenario affect overall diversity at the end of simulation? ----
## (i.e. HV size)
modelData <- allHVData[year == end(preSimList) & HVtype == "vegHV",
                       .(HVtype, HV_noPM, HV_PM, Intersection, MinDist, CentroidDist,
                         rep, repHV, vegType)]
modelData <- melt.data.table(modelData, measure.vars = c("HV_noPM", "HV_PM"),
                             variable.name = "scenario", value.name = "Volume")

modelData[, scenario := relevel(scenario, "HV_noPM")]
modelData[, vegType := relevel(vegType, "No veg.")]

## model landscape separately
vegHVVolumeLandscape.lm <- lm(Volume ~ scenario, modelData[vegType == "landscape"])

tiff(file.path(figOutputPath, "vegHVVolumeLandscapelmRESIDUALS.tiff"))
sets <- par(mfrow = c(2,2))
plot(vegHVVolumeLandscape.lm) ## not too bad for an anova
par(sets)
dev.off()

sink(file.path(statsOutputPath, "vegHVVolumeLandscapelmSUMMARY.txt"))
summary(vegHVVolumeLandscape.lm)
sink()

## model vegTypes
vegHVVolumeVegTypes.lm <- lm(Volume ~ scenario * vegType, modelData[vegType != "landscape"])

## some heteroscedasticity
sets <- par(mfrow = c(2,2))
plot(vegHVVolumeVegTypes.lm)
par(sets)

modelData2 <- as.data.frame(modelData[vegType != "landscape"])
modelData2$scenario <- factor(modelData2$scenario)
modelData2$vegType <- factor(modelData2$vegType)
vegHVVolumeVegTypes.gls <- gls(Volume ~ scenario * vegType,
                               weights = varIdent(form = ~ 1 | scenario * vegType),
                               modelData2)

tiff(file.path(figOutputPath, "vegHVVolumeVegTypesglsRESIDUALS.tiff"))
plot(vegHVVolumeVegTypes.gls) ## better
dev.off()

sink(file.path(statsOutputPath, "vegHVVolumeVegTypesglsSUMMARY.txt"))
anova(vegHVVolumeVegTypes.gls)
cat("\n*********************\n")
summary(vegHVVolumeVegTypes.gls)
cat("\n*********************\n")
emmeans(vegHVVolumeVegTypes.gls, pairwise ~ scenario | vegType, adjust = "tukey")    ## tukey contrasts for scenarios within vegTypes
cat("\n*********************\n")
sink()

## ------------------------------------------------------------------------
## How much did communities change in time? ---------------------------------------
## as in Barros et al 2016, calculate the prop overlap
modelData <- allHVData[HVtype == "vegHV" & compare == "2011_2111",
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

tiff(file.path(figOutputPath, "vegHVOverlapLandscapelmRESIDUALS.tiff"))
sets <- par(mfrow = c(2,2))
plot(vegHVOverlapLandscape.lm) ## still not great
par(sets)
dev.off()

sink(file.path(statsOutputPath, "vegHVOverlapLandscapelmSUMMARY.txt"))
summary(vegHVOverlapLandscape.lm)
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

## can't fix de variance on vegtype because some levels have 0 variance
# vegHVOverlapVegTypes.gls <- gls(overlap ~ scenario * vegType,
#                                 weights = varIdent(form = ~ 1 | scenario * vegType),
#                                 modelData2)   ## can't converge

vegHVOverlapVegTypes.gls <- gls(overlap ~ scenario * vegType,
                                weights = varIdent(form = ~ 1 | scenario),
                                modelData2)

tiff(file.path(figOutputPath, "vegHVOvelapVegTypesglsRESIDUALS.tiff"))
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


## ------------------------------------------------------------------------
## STATISTICS: EFFECT OF SCENARIO ON PYRO-BIODIVERSITY RELATIONSHIP  ------

modelData <- allHVData[year == end(preSimList),
                       .(HVtype, HV_noPM, HV_PM, rep, repHV, vegType)]
modelData <- melt.data.table(modelData, measure.vars = c("HV_noPM", "HV_PM"),
                             variable.name = "scenario", value.name = "Volume")

modelData <- dcast.data.table(modelData, ... ~ HVtype, value.var = "Volume")

## model landscape separately
pyroVSbiodiversityLandscape.lm <- lm(vegHV ~ fireHV*scenario, data = modelData[vegType == "landscape"])
## the data seems very dispersed for fire HVs - some reps are extreme outliers - logging helps
hist(modelData[vegType == "landscape", fireHV], breaks = 1000)
hist(modelData[vegType == "landscape" & fireHV > 500, fireHV], breaks = 1000)
pyroVSbiodiversityLandscape.lm2 <- lm(log(vegHV) ~ log(fireHV)*scenario, data = modelData[vegType == "landscape"])

tiff(file.path(figOutputPath, "pyroVSbiodiversityLandscapelmRESIDUALS.tiff"))
sets <- par(mfrow = c(2,2))
plot(pyroVSbiodiversityLandscape.lm2)
par(sets)
dev.off()

sink(file.path(statsOutputPath, "pyroVSbiodiversityLandscapelmSUMMARY.txt"))
anova(pyroVSbiodiversityLandscape.lm2)
cat("\n*********************\n")
summary(pyroVSbiodiversityLandscape.lm2)
sink()

## by vetType
pyroVSbiodiversityVegTypes.lm <- lm(vegHV ~ fireHV*scenario*vegType, data = modelData[vegType != "landscape"])
## the data seems very dispersed for fire HVs - some reps are extreme outliers - logging helps
hist(modelData[vegType != "landscape", fireHV], breaks = 1000)
hist(modelData[vegType != "landscape" & fireHV > 500, fireHV], breaks = 1000)
hist(modelData[vegType != "landscape" & fireHV < 500, fireHV], breaks = 1000)
pyroVSbiodiversityVegTypes.lm2 <- lm(log(vegHV) ~ log(fireHV)*scenario*vegType, data = modelData[vegType != "landscape"])

tiff(file.path(figOutputPath, "pyroVSbiodiversityVegTypeslmRESIDUALS.tiff"))
sets <- par(mfrow = c(2,2))
plot(pyroVSbiodiversityVegTypes.lm2)
par(sets)
dev.off()

sink(file.path(statsOutputPath, "pyroVSbiodiversityVegTypeslmSUMMARY.txt"))
anova(pyroVSbiodiversityVegTypes.lm2)
cat("\n*********************\n")
summary(pyroVSbiodiversityVegTypes.lm2)
cat("\n*********************\n")
emtrends(pyroVSbiodiversityVegTypes.lm2, pairwise ~ scenario | vegType,
         var = "fireHV", mode = "df.error", data = modelData[vegType != "landscape"])
emtrends(pyroVSbiodiversityVegTypes.lm2, pairwise ~ vegType | scenario,
         var = "fireHV", mode = "df.error", data = modelData[vegType != "landscape"])
sink()


## PLOTS: EFFECT OF SCENARIO ON PYRODIVERSITY AND BIODIVERSITY ----------------------
## how does scenario affect hypervolume sizes at the end of the simulation?
## general plot

## reorder vegType for plotting
allHVData[, vegType := factor(vegType, levels = names(vegTypeCNLabels))]

## melt volume data
plotData <- allHVData[year == end(preSimList),
                      .(HVtype, HV_noPM, HV_PM, Intersection, MinDist, CentroidDist,
                        rep, repHV, vegType)]
plotData <- melt.data.table(plotData, measure.vars = c("HV_noPM", "HV_PM"),
                            variable.name = "scenario", value.name = "Volume")

## log fire volumes
plotData[HVtype == "fireHV", Volume := log(Volume)]

HVvolumeVegTypesPlot <- ggplot(plotData[vegType != "landscape"],
                               aes(x = vegType, y = Volume, alpha = scenario, fill = vegType)) +
  geom_boxplot() +
  scale_x_discrete(labels = vegTypeCNLabels) +
  scale_fill_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
  scale_alpha_manual(values = c("HV_noPM" = 0.4, "HV_PM" = 1.0),
                     labels = c("HV_noPM" = "no PM", "HV_PM" = "PM")) +
  theme_pubr(base_size = 12, x.text.angle = 30, margin = FALSE) +
  theme(legend.box = "vertical",
        strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size = 11/22, linetype = "dotted")) +
  labs(x = "", fill = "", alpha = "") +
  guides(alpha = guide_legend(override.aes = list(fill = "grey50"))) +
  facet_wrap(~ HVtype, nrow = 2, scales = "free_y",
             labeller = labeller(HVtype = c("vegHV" = "forest diversity",
                                            "fireHV" = "pyrodiversity")))

HVvolumeLandscapePlot <- ggplot(plotData[vegType == "landscape"],
                                aes(x = vegType, y = Volume, fill = scenario)) +
  geom_boxplot() +
  scale_x_discrete(labels = vegTypeCNLabels) +
  scale_fill_brewer(palette = "Greys",
                    labels = c("HV_noPM" = "no PM", "HV_PM" = "PM")) +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size = 11/22, linetype = "dotted")) +
  labs(x = "", fill = "") +
  guides(alpha = "none", fill = "none") +
  facet_wrap(~ HVtype, nrow = 2, scales = "free_y",
             labeller = labeller(HVtype = c("vegHV" = "forest diversity",
                                            "fireHV" = "pyrodiversity")))

plotSave <- ggarrange(HVvolumeVegTypesPlot, HVvolumeLandscapePlot + labs(y = ""),
                      ncol = 2, nrow = 1, align = "h", widths = c(1, 0.5),
                      common.legend = TRUE, legend = "bottom")
ggsave(plot = plotSave, filename = file.path(figOutputPath, "HVVolumes.tiff"),
       width = 10, height = 8)

## were biodiversity volumes different at the start?
plotData <- allHVData[year == start(preSimList),
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
                     labels = c("HV_noPM" = "no PM", "HV_PM" = "PM")) +
  theme_pubr(base_size = 12, x.text.angle = 30,
             margin = FALSE, legend = "right") +
  theme(legend.box = "vertical",
        strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size = 11/22, linetype = "dotted")) +
  labs(x = "", fill = "", alpha = "") +
  guides(alpha = guide_legend(override.aes = list(fill = "grey50")))

ggsave(plot = HVvolumeVegTypesStartPlot, filename = file.path(figOutputPath, "HVVolumesVegStart.tiff"),
       width = 7, height = 5)


## How much did communities change in time? ---------------------------------------
plotData <- allHVData[HVtype == "vegHV" & compare == "2011_2111",
                      .(overlap, MinDist, CentroidDist, rep, repHV, vegType, scenario)]

HVOverlapVegTypesPlot <- ggplot(plotData[vegType != "landscape"],
                                aes(x = vegType, y = overlap,
                                    alpha = scenario, fill = vegType)) +
  geom_boxplot() +
  scale_x_discrete(labels = vegTypeCNLabels) +
  scale_fill_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
  scale_alpha_manual(values = c("noPM" = 0.4, "PM" = 1.0),
                     labels = c("noPM" = "no PM", "PM" = "PM")) +
  theme_pubr(base_size = 12, x.text.angle = 30, margin = FALSE) +
  theme(legend.box = "vertical",
        strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size = 11/22, linetype = "dotted")) +
  labs(x = "", y = "Overlap", fill = "", alpha = "") +
  guides(alpha = guide_legend(override.aes = list(fill = "grey50")))

HVOverlapLandscapePlot <- ggplot(plotData[vegType == "landscape"],
                                 aes(x = vegType, y = overlap, fill = scenario)) +
  geom_boxplot() +
  scale_x_discrete(labels = vegTypeCNLabels) +
  scale_fill_brewer(palette = "Greys",
                    labels = c("noPM" = "no PM", "PM" = "PM")) +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size = 11/22, linetype = "dotted")) +
  labs(x = "", y = "Overlap", fill = "") +
  guides(alpha = "none", fill = "none")

plotSave <- ggarrange(HVOverlapVegTypesPlot, HVOverlapLandscapePlot + labs(y = ""),
                      ncol = 2, nrow = 1, align = "h", widths = c(1, 0.5),
                      common.legend = TRUE, legend = "bottom")
ggsave(plot = plotSave, filename = file.path(figOutputPath, "HVOverlap.tiff"),
       width = 10, height = 8)


## relationship between pyrodiversity and biodiversity --------------------
## melt volume data and dcast by volume type to relate the two
plotData <- allHVData[year == end(preSimList),
                      .(HVtype, HV_noPM, HV_PM, rep, repHV, vegType)]
plotData <- melt.data.table(plotData, measure.vars = c("HV_noPM", "HV_PM"),
                            variable.name = "scenario", value.name = "Volume")

plotData <- dcast(plotData, ... ~ HVtype, value.var = "Volume")

plotData2 <- as.data.frame(plotData[vegType != "landscape"])
plotData2$scenario <- factor(plotData2$scenario)
plotData2$vegType <- factor(plotData2$vegType)

## refit models for plot
pyroVSbiodiversityLandscape.lm2 <- lm(log(vegHV) ~ log(fireHV)*scenario, data = plotData[vegType == "landscape"])
pyroVSbiodiversityVegTypes.lm2 <- lm(log(vegHV) ~ log(fireHV)*scenario*vegType, data = plotData[vegType != "landscape"])

plotData[vegType == "landscape", pred := predict(pyroVSbiodiversityLandscape.lm2)]
plotData[vegType != "landscape", pred := predict(pyroVSbiodiversityVegTypes.lm2)]

pyroVSbioDivVegTypesPlot <- ggplot(plotData[vegType != "landscape" & fireHV < 1000],
                                   aes(x = log(fireHV), y = log(vegHV),
                                       linetype = scenario, shape = scenario,
                                       colour = vegType)) +
  geom_point() +
  geom_line(aes(y = pred)) +
  scale_colour_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
  scale_linetype_discrete(labels = c("HV_noPM" = "no PM", "HV_PM" = "PM")) +
  scale_shape_discrete(labels = c("HV_noPM" = "no PM", "HV_PM" = "PM")) +
  theme_pubr(base_size = 12, margin = TRUE) +
  theme(legend.box = "vertical",
        strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size = 11/22, linetype = "dotted")) +
  labs(x = "pyrodiversity", y = "forest diversity", colour = "", linetype = "", shape = "") +
  facet_wrap( ~ vegType, labeller = labeller(vegType = vegTypeCNLabels),
              scales = "free")

pyroVSbioDivVegLandscapePlot <- ggplot(plotData[vegType == "landscape"],
                                       aes(x = log(fireHV), y = log(vegHV), shape = scenario,
                                           linetype = scenario, colour = vegType)) +
  geom_point() +
  geom_line(aes(y = pred)) +
  scale_colour_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
  scale_linetype_discrete(labels = c("HV_noPM" = "no PM", "HV_PM" = "PM")) +
  scale_shape_discrete(labels = c("HV_noPM" = "no PM", "HV_PM" = "PM")) +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(legend.box = "vertical",
        strip.background = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5),
        panel.grid.major.y = element_line(colour = "grey", size = 11/22, linetype = "dotted")) +
  labs(x = "pyrodiversity", y = "forest diversity", title = "landscape",
       colour = "", linetype = "", shape = "")

plotSave <- ggarrange(pyroVSbioDivVegTypesPlot, pyroVSbioDivVegLandscapePlot + labs(y = ""),
                      ncol = 2, nrow = 1, widths = c(1.1, 0.6),
                      common.legend = TRUE, legend = "bottom")

ggsave(plot = plotSave, filename = file.path(figOutputPath, "pyroVsbiodiversity.tiff"),
       width = 12, height = 7)









