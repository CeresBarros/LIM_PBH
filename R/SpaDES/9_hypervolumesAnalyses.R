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
mergeDMCPSME <- TRUE
if (mergeDMCPSME) {
  HVoutputPathMergedVegType <- file.path(simPaths$outputPath, "hypervolumes", "mergeDMCPSME")
  figOutputPath <- file.path(simPaths$outputPath, "figuresAnalysis/Montane/mergeDMCPSME")
  statsOutputPath <- file.path(simPaths$outputPath, "statsAnalysis/mergeDMCPSME")
}

dir.create(figOutputPath, recursive = TRUE)
dir.create(statsOutputPath, recursive = TRUE)

## LOAD SIM LIST  ---------------------------------
preSimList <- loadSimList(file.path(simPaths$outputPath, "noPM", "LIM_simInit_noPM"))

## LOAD HYPERVOLUMES RESULTS ----------------------
allFiles <- list.files(HVoutputPath, "Intersection.*.rds", full.names = TRUE)
if (mergeDMCPSME) {
  allFiles <- grep("_DMCPSME_|_PSME_|_dryPSME_", allFiles, value = TRUE, invert = TRUE) ## remove HV for vegTypes that were merged
  allFiles <- c(allFiles, list.files(HVoutputPathMergedVegType, "Intersection.*.rds", full.names = TRUE)) ## add HV for merged vegType
}

fireHVData <- loadHVResultsFromRDS("fireHVs", allFiles)
if ("HVid" %in% names(fireHVData)) {
  fireHVData[, scenario := HVid]   ## no HVid with intersecion results
}
fireHVData[, year := as.integer(end(preSimList))]   ## fire HV are from last year, but integrate the whole simulation period
## drop unique components
set(fireHVData, NULL, grep("Unique", names(fireHVData)), NULL)

## add comparison type
comp <- sub(".*_", "", grep("Volume", names(fireHVData), value = TRUE))
comp <- paste(comp, collapse = "_")
fireHVData[, compare := comp]
setnames(fireHVData, new = sub("Volume_(HV)[0-9]_(.*)", "\\1_\\2", names(fireHVData)))

## check reps/years/scenario per veg type
if (getOption("LandR.assertions")) {
  temp <- split(fireHVData[, .(scenario, year, rep, repHV, vegType)], by = "vegType", keep.by = FALSE)
  temp <- lapply(temp, FUN = function(x) paste0(x$scenario, x$year, x$rep))
  temp <- lapply(temp, FUN = unique)
  test <- lapply(1:length(temp), function(n) setdiff(temp[[n]], unlist(temp[-n])))
  test <- sapply(test, length)

  if (any(test)) {
    stop("fireHVData has different combinations of scenario, year, rep across vegTypes")
  }
}

## get between and within year (between scenarios) comparisons separately
withinYearFiles <- grep("yr", allFiles, value = TRUE)
vegHVDataWYrComparisons <- loadHVResultsFromRDS("vegHVs", withinYearFiles)
if ("HVid" %in% names(vegHVDataWYrComparisons)) {
  vegHVDataWYrComparisons[, tempHVid := as.numeric(HVid)]
  vegHVDataWYrComparisons[is.na(tempHVid), scenario := HVid]
  vegHVDataWYrComparisons[!is.na(tempHVid), year := HVid]
  vegHVDataWYrComparisons[, tempHVid := NULL]
}
## drop unique components
set(vegHVDataWYrComparisons, NULL, grep("Unique", names(vegHVDataWYrComparisons)), NULL)

## add comparison type
comp <- sub(".*_", "", grep("Volume", names(vegHVDataWYrComparisons), value = TRUE))
comp <- paste(comp, collapse = "_")
vegHVDataWYrComparisons[, compare := comp]
setnames(vegHVDataWYrComparisons, new = sub("Volume_(HV)[0-9]_(.*)", "\\1_\\2", names(vegHVDataWYrComparisons)))

## between year comparisons data
betweenYearFiles <- grep("yr", allFiles, value = TRUE, invert = TRUE)
vegHVDataBYrComparisons <- loadHVResultsFromRDS("vegHVs", betweenYearFiles)
if ("HVid" %in% names(vegHVDataBYrComparisons)) {
  vegHVDataBYrComparisons[, tempHVid := as.numeric(HVid)]
  vegHVDataBYrComparisons[is.na(tempHVid), scenario := HVid]
  vegHVDataBYrComparisons[!is.na(tempHVid), year := HVid]
  vegHVDataBYrComparisons[, tempHVid := NULL]
}
## drop unique components
set(vegHVDataBYrComparisons, NULL, grep("Unique", names(vegHVDataBYrComparisons)), NULL)

## sometimes HV1 is 2011, others is 2111 (probably because data wasn't sorted) -
## make sure there aren't duplicate comparisons
if (getOption("LandR.assertions")) {
  test <- unique(vegHVDataBYrComparisons[is.na(Volume_HV1_2011), .(scenario, rep, vegType)])
  test2 <- unique(vegHVDataBYrComparisons[!is.na(Volume_HV1_2011), .(scenario, rep, vegType)])
  test[, test := paste(scenario, rep, vegType)]
  test2[, test := paste(scenario, rep, vegType)]
  if (length(intersect(test$test, test2$test)) |
      length(intersect(test2$test, test$test))) {
    stop("There are duplicated HV comparisions between years\n",
         "(per scenario/rep/vegType combination)")
  }

  test <- any(is.na(vegHVDataBYrComparisons[is.na(Volume_HV1_2011), Volume_HV1_2111]))
  test2 <- any(is.na(vegHVDataBYrComparisons[!is.na(Volume_HV1_2011), Volume_HV1_2111]))
  test3 <- any(is.na(vegHVDataBYrComparisons[is.na(Volume_HV2_2011), Volume_HV2_2111]))
  test4 <- any(is.na(vegHVDataBYrComparisons[!is.na(Volume_HV2_2011), Volume_HV2_2111]))

  if ((isTRUE(test) | isFALSE(test2)) |
      (isTRUE(test3) | isFALSE(test4))) {
    stop("There are either duplicated or missing HV comparisions between years\n",
         "(per scenario/rep/vegType combination)")
  }
}

## add comparison type - because sometimes HV1 is 2011, others is 2111
## we need to use `sort(unique())` bellow
comp <- sort(unique(sub(".*_", "", grep("Volume", names(vegHVDataBYrComparisons), value = TRUE))))
comp <- paste(comp, collapse = "_")
vegHVDataBYrComparisons[, compare := comp]

## break into 2 tables, remove empty volluem columns,
## make column of comparison ID, change names and re-rbind
tempData  <- vegHVDataBYrComparisons[is.na(Volume_HV1_2011),]
tempData2  <- vegHVDataBYrComparisons[!is.na(Volume_HV1_2011),]

cols <- grep("Volume", names(tempData), value = TRUE)
set(tempData, NULL, cols[which(is.na(colSums(tempData[, ..cols])))], NULL)
set(tempData2, NULL, cols[which(is.na(colSums(tempData2[, ..cols])))], NULL)

setnames(tempData, new = sub("Volume_(HV)[0-9]_(.*)", "\\1_\\2", names(tempData)))
setnames(tempData2, new = sub("Volume_(HV)[0-9]_(.*)", "\\1_\\2", names(tempData2)))

vegHVDataBYrComparisons <- rbind(tempData, tempData2, use.names = TRUE)

vegHVData <- rbind(vegHVDataWYrComparisons, vegHVDataBYrComparisons, use.names = TRUE,
                   fill = TRUE)

## check reps/years/scenario per veg type
if (getOption("LandR.assertions")) {
  temp <- split(vegHVData[, .(scenario, year, rep, vegType)], by = "vegType", keep.by = FALSE)
  temp <- lapply(temp, FUN = function(x) paste0(x$scenario, x$year, x$rep))
  temp <- lapply(temp, FUN = unique)
  test <- lapply(1:length(temp), function(n) setdiff(temp[[n]], unlist(temp[-n])))
  test <- sapply(test, length)

  if (any(test)) {
    stop("vegHVData has different combinations of scenario, year, rep across vegTypes")
  }
}

## bind the two tables
allHVData <- rbindlist(list("fireHV" = fireHVData, "vegHV" = vegHVData),
                       use.names = TRUE, fill = TRUE, idcol = "HVtype")

## calculate overlap following Barros et al 2016
allHVData[, overlap := Intersection/Union]

## LABELS AND COLOURS FOR PLOTTING ----------------
vegTypeCNLabels <- unique(as.character(allHVData$vegType))
names(vegTypeCNLabels) <- vegTypeCNLabels
vegTypeCNLabels <- sub("PIEN", "Spruce", vegTypeCNLabels)
vegTypeCNLabels <- sub("MMC", "Moist conif.", vegTypeCNLabels)
vegTypeCNLabels <- sub("mixedwood", "Mixed", vegTypeCNLabels)
vegTypeCNLabels <- sub("broadleaf", "Deciduous", vegTypeCNLabels)
vegTypeCNLabels <- sub("PICO", "Pine", vegTypeCNLabels)
vegTypeCNLabels <- sub("DMCPSME", "Dry conif./Douglas-fir", vegTypeCNLabels)
vegTypeCNLabels <- sub("^PSME$", "Douglas-fir", vegTypeCNLabels)
vegTypeCNLabels <- sub("dryPSME", "Dry Douglas-fir", vegTypeCNLabels)

## reorder
vegTypeCNLabels <- vegTypeCNLabels[c(grep("No veg.|landscape", vegTypeCNLabels, invert = TRUE),
                                     grep("No veg.", vegTypeCNLabels),
                                     grep("landscape", vegTypeCNLabels))]

## landscape gets a different colour
vegTypeCNColours <- vegTypeCNLabels
vegTypeCNColours[1:length(vegTypeCNColours)-1] <- RColorBrewer::brewer.pal(length(vegTypeCNColours)-1, name = "Set1")
vegTypeCNColours["landscape"] <- "darkgreen"

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
summary(fireHVVolumeLandscape.lm)

## a bit heteroskedastic
sets <- par(mfrow = c(2,2))
plot(fireHVVolumeLandscape.lm)
par(sets)

fireHVVolumeLandscape.gls <- gls(Volume ~ scenario, weights = varIdent(form = ~ 1 | scenario),
                                 modelData[vegType == "landscape"])

tiff(file.path(figOutputPath, "fireHVVolumeLandscapeglsRESIDUALS.tiff"))
plot(fireHVVolumeLandscape.gls) ## much better :)
dev.off()

sink(file.path(statsOutputPath, "fireHVVolumeLandscapeglsSUMMARY.txt"))
summary(fireHVVolumeLandscape.gls)   ## same p-value estimate.
sink()

## model vegTypes - go directly to gls
fireHVVolumeVegTypes.gls <- gls(Volume ~ scenario * vegType, weights = varIdent(form = ~ 1 | scenario * vegType),
                                modelData[vegType != "landscape"])
anova(fireHVVolumeVegTypes.gls)   ## no significant interaction between scenario and vegtype. let's drop it
# Denom. DF: 522
# numDF  F-value p-value
# (Intercept)          1 52.52881  <.0001
# scenario             1 39.63523  <.0001
# vegType              8 78.44978  <.0001
# scenario:vegType     8  1.21105  0.2901

## glht can't see the missing levels, so we fit the model on a separate data table
modelData2 <- as.data.frame(modelData[vegType != "landscape"])
modelData2$scenario <- factor(modelData2$scenario)
modelData2$vegType <- factor(modelData2$vegType)
fireHVVolumeVegTypes.gls2 <- gls(Volume ~ scenario + vegType,
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
summary(vegHVVolumeLandscape.lm)


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

## evidence of heteroscedasticity
sets <- par(mfrow = c(2,2))
plot(vegHVVolumeVegTypes.lm)
par(sets)

modelData2 <- as.data.frame(modelData[vegType != "landscape"])
modelData2$scenario <- factor(modelData2$scenario)
modelData2$vegType <- factor(modelData2$vegType)
vegHVVolumeVegTypes.gls <- gls(Volume ~ scenario * vegType,
                               weights = varIdent(form = ~ 1 | scenario * vegType),
                               modelData2)

summary(glht(cell, linfct = K))

tiff(file.path(figOutputPath, "vegHVVolumeVegTypesglsRESIDUALS.tiff"))
plot(vegHVVolumeVegTypes.gls)
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

## heteroscedasticity
plot(vegHVOverlapLandscape.lm)


vegHVOverlapLandscape.gls <- gls(overlap ~ scenario, weights = varIdent(form = ~1 | scenario),
                                 data = modelData[vegType == "landscape"])


tiff(file.path(figOutputPath, "vegHVOverlapLandscapeglsRESIDUALS.tiff"))
plot(vegHVOverlapLandscape.gls) ## still not great
dev.off()


sink(file.path(statsOutputPath, "vegHVOverlapLandscapeglsSUMMARY.txt"))
summary(vegHVOverlapLandscape.gls)
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
plot(vegHVOverlapVegTypes.gls) ## not too bad
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


tiff(file.path(figOutputPath, "pyroVSbiodiversityLandscapelmRESIDUALS.tiff"))
sets <- par(mfrow = c(2,2))
plot(pyroVSbiodiversityLandscape.lm)  ## not too bad
par(sets)
dev.off()

sink(file.path(statsOutputPath, "pyroVSbiodiversityLandscapelmSUMMARY.txt"))
anova(pyroVSbiodiversityLandscape.lm)
cat("\n*********************\n")
summary(pyroVSbiodiversityLandscape.lm)
sink()

## by vetType
pyroVSbiodiversityVegTypes.lm <- lm(vegHV ~ fireHV*scenario*vegType, data = modelData[vegType != "landscape"])

## very heteroscedastic
sets <- par(mfrow = c(2,2))
plot(pyroVSbiodiversityVegTypes.lm)  ## not too bad
par(sets)


modelData2 <- as.data.frame(modelData[vegType != "landscape"])
modelData2$scenario <- factor(modelData2$scenario)
modelData2$vegType <- factor(modelData2$vegType)

pyroVSbiodiversityVegTypes.gls <- gls(vegHV ~ fireHV*scenario*vegType,
                                      weights = varIdent(form = ~ 1 | scenario * vegType),
                                      data = modelData2)
tiff(file.path(figOutputPath, "pyroVSbiodiversityVegTypeslmRESIDUALS.tiff"))
plot(pyroVSbiodiversityVegTypes.gls) ## not great but best structure possible (tried all combos in varIdent)
dev.off()

sink(file.path(statsOutputPath, "pyroVSbiodiversityVegTypeslmSUMMARY.txt"))
anova(pyroVSbiodiversityVegTypes.gls)
cat("\n*********************\n")
summary(pyroVSbiodiversityVegTypes.gls)
cat("\n*********************\n")
emtrends(pyroVSbiodiversityVegTypes.gls, pairwise ~ scenario | vegType,
         var = "fireHV", mode = "df.error", data = modelData2)
emtrends(pyroVSbiodiversityVegTypes.gls, pairwise ~ vegType | scenario,
         var = "fireHV", mode = "df.error", data = modelData2)
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

HVvolumeVegTypesPlot <- ggplot(plotData[vegType != "landscape"],
                               aes(x = vegType, y = Volume, alpha = scenario, fill = vegType)) +
  geom_boxplot() +
  scale_x_discrete(labels = vegTypeCNLabels) +
  scale_fill_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
  scale_alpha_manual(values = c("HV_noPM" = 0.4, "HV_PM" = 1.0),
                     labels = c("HV_noPM" = "no PM", "HV_PM" = "PM")) +
  theme_pubr(base_size = 12, x.text.angle = 30, margin = FALSE) +
  theme(legend.box = "vertical",
        strip.background = element_blank()) +
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
  theme(strip.background = element_blank()) +
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
        strip.background = element_blank()) +
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
        strip.background = element_blank()) +
  labs(x = "", y = "Overlap", fill = "", alpha = "") +
  guides(alpha = guide_legend(override.aes = list(fill = "grey50")))

HVOverlapLandscapePlot <- ggplot(plotData[vegType == "landscape"],
                                 aes(x = vegType, y = overlap, fill = scenario)) +
  geom_boxplot() +
  scale_x_discrete(labels = vegTypeCNLabels) +
  scale_fill_brewer(palette = "Greys",
                    labels = c("noPM" = "no PM", "PM" = "PM")) +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(strip.background = element_blank()) +
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
pyroVSbiodiversityLandscape.lm <- lm(vegHV ~ fireHV*scenario, data = plotData[vegType == "landscape"])
pyroVSbiodiversityVegTypes.gls <- gls(vegHV ~ fireHV*scenario*vegType,
                                      weights = varIdent(form = ~ 1 | scenario * vegType),
                                      data = plotData2)

plotData[vegType == "landscape", pred := predict(pyroVSbiodiversityLandscape.lm)]
plotData[vegType != "landscape", pred := predict(pyroVSbiodiversityVegTypes.gls)]

pyroVSbioDivVegTypesPlot <- ggplot(plotData[vegType != "landscape"],
                                   aes(x = fireHV, y = vegHV,
                                       linetype = scenario, shape = scenario,
                                       colour = vegType)) +
  geom_point() +
  geom_line(aes(y = pred)) +
  scale_colour_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
  scale_linetype_discrete(labels = c("HV_noPM" = "no PM", "HV_PM" = "PM")) +
  scale_shape_discrete(labels = c("HV_noPM" = "no PM", "HV_PM" = "PM")) +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(legend.box = "vertical",
        strip.background = element_blank()) +
  labs(x = "pyrodiversity", y = "forest diversity", colour = "", linetype = "", shape = "") +
  facet_wrap( ~ vegType, labeller = labeller(vegType = vegTypeCNLabels),
              scales = "free")

pyroVSbioDivVegLandscapePlot <- ggplot(plotData[vegType == "landscape"],
                                       aes(x = fireHV, y = vegHV, shape = scenario,
                                           linetype = scenario, colour = vegType)) +
  geom_point() +
  geom_line(aes(y = pred)) +
  scale_colour_manual(labels = vegTypeCNLabels, values = vegTypeCNColours) +
  scale_linetype_discrete(labels = c("HV_noPM" = "no PM", "HV_PM" = "PM")) +
  scale_shape_discrete(labels = c("HV_noPM" = "no PM", "HV_PM" = "PM")) +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(legend.box = "vertical",
        strip.background = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5)) +
  labs(x = "pyrodiversity", y = "forest diversity", title = "landscape",
       colour = "", linetype = "", shape = "")

plotSave <- ggarrange(pyroVSbioDivVegTypesPlot, pyroVSbioDivVegLandscapePlot + labs(y = ""),
                      ncol = 2, nrow = 1, widths = c(1, 0.6),
                      common.legend = TRUE, legend = "bottom")

ggsave(plot = plotSave, filename = file.path(figOutputPath, "pyroVsbiodiversity.tiff"),
       width = 12, height = 7)
