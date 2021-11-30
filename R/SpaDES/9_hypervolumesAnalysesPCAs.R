## --------------------------------------------------
##  HYPERVOLUMES OF PYRODIVERSITY AND BIODIVERSITY
##  ANALYSING PCAs
## --------------------------------------------------

library(SpaDES)
library(ToolsCB)
library(data.table)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(FD)
library(vegan)

source("R/R_tools/convertToCNVegType.R")
source("R/R_tools/Useful_functions.R")
source("R/R_tools/hypervolumesHelpers.R")
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
figOutputPath <- file.path(simPaths$outputPath, "figuresAnalysis")
HVoutputPath <- file.path(simPaths$outputPath, "hypervolumes")
# bw.outputPath <- file.path(HVoutputPath, "bwTest")

## LOAD DATA (RESULTS)  ---------------------
yearSubset <- c(seq(2011, 2111, 5), 2111)
source("R/SpaDES/6_resultsDataPrep.R")


## PREP FIRE AND VEG DATA -------------------
## Fire properties (fire patch size in pixels, fire frequency, fire severity as biomass loss)
opts <- options("LandR.assertions" = FALSE)
source("R/R_tools/prepFireData4HVs.R")
source("R/R_tools/prepVegData4HVs.R")
options(opts)
## get labels and colours
source("R/R_tools/plotLabels&Cols.R")

## CALCULATE CWM TRAIT VALUES ----------------
## averages across reps (species abundances were averaged across rep per pixel)
## for categorical/ordinal traits the dominant trait value (highest abundance) was taken

traitsTable <- preSimList$species
traitsTable[, firetolerance := as.ordered(firetolerance)]
traitsTable <- data.frame(traitsTable[, .(longevity, shadetolerance, firetolerance, postfireregen)],
                          row.names = traitsTable$speciesCode)

vegData <- allPixelCohortDataMnt[year %in% c(start(preSimList), end(preSimList))]
cols <- c("speciesCode", "scenario", "rep", "year", "pixelIndex", "B")   ## keep rep for wrapper.
vegData <- vegData[, ..cols]
## sum B across cohorts
vegData <- vegData[, list(B = sum(B)), by = c("speciesCode", "scenario", "rep", "year", "pixelIndex")]
vegData <- dcast.data.table(vegData, as.formula("... ~ speciesCode"), value.var = "B")

spp <- rownames(traitsTable)
combos <- unique(vegData[, .(scenario, rep, year)])

## faster to do it in chunks
traitCWMs <- Cache(Map,
                   f = function(scen, rep, yr, vegData, traitsTable, spp) {
                     vegData2 <- vegData[rowSums(vegData[, ..spp]) != 0,]
                     vegData2 <- vegData2[scenario == scen & year == yr & rep == rep]
                     vegData2 <- vegData2[, lapply(.SD, mean), .SDcols = spp,
                                        by = .(pixelIndex)]

                     traitCWMs <- functcomp(traitsTable,
                                            as.matrix(data.frame(vegData2[, ..spp])),
                                            CWM.type = "dom")
                     traitCWMs$scenario <- scen
                     traitCWMs$rep <- rep
                     traitCWMs$year <- yr
                     traitCWMs$pixelIndex <- vegData2$pixelIndex

                     return(as.data.table(traitCWMs))
                   },
                   scen = combos$scenario,
                   rep = combos$rep,
                   yr = combos$year,
                   MoreArgs = list(vegData = vegData,
                                   traitsTable = traitsTable,
                                   spp = spp),
                   cacheRepo = simPaths$cachePath,
                   userTags = c("traitCWMs"),
                   omitArgs = c("userTags"))

traitCWMs <- rbindlist(traitCWMs)

## add the 0 rows again, and replace NAs with 0s
traitCWMs <- traitCWMs[vegData[, .(scenario, rep, year, pixelIndex)],
                       on = .(scenario, rep, year, pixelIndex)]
traitCWMs <- replaceNAs(traitCWMs, val = 0)

## LOAD PCAs ---------------------------
vegHVPCA <- grep("vegHVs", list.files(HVoutputPath, "OrdinationObj", full.names = TRUE), value = TRUE)
vegHVPCA <- readRDS(vegHVPCA)
vegHVPCAscores <- cbind(as.data.table(vegHVPCA$x), vegDataForHVs)

## clean wd
rm(summaryFireAttributes)
gc()

## FIT TRAIT VECTORS TO PCA ------------
trait.fit <- Cache(
  envfit,
  ord = vegHVPCA,
  env = traitCWMs,
  choices = c(1:3),
  cacheRepo = simPaths$cachePath,
  userTags = c("traitCWMs", "envfitPCA"),
  omitArgs = "userTags")


## get vector coordinates (adding the origin to each vector value)
trait_coords <- NULL
traits_loads <- as.data.frame(trait.fit$vectors$arrows*sqrt(trait.fit$vectors$r))
for (i in 1:nrow(traits_loads)) {
  trait_coords <- rbind(trait_coords, rbind(c(0,0,0,0), traits_loads[i,1:4]))
}
trait_coords <- cbind(as.data.frame(trait_coords), rep(rownames(traits_loads), each = 2))
colnames(trait_coords)[4] <- "Trait"

## get highest environmental loadings from PCA
trts <- as.character(trait_coords$Trait[abs(trait_coords$PC1) >= 0.8])

## -----------------------------------------------
## HYPERVOLUMES 3D PLOTS WITH LOADINGS AND TRAITS
colsHV <- c("PC1", "PC2", "PC3", "PC4")
userTags <- c("hypervolume", "PSME", "averagePlots", "2011", "noPM")
HV2011_noPM  <- Cache(hypervolume::hypervolume,
                      data = plotData[grepl("^dryPSME|^PSME", vegTypeCN) & year == "2011" & scenario == "noPM", ..colsHV],
                      method = "svm",
                      .cacheExtra = paste(userTags, collapse = "_"),
                      cacheRepo = simPaths$cachePath,
                      userTags = userTags,
                      omitArgs = c("userTags"))

userTags <- c("hypervolume", "PSME", "averagePlots", "2111", "noPM")
HV2111_noPM  <- Cache(hypervolume::hypervolume,
                      data = plotData[grepl("^dryPSME|^PSME", vegTypeCN) & year == "2111" & scenario == "noPM", ..colsHV],
                      method = "svm",
                      .cacheExtra = paste(userTags, collapse = "_"),
                      cacheRepo = simPaths$cachePath,
                      userTags = userTags,
                      omitArgs = c("userTags"))

userTags <- c("hypervolume", "PSME", "averagePlots", "2011", "PM")
HV2011_PM  <- Cache(hypervolume::hypervolume,
                    data = plotData[grepl("^dryPSME|^PSME", vegTypeCN) & year == "2011" & scenario == "PM", ..colsHV],
                    method = "svm",
                    .cacheExtra = paste(userTags, collapse = "_"),
                    cacheRepo = simPaths$cachePath,
                    userTags = userTags,
                    omitArgs = c("userTags"))

userTags <- c("hypervolume", "PSME", "averagePlots", "2111", "PM")
HV2111_PM  <- Cache(hypervolume::hypervolume,
                    data = plotData[grepl("^dryPSME|^PSME", vegTypeCN) & year == "2111" & scenario == "PM", ..colsHV],
                    method = "svm",
                    .cacheExtra = paste(userTags, collapse = "_"),
                    cacheRepo = simPaths$cachePath,
                    userTags = userTags,
                    omitArgs = c("userTags"))

HVls_noPM <- hypervolume::hypervolume_join(HV2011_noPM, HV2111_noPM)
HVls_PM <- hypervolume::hypervolume_join(HV2011_PM, HV2111_PM)

plotHypervolumes3D(HVls_noPM, loadings_coords = , PHvect_coords = ,
                   show.random = TRUE, show.data = FALSE,
                   show.legend = FALSE, cex.axis = 1, cex.lab = 1.5, cex.random = 0.5, cex.centroid = 1,
                   show.contour = TRUE,
                   colors = c("black", scales::hue_pal()(2)[1]), centroid.cols = rep("blue", 3), grid = FALSE, box = TRUE,
                   names = c("PC1\n", "PC2\n", "\nPC3"), limits = c(-1, 1), y.margin.add = 0.6,
                   angle = 40, pch = 16)

randomPoints <- rbindlist(list(
  HV2011_noPM = as.data.table(HV2011_noPM@RandomPoints),
  HV2111_noPM = as.data.table(HV2111_noPM@RandomPoints),
  HV2011_PM = as.data.table(HV2011_PM@RandomPoints),
  HV2111_PM = as.data.table(HV2111_PM@RandomPoints)), idcol = "yrScen")

dataPoints <- rbindlist(list(
  HV2011_noPM = as.data.table(HV2011_noPM@Data),
  HV2111_noPM = as.data.table(HV2111_noPM@Data),
  HV2011_PM = as.data.table(HV2011_PM@Data),
  HV2111_PM = as.data.table(HV2111_PM@Data)), idcol = "yrScen")

allPoints <- rbindlist(list(randomPoints = randomPoints, dataPoints = dataPoints), idcol = "type")
allPoints <- cbind(allPoints, as.data.table(do.call(rbind, strsplit(allPoints$yrScen, split = "_"))))
setnames(allPoints, c("V1", "V2"), c("year", "scen"))
set(allPoints, j = "yrScen", value = NULL)



## HOW SIMILAR ARE DOUG-FIR FOREST TYPES?
## If we look at the first year (top row), there seems to be an indication that they
## could all overlap a similar region in space, and that the fact that some have a smaller
## sample size may be not showing how much these forest types are overlapping in terms of forest composition.
## However, if we analyse this across years, it seems that DMCPSME may have shifted less than the other two.
## If we were to group all three the between-year shift would be overall smaller (and overlap larger),
## because both the 2011 and 2111 HVs would be larger to start with
## based on these I think we can merge the pure PSME stands and leave the dry mixed conifer/PSME separate

## average scores across reps
cols <- grep("^PC", names(vegHVPCAscores), value = TRUE)
plotData <- vegHVPCAscores[, lapply(.SD, mean), .SDcols = cols,
                           by = .(scenario, rep, pixelIndex, vegTypeCN, year)]
ggplot(plotData[grepl("PSME", vegTypeCN)], aes(PC1, PC2, colour = vegTypeCN, shape = scenario)) +
  geom_point(alpha = 0.5) +
  facet_wrap(scenario ~ year)

ggplot(plotData[grepl("PSME", vegTypeCN) & scenario == "PM"], aes(PC1, PC3, colour = as.factor(rep), shape = scenario)) +
  geom_point(alpha = 0.5) +
  facet_wrap(year ~ vegTypeCN)

ggplot(plotData[grepl("PSME", vegTypeCN) & scenario == "PM"], aes(PC1, PC4, colour = as.factor(rep), shape = scenario)) +
  geom_point(alpha = 0.5) +
  facet_wrap(year ~ vegTypeCN)

## DOUGLAS-FIR

plotData$year <- as.character(plotData$year)


plotDataPCA <- vegHVPCA
plotDataPCA$x <- plotDataPCA$x[grepl("PSME", plotData$vegTypeCN),]
PSMEPCAplot <- autoplot(vegHVPCA2, data = plotData[grepl("PSME", vegTypeCN)],
                        colour = "year", shape = "year", label.colour = "blue",
                        loadings = TRUE, loadings.colour = "blue", scale = 0,
                        loadings.label = TRUE, loadings.label.size = 3, alpha = 0.3) +
  scale_x_continuous(expand = expansion(add = 0.1)) +
  theme_pubr(base_size = 12,  margin = FALSE) +
  theme(legend.position = "right",
        panel.grid.major.y = element_line(colour = "grey", size = 11/22, linetype = "dotted")) +
  labs(colour = "", shape = "") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  facet_wrap(~ scenario)


