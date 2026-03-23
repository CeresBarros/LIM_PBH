## --------------------------------------------------
##  HYPERVOLUMES OF PYRODIVERSITY AND BIODIVERSITY
##  ANALYSING PCAs
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

Require(c("data.table",
          "FD",
          "ggfortify",
          "ggplot2",
          "ggpubr",
          "ggvegan",
          # "LandR (==1.0.7.9026)", ## needed but don't load -- here for install.
          "reproducible",
          "SpaDES",
          "ToolsCB",
          "vegan"
), install = FALSE)

source("R/R_tools/convertToCNVegType.R")
source("R/R_tools/Useful_functions.R")
source("R/R_tools/hypervolumesHelpers.R")
source("R/R_tools/utilsForResultsAnalyses.R")
options("reproducible.useNewDigestAlgorithm" = 2)
options("spades.moduleCodeChecks" = FALSE)
options("reproducible.useCache" = TRUE)
options("reproducible.destinationPath" = normPath("R/SpaDES/inputs"))
options("reproducible.useGDAL" = FALSE)


## general paths
# simDirName <- "jun2021Runs"
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
figOutputPath <- file.path(simPaths$outputPath, "figuresAnalysis", "Montane")
HVoutputPath <- file.path(simPaths$outputPath, "hypervolumes")
# bw.outputPath <- file.path(HVoutputPath, "bwTest")

## merge Douglas-fir/dry-conifer stands?
mergeDMCPSME <- FALSE  ## merge DMCPSME PSME dryPSME
mergePSME <- TRUE ## merge PSME dryPSME
if (mergeDMCPSME) {
  figOutputPath <- file.path(simPaths$outputPath, "figuresAnalysis/Montane/mergeDMCPSME")
}
if (mergePSME) {
  figOutputPath <- file.path(simPaths$outputPath, "figuresAnalysis/Montane/mergePSME")
}
dir.create(figOutputPath, recursive = TRUE)

## LOAD DATA (RESULTS)  ---------------------
yearSubset <- unique(as.integer(c(seq(3511, 4011, 5), 4011)))
runPrepResultsModule <- FALSE
source("R/SpaDES/simResultsDataPrep.R")

## PREP FIRE AND VEG DATA -------------------
## Fire properties (fire patch size in pixels, fire frequency, fire severity as biomass loss)
opts <- options("LandR.assertions" = FALSE)

opts <- options(reproducible.cachePath = simPaths$cachePath)
source("R/R_tools/prepFireData4HVs.R")

yearSamples <- sample5SimYears(allPixelCohortDataMnt[, .(year, rep)])   ## seed ensures same years are drawn
useFirstLastYear <- FALSE
source("R/R_tools/prepVegData4HVs.R")
options(opts)

## get labels and colours
source("R/R_tools/plotLabels&Cols.R")

## don't need these
rm(allPixelBurnData, allPixelCohortData)
gc(reset = TRUE)

## CALCULATE CWM TRAIT VALUES ----------------
## averages across reps (species abundances were averaged across rep per pixel)
## for categorical/ordinal traits the dominant trait value (highest abundance) was taken
traitsTable <- preSimList$species
traitsTable[, firetolerance := as.ordered(firetolerance)]
traitsTable <- data.frame(traitsTable[, .(longevity, shadetolerance, firetolerance, postfireregen)],
                          row.names = traitsTable$speciesCode)
rm(preSimList)  ## not needed anymore
gc(reset = TRUE)

if (useFirstLastYear) {
  vegData <- allPixelCohortDataMnt[year %in% c(min(yearSubset), max(yearSubset))]
} else {
  if (exists("yearSamples")) {
    vegData <- allPixelCohortDataMnt[yearSamples, on = .(year, rep)]
  } else {
    vegData <- allPixelCohortDataMnt[year %in% yearSubset]
  }
}

cols <- c("speciesCode", "scenario", "rep", "year", "pixelIndex", "B")   ## keep rep for wrapper.
vegData <- vegData[, ..cols]
## sum B across cohorts
vegData <- vegData[, list(B = sum(B)), by = c("speciesCode", "scenario", "rep", "year", "pixelIndex")]
vegData <- dcast.data.table(vegData, as.formula("... ~ speciesCode"), value.var = "B")

if (!useFirstLastYear) {
  ## if not comparing years, make all years equal to last
  vegData[, year := max(year)]
}

spp <- rownames(traitsTable)
combos <- unique(vegData[, .(scenario, rep, year)])

cacheExtra <- CacheDigest(list(vegData = vegData,
                               traitsTable = traitsTable,
                               spp = spp))
## faster to do it in chunks
traitCWMs <- Cache(Map,
                   f = calcTraitCWMs,
                   scen = combos$scenario,
                   rep = combos$rep,
                   yr = combos$year,
                   MoreArgs = list(vegData = vegData,
                                   traitsTable = traitsTable,
                                   spp = spp),
                   cacheRepo = simPaths$cachePath,
                   userTags = c("traitCWMs"),
                   omitArgs = c("userTags", "MoreArgs"),
                   .cacheExtra = cacheExtra)

traitCWMs <- rbindlist(traitCWMs)

## add the 0 rows again, and replace NAs with 0s
traitCWMs <- traitCWMs[vegData[, .(scenario, rep, year, pixelIndex)],
                       on = .(scenario, rep, year, pixelIndex)]
traitCWMs <- replaceNAs(traitCWMs, val = 0)

## LOAD PCAs ---------------------------
vegHVPCA <- grep("vegHVs", list.files(HVoutputPath, "OrdinationObj", full.names = TRUE), value = TRUE)
vegHVPCA <- readRDS(vegHVPCA)
vegHVPCAscores <- cbind(as.data.table(vegHVPCA$x), vegDataForHVs)

## get highest/lowest 3 factor loadings from PCA
loadings_coords <- NULL
for (i in 1:nrow(vegHVPCA$rotation)) {
  loadings_coords <- rbind(loadings_coords, rbind(c(0,0,0), vegHVPCA$rotation[i, 1:3]))
}
loadings_coords <- as.data.table(cbind(as.data.frame(loadings_coords), rep(rownames(vegHVPCA$rotation), each = 2)))
names(loadings_coords)[4] <- "Var"

## rename some of the variables for plotting
loadings_coords$Var <- sub("meanStandAge", "mean_age", loadings_coords$Var)
loadings_coords$Var <- sub("sdStandAge", "sd_age", loadings_coords$Var)

## most influential variables across PCs
Vars <- unique(unlist(loadings_coords[, lapply(.SD, function(x) which(abs(x) >= 0.5)), .SDcols = c("PC1", "PC2", "PC3")]))
Vars <- loadings_coords$Var[Vars]

## clean wd
rm(summaryFireAttributes)
gc(reset = TRUE)

## FIT TRAIT VECTORS TO PCA ------------
## treat all ordered vectors as continuous
traitCWMs[, firetoleranceCont := as.numeric(firetolerance)]
trait.fit <- Cache(
  envfit,
  ord = vegHVPCA,
  env = traitCWMs[, .(longevity, shadetolerance, firetoleranceCont, postfireregen)],
  choices = c(1:3),
  cacheRepo = simPaths$cachePath,
  userTags = c("traitCWMs", "envfitPCA"),
  omitArgs = "userTags")


## get vector coordinates (adding the origin to each vector value)
traits_coords <- NULL  ## automatically scales by correlation
for (i in 1:nrow(fortify(trait.fit))) {
  traits_coords <- rbind(traits_coords, rbind(c(0,0,0), as.data.frame(fortify(trait.fit))[i, c("PC1", "PC2", "PC3")]))
}
traits_coords <- as.data.table(traits_coords)
## rename traits
traits_coords$Label <- rep(fortify(trait.fit)$Label, each = 2)
traits_coords[, Label := sub("postfireregen", "regen_", Label)]
traits_coords[, Label := sub("Cont", "", Label)]
traits_coords[, Label := sub("tolerance", "_tol.", Label)]

## all have very low correlations, so use all
trts <- unique(unlist(traits_coords[, lapply(.SD, function(x) which(abs(x) >= 0)), .SDcols = c("PC1", "PC2", "PC3")]))
trts <- na.omit(traits_coords$Label[trts])
trts <- trts[trts != "regen_0"]  ## exclude this - comes from postfireregen in pixels where there is no B

## -----------------------------------------------
## HYPERVOLUMES 3D PLOTS WITH LOADINGS AND TRAITS

## loop through all vegTypes
vegTypes <- as.character(unique(vegHVPCAscores$vegTypeCN))

if (mergeDMCPSME) {
  vegTypes[grep("PSME", vegTypes)] <- "PSME"
  vegTypes <- unique(vegTypes)
}
if (mergePSME) {
  vegTypes[grep("^dryPSME|^PSME", vegTypes)] <- "PSME"
  vegTypes <- unique(vegTypes)
}

vegTypes <- vegTypes[vegTypes != "No veg."]

pixelIndexDT <- vegHVPCAscores[year == max(yearSubset), .(rep, scenario, vegTypeCN, pixelIndex)]
pixelIndexDT <- unique(pixelIndexDT)

if (useFirstLastYear) {
  startYear <- min(yearSubset)
  endYear <- max(yearSubset)
  cols <- c("black", "black", scales::hue_pal()(2)[1], scales::hue_pal()(2)[2])
  centroid.cols <- c("grey", "grey", "red", "blue")
} else {
  startYear <- NULL
  endYear <- NULL
  cols <- c(scales::hue_pal()(2)[1], scales::hue_pal()(2)[2])
  centroid.cols <- c("red", "blue")
}

rep2Plot <- unique(pixelIndexDT$rep)[1]
vegHVPCAscores2 <- copy(vegHVPCAscores)
pixelIndexDT2 <- copy(pixelIndexDT)

## plot all reps by assigning the same rep to all
plotAllreps <- TRUE
if (plotAllreps) {
  vegHVPCAscores2[, rep := rep2Plot]
  pixelIndexDT2[, rep := rep2Plot]
  figOutputPath2 <- file.path(figOutputPath, "HVplotsAllreps")
  dir.create(figOutputPath2, showWarnings = FALSE)
} else {
  figOutputPath2 <- file.path(figOutputPath, paste0("HVplots_rep", rep2Plot))
  dir.create(figOutputPath2, showWarnings = FALSE)
}

lapply(vegTypes, FUN = plotHVs3DWrapper,
       vegHVPCAscores = vegHVPCAscores2[rep == rep2Plot],
       pixelIndexDT = pixelIndexDT2[rep == rep2Plot],
       vegTypeCNLabels = vegTypeCNLabels,
       figOutputPath = figOutputPath2,
       cacheRepo = simPaths$cachePath,
       mergeVegType = "mergePSME",
       startYear = startYear,
       endYear = endYear,
       colsHV = c("PC1", "PC2", "PC3", "PC4"),
       ## plotHypervolumes3D args:
       loadings_coords = as.data.frame(loadings_coords[Var %in% Vars, .(PC1, PC2, PC3)]) * ordiArrowMul(vegHVPCA, fill = 2, choices = 1:3, display = "species"),
       PHvect_coords = as.data.frame(traits_coords[Label %in% trts, .(PC1, PC2, PC3)] * ordiArrowMul(trait.fit, choices = 1:3)), ## ordiArrowMul finds the appropriate multiplifer to plot axes.
       loadings_labels = loadings_coords[Var %in% Vars, Var],
       PHvect_labels = traits_coords[Label %in% trts, Label],
       show.random = TRUE,
       show.data = FALSE,
       show.legend = TRUE,
       cex.axis = 1,
       cex.lab = 1,
       cex.random = 0.5,
       cex.centroid = 1.7,
       lwd = 2,
       colors = cols,
       centroid.cols = centroid.cols,
       grid = FALSE,
       box = TRUE,
       names = c("PC1\n", "PC2\n", "\nPC3"),
       # limits = c(round(min(vegHVPCAscores2[, c("PC1", "PC2", "PC3", "PC4")]), 2) - 0.1,
       #            round(max(vegHVPCAscores2[, c("PC1", "PC2", "PC3", "PC4")]), 2) + 0.1),
       # limits = c(-6, 6),
       y.margin.add = 0.6,
       angle = 50,
       pch = 16)


## HOW SIMILAR ARE DOUG-FIR FOREST TYPES?

## If we look at the first year (top row), there seems to be an indication that they
## could all overlap a similar region in space, and that the fact that some have a smaller
## sample size may be not showing how much these forest types are overlapping in terms of forest composition.
## However, if we analyse this across years, it seems that DMCPSME may have shifted less than the other two.
## If we were to group all three the between-year shift would be overall smaller (and overlap larger),
## because both the start and end HVs would be larger to start with
## based on these I think we can merge the pure PSME stands and leave the dry mixed conifer/PSME separate

## average scores across reps
cols <- grep("^PC", names(vegHVPCAscores), value = TRUE)

colsBy <- c("scenario", "rep", "pixelIndex", "vegTypeCN")
if (useFirstLastYear) {
  colsBy <- c(colsBy, "year")
}
plotData <- vegHVPCAscores[, lapply(.SD, mean), .SDcols = cols,
                           by = colsBy]
ggplot(plotData[grepl("PSME", vegTypeCN)], aes(PC1, PC2, colour = vegTypeCN, shape = scenario)) +
  geom_point(alpha = 0.5) +
  if (useFirstLastYear) {
    facet_wrap(scenario ~ year)
  } else {
    facet_wrap(~ scenario)
  }


ggplot(plotData[grepl("PSME", vegTypeCN) & scenario == "PM"], aes(PC1, PC3, colour = as.factor(rep), shape = scenario)) +
  geom_point(alpha = 0.5) +
  if (useFirstLastYear) {
    facet_wrap(year ~ vegTypeCN)
  } else {
    facet_wrap(~ vegTypeCN)
  }


ggplot(plotData[grepl("PSME", vegTypeCN) & scenario == "PM"], aes(PC1, PC2, colour = as.factor(rep), shape = scenario)) +
  geom_point(alpha = 0.5) +
  if (useFirstLastYear) {
    facet_wrap(year ~ vegTypeCN)
  } else {
    facet_wrap(~ vegTypeCN)
  }

## DOUGLAS-FIR
plotData$year <- as.character(plotData$year)
plotDataPCA <- vegHVPCA
plotDataPCA$x <- plotDataPCA$x[grepl("PSME", plotData$vegTypeCN),]
if (useFirstLastYear) {
  col <- "year"
  shape <- "year"
} else {
  col <- "vegTypeCN"
  shape <- "vegTypeCN"
}

PSMEPCAplot <- autoplot(plotDataPCA, data = plotData[grepl("PSME", vegTypeCN)],
                        colour = col, shape = shape,
                        label.colour = "blue",
                        loadings = TRUE, loadings.colour = "blue", scale = 0,
                        loadings.label = TRUE, loadings.label.size = 3, alpha = 0.3) +
  # scale_x_continuous(expand = expansion(add = 0.1)) +
  theme_pubr(base_size = 12,  margin = FALSE) +
  theme(legend.position = "right",
        panel.grid.major.y = element_line(colour = "grey", size = 11/22, linetype = "dotted")) +
  labs(colour = "", shape = "") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  facet_wrap(~ scenario)


