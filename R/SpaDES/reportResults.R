## clean workspace
rm(list=ls()); amc::.gc()

library(data.table)
library(ggplot2)
library(ggpubr)
library(raster)
library(quickPlot)
library(SpaDES)
library(qs)
library(purrr)
library(magick)
library(LandR)
library(future)
library(future.apply)

source("R/R_tools/convertToCNVegType.R")

simPaths <- list(cachePath = file.path("R/SpaDES/cache/AI_report")
                 , modulePath = file.path("R/SpaDES/m")
                 , inputPath = file.path("R/SpaDES/inputs")
                 , outputPath = file.path("R/SpaDES/outputs"))

## path to figure folder and cache folder
figOutputPath <- file.path(simPaths$outputPath, "figuresAnalysis")
dir.create(figOutputPath)
cPath <- file.path(simPaths$cachePath, "postSimAnalyses")

## GET CAMERON'S AGE DATA AND STAND VEG TYPES
ageDataCN <- fread("data/CameronsAgeData/treelist_outputs_for Ceres.csv")
patchVegTypeCN <- fread("data/CameronsAgeData/patch outputs_for Ceres.csv")

ageDataCN <- patchVegTypeCN[ageDataCN, on = .(Patch.ID)]
ageDataCN$Cover.dendro <- sub("Mixedwood", "mixedwood", ageDataCN$Cover.dendro)
ageDataCN$Cover.dendro <- sub("Broadleaf", "broadleaf", ageDataCN$Cover.dendro)
ageDataCN$Cover.dendro <- sub("-", "", ageDataCN$Cover.dendro)
rm(patchVegTypeCN)

## remove a record that seems funky (maybe it's a new cohort?)
ageDataCN <- ageDataCN[Reconstructed.age != 2018]

## calculate no. cohorts:
ageDataCN[, noCohorts := length(unique(Est.bin)) , by = .(Patch.ID)]

## GET OUTPUTS FOLDERS FOR EACH SCENARIO
outputs_PM <- list.dirs(simPaths$outputPath, full.names = TRUE, recursive = FALSE) %>%
  grep("/PM_rep", ., value = TRUE)
outputs_noPM <- list.dirs(simPaths$outputPath, full.names = TRUE, recursive = FALSE) %>%
  grep("/noPM_rep", ., value = TRUE)

## GET FILE NAMES
cohortDataFiles <- c(sapply(outputs_PM, FUN = function(x) list.files(x, pattern = "cohortData", full.names = TRUE)),
                     sapply(outputs_noPM, FUN = function(x) list.files(x, pattern = "cohortData", full.names = TRUE))) %>%
  unique(.)

rstCurrentBurnFiles <- c(sapply(outputs_PM, FUN = function(x) list.files(x, pattern = "rstCurrentBurn", full.names = TRUE)),
                         sapply(outputs_noPM, FUN = function(x) list.files(x, pattern = "rstCurrentBurn", full.names = TRUE))) %>%
  unique(.)

pixelGroupMapFiles <- c(sapply(outputs_PM, FUN = function(x) list.files(x, pattern = "pixelGroupMap", full.names = TRUE)),
                        sapply(outputs_noPM, FUN = function(x) list.files(x, pattern = "pixelGroupMap", full.names = TRUE))) %>%
  unique(.)

vegTypeMapFiles <- c(sapply(outputs_PM, FUN = function(x) list.files(x, pattern = "vegTypeMap", full.names = TRUE)),
                     sapply(outputs_noPM, FUN = function(x) list.files(x, pattern = "vegTypeMap", full.names = TRUE))) %>%
  unique(.)

## load rasters as stacks
rstCurrentBurnStk_noPM <- lapply(grep("noPM", rstCurrentBurnFiles, value = TRUE), readRDS) %>%
  stack(.)
rstCurrentBurnStk_PM <- lapply(grep("noPM", rstCurrentBurnFiles, value = TRUE, invert = TRUE), readRDS) %>%
  stack(.)
pixelGroupMapStk_noPM <- lapply(grep("noPM", pixelGroupMapFiles, value = TRUE), readRDS) %>%
  stack(.)
pixelGroupMapStk_PM <- lapply(grep("noPM", pixelGroupMapFiles, value = TRUE, invert = TRUE), readRDS) %>%
  stack(.)
vegTypeMapStk_noPM <- lapply(grep("noPM", vegTypeMapFiles, value = TRUE), readRDS) %>%
  stack(.)
vegTypeMapStk_PM <- lapply(grep("noPM", vegTypeMapFiles, value = TRUE, invert = TRUE), readRDS) %>%
  stack(.)

names(rstCurrentBurnStk_noPM) <- paste(sub(".*year", "year", sub("\\.rds", "", grep("noPM", rstCurrentBurnFiles, value = TRUE))),
                                       sub(".*(rep)([0-9]+)/.*", "\\1\\2", grep("noPM", rstCurrentBurnFiles, value = TRUE)), sep = "_")
names(rstCurrentBurnStk_PM) <- paste(sub(".*year", "year", sub("\\.rds", "", grep("noPM", rstCurrentBurnFiles, value = TRUE, invert = TRUE))),
                                     sub(".*(rep)([0-9]+)/.*", "\\1\\2", grep("noPM", rstCurrentBurnFiles, value = TRUE, invert = TRUE)), sep = "_")

names(pixelGroupMapStk_noPM) <- paste(sub(".*year", "year", sub("\\.rds", "", grep("noPM", pixelGroupMapFiles, value = TRUE))),
                                      sub(".*(rep)([0-9]+)/.*", "\\1\\2", grep("noPM", pixelGroupMapFiles, value = TRUE)), sep = "_")
names(pixelGroupMapStk_PM) <- paste(sub(".*year", "year", sub("\\.rds", "", grep("noPM", pixelGroupMapFiles, value = TRUE, invert = TRUE))),
                                    sub(".*(rep)([0-9]+)/.*", "\\1\\2", grep("noPM", pixelGroupMapFiles, value = TRUE, invert = TRUE)), sep = "_")

names(vegTypeMapStk_noPM) <- paste(sub(".*year", "year", sub("\\.rds", "", grep("noPM", vegTypeMapFiles, value = TRUE))),
                                   sub(".*(rep)([0-9]+)/.*", "\\1\\2", grep("noPM", vegTypeMapFiles, value = TRUE)), sep = "_")
names(vegTypeMapStk_PM) <- paste(sub(".*year", "year", sub("\\.rds", "", grep("noPM", vegTypeMapFiles, value = TRUE, invert = TRUE))),
                                 sub(".*(rep)([0-9]+)/.*", "\\1\\2", grep("noPM", vegTypeMapFiles, value = TRUE, invert = TRUE)), sep = "_")

## BUILD TABLES OF RESULTS ------------------------
## pixelCohortData tables
files <- grep("noPM", cohortDataFiles, value = TRUE)
pixelCohortData_noPM <- lapply(files, FUN = function(ff, pixelGroupMapStk) {
  cohortData <- readRDS(ff)
  yr <- sub(".*year", "",  sub(".rds", "", ff))
  r <- sub(".*(rep)([0-9]+)/.*", "\\2", ff)
  pixelGroupMap <- pixelGroupMapStk[[paste0("year", yr, "_rep", r)]]
  cohortData <- addPixels2CohortData(cohortData, pixelGroupMap)
  cohortData[, year := as.integer(yr)]
  cohortData[, rep := as.integer(r)]
  return(cohortData)
}, pixelGroupMapStk = pixelGroupMapStk_noPM) %>%
  rbindlist(fill = TRUE, l = .)

files <- grep("noPM", cohortDataFiles, value = TRUE, invert = TRUE)
pixelCohortData_PM <- lapply(files, FUN = function(ff, pixelGroupMapStk) {
  cohortData <- readRDS(ff)
  yr <- sub(".*year", "",  sub(".rds", "", ff))
  r <- sub(".*(rep)([0-9]+)/.*", "\\2", ff)
  pixelGroupMap <- pixelGroupMapStk[[paste0("year", yr, "_rep", r)]]
  cohortData <- addPixels2CohortData(cohortData, pixelGroupMap)
  cohortData[, year := as.integer(yr)]
  cohortData[, rep := as.integer(r)]
  return(cohortData)
}, pixelGroupMapStk = pixelGroupMapStk_PM) %>%
  rbindlist(fill = TRUE, l = .)

## vegTypeData tables
vegTypeSubset <- intersect(names(vegTypeMapStk_noPM), names(pixelGroupMapStk_noPM))
vegTypeData_noPM <- lapply(vegTypeSubset, FUN = function(x) {
  yr <- grep("year", unlist(strsplit(x, split = "_")), value = TRUE)
  r <- grep("rep", unlist(strsplit(x, split = "_")), value = TRUE)
  data.table(pixelIndex = seq_len(ncell(vegTypeMapStk_noPM[[x]])),
             pixelGroup = getValues(pixelGroupMapStk_noPM[[x]]),
             vegType = vegTypeMapStk_noPM[[x]][],
             year = as.integer(sub("year", "", yr)),
             rep = as.integer(sub("rep", "", r)))
}) %>%
  rbindlist(.)
vegTypeData_noPM <- vegTypeData_noPM[!is.na(pixelGroup)]

vegTypeSubset <- intersect(names(vegTypeMapStk_PM), names(pixelGroupMapStk_PM))
vegTypeData_PM <- lapply(vegTypeSubset, FUN = function(x) {
  yr <- grep("year", unlist(strsplit(x, split = "_")), value = TRUE)
  r <- grep("rep", unlist(strsplit(x, split = "_")), value = TRUE)
  data.table(pixelIndex = seq_len(ncell(vegTypeMapStk_PM[[x]])),
             pixelGroup = getValues(pixelGroupMapStk_PM[[x]]),
             vegType = vegTypeMapStk_PM[[x]][],
             year = as.integer(sub("year", "", yr)),
             rep = as.integer(sub("rep", "", r)))
})  %>%
  rbindlist(.)
vegTypeData_PM <- vegTypeData_PM[!is.na(pixelGroup)]

## pixelBurnData tables - all rasters
pixelBurnData_noPM <- lapply(unstack(rstCurrentBurnStk_noPM), FUN = function(ras) {
  yr <- grep("year", unlist(strsplit(names(ras), split = "_")), value = TRUE)
  r <- grep("rep", unlist(strsplit(names(ras), split = "_")), value = TRUE)
  data.table(pixelIndex = seq_len(ncell(ras)),
             burnt = as.integer(!is.na(getValues(ras))),
             year = as.integer(sub("year", "", yr)),
             rep = as.integer(sub("rep", "", r)))
}) %>%
  rbindlist(.)

pixelBurnData_PM <- lapply(unstack(rstCurrentBurnStk_PM), FUN = function(ras) {
  yr <- grep("year", unlist(strsplit(names(ras), split = "_")), value = TRUE)
  r <- grep("rep", unlist(strsplit(names(ras), split = "_")), value = TRUE)
  data.table(pixelIndex = seq_len(ncell(ras)),
             burnt = as.integer(!is.na(getValues(ras))),
             year = as.integer(sub("year", "", yr)),
             rep = as.integer(sub("rep", "", r)))
})  %>%
  rbindlist(.)

rm(files,vegTypeSubset)

## join tables, add scenario col and rbind.
## note that vegTypeData and pixelBurnData have more pixels because PGs of 0 are there, but not on cohortData
## so join by keeping all pixels
pixelCohortData_noPM <- merge(vegTypeData_noPM, pixelCohortData_noPM,
                              by = c("pixelIndex", "pixelGroup", "year", "rep"), all = TRUE)
pixelCohortData_noPM <- merge(pixelCohortData_noPM, pixelBurnData_noPM,
                              by = c("pixelIndex", "year", "rep"), all = TRUE)
pixelCohortData_noPM[is.na(burnt), burnt := 0]
pixelCohortData_noPM[, scenario := "noPM"]

pixelCohortData_PM <- merge(vegTypeData_PM, pixelCohortData_PM,
                            by = c("pixelIndex", "pixelGroup", "year", "rep"), all = TRUE)
pixelCohortData_PM <- merge(pixelCohortData_PM, pixelBurnData_PM,
                            by = c("pixelIndex", "year", "rep"), all = TRUE)
pixelCohortData_PM[is.na(burnt), burnt := 0]
pixelCohortData_PM[, scenario := "PM"]

allPixelCohortData <- rbind(pixelCohortData_noPM, pixelCohortData_PM, use.names = TRUE)
rm(list = grep("^[pixel|vegType].*Data", ls(), value = TRUE))
amc::.gc()

## ECOLOGICAL ZONATION -----------------------------
preSimList <- qread(file.path(simPaths$outputPath, "preSim", "preSimList.qs"))
ecoregionLayerRas <- rasterize(preSimList$ecoregionLayer, preSimList$rasterToMatch, field = "ecozoneCode")
ecoregionLayerDT <- data.table(ecozoneCode = getValues(ecoregionLayerRas),
                               pixelIndex = seq_len(ncell(ecoregionLayerRas)))
ecoregionLayerDT <- ecoregionLayerDT[!is.na(ecozoneCode)]

ecoregionLayerLabels <- data.table(ecozoneCode = preSimList$ecoregionLayer$ecozoneCode,
                                   ecozoneName = paste(preSimList$ecoregionLayer$NRNAME,
                                                       preSimList$ecoregionLayer$NSRNAME, sep = " - ")) %>%
  unique(.)

ecoregionLayerDT <- ecoregionLayerLabels[ecoregionLayerDT, on = .(ecozoneCode)]
allPixelCohortData <- ecoregionLayerDT[allPixelCohortData, on = .(pixelIndex)]
amc::.gc()

## NO. FIRES ---------------------------------------
## add no. fires per pixel, then remove lines with NA pixelGroup (not recorded)
## how many times did each pixel burn?
noFiresPixel <- unique(allPixelCohortData[, .(pixelIndex, scenario, year, rep, burnt)])
noFiresPixel[, noFires := sum(burnt), by = .(pixelIndex, scenario, rep)]
noFiresPixel <- unique(noFiresPixel[, .(pixelIndex, scenario, rep, noFires)])

test <- sapply(split(noFiresPixel, by = c("scenario", "rep")), FUN = function(x){
  any(duplicated(x[, pixelIndex]))
})
if (any(test))
  stop("Each pixel should only have one record of no. fires per scenario")

allPixelCohortData <- noFiresPixel[allPixelCohortData, on = .(pixelIndex, scenario, rep)]
allPixelCohortData[, burnt := NULL] ## no longer necessary
allPixelCohortData <- allPixelCohortData[!is.na(pixelGroup),]
rm(noFiresPixel, test)
amc::.gc()

## ADD MISSING SPECIES IN YEAR/SCENARIO/PIXEL COMBINATION
## cohortData doens't track absent cohorts, so they need to ba added back
## for now pixels from the 0s pixelGroup  have one entry with NAs for speciesCode
## they will be ignored for now and removed later, after adding one species entry for each of these pixels.
## for reporting consistency add to show losses in B
combinations <- unique(allPixelCohortData[, .(scenario, rep, year, pixelIndex, pixelGroup)])
spp <- as.character(na.omit(unique(allPixelCohortData$speciesCode)))
combinations <- lapply(spp, FUN = function(x) {
  data.table(combinations,
             speciesCode = x)
}) %>%
  rbindlist(., use.names = TRUE)

## join while keeping all combos, NA species will now disappear.
allPixelCohortData <- allPixelCohortData[combinations,
                                         on = .(scenario, rep, year, pixelIndex,
                                                pixelGroup, speciesCode)]
rm(spp, combinations)
amc::.gc()

## checks
test <- length(unique(allPixelCohortData[, length(unique(pixelIndex)), by = .(scenario, rep, year)]$V1)) == 1
test2 <- length(unique(allPixelCohortData[, length(unique(speciesCode)), by = .(scenario, rep, year, pixelIndex)]$V1)) == 1
test3 <- any(is.na(allPixelCohortData$speciesCode))

if (!test)
  stop("No. pixels should be the same across years, for a given scenario/rep")
if (!test2)
  stop("No. species per pixel should be the same across pixels, for a given scenario/rep/year")
if (test3)
  stop("There are NA speciesCodes")

## add ecoregion group where it's missing
## add vegType where it's missing, but it's a pixel with some veg
allPixelCohortData[, `:=`(ecoregionGroup = unique(na.omit(ecoregionGroup)),
                          vegType = max(vegType, na.rm = TRUE)),
                   by = .(scenario, rep, year, pixelGroup)]
amc::.gc()

## add noFires where it's missing
allPixelCohortData[, noFires := max(noFires, na.rm = TRUE),
                   by = .(scenario, rep, pixelIndex)]
amc::.gc()

## replace NAs of cohortData by 0s
replaceNAs <- function(x, val = 0) {
  x[is.na(x)] <- val
  x
}

cols <- c("age", "B", "mortality", "aNPPAct", "vegType", "noFires")
allPixelCohortData[, (cols) := lapply(.SD, replaceNAs), .SDcols = cols]
amc::.gc()

## add presence/absence of fire across simulation per pixel/scenario
allPixelCohortData[, firePresAbs := as.integer(any(noFires > 0)), by = .(scenario, rep, pixelIndex)]
amc::.gc()

## USING CAMERON'S CLASSIFICATION/SUMMARY ---------------------
## Cameron uses relative basal area to classify stand structure, we can use relative Biomass.
## we subset to the montane ecological zone, from where Cameron's data comes from
allPixelCohortDataMnt <- allPixelCohortData[grep("Montane", ecozoneName)]
allPixelCohortDataMnt[, sumB := sum(B), by = .(scenario, year, rep, pixelGroup)]
allPixelCohortDataMnt[, relB := sum(B)/sumB, by = .(scenario, rep, year, pixelGroup, speciesCode)]
allPixelCohortDataMnt[is.na(relB) & sumB == 0, relB := 0]

if (any(is.na(allPixelCohortDataMnt$relB)))
  stop("Missing values in relative biomass")

## subset to a smaller DT and join Cameron's species names
vegTypesCN <- unique(allPixelCohortDataMnt[, .(scenario, rep, year, pixelGroup, speciesCode, relB)])
vegTypesCN <- unique(na.omit(preSimList$sppEquiv[, .(Cameron, LIM)]))[vegTypesCN, on = "LIM==speciesCode",
                                                                      allow.cartesian = TRUE]
setnames(vegTypesCN, "LIM", "speciesCode")

parallelFUN <- function(DT) {
  set.seed(123)
  tempArg <- sample(1:nrow(DT), 100, replace = FALSE)
  tempArg <- DT[tempArg,]
  setkey(DT, scenario, rep, year, pixelGroup)
  out <- Cache(convertToCNVegType,
               DT = DT,
               groupingCols = c("scenario", "rep", "year", "pixelGroup"),
               cachingArg = tempArg,
               omitArgs = c("DT"),
               cacheRepo = cPath,
               userTags = c("reportResults"))
  out
}

amc::.gc()
plan("multiprocess", workers = 10)
vegTypesCN <- future_lapply(split(vegTypesCN, by = c("scenario", "rep")),
                            FUN = parallelFUN)
future:::ClusterRegistry("stop")
amc::.gc()

vegTypesCN <- rbindlist(vegTypesCN, use.names = TRUE)

## add Cameron's veg types and get rid of useless columns
cols <- c("scenario", "rep", "year", "pixelGroup", "vegTypeCN")
allPixelCohortDataMnt <- unique(vegTypesCN[, ..cols])[allPixelCohortDataMnt,
                                                      on = c("scenario", "rep", "year", "pixelGroup")]
allPixelCohortDataMnt[, `:=`(sumB = NULL,
                             relB = NULL,
                             vegType = NULL)]
amc::.gc()

## SUMMARY ACROSS LANDSCAPE -----------------------------------
## BY SPECIES
## remember that biomass is multiplied by 100 in *boreal*, this will revert the units to tonnes/ha
summaryBurnCohortDataSpp <- allPixelCohortData[, list(BiomassBySpecies = as.numeric(sum((B/100), na.rm = TRUE)),
                                                      MortalityBySpecies = as.numeric(sum((mortality/100), na.rm = TRUE)),
                                                      aNPPBySpecies = as.numeric(sum((aNPPAct/100), na.rm = TRUE)),
                                                      AgeBySppWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                                      sum((B/100), na.rm = TRUE))),
                                               by = .(scenario, rep, year, noFires, speciesCode)]
summaryBurnCohortDataSpp[is.nan(AgeBySppWeighted), AgeBySppWeighted := 0]
summaryBurnCohortDataSpp[, firePresAbs := as.integer(noFires != 0)]
amc::.gc()

## BY FOREST TYPE (def. as dominant species)
## remember that biomass is multiplied by 100 in *boreal*, this will revert the units to tonnes/ha
summaryBurnCohortDataVegType <- allPixelCohortData[, list(BiomassByVegType = as.numeric(sum((B/100), na.rm = TRUE)),
                                                          MortalityByVegType = as.numeric(sum((mortality/100), na.rm = TRUE)),
                                                          aNPPByVegType = as.numeric(sum((aNPPAct/100), na.rm = TRUE)),
                                                          AgeBySppWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                                          sum((B/100), na.rm = TRUE))),
                                                   by = .(scenario, rep, year, noFires, vegType)]
summaryBurnCohortDataVegType[is.nan(AgeBySppWeighted), AgeBySppWeighted := 0]
summaryBurnCohortDataVegType[, firePresAbs := as.integer(noFires != 0)]
amc::.gc()

## BY CAMERON'S VEG TYPE
summaryBurnCohortDataVegTypeCN <- allPixelCohortDataMnt[, list(BiomassByVegType = as.numeric(sum((B/100), na.rm = TRUE)),
                                                               MortalityByVegType = as.numeric(sum((mortality/100), na.rm = TRUE)),
                                                               aNPPByVegType = as.numeric(sum((aNPPAct/100), na.rm = TRUE)),
                                                               AgeBySppWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                                               sum((B/100), na.rm = TRUE))),
                                                        by = .(scenario, rep, year, noFires, vegTypeCN)]
summaryBurnCohortDataVegTypeCN[is.nan(AgeBySppWeighted), AgeBySppWeighted := 0]
summaryBurnCohortDataVegTypeCN[, firePresAbs := as.integer(noFires != 0)]
amc::.gc()

## make species labels/colours
speciesLabels <- LandR::equivalentName(value = unique(summaryBurnCohortDataSpp$speciesCode), column = "EN_generic_full",
                                       df = preSimList$sppEquiv)
names(speciesLabels) <- unique(summaryBurnCohortDataSpp$speciesCode)

speciesColours <- levels(vegTypeMapStk_noPM[[1]])[[1]]$colors
names(speciesColours) <- levels(vegTypeMapStk_noPM[[1]])[[1]]$VALUE

## make vegType labels/colours
vegTypeLabels <- as.character(levels(vegTypeMapStk_noPM[[1]])[[1]]$VALUE)
vegTypeLabels <- LandR::equivalentName(value = vegTypeLabels, column = "EN_generic_full",
                                       df = preSimList$sppEquiv)
vegTypeLabels[length(vegTypeLabels) + 1] <- "No veg."
names(vegTypeLabels) <- c(levels(vegTypeMapStk_noPM[[1]])[[1]]$ID, "0")

vegTypeColours <- as.character(levels(vegTypeMapStk_noPM[[1]])[[1]]$colors)
vegTypeColours[length(vegTypeColours) + 1] <- "grey40"
names(vegTypeColours) <- c(levels(vegTypeMapStk_noPM[[1]])[[1]]$ID, "0")

vegTypeCNLabels <- unique(summaryBurnCohortDataVegTypeCN$vegTypeCN)
names(vegTypeCNLabels) <- vegTypeCNLabels
vegTypeCNLabels <- sub("PIEN", "Spruce", vegTypeCNLabels)
vegTypeCNLabels <- sub("MMC", "Moist conif.", vegTypeCNLabels)
vegTypeCNLabels <- sub("mixedwood", "Mixed", vegTypeCNLabels)
vegTypeCNLabels <- sub("broadleaf", "Deciduous", vegTypeCNLabels)
vegTypeCNLabels <- sub("PICO", "Pine", vegTypeCNLabels)
vegTypeCNLabels <- sub("DMCPSME", "Dry conif./Douglas-fir", vegTypeCNLabels)
vegTypeCNLabels <- sub("^PSME$", "Douglas-fir", vegTypeCNLabels)
vegTypeCNLabels <- sub("dryPSME", "Dry Douglas-fir", vegTypeCNLabels)

vegTypeCNColours <- vegTypeCNLabels
vegTypeCNColours <- RColorBrewer::brewer.pal(length(vegTypeCNColours), name = "Accent")

## remove rasters to save memory
rm(list = grep("Stk", ls(), value = TRUE))
amc::.gc()


## STAT TESTS ------------------------------------------
## AGE SIM VS OBS ----
## difference between simulated and observed ages - end of sim
plotData <- allPixelCohortDataMnt[year == max(year) & vegTypeCN %in% unique(ageDataCN$Cover.dendro),
                                  list(AgeBySppWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                       sum((B/100), na.rm = TRUE))),
                                  by = .(scenario, rep, pixelIndex, firePresAbs, vegTypeCN)]

plotData2 <- ageDataCN[, list(avgAgeBySppWeightedObs = mean(Reconstructed.age)),
                       by = .(Cover.dendro)]
plotData <- plotData2[plotData, on = "Cover.dendro==vegTypeCN"]
setnames(plotData, "Cover.dendro", "vegTypeCN")

plotData[, ageDiffSimObs := AgeBySppWeighted - avgAgeBySppWeightedObs]
plotData[, firePresAbs := as.factor(firePresAbs)]
plotData[, vegTypeCN := as.factor(vegTypeCN)]
plotData[, vegTypeCN := relevel(vegTypeCN, "PICO")]

ageDiffLMList <- lapply(split(plotData, by = "firePresAbs"), function(DT) {
  lme4::lmer(abs(ageDiffSimObs) ~ scenario + (scenario | vegTypeCN), data = DT)
})
summary(ageDiffLMList$`0`)

plotData[firePresAbs == 0, predVals := predict(ageDiffLMList$`0`, type = "response")]
plotData[firePresAbs == 1, predVals := predict(ageDiffLMList$`1`, type = "response")]

plotTest1 <- ggplot(plotData) +
  geom_boxplot(aes(x = firePresAbs, y = abs(ageDiffSimObs), fill = scenario)) +
  geom_point(aes(x = firePresAbs, y = predVals, colour = scenario),
             position = position_dodge(width = 0.75), size = 2, show.legend = FALSE) +
  geom_hline(yintercept = 0, colour = "grey", linetype = "dashed") +
  scale_color_manual(values = c("red", "red")) +
  scale_x_discrete(labels = c("0" = "no fire", "1" = "fire")) +
  theme_pubr(base_size = 16) +
  theme(legend.title = element_blank()) +
  labs(title = "Age differences between sim. and obs. data",
       y = expression(bgroup("|", "sim." - bar("obs."), "|")),
       x = "") +
  facet_wrap(~ vegTypeCN, labeller = labeller(vegTypeCN = vegTypeCNLabels))

## SCENARIO EFFECTS ON STAND AGE------
## scenario effects on stand age across landscape
plotData <- allPixelCohortData[year == max(year),
                               list(AgeBySppWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                    sum((B/100), na.rm = TRUE))),
                               by = .(scenario, rep, pixelIndex, firePresAbs, ecoregionGroup)]

plotData[, firePresAbs := as.factor(firePresAbs)]
plotData[, ecoregionGroup := as.factor(ecoregionGroup)]
plotData <- plotData[!is.na(ecoregionGroup) & !is.na(AgeBySppWeighted)]

## some ecoregionGroups only have one observation, leading to signular fit
ageEffectsLMList <- lapply(split(plotData, by = "firePresAbs"), function(DT) {
  lme4::lmer(AgeBySppWeighted ~ scenario + (scenario | ecoregionGroup), data = DT)
})
summary(ageEffectsLMList$`1`)

## fitted values vector is shorter than nrow, prob. due to signularity. predict using the same data.
plotData[firePresAbs == 0, predVals := predict(ageEffectsLMList$`0`, type = "response",
                                               re.form = NA,
                                               newdata = plotData[firePresAbs == 0])]
plotData[firePresAbs == 1, predVals := predict(ageEffectsLMList$`1`, type = "response",
                                               re.form = NA,
                                               newdata = plotData[firePresAbs == 1])]

plotTest2 <- ggplot(plotData) +
  geom_boxplot(aes(x = firePresAbs, y = AgeBySppWeighted, fill = scenario)) +
  geom_point(aes(x = firePresAbs, y = predVals, colour = scenario),
             position = position_dodge(width = 0.75), size = 2, show.legend = FALSE) +
  geom_hline(yintercept = 0, colour = "grey", linetype = "dashed") +
  scale_color_manual(values = c("red", "red")) +
  scale_x_discrete(labels = c("0" = "no fire", "1" = "fire")) +
  theme_pubr(base_size = 16) +
  theme(legend.title = element_blank()) +
  labs(y = "biomass-weighted age (years)", x = "")

## ALPHA DIVERSITY ------------------
## scenario effects on alpha-div across landscape
## inverse Simpson concentration (Whittaker 1972)
plotData <- allPixelCohortData[year == max(year)]
plotData[, sumB := sum(B, na.rm = TRUE),
         by = .(scenario, rep, pixelIndex, firePresAbs, ecoregionGroup)]
plotData <- plotData[, list(relB = sum(B, na.rm = TRUE)/sumB),
                     by = .(scenario, rep, pixelIndex, speciesCode, firePresAbs, ecoregionGroup)]
plotData <- plotData[, list(alphaDiv = 1/sum(relB^2)),
                     by = .(scenario, rep, pixelIndex, firePresAbs, ecoregionGroup)]

plotData[, firePresAbs := as.factor(firePresAbs)]
plotData[, ecoregionGroup := as.factor(ecoregionGroup)]
plotData <- plotData[!is.na(ecoregionGroup) & !is.na(alphaDiv)]

## some ecoregionGroups only have one observation, leading to signular fit
alphaEffectsLMList <- lapply(split(plotData, by = "firePresAbs"), function(DT) {
  lme4::lmer(alphaDiv ~ scenario + (scenario | ecoregionGroup), data = DT)
})
summary(alphaEffectsLMList$`1`)

## fitted values vector is shorter than nrow, prob. due to signularity. predict using the same data.
plotData[firePresAbs == 0, predVals := predict(alphaEffectsLMList$`0`, type = "response",
                                               re.form = NA,
                                               newdata = plotData[firePresAbs == 0])]
plotData[firePresAbs == 1, predVals := predict(alphaEffectsLMList$`1`, type = "response",
                                               re.form = NA,
                                               newdata = plotData[firePresAbs == 1])]

plotTest3 <- ggplot(plotData) +
  geom_boxplot(aes(x = firePresAbs, y = alphaDiv, fill = scenario)) +
  geom_point(aes(x = firePresAbs, y = predVals, colour = scenario),
             position = position_dodge(width = 0.75), size = 2, show.legend = FALSE) +
  geom_hline(yintercept = 0, colour = "grey", linetype = "dashed") +
  scale_color_manual(values = c("red", "red")) +
  scale_x_discrete(labels = c("0" = "no fire", "1" = "fire")) +
  theme_pubr(base_size = 16) +
  theme(legend.title = element_blank()) +
  labs(y = "alpha-div.", x = "")

## BETA DIVERSITY ------------------
## scenario effects on beta-div across landscape
## from multiplicative decomposition where both alpha
## and gamma diversity are inverse Simpson concentration (Whittaker 1972)
# alpha-div
plotData <- allPixelCohortData[year == max(year)]
plotData[, sumB := sum(B, na.rm = TRUE),
         by = .(scenario, rep, pixelIndex, firePresAbs, ecoregionGroup)]
plotData <- plotData[, list(relB = sum(B, na.rm = TRUE)/sumB),
                     by = .(scenario, rep, pixelIndex, speciesCode, firePresAbs, ecoregionGroup)]
plotData <- plotData[, list(alphaDiv = 1/sum(relB^2)),
                     by = .(scenario, rep, pixelIndex, firePresAbs, ecoregionGroup)]

## gamma-div
plotData2 <- allPixelCohortData[year == max(year),
                                list(BiomassBySpecies = sum(B, na.rm = TRUE)),
                                by = .(scenario, rep, firePresAbs, speciesCode, ecoregionGroup)]
plotData2[, sumB := sum(BiomassBySpecies, na.rm = TRUE),
          .(scenario, rep, firePresAbs, ecoregionGroup)]
plotData2 <- plotData2[, list(relB = sum(BiomassBySpecies, na.rm = TRUE)/sumB),
                       by = .(scenario, rep, firePresAbs, speciesCode, ecoregionGroup)]
plotData2 <- plotData2[, list(gammaDiv = 1/sum(relB^2)),
                       by = .(scenario, rep, firePresAbs, ecoregionGroup)]

## beta-div
plotData2 <- plotData2[plotData, on = .(scenario, rep, firePresAbs, ecoregionGroup)]
plotData2 <- plotData2[, list(betaDiv = gammaDiv * sum((1/.N) * (1/alphaDiv), na.rm = TRUE)),
                       by = .(scenario, rep, firePresAbs, ecoregionGroup)]

plotData2[, firePresAbs := as.factor(firePresAbs)]
plotData2[, ecoregionGroup := as.factor(ecoregionGroup)]
plotData2 <- unique(plotData2[!is.na(ecoregionGroup) & !is.na(betaDiv)])

## some ecoregionGroups only have one observation, leading to signular fit
betaEffectsLMList <- lapply(split(plotData2, by = "firePresAbs"), function(DT) {
  lme4::lmer(betaDiv ~ scenario + (1|ecoregionGroup), data = DT)
})
summary(betaEffectsLMList$`1`)

## fitted values vector is shorter than nrow, prob. due to signularity. predict using the same data.
plotData2[firePresAbs == 0, predVals := predict(betaEffectsLMList$`0`, type = "response",
                                                re.form = NA,
                                                newdata = plotData2[firePresAbs == 0])]
plotData2[firePresAbs == 1, predVals := predict(betaEffectsLMList$`1`, type = "response",
                                                re.form = NA,
                                                newdata = plotData2[firePresAbs == 1])]

plotTest4 <- ggplot(plotData2) +
  geom_boxplot(aes(x = firePresAbs, y = betaDiv, fill = scenario)) +
  geom_point(aes(x = firePresAbs, y = predVals, colour = scenario),
             position = position_dodge(width = 0.75), size = 2, show.legend = FALSE) +
  geom_hline(yintercept = 0, colour = "grey", linetype = "dashed") +
  scale_color_manual(values = c("red", "red")) +
  scale_x_discrete(labels = c("0" = "no fire", "1" = "fire")) +
  theme_pubr(base_size = 16) +
  theme(legend.title = element_blank()) +
  labs(y = "beta-div.", x = "")

saveRDS(ageDiffLMList, file = file.path(simPaths$outputPath, "ageDiffModel.rds"))
saveRDS(ageEffectsLMList, file = file.path(simPaths$outputPath, "ageEffectsModel.rds"))
saveRDS(alphaEffectsLMList, file = file.path(simPaths$outputPath, "alphaEffectsModel.rds"))
saveRDS(betaEffectsLMList, file = file.path(simPaths$outputPath, "betaEffectsModel.rds"))

## PLOTS -----------------------------------------------
## MODELLED PROPERTIES ---------------------------------
## BY SPECIES ------------
## total biomass
plotData <- summaryBurnCohortDataSpp[, list(BiomassBySpecies = sum(BiomassBySpecies)),
                                     by = .(scenario, year, rep, speciesCode, firePresAbs)]
plot1 <- ggplot(data = plotData,
                aes(x = year, y = log(BiomassBySpecies + 0.000001), colour = speciesCode)) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours, labels = speciesLabels) +
  labs(title = "Total landscape biomass", y = "log-biomass (ton/ha)",
       subtitle = "by species") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot2 <- ggplot(data = summaryBurnCohortDataSpp,
                aes(x = year, y = log(BiomassBySpecies + 0.000001), colour = speciesCode)) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours, labels = speciesLabels) +
  labs(title = "Total landscape biomass", y = "log-biomass (ton/ha)",
       subtitle = "by species") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ noFires)

## total mortality
plotData <- summaryBurnCohortDataSpp[, list(MortalityBySpecies = sum(MortalityBySpecies)),
                                     by = .(scenario, year, rep, speciesCode, firePresAbs)]
plot3 <- ggplot(data = plotData,
                aes(x = year, y = log(MortalityBySpecies + 0.000001), colour = speciesCode)) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours, labels = speciesLabels) +
  labs(title = "Total landscape mortality", y = "log-biomass (ton/ha)",
       subtitle = "by species") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot4 <- ggplot(data = summaryBurnCohortDataSpp,
                aes(x = year, y = log(MortalityBySpecies + 0.000001), colour = speciesCode)) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours, labels = speciesLabels) +
  labs(title = "Total landscape mortality", y = "log-biomass (ton/ha)",
       subtitle = "by species") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ noFires)
amc::.gc()

## average age
plotData <- allPixelCohortData[, list(AgeBySppWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                      sum((B/100), na.rm = TRUE))),
                               by = .(scenario, year, rep, firePresAbs, speciesCode)]
plot5 <- ggplot(data = plotData,
                aes(x = year, y = AgeBySppWeighted, colour = speciesCode)) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours, labels = speciesLabels) +
  labs(title = "Avg. age across landscape", y = "biomass-weighted age",
       subtitle = "by species") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot6 <- ggplot(data = summaryBurnCohortDataSpp,
                aes(x = year, y = AgeBySppWeighted, colour = speciesCode)) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours, labels = speciesLabels) +
  labs(title = "Avg. age across landscape", y = "biomass-weighted age",
       subtitle = "by species") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ noFires)

## no. cohorts
plotData <- allPixelCohortData[B > 0, list(noCohorts = length(unique(age))),
                               by = .(scenario, year, rep, pixelIndex, firePresAbs, noFires, speciesCode)]
plot7 <- ggplot(data = plotData,
                aes(x = year, y = noCohorts, colour = speciesCode)) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours, labels = speciesLabels) +
  labs(title = "Avg. no. cohorts", y = "no. cohorts",
       subtitle = "by species") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot8 <- ggplot(data = plotData,
                aes(x = year, y = noCohorts, colour = speciesCode)) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours, labels = speciesLabels) +
  labs(title = "Avg. no. cohorts by species", y = "no. cohorts",
       subtitle = "by species") +
  facet_grid(scenario ~ noFires)
amc::.gc()


## BY FOREST TYPE (DOMINANT SPECIES) ------------
## no. pixels
plotData <- allPixelCohortData[, list(noPixelsVeg = length(unique(pixelIndex))),
                               by = .(scenario, year, rep, firePresAbs, vegType)]
plot9 <- ggplot(data = plotData,
                aes(x = year, y = noPixelsVeg, colour = as.factor(vegType))) +
  stat_summary(fun = "mean", geom = "area", position = "stack") +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Forest type (dom. species)", y = "no. pixels") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot9.2 <- ggplot(data = plotData,
                  aes(x = year, y = noPixelsVeg, colour = as.factor(vegType))) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Forest type (dom. species)", y = "prop. pixels") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))
amc::.gc()

plotData <- allPixelCohortData[, list(noPixelsVeg = length(unique(pixelIndex))),
                               by = .(scenario, year, rep, noFires, vegType)]
plot10 <- ggplot(data = plotData,
                 aes(x = year, y = noPixelsVeg, colour = as.factor(vegType))) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Forest type (dom. species)", y = "no. pixels") +
  facet_grid(scenario ~ noFires)

plot10.2 <- ggplot(data = plotData,
                   aes(x = year, y = noPixelsVeg, colour = as.factor(vegType))) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Forest type (dom. species)", y = "prop. pixels") +
  facet_grid(scenario ~ noFires)

## no. cohorts
plotData <- allPixelCohortData[B > 0, list(noCohorts = length(unique(age))),
                               by = .(scenario, year, rep, pixelIndex, firePresAbs, noFires, vegType)]
plot11 <- ggplot(data = plotData,
                 aes(x = year, y = noCohorts, colour = as.factor(vegType))) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Avg. no. cohorts", y = "no. cohorts",
       subtitle = "by forest type (dom. species)") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot12 <- ggplot(data = plotData,
                 aes(x = year, y = noCohorts, colour = as.factor(vegType))) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Avg. no. cohorts", y = "no. cohorts",
       subtitle = "by forest type (dom. species)") +
  facet_grid(scenario ~ noFires)

## total biomass
plotData <- summaryBurnCohortDataVegType[, list(BiomassByVegType = sum(BiomassByVegType)),
                                         by = .(scenario, year, rep, vegType, firePresAbs)]
plot13 <- ggplot(data = plotData,
                 aes(x = year, y = log(BiomassByVegType + 0.000001), colour = as.factor(vegType))) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Total landscape biomass", y = "log-biomass (ton/ha)",
       subtitle = "by forest type (dom. species)") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot14 <- ggplot(data = summaryBurnCohortDataVegType,
                 aes(x = year, y = log(BiomassByVegType + 0.000001), colour = as.factor(vegType))) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Total landscape biomass", y = "log-biomass (ton/ha)",
       subtitle = "by forest type (dom. species)") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ noFires)

## total mortality
plotData <- summaryBurnCohortDataVegType[, list(MortalityByVegType = sum(MortalityByVegType)),
                                         by = .(scenario, year, rep, vegType, firePresAbs)]
plot15 <- ggplot(data = plotData,
                 aes(x = year, y = log(MortalityByVegType + 0.000001), colour = as.factor(vegType))) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Total landscape mortality", y = "log-biomass (ton/ha)",
       subtitle = "by forest type (dom. species)") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot16 <- ggplot(data = summaryBurnCohortDataVegType,
                 aes(x = year, y = log(MortalityByVegType + 0.000001), colour = as.factor(vegType))) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Total landscape mortality", y = "log-biomass (ton/ha)",
       subtitle = "by forest type (dom. species)") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ noFires)

## average age
plotData <- allPixelCohortData[, list(AgeBySppWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                      sum((B/100), na.rm = TRUE))),
                               by = .(scenario, year, rep, firePresAbs, vegType)]

plot17 <- ggplot(data = plotData,
                 aes(x = year, y = AgeBySppWeighted, colour = as.factor(vegType))) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Avg. forest type age (dom. species)", y = "biomass-weighted age",
       subtitle = "presence/absence of fire") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot18 <- ggplot(data = summaryBurnCohortDataVegType,
                 aes(x = year, y = AgeBySppWeighted, colour = as.factor(vegType))) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Avg. forest type age (dom. species)", y = "biomass-weighted age",
       subtitle = "no. fires") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ noFires)
amc::.gc()


## BY CAMERON'S VEG TYPE ------------
## no. pixels
plotData <- allPixelCohortDataMnt[, list(noPixelsVeg = length(unique(pixelIndex))),
                                  by = .(scenario, year, rep, firePresAbs, vegTypeCN)]
plot19 <- ggplot(data = plotData,
                 aes(x = year, y = noPixelsVeg, fill = as.factor(vegTypeCN))) +
  stat_summary(fun = "mean", geom = "area", position = "stack") +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  labs(title = "Forest type", y = "no. pixels",
       subtitle = "presence/absence of fire") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot19.2 <- ggplot(data = plotData,
                   aes(x = year, y = noPixelsVeg, fill = as.factor(vegTypeCN))) +
  stat_summary(fun = "mean", geom = "area", position = "fill") +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  labs(title = "Forest type", y = "prop. pixels",
       subtitle = "presence/absence of fire") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))
amc::.gc()

plotData <- allPixelCohortDataMnt[, list(noPixelsVeg = length(unique(pixelIndex))),
                                  by = .(scenario, year, rep, noFires, vegTypeCN)]
plot20 <- ggplot(data = plotData,
                 aes(x = year, y = noPixelsVeg, fill = as.factor(vegTypeCN))) +
  stat_summary(fun = "mean", geom = "area", position = "stack") +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  labs(title = "Forest type", y = "no. pixels",
       subtitle = "no. fires") +
  facet_grid(scenario ~ noFires)

plot20.2 <- ggplot(data = plotData,
                   aes(x = year, y = noPixelsVeg, fill = as.factor(vegTypeCN))) +
  stat_summary(fun = "mean", geom = "area", position = "fill") +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  labs(title = "Forest type", y = "prop. pixels",
       subtitle = "no. fires") +
  facet_grid(scenario ~ noFires)

## no. cohorts
plotData <- allPixelCohortDataMnt[B > 0, list(noCohorts = length(unique(age))),
                                  by = .(scenario, year, rep, pixelIndex, firePresAbs, noFires, vegTypeCN)]
plot21 <- ggplot(plotData,
                 aes(x = year, y = noCohorts, colour = as.factor(vegTypeCN))) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  labs(title = "Avg. no. cohorts", y = "no. cohorts",
       subtitle = "by forest type") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot22 <- ggplot(data = plotData,
                 aes(x = year, y = noCohorts, colour = as.factor(vegTypeCN))) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  labs(title = "Avg. no. cohorts", y = "no. cohorts",
       subtitle = "by forest type") +
  facet_grid(scenario ~ noFires)

## total biomass
plotData <- summaryBurnCohortDataVegTypeCN[, list(BiomassByVegType = sum(BiomassByVegType)),
                                           by = .(scenario, year, rep, vegTypeCN, firePresAbs)]
plot23 <- ggplot(data = plotData,
                 aes(x = year, y = log(BiomassByVegType + 0.000001), colour = as.factor(vegTypeCN))) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  labs(title = "Total landscape biomass", y = "log-biomass (ton/ha)",
       subtitle = "by forest type") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot24 <- ggplot(data = summaryBurnCohortDataVegTypeCN,
                 aes(x = year, y = log(BiomassByVegType + 0.000001), colour = as.factor(vegTypeCN))) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  labs(title = "Total landscape biomass", y = "log-biomass (ton/ha)",
       subtitle = "by forest type") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ noFires)

## total mortality
plotData <- summaryBurnCohortDataVegTypeCN[, list(MortalityByVegType = sum(MortalityByVegType)),
                                           by = .(scenario, year, rep, vegTypeCN, firePresAbs)]
plot25 <- ggplot(data = plotData,
                 aes(x = year, y = log(MortalityByVegType + 0.000001), colour = as.factor(vegTypeCN))) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  labs(title = "Total landscape mortality", y = "log-biomass (ton/ha)",
       subtitle = "by forest type") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot26 <- ggplot(data = summaryBurnCohortDataVegTypeCN,
                 aes(x = year, y = log(MortalityByVegType + 0.000001), colour = as.factor(vegTypeCN))) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  labs(title = "Total landscape mortality", y = "log-biomass (ton/ha)",
       subtitle = "by forest type") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ noFires)

## average age
plotData <- allPixelCohortDataMnt[, list(AgeBySppWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                         sum((B/100), na.rm = TRUE))),
                                  by = .(scenario, year, rep, firePresAbs, vegTypeCN)]
plot27 <- ggplot(data = plotData,
                 aes(x = year, y = AgeBySppWeighted, colour = as.factor(vegTypeCN))) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  labs(title = "Avg. forest type age", y = "biomass-weighted age",
       subtitle = "presence/absence of fire") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot28 <- ggplot(data = summaryBurnCohortDataVegTypeCN,
                 aes(x = year, y = AgeBySppWeighted, colour = as.factor(vegTypeCN))) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  labs(title = "Avg. forest type age", y = "biomass-weighted age",
       subtitle = "no. fires") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ noFires)

## SIMULATED VS. OBSERVED DATA ------------
## age distribution at end of simulation vs cameron's - subset to common veg types
plotData <- allPixelCohortDataMnt[year == max(year) & vegTypeCN %in% unique(ageDataCN$Cover.dendro),
                                  list(AgeBySppWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                       sum((B/100), na.rm = TRUE))),
                                  by = .(scenario, rep, pixelIndex, firePresAbs, vegTypeCN)]
plot29 <- ggplot(plotData,
                 aes(x = vegTypeCN, y = AgeBySppWeighted, fill = vegTypeCN)) +
  geom_violin(position = position_nudge(x = -0.15),
              alpha = 0.7, show.legend = FALSE) +
  stat_summary(fun = "mean", geom = "point", colour = "black",
               position = position_nudge(x = -0.15), show.legend = FALSE) +
  stat_summary(data = ageDataCN,
               aes(x = Cover.dendro, y = Reconstructed.age),
               fun.data = "mean_sd", fill = "black", colour = "red",
               position = position_nudge(x = 0.15), show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_x_discrete(labels = vegTypeCNLabels) +
  scale_fill_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  scale_colour_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  labs(title = "Forest type age", y = "biomass-weighted age", x = "",
       subtitle = "presence/absence of fire") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot29.2 <- ggplot(plotData,
                   aes(x = vegTypeCN, y = AgeBySppWeighted, fill = vegTypeCN)) +
  geom_violin(position = position_nudge(x = -0.15),
              alpha = 0.7, show.legend = FALSE) +
  stat_summary(fun = "mean", geom = "point", colour = "black",
               position = position_nudge(x = -0.15), show.legend = FALSE) +
  stat_summary(data = ageDataCN,
               aes(x = Cover.dendro, y = Reconstructed.age),
               fun.data = "mean_sd", fill = "black", colour = "red",
               position = position_nudge(x = 0.15), show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_x_discrete(labels = vegTypeCNLabels) +
  scale_fill_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  scale_colour_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  labs(title = "Forest type age", y = "biomass-weighted age", x = "") +
  facet_grid(~ scenario)


plotData <- allPixelCohortDataMnt[year == max(year) & vegTypeCN %in% unique(ageDataCN$Cover.dendro),
                                  list(avgAge = mean(age, na.rm = TRUE)),
                                  by = .(scenario, rep, pixelIndex, firePresAbs, vegTypeCN)]
plot30 <- ggplot(plotData,
                 aes(x = vegTypeCN, y = avgAge, fill = vegTypeCN)) +
  geom_violin(position = position_nudge(x = -0.15),
              alpha = 0.7, show.legend = FALSE) +
  stat_summary(fun = "mean", geom = "point", colour = "black",
               position = position_nudge(x = -0.15), show.legend = FALSE) +
  stat_summary(data = ageDataCN,
               aes(x = Cover.dendro, y = Reconstructed.age),
               fun.data = "mean_sd", fill = "black", colour = "red",
               position = position_nudge(x = 0.15), show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_x_discrete(labels = vegTypeCNLabels) +
  scale_fill_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  scale_colour_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  labs(title = "Forest type age", y = "average age", x = "",
       subtitle = "presence/absence of fire") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot30.2 <- ggplot(plotData,
                   aes(x = vegTypeCN, y = avgAge, fill = vegTypeCN)) +
  geom_violin(position = position_nudge(x = -0.15) , show.legend = FALSE,
              alpha = 0.7) +
  stat_summary(fun = "mean", geom = "point", colour = "black",
               position = position_nudge(x = -0.15), show.legend = FALSE) +
  stat_summary(data = ageDataCN,
               aes(x = Cover.dendro, y = Reconstructed.age),
               fun.data = "mean_sd", fill = "black", colour = "red",
               position = position_nudge(x = 0.15), show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_x_discrete(labels = vegTypeCNLabels) +
  scale_fill_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  scale_colour_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  labs(title = "Forest type age", y = "average age", x = "") +
  facet_grid(~ scenario)

## no. cohorts at end of simulation vs cameron's - subset to common veg types
plotData <- allPixelCohortDataMnt[year == max(year) & B > 0 & vegTypeCN %in% unique(ageDataCN$Cover.dendro),
                                  list(noCohorts = length(unique(age))),
                                  by = .(scenario, rep, pixelIndex, firePresAbs, vegTypeCN)]
plot31 <- ggplot(data = plotData,
                 aes(x = vegTypeCN, y = noCohorts, fill = vegTypeCN)) +
  geom_boxplot(alpha = 0.7, show.legend = FALSE) +
  stat_summary(data = ageDataCN,
               aes(x = Cover.dendro, y = noCohorts),
               fun.data = "mean_sd", fill = "black", colour = "red",
               size = 0.7, show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_x_discrete(labels = vegTypeCNLabels) +
  scale_fill_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  scale_colour_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  labs(title = "Avg. no. cohorts by forest type", y = "no. cohorts", x = "",
       subtitle = "presence/absence of fire") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

## mean age differences
plotData <- allPixelCohortDataMnt[year == max(year) & vegTypeCN %in% unique(ageDataCN$Cover.dendro),
                                  list(AgeBySppWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                       sum((B/100), na.rm = TRUE))),
                                  by = .(scenario, rep, pixelIndex, firePresAbs, vegTypeCN)]

plotData2 <- ageDataCN[, list(avgAgeBySppWeightedObs = mean(Reconstructed.age)),
                       by = .(Cover.dendro)]
plotData <- plotData2[plotData, on = "Cover.dendro==vegTypeCN"]
setnames(plotData, "Cover.dendro", "vegTypeCN")

plotData2 <- plotData[, list(meanAbsDevSimObs = 1/.N * sum(AgeBySppWeighted - mean(avgAgeBySppWeightedObs))),
                      by = .(scenario, firePresAbs)]
plotData <- plotData2[plotData, on = .(scenario, firePresAbs)]
rm(plotData2)

plot32 <- ggplot(plotData,
                 aes(x = vegTypeCN,
                     y = AgeBySppWeighted - avgAgeBySppWeightedObs,
                     fill = vegTypeCN)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") +
  geom_hline(mapping = aes(yintercept = meanAbsDevSimObs), size = 1, colour = "red") +
  geom_boxplot(show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_x_discrete(labels = vegTypeCNLabels) +
  scale_fill_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  labs(title = "Age differences between sim. and obs. data",
       y = "mean age diff.", x = "",
       subtitle = "presence/absence of fire") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

## BIODIVERSITY METRICS --------------------------------
## ALPHA DIVERSITY --------------
## inverse Simpson concentration (Whittaker 1972)
plotData <- copy(allPixelCohortData)
plotData[, sumB := sum(B, na.rm = TRUE),
         by = .(scenario, rep, year, pixelIndex, firePresAbs, noFires)]
plotData <- plotData[, list(relB = sum(B, na.rm = TRUE)/sumB),
                     by = .(scenario, rep, year, pixelIndex, speciesCode, firePresAbs, noFires)]
plotData <- plotData[, list(alphaDiv = 1/sum(relB^2)),
                     by = .(scenario, rep, year, pixelIndex, firePresAbs, noFires)]

plot33 <- ggplot(data = plotData,
                 aes(x = year, y = alphaDiv, colour = scenario)) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  labs(y = "mean alpha-div.") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot34 <- ggplot(data = plotData,
                 aes(x = year, y = alphaDiv, colour = scenario)) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  labs(y = "mean alpha-div.") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(~ noFires)

## map of alpha - avg across reps
plotData2 <- plotData[year == 100, list(alphaDiv = mean(alphaDiv)),
                      by = .(scenario, pixelIndex)]
alphaStk <- lapply(split(plotData2, by = "scenario"), function(DT, RTM) {
  rasAlpha <- RTM
  rasAlpha[] <- NA
  rasAlpha[DT$pixelIndex] <- DT$alphaDiv
  rasAlpha
}, RTM = preSimList$rasterToMatch)
alphaStk <- stack(alphaStk)

rasData <- data.table(vals = getValues(alphaStk), coordinates(alphaStk))
rasData <- melt(rasData, id.vars = c("x", "y"))
rasData$variable <- sub("vals.", "", rasData$variable)

plotMap1 <- ggplot() +
  layer_spatial(data = preSimList$studyArea, col = "black",
                fill = "transparent") +
  geom_raster(data = rasData,
              mapping = aes(x, y, fill = value)) +
  annotation_north_arrow(style = north_arrow_minimal,
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         location = "tr", which_north = "true") +
  theme_pubr(base_size = 16) +
  scale_fill_distiller(palette = "YlGnBu", direction = 1,
                       breaks = c(0:length(unique(allPixelCohortData$speciesCode))),
                       labels = c(0:length(unique(allPixelCohortData$speciesCode))),
                       limits = c(0, length(unique(allPixelCohortData$speciesCode))),
                       na.value = "transparent") +
  labs(x = "longitude", y = "latitude", fill = "alpha-div.") +
  facet_wrap(~ variable)

plotMap1hist <- ggplot(rasData[!is.na(value)]) +
  geom_density(aes(x = value, fill = variable, alpha = variable)) +
  theme_pubr(base_size = 16) +
  theme(legend.title = element_blank()) +
  scale_alpha_manual(values = c("PM" = 0.6, "noPM" = 1)) +
  labs(y = "density", x = "mean alpha-div (years)")


## BETA DIVERSITY --------------
## from multiplicative decomposition where both alpha
## and gamma diversity are inverse Simpson concentration (Whittaker 1972)
## alpha-div
plotData <- copy(allPixelCohortData)
plotData[, sumB := sum(B, na.rm = TRUE),
         by = .(scenario, rep, year, pixelIndex, firePresAbs, noFires)]
plotData <- plotData[, list(relB = sum(B, na.rm = TRUE)/sumB),
                     by = .(scenario, rep, year, pixelIndex, speciesCode, firePresAbs, noFires)]
plotData <- plotData[, list(alphaDiv = 1/sum(relB^2)),
                     by = .(scenario, rep, year, pixelIndex, firePresAbs, noFires)]

## gamma-div - firePresAbs
plotData2 <- copy(summaryBurnCohortDataSpp)
plotData2[, sumB := sum(BiomassBySpecies, na.rm = TRUE),
          .(scenario, rep, year, firePresAbs)]
plotData2 <- plotData2[, list(relB = sum(BiomassBySpecies, na.rm = TRUE)/sumB),
                       by = .(scenario, rep, year, firePresAbs, speciesCode)]
plotData2 <- plotData2[, list(gammaDiv = 1/sum(relB^2)),
                       by = .(scenario, rep, year, firePresAbs)]

## beta-div  - firePresAbs
plotData2 <- plotData2[plotData, on = .(scenario, rep, year, firePresAbs)]
plotData2 <- plotData2[, list(betaDiv = gammaDiv * sum((1/.N) * (1/alphaDiv), na.rm = TRUE)),
                       by = .(scenario, rep, year, firePresAbs)]

plot35 <- ggplot(data = plotData2,
                 aes(x = year, y = betaDiv, colour = scenario)) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  labs(y = "mean beta-div.") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))


## gamma-div - no fires
plotData2 <- copy(summaryBurnCohortDataSpp)
plotData2[, sumB := sum(BiomassBySpecies, na.rm = TRUE),
          .(scenario, rep, year, noFires)]
plotData2 <- plotData2[, list(relB = sum(BiomassBySpecies, na.rm = TRUE)/sumB),
                       by = .(scenario, rep, year, noFires, speciesCode)]
plotData2 <- plotData2[, list(gammaDiv = 1/sum(relB^2)),
                       by = .(scenario, rep, year, noFires)]
## beta-div
plotData2 <- plotData2[plotData, on = .(scenario, rep, year, noFires)]
plotData2 <- plotData2[, list(betaDiv = gammaDiv * sum((1/.N) * (1/alphaDiv), na.rm = TRUE)),
                       by = .(scenario, rep, year, noFires)]

plot36 <- ggplot(data = plotData2,
                 aes(x = year, y = betaDiv, colour = scenario)) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  labs(y = "mean beta-div.") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(~ noFires)


## VARIABILITY ACROSS REPS -----------------------------
## BY SPECIES ---------------------
## total biomass across landscape
plotData <- summaryBurnCohortDataSpp[year == 100, list(BiomassBySpecies = sum(BiomassBySpecies)),
                                     by = .(scenario, year, rep, speciesCode, firePresAbs)]
plot1var <- ggplot(data = plotData,
                   aes(x = speciesCode, y = log(BiomassBySpecies + 0.000001),
                       fill = speciesCode)) +
  geom_boxplot(show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = speciesColours, labels = speciesLabels) +
  scale_x_discrete(labels = speciesLabels) +
  labs(title = "Total biomass (year 100)", y = "log-biomass (ton/ha)",
       x = "", subtitle = "by species") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

## average stand biomass
plotData <- allPixelCohortData[year == 100, list(B = sum(B/100, na.rm = TRUE)),   ## copllapse cohorts first
                               by = .(scenario, rep, firePresAbs, pixelIndex, speciesCode)]
plotData <- plotData[, list(avgBspp = mean(B, na.rm = TRUE)),    ## now average across landscape
                     by = .(scenario, rep, firePresAbs, speciesCode)]
plot2var <- ggplot(data = plotData,
                   aes(x = speciesCode, y = avgBspp, fill = speciesCode)) +
  geom_boxplot(show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = speciesColours, labels = speciesLabels) +
  scale_x_discrete(labels = speciesLabels) +
  labs(title = "Avg. biomass in a stand (year 100)", y = "biomass (ton/ha)",
       x = "" , subtitle = "by species") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

## average age across landscape
plotData <- allPixelCohortData[year == 100, list(AgeBySppWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                                 sum((B/100), na.rm = TRUE))),
                               by = .(scenario, year, rep, firePresAbs, speciesCode)]
plot3var <- ggplot(data = plotData,
                   aes(x = speciesCode, y = AgeBySppWeighted, fill = speciesCode)) +
  geom_boxplot(show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = speciesColours, labels = speciesLabels) +
  scale_x_discrete(labels = speciesLabels) +
  labs(title = "Avg. biomass-weighted age (year 100)", y = "years",
       x = "" , subtitle = "by species") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

## average stand age
plotData <- allPixelCohortData[year == 100,
                               list(AgeBySppWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                    sum((B/100), na.rm = TRUE))),   ## calculate per stand first
                               by = .(scenario, rep, firePresAbs, pixelIndex, speciesCode)]
plotData <- plotData[, list(AgeBySppWeighted = mean(AgeBySppWeighted, na.rm = TRUE)),    ## now average across landscape
                     by = .(scenario, rep, firePresAbs, speciesCode)]
plot4var <- ggplot(data = plotData,
                   aes(x = speciesCode, y = AgeBySppWeighted, fill = speciesCode)) +
  geom_boxplot(show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = speciesColours, labels = speciesLabels) +
  scale_x_discrete(labels = speciesLabels) +
  labs(title = "Avg. biomass-weighted age in a stand (year 100)", y = "years",
       x = "" , subtitle = "by species") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

## avg no. cohorts in stand
plotData <- allPixelCohortData[year == 100 & B > 0,
                               list(noCohorts = length(unique(age))),   ## calculate per stand first
                               by = .(scenario, rep, firePresAbs, pixelIndex, speciesCode)]
plotData <- plotData[, list(avgNoCohorts = mean(noCohorts, na.rm = TRUE)),    ## now average across landscape
                     by = .(scenario, rep, firePresAbs, speciesCode)]
plot5var <- ggplot(data = plotData,
                   aes(x = speciesCode, y = avgNoCohorts, fill = speciesCode)) +
  geom_boxplot(show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = speciesColours, labels = speciesLabels) +
  scale_x_discrete(labels = speciesLabels) +
  labs(title = "Avg. no. cohorts in a stand (year 100)", y = "no. cohorts",
       x = "" , subtitle = "by species") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

## (DOMINANT SPECIES) ------------
## total biomass across landscape
plotData <- summaryBurnCohortDataVegType[year == 100, list(BiomassBySpecies = sum(BiomassByVegType)),
                                         by = .(scenario, rep, vegType, firePresAbs)]
plot6var <- ggplot(data = plotData,
                   aes(x = as.factor(vegType), y = log(BiomassBySpecies + 0.000001),
                       fill = as.factor(vegType))) +
  geom_boxplot(show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeColours, labels = vegTypeLabels) +
  scale_x_discrete(labels = vegTypeLabels) +
  labs(title = "Total biomass (year 100)", y = "log-biomass (ton/ha)",
       x = "" , subtitle = "by forest type (dom. species)") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

## average stand biomass
plotData <- allPixelCohortData[year == 100, list(B = sum(B/100, na.rm = TRUE)),   ## copllapse cohorts first
                               by = .(scenario, rep, firePresAbs, pixelIndex, vegType)]
plotData <- plotData[, list(avgBvegType = mean(B, na.rm = TRUE)),    ## now average across landscape
                     by = .(scenario, rep, firePresAbs, vegType)]
plot7var <- ggplot(data = plotData,
                   aes(x = as.factor(vegType), y = avgBvegType, fill = as.factor(vegType))) +
  geom_boxplot(show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeColours, labels = vegTypeLabels) +
  scale_x_discrete(labels = vegTypeLabels) +
  labs(title = "Avg. biomass in a stand (year 100)", y = "biomass (ton/ha)",
       x = "" , subtitle = "by forest type (dom. species)") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

## average age across landscape
plotData <- allPixelCohortData[year == 100, list(AgeBySppWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                                 sum((B/100), na.rm = TRUE))),
                               by = .(scenario, year, rep, firePresAbs, vegType)]
plot8var <- ggplot(data = plotData,
                   aes(x = as.factor(vegType), y = AgeBySppWeighted, fill = as.factor(vegType))) +
  geom_boxplot(show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeColours, labels = vegTypeLabels) +
  scale_x_discrete(labels = vegTypeLabels) +
  labs(title = "Avg. biomass-weighted age (year 100)", y = "years",
       x = "" , subtitle = "by forest type (dom. species)") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

## average stand age
plotData <- allPixelCohortData[year == 100,
                               list(AgeBySppWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                    sum((B/100), na.rm = TRUE))),   ## calculate per stand first
                               by = .(scenario, rep, firePresAbs, pixelIndex, vegType)]
plotData <- plotData[, list(AgeBySppWeighted = mean(AgeBySppWeighted, na.rm = TRUE)),    ## now average across landscape
                     by = .(scenario, rep, firePresAbs, vegType)]
plot9var <- ggplot(data = plotData,
                   aes(x = as.factor(vegType), y = AgeBySppWeighted, fill = as.factor(vegType))) +
  geom_boxplot(show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeColours, labels = vegTypeLabels) +
  scale_x_discrete(labels = vegTypeLabels) +
  labs(title = "Avg. biomass-weighted age in a stand (year 100)", y = "years",
       x = "" , subtitle = "by forest type (dom. species)") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

## avg no. cohorts in stand
plotData <- allPixelCohortData[year == 100 & B > 0,
                               list(noCohorts = length(unique(age))),   ## calculate per stand first
                               by = .(scenario, rep, firePresAbs, pixelIndex, vegType)]
plotData <- plotData[, list(avgNoCohorts = mean(noCohorts, na.rm = TRUE)),    ## now average across landscape
                     by = .(scenario, rep, firePresAbs, vegType)]
plotData[, vegType := factor(vegType, levels = sort(names(vegTypeLabels)),
                             labels = vegTypeLabels[order(names(vegTypeLabels))])]

plot10var <- ggplot(data = plotData,
                    aes(x = as.factor(vegType), y = avgNoCohorts, fill = as.factor(vegType))) +
  geom_boxplot(show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeColours, labels = vegTypeLabels) +
  scale_x_discrete(labels = vegTypeLabels, drop = FALSE) +
  labs(title = "Avg. no. cohorts in a stand (year 100)", y = "no. cohorts",
       x = "" , subtitle = "by forest type (dom. species)") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

## no. pixels
plotData <- allPixelCohortData[year == 100, list(noPixelsVeg = length(unique(pixelIndex))),
                               by = .(scenario, year, rep, firePresAbs, vegType)]
plot11var <- ggplot(data = plotData,
                    aes(x = as.factor(vegType), y = noPixelsVeg, fill = as.factor(vegType))) +
  geom_boxplot(show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeColours, labels = vegTypeLabels) +
  scale_x_discrete(labels = vegTypeLabels) +
  labs(title = "Forest type (year 100)", y = "no. pixels",
       x = "" , subtitle = "by forest type (dom. species)") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

## BY CAMERON'S VEG TYPE ------------------------
## total biomass across landscape
plotData <- summaryBurnCohortDataVegTypeCN[year == 100, list(BiomassBySpecies = sum(BiomassByVegType)),
                                           by = .(scenario, rep, vegTypeCN, firePresAbs)]
plot12var <- ggplot(data = plotData,
                    aes(x = as.factor(vegTypeCN), y = log(BiomassBySpecies + 0.000001),
                        fill = as.factor(vegTypeCN))) +
  geom_boxplot(show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  scale_x_discrete(labels = vegTypeCNLabels) +
  labs(title = "Total biomass (year 100)", y = "log-biomass (ton/ha)",
       x = "" , subtitle = "by forest type") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

## average stand biomass
plotData <- allPixelCohortDataMnt[year == 100, list(B = sum(B/100, na.rm = TRUE)),   ## copllapse cohorts first
                                  by = .(scenario, rep, firePresAbs, pixelIndex, vegTypeCN)]
plotData <- plotData[, list(avgBvegType = mean(B, na.rm = TRUE)),    ## now average across landscape
                     by = .(scenario, rep, firePresAbs, vegTypeCN)]
plot13var <- ggplot(data = plotData,
                    aes(x = vegTypeCN, y = avgBvegType, fill = vegTypeCN)) +
  geom_boxplot(show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  scale_x_discrete(labels = vegTypeCNLabels) +
  labs(title = "Avg. biomass in a stand (year 100)", y = "biomass (ton/ha)",
       x = "" , subtitle = "by forest type") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

## average age across landscape
plotData <- allPixelCohortDataMnt[year == 100, list(AgeBySppWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                                    sum((B/100), na.rm = TRUE))),
                                  by = .(scenario, year, rep, firePresAbs, vegTypeCN)]
plot14var <- ggplot(data = plotData,
                    aes(x = vegTypeCN, y = AgeBySppWeighted, fill = vegTypeCN)) +
  geom_boxplot(show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  scale_x_discrete(labels = vegTypeCNLabels) +
  labs(title = "Avg. biomass-weighted age (year 100)", y = "years",
       x = "" , subtitle = "by forest type") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

## average stand age
plotData <- allPixelCohortDataMnt[year == 100,
                                  list(AgeBySppWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                       sum((B/100), na.rm = TRUE))),   ## calculate per stand first
                                  by = .(scenario, rep, firePresAbs, pixelIndex, vegTypeCN)]
plotData <- plotData[, list(AgeBySppWeighted = mean(AgeBySppWeighted, na.rm = TRUE)),    ## now average across landscape
                     by = .(scenario, rep, firePresAbs, vegTypeCN)]
plot15var <- ggplot(data = plotData,
                    aes(x = vegTypeCN, y = AgeBySppWeighted, fill = vegTypeCN)) +
  geom_boxplot(show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  scale_x_discrete(labels = vegTypeCNLabels) +
  labs(title = "Avg. biomass-weighted age in a stand (year 100)", y = "years",
       x = "" , subtitle = "by forest type") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

## avg no. cohorts in stand
plotData <- allPixelCohortDataMnt[year == 100 & B > 0,
                                  list(noCohorts = length(unique(age))),   ## calculate per stand first
                                  by = .(scenario, rep, firePresAbs, pixelIndex, vegTypeCN)]
plotData <- plotData[, list(avgNoCohorts = mean(noCohorts, na.rm = TRUE)),    ## now average across landscape
                     by = .(scenario, rep, firePresAbs, vegTypeCN)]
plot16var <- ggplot(data = plotData,
                    aes(x = vegTypeCN, y = avgNoCohorts, fill = vegTypeCN)) +
  geom_boxplot(show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  scale_x_discrete(labels = vegTypeCNLabels) +
  labs(title = "Avg. no. cohorts in a stand (year 100)", y = "no. cohorts",
       x = "" , subtitle = "by forest type") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

## no. pixels
plotData <- allPixelCohortDataMnt[year == 100, list(noPixelsVeg = length(unique(pixelIndex))),
                                  by = .(scenario, year, rep, firePresAbs, vegTypeCN)]
plot17var <- ggplot(data = plotData,
                    aes(x = vegTypeCN, y = noPixelsVeg, fill = vegTypeCN)) +
  geom_boxplot(show.legend = FALSE) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeCNColours, labels = vegTypeCNLabels) +
  scale_x_discrete(labels = vegTypeCNLabels) +
  labs(title = "Forest type (year 100)", y = "no. pixels",
       x = "" , subtitle = "by forest type") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))


## MAPS ------------------------
## dominant species - to deal with reps we need to reclassify the stands after averaging across reps.
plotData <- allPixelCohortData[year == 100, list(B = mean(B, na.rm = TRUE),
                                                 age =  as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                     sum((B/100), na.rm = TRUE))),
                               by = .(scenario, pixelIndex, speciesCode, firePresAbs)]

plotData <- lapply(split(plotData[B > 0], by = "scenario"), vegTypeGenerator,
                   vegLeadingProportion = 0, sppEquiv = preSimList$sppEquiv,
                   sppEquivCol = P(preSimList)$Biomass_borealDataPrep$sppEquivCol,
                   pixelGroupColName = "pixelIndex", doAssertion = FALSE)  ## turn assertions off so that the algorythm choice is based on table size
plotData <- rbindlist(plotData, use.names = TRUE)
setnames(plotData, "leading", "vegType")

plotData[, vegType := as.numeric(factor(vegType))]

vegTypeStk <- lapply(split(plotData, by = "scenario"), function(DT, RTM) {
  rasVeg <- RTM
  rasVeg[] <- NA
  rasVeg[DT$pixelIndex] <- DT$vegType
  rasVeg
}, RTM = preSimList$rasterToMatch)
vegTypeStk <- stack(vegTypeStk)

ageStk <- lapply(split(plotData, by = "scenario"), function(DT, RTM) {
  rasVeg <- RTM
  rasVeg[] <- NA
  rasVeg[DT$pixelIndex] <- DT$age
  rasVeg
}, RTM = preSimList$rasterToMatch)
ageStk <- stack(ageStk)

rasData <- data.table(vals = getValues(vegTypeStk), coordinates(vegTypeStk))
rasData <- melt(rasData, id.vars = c("x", "y"))
rasData$variable <- sub("vals.", "", rasData$variable)

plotMap2 <- ggplot() +
  layer_spatial(data = preSimList$studyArea, col = "black",
                fill = "transparent") +
  geom_raster(data = rasData,
              mapping = aes(x, y, fill = as.factor(value))) +
  annotation_north_arrow(style = north_arrow_minimal,
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         location = "tr", which_north = "true") +
  theme_pubr(base_size = 16) +
  theme(plot.margin = unit(c(0,0,0,0), units = "mm"), legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeColours, labels = vegTypeLabels,
                    na.translate = FALSE) +
  labs(x = "longitude", y = "latitude") +
  facet_wrap(~ variable)

plotMap2hist <- ggplot(rasData[!is.na(value)]) +
  geom_bar(aes(x = as.factor(value), fill = as.factor(value),
               alpha = variable == "PM"), position = "dodge") +
  theme_pubr(base_size = 16, x.text.angle = 45) +
  theme(plot.margin = unit(c(0,0,0,0), units = "mm"), legend.title = element_blank(),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = vegTypeColours, labels = vegTypeLabels,
                    na.translate = FALSE) +
  scale_alpha_manual(labels = c("TRUE" = "PM", "FALSE" = "noPM"),
                     values = c("TRUE" = 1, "FALSE" = 0.6)) +
  scale_x_discrete(labels = vegTypeLabels) +
  guides(alpha = guide_legend(nrow = 2)) +
  labs(y = "no. pixels")

rasData <- data.table(vals = getValues(ageStk), coordinates(ageStk))
rasData <- melt(rasData, id.vars = c("x", "y"))
rasData$variable <- sub("vals.", "", rasData$variable)

plotMap3 <- ggplot() +
  layer_spatial(data = preSimList$studyArea, col = "black",
                fill = "transparent") +
  geom_raster(data = rasData,
              mapping = aes(x, y, fill = value)) +
  annotation_north_arrow(style = north_arrow_minimal,
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         location = "tr", which_north = "true") +
  theme_pubr(base_size = 16) +
  theme(plot.margin = unit(c(0,0,0,0), units = "mm")) +
  scale_fill_distiller(palette = "Greens", na.value = "transparent",
                       direction = 1) +
  labs(x = "longitude", y = "latitude", fill = "age") +
  facet_wrap(~ variable)

plotMap3hist <- ggplot(rasData[!is.na(value)]) +
  geom_density(aes(x = value, fill = variable, alpha = variable)) +
  theme_pubr(base_size = 16) +
  theme(plot.margin = unit(c(0,0,0,0), units = "mm"), legend.title = element_blank()) +
  scale_alpha_manual(values = c("PM" = 0.6, "noPM" = 1)) +
  labs(y = "density", x = "age (years)")

## SAVE PLOTS ------------------------------------------
amc::.gc()
plotSave <- ggarrange(plot1,
                      plot2 +
                        theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
                        labs(title = "", subtitle = ""),
                      widths = c(0.4, 0.7),
                      legend = "bottom", common.legend = TRUE)
ggsave(plot = plotSave, filename = file.path(figOutputPath, "results_landscapeB.tiff"),
       width = 14, height = 7)

plotSave <- ggarrange(plot13 + theme(plot.title = element_text(size = 17)),
                      plot14 +
                        theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
                        labs(title = "", subtitle = ""),
                      widths = c(0.4, 0.7),
                      legend = "bottom", common.legend = TRUE)
ggsave(plot = plotSave, filename = file.path(figOutputPath, "results_landscapeBVegType.tiff"),
       width = 14, height = 7)

plotSave <- ggarrange(plot23,
                      plot24 +
                        theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
                        labs(title = "", subtitle = ""),
                      widths = c(0.4, 0.7),
                      legend = "bottom", common.legend = TRUE)
ggsave(plot = plotSave, filename = file.path(figOutputPath, "results_landscapeBVegTypeCN.tiff"),
       width = 14, height = 7)

# plotSave <- ggarrange(plot3,
#                   plot4 +
#                     theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
#                     labs(title = "", subtitle = ""),
#                   widths = c(0.4, 0.6),
#                   legend = "bottom", common.legend = TRUE)
# ggsave(plot = plotSave, filename = file.path(figOutputPath, "results_landscapeMort.tiff"),
#        width = 14, height = 7)

plotSave <- ggarrange(plot5,
                      plot6 +
                        theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
                        labs(title = "", subtitle = ""),
                      widths = c(0.4, 0.7),
                      legend = "bottom", common.legend = TRUE)
ggsave(plot = plotSave, filename = file.path(figOutputPath, "results_landscapeAge.tiff"),
       width = 14, height = 7)

plotSave <- ggarrange(plot17,
                      plot18 +
                        theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
                        labs(title = "", subtitle = ""),
                      widths = c(0.4, 0.7),
                      legend = "bottom", common.legend = TRUE)
ggsave(plot = plotSave, filename = file.path(figOutputPath, "results_landscapeAgeVegType.tiff"),
       width = 14, height = 7)

plotSave <- ggarrange(plot27,
                      plot28 +
                        theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
                        labs(title = "", subtitle = ""),
                      widths = c(0.4, 0.7),
                      legend = "bottom", common.legend = TRUE)
ggsave(plot = plotSave, filename = file.path(figOutputPath, "results_landscapeAgeVegTypeCN.tiff"),
       width = 14, height = 7)

plotSave <- ggarrange(plot9.2,
                      plot10.2 +
                        theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
                        labs(title = ""),
                      widths = c(0.4, 0.7),
                      legend = "bottom", common.legend = TRUE)
ggsave(plot = plotSave, filename = file.path(figOutputPath, "results_landscapeVegTypes.tiff"),
       width = 14, height = 7)

plotSave <- ggarrange(plot19.2 + theme(plot.subtitle = element_blank()),
                      plot20.2 +
                        theme(plot.subtitle = element_blank(),
                              axis.title.y = element_blank(), axis.text.y = element_blank()) +
                        labs(title = ""),
                      widths = c(0.4, 0.7),
                      legend = "bottom", common.legend = TRUE)
ggsave(plot = plotSave, filename = file.path(figOutputPath, "results_landscapeVegTypesCN.tiff"),
       width = 14, height = 7)

plotSave <- ggarrange(plot7 + theme(plot.subtitle = element_blank()),
                      plot8 +
                        theme(plot.subtitle = element_blank(),
                              axis.title.y = element_blank(), axis.text.y = element_blank()) +
                        labs(title = "", subtitle = ""),
                      widths = c(0.4, 0.7),
                      legend = "bottom", common.legend = TRUE)
ggsave(plot = plotSave, filename = file.path(figOutputPath, "results_noCohorts.tiff"),
       width = 14, height = 7)

plotSave <- ggarrange(plot11,
                      plot12 +
                        theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
                        labs(title = "", subtitle = ""),
                      widths = c(0.4, 0.7),
                      legend = "bottom", common.legend = TRUE)
ggsave(plot = plotSave, filename = file.path(figOutputPath, "results_noCohortsVegType.tiff"),
       width = 14, height = 7)

plotSave <- ggarrange(plot21,
                      plot22 +
                        theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
                        labs(title = "", subtitle = ""),
                      widths = c(0.4, 0.7),
                      legend = "bottom", common.legend = TRUE)
ggsave(plot = plotSave, filename = file.path(figOutputPath, "results_noCohortsVegTypeCN.tiff"),
       width = 14, height = 7)

ggsave(plot = plot29, filename = file.path(figOutputPath, "results_ageSimVsObs.tiff"),
       width = 14, height = 7)

ggsave(plot = plot32, filename = file.path(figOutputPath, "results_meanAgeDiff.tiff"),
       width = 14, height = 7)

ggsave(plot = plot31, filename = file.path(figOutputPath, "results_noCohortsSimVsObs.tiff"),
       width = 14, height = 7)

ggsave(plot = plot31.2, filename = file.path(figOutputPath, "results_meanCohortDiff.tiff"),
       width = 14, height = 7)

plotSave <- ggarrange(plot1var +
                        theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                              plot.margin = margin(0,0,0,0)),
                      plot3var +
                        theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                              plot.title = element_blank(), plot.margin = margin(0,0,0,0)) +
                        labs(title = "", subtitle = ""),
                      plot5var + theme(plot.title = element_blank(),
                                       plot.margin = margin(0,0,0,5, unit = "mm")) +
                        labs(title = "", subtitle = ""), ncol = 1,
                      heights = c(1, 0.9, 1.3),
                      legend = "bottom", common.legend = TRUE)
ggsave(plot = plotSave, filename = file.path(figOutputPath, "results_landscapeBAgeCohortsVar.tiff"),
       width = 7, height = 12)

plotSave <- ggarrange(plot6var + theme(plot.subtitle = element_blank(),
                                       axis.title.x = element_blank(), axis.text.x = element_blank()),
                      plot8var +
                        theme(plot.margin = margin(l = 10, unit = "mm"), plot.subtitle = element_blank(),
                              axis.title.x = element_blank(), axis.text.x = element_blank()),
                      plot10var +
                        theme(plot.margin = margin(l = 9, unit = "mm"), plot.subtitle = element_blank()),
                      plot11var + theme(plot.subtitle = element_blank()),
                      ncol = 2, nrow = 2,
                      heights = c(1, 1.6),
                      legend = "bottom", common.legend = TRUE)
ggsave(plot = plotSave, filename = file.path(figOutputPath, "results_landscapeBAgeCohortsVarVegType.tiff"),
       width = 14, height = 9)

plotSave <- ggarrange(plot12var +
                        theme(plot.subtitle = element_blank(), plot.margin = margin(unit = "mm"),
                              axis.title.x = element_blank(), axis.text.x = element_blank()),
                      plot14var +
                        theme(plot.margin = margin(l = 10, unit = "mm"), plot.subtitle = element_blank(),
                              axis.title.x = element_blank(), axis.text.x = element_blank()),
                      plot16var +
                        theme(plot.margin = margin(t = 6, l = 6, unit = "mm"), plot.subtitle = element_blank()),
                      plot17var +
                        theme(plot.margin = margin(t = 6, unit = "mm"), plot.subtitle = element_blank()),
                      ncol = 2, nrow = 2,
                      heights = c(1, 1.6),
                      legend = "bottom", common.legend = TRUE)
ggsave(plot = plotSave, filename = file.path(figOutputPath, "results_landscapeBAgeCohortsVarVegTypeCN.tiff"),
       width = 14, height = 9)

plotSave <- ggarrange(plotMap1 + theme(legend.position = "right"),
                      plotMap1hist +
                        theme(legend.position = c(0.7, 0.95), legend.justification = c("left", "top")),
                      ncol = 1, heights = c(2,1), widths = c(0.98, 1))
ggsave(plot = plotSave, filename = file.path(figOutputPath, "results_avgAlphaDivMap.tiff"),
       width = 8, height = 12)

plotSave <- ggarrange(plotMap2 + theme(legend.position = "none"),
                      plotMap3 +
                        theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                              legend.position = "right"),
                      plotMap2hist +
                        theme(plot.margin = margin(b = 12, unit = "mm"),
                              axis.text.x = element_blank(),
                              legend.position = "none") +
                        guides(fill = guide_legend(nrow = 2), alpha = guide_legend(nrow = 2)),
                      plotMap3hist +
                        theme(legend.position = c(0.95, 0.95),
                              legend.justification = c("right", "top")),
                      get_legend(plotMap2hist),
                      ncol = 2, nrow = 3, heights = c(2,1,0.2), widths = c(0.98, 1))
ggsave(plot = plotSave, filename = file.path(figOutputPath, "results_avgVegTypeAgeMap.tiff"),
       width = 13, height = 12)

plotSave <- ggarrange(plotTest1 +
                        theme(legend.position = c(0.6, 0), legend.justification = c("left", "bottom"),
                              legend.key.size = unit(1.5, "cm")))
ggsave(plot = plotSave, filename = file.path(figOutputPath, "results_ageSimVsObsTestBP.tiff"),
       width = 12, height = 7)

plotSave <- ggarrange(plotTest2 + theme(legend.position = "none"),
                      plotTest3 + theme(legend.position = "none"),
                      plotTest4 + theme(legend.position = "none"),
                      get_legend(plotTest2 + theme(legend.key.size = unit(1.5, "cm"))),
                      ncol = 2, nrow = 2)
ggsave(plot = plotSave, filename = file.path(figOutputPath, "results_ageAlphaBetaDivEffectsTestBP.tiff"),
       width = 12, height = 7)

