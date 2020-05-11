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

source("R/R_tools/convertToCNVegType.R")

simPaths <- list(cachePath = file.path("R/SpaDES/cache/AI_report")
                 , modulePath = file.path("R/SpaDES/m")
                 , inputPath = file.path("R/SpaDES/inputs")
                 , outputPath = file.path("R/SpaDES/outputs"))

## path to figure folder and cache folder
figOutputPath <- "C:/Users/Ceres Barros/Google Drive/Shared/McIntire-lab/Manuscripts_inPrep/LIMmodel_paper"
cPath <- file.path(simPaths$cachePath, "postSimAnalyses")

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

test <- sapply(unique(noFiresPixel$scenario), FUN = function(x){
  any(duplicated(noFiresPixel[scenario == x, pixelIndex]))
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
allPixelCohortDataMnt[, sumB := sum(B), by = .(scenario, year, rep, pixelIndex)]
allPixelCohortDataMnt[, relB := sum(B)/sumB, by = .(scenario, rep, year, pixelIndex, speciesCode)]
allPixelCohortDataMnt[is.na(relB) & sumB == 0, relB := 0]

if (any(is.na(allPixelCohortDataMnt$relB)))
  stop("Missing values in relative biomass")

## subset to a smaller DT and join Cameron's species names
vegTypesCN <- unique(allPixelCohortDataMnt[, .(scenario, rep, year, pixelGroup, speciesCode, relB)])
vegTypesCN <- unique(na.omit(preSimList$sppEquiv[, .(Cameron, LIM)]))[vegTypesCN, on = "LIM==speciesCode",
                                                                      allow.cartesian = TRUE]
setnames(vegTypesCN, "LIM", "speciesCode")

set.seed(123)
tempArg <- sample(1:nrow(vegTypesCN), 100, replace = FALSE)
tempArg <- vegTypesCN[tempArg,]
vegTypesCN <- Cache(convertToCNVegType,
                    DT = vegTypesCN,
                    groupingCols = c("scenario", "rep", "year", "pixelGroup"),
                    cachingArg = tempArg,
                    omitArgs = c("DT"),
                    cacheRepo = cPath)
rm(tempArg)
amc::.gc()
## SUMMARY ACROSS LANDSCAPE -----------------------------------
## BY SPECIES
## remember that biomass is multiplied by 100 in *boreal*, this will revert the units to tonnes/ha
summaryBurnCohortDataSpp <- allPixelCohortData[, list(BiomassBySpecies = as.numeric(sum((B/100), na.rm = TRUE)),
                                                      MortalityBySpecies = as.numeric(sum((mortality/100), na.rm = TRUE)),
                                                      aNPPBySpecies = as.numeric(sum((aNPPAct/100), na.rm = TRUE)),
                                                      AgeBySppWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                                      sum((B/100), na.rm = TRUE))),
                                               by = .(scenario, year, noFires, speciesCode)]
summaryBurnCohortDataSpp[is.nan(AgeBySppWeighted), AgeBySppWeighted := 0]
summaryBurnCohortDataSpp[, firePresAbs := as.integer(noFires != 0)]
amc::.gc()

## BY FOREST TYPE (def. as dominant species)
## remember that biomass is multiplied by 100 in *boreal*, this will revert the units to tonnes/ha
summaryBurnCohortDataVegType <- allPixelCohortData[, list(BiomassBySpecies = as.numeric(sum((B/100), na.rm = TRUE)),
                                                          MortalityBySpecies = as.numeric(sum((mortality/100), na.rm = TRUE)),
                                                          aNPPBySpecies = as.numeric(sum((aNPPAct/100), na.rm = TRUE)),
                                                          AgeBySppWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                                          sum((B/100), na.rm = TRUE))),
                                                   by = .(scenario, year, noFires, vegType)]
summaryBurnCohortDataVegType[is.nan(AgeBySppWeighted), AgeBySppWeighted := 0]
summaryBurnCohortDataVegType[, firePresAbs := as.integer(noFires != 0)]
amc::.gc()

## make species labels/colours
speciesLabels <- LandR::equivalentName(value = unique(summaryBurnCohortDataSpp$speciesCode), column = "EN_generic_full",
                                       df = simList_noPM$sppEquiv)
names(speciesLabels) <- unique(summaryBurnCohortDataSpp$speciesCode)

speciesColours <- levels(vegTypeMapStk_noPM[[1]])[[1]]$colors
names(speciesColours) <- levels(vegTypeMapStk_noPM[[1]])[[1]]$VALUE

## make vegType labels/colours
vegTypeLabels <- as.character(levels(vegTypeMapStk_noPM[[1]])[[1]]$VALUE)
vegTypeLabels <- LandR::equivalentName(value = vegTypeLabels, column = "EN_generic_full",
                                       df = simList_noPM$sppEquiv)
vegTypeLabels[length(vegTypeLabels) + 1] <- "No veg."
names(vegTypeLabels) <- c(levels(vegTypeMapStk_noPM[[1]])[[1]]$ID, "0")

vegTypeColours <- as.character(levels(vegTypeMapStk_noPM[[1]])[[1]]$colors)
vegTypeColours[length(vegTypeColours) + 1] <- "grey40"
names(vegTypeColours) <- c(levels(vegTypeMapStk_noPM[[1]])[[1]]$ID, "0")

## remove rasters to save memory
rm(list = grep("Stk", ls(), value = TRUE))

## PLOTS -----------------------------------------------

## BY SPECIES ------------
fireYears <- as.numeric(sub("year", "", names(rstCurrentBurnStk_PM)))

plotData <- summaryBurnCohortDataSpp[, list(BiomassBySpecies = sum(BiomassBySpecies)),
                                     by = .(scenario, year, speciesCode, firePresAbs)]
plot1 <- ggplot(data = plotData,
                aes(x = year, y = log(BiomassBySpecies + 0.000001), colour = speciesCode)) +
  # geom_vline(xintercept = fireYears, size = 1, linetype = "dashed", colour = "grey") +
  geom_line(size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours, labels = speciesLabels) +
  labs(title = "Total landscape biomass by species", y = expression(paste("log-Biomass"~~"(g/m"^2, ")")),
       subtitle = "presence/absence of fire") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot2 <- ggplot(data = summaryBurnCohortDataSpp,
                aes(x = year, y = log(BiomassBySpecies + 0.000001), colour = speciesCode)) +
  # geom_vline(xintercept = fireYears, size = 1, linetype = "dashed", colour = "grey") +
  geom_line(size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours, labels = speciesLabels) +
  labs(title = "Total landscape biomass by species", y = expression(paste("log-Biomass"~~"(g/m"^2, ")")),
       subtitle = "no. fires") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ noFires)

plotData <- summaryBurnCohortDataSpp[, list(MortalityBySpecies = sum(MortalityBySpecies)),
                                     by = .(scenario, year, speciesCode, firePresAbs)]
plot3 <- ggplot(data = plotData,
                aes(x = year, y = log(MortalityBySpecies + 0.000001), colour = speciesCode)) +
  geom_line(size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours, labels = speciesLabels) +
  labs(title = "Total landscape mortality by species", y = expression(paste("log-Biomass"~~"(g/m"^2, ")")),
       subtitle = "presence/absence of fire") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot4 <- ggplot(data = summaryBurnCohortDataSpp,
                aes(x = year, y = log(MortalityBySpecies + 0.000001), colour = speciesCode)) +
  # geom_vline(xintercept = fireYears, size = 1, linetype = "dashed", colour = "grey") +
  geom_line(size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours, labels = speciesLabels) +
  labs(title = "Total landscape mortality by species", y = expression(paste("log-Biomass"~~"(g/m"^2, ")")),
       subtitle = "no. fires") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ noFires)

amc::.gc()
plotData <- allPixelCohortData[, list(AgeBySppWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                      sum((B/100), na.rm = TRUE))),
                               by = .(scenario, year, firePresAbs, speciesCode)]

plot5 <- ggplot(data = plotData,
                aes(x = year, y = AgeBySppWeighted, colour = speciesCode)) +
  geom_line(size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours, labels = speciesLabels) +
  labs(title = "Avg. species age across landscape", y = "years",
       subtitle = "presence/absence of fire") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot6 <- ggplot(data = summaryBurnCohortDataSpp,
                aes(x = year, y = AgeBySppWeighted, colour = speciesCode)) +
  # geom_vline(xintercept = fireYears, size = 1, linetype = "dashed", colour = "grey") +
  geom_line(size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours, labels = speciesLabels) +
  labs(title = "Avg. species age across landscape", y = "years",
       subtitle = "no. fires") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ noFires)


plotData <- allPixelCohortData[B > 0, list(noCohorts = length(unique(age))),
                               by = .(scenario, year, pixelIndex, firePresAbs, noFires, speciesCode)]
plot7 <- ggplot(data = plotData[, list(noCohorts = mean(noCohorts)),
                                by = .(scenario, year, firePresAbs, speciesCode)],
                aes(x = year, y = noCohorts, colour = speciesCode)) +
  # geom_vline(xintercept = fireYears, size = 1, linetype = "dashed", colour = "grey") +
  geom_line(size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours, labels = speciesLabels) +
  labs(title = "Avg. no. cohorts by species", y = "no. cohorts",
       subtitle = "presence/absence of fire") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot8 <- ggplot(data = plotData[, list(noCohorts = mean(noCohorts)),
                                by = .(scenario, year, noFires, speciesCode)],
                aes(x = year, y = noCohorts, colour = speciesCode)) +
  # geom_vline(xintercept = fireYears, size = 1, linetype = "dashed", colour = "grey") +
  geom_line(size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours, labels = speciesLabels) +
  labs(title = "Avg. no. cohorts by species", y = "no. cohorts",
       subtitle = "no. fires") +
  facet_grid(scenario ~ noFires)


## BY FOREST TYPE (DOMINANT SPECIES) ------------
amc::.gc()
plotData <- allPixelCohortData[, list(noPixelsVeg = length(unique(pixelIndex))),
                               by = .(scenario, year, firePresAbs, vegType)]
plot9 <- ggplot(data = plotData,
                aes(x = year, y = noPixelsVeg, fill = as.factor(vegType))) +
  geom_area(stat = "identity", position = "stack") +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Forest type (dominant species)", y = "no. pixels",
       subtitle = "presence/absence of fire") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot9.2 <- ggplot(data = plotData,
                  aes(x = year, y = noPixelsVeg, fill = as.factor(vegType))) +
  geom_area(stat = "identity", position = "fill") +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Forest type (dominant species)", y = "prop. pixels",
       subtitle = "presence/absence of fire") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

amc::.gc()
plotData <- allPixelCohortData[, list(noPixelsVeg = length(unique(pixelIndex))),
                               by = .(scenario, year, noFires, vegType)]
plot10 <- ggplot(data = plotData,
                 aes(x = year, y = noPixelsVeg, fill = as.factor(vegType))) +
  # geom_vline(xintercept = fireYears, size = 1, linetype = "dashed", colour = "grey") +
  geom_area(stat = "identity", position = "stack") +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Forest type (dominant species)", y = "no. pixels",
       subtitle = "no. fires") +
  facet_grid(scenario ~ noFires)

plot10.2 <- ggplot(data = plotData,
                   aes(x = year, y = noPixelsVeg, fill = as.factor(vegType))) +
  # geom_vline(xintercept = fireYears, size = 1, linetype = "dashed", colour = "grey") +
  geom_area(stat = "identity", position = "fill") +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Forest type (dominant species)", y = "prop. pixels",
       subtitle = "no. fires") +
  facet_grid(scenario ~ noFires)


plotData <- allPixelCohortData[B > 0, list(noCohorts = length(unique(age))),
                               by = .(scenario, year, pixelIndex, firePresAbs, noFires, vegType)]
plot11 <- ggplot(data = plotData[, list(noCohorts = mean(noCohorts)),
                                 by = .(scenario, year, firePresAbs, vegType)],
                 aes(x = year, y = noCohorts, colour = as.factor(vegType))) +
  # geom_vline(xintercept = fireYears, size = 1, linetype = "dashed", colour = "grey") +
  geom_line(size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Avg. no. cohorts by forest type", y = "no. cohorts",
       subtitle = "presence/absence of fire") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot12 <- ggplot(data = plotData[, list(noCohorts = mean(noCohorts)),
                                 by = .(scenario, year, noFires, vegType)],
                 aes(x = year, y = noCohorts, colour = as.factor(vegType))) +
  # geom_vline(xintercept = fireYears, size = 1, linetype = "dashed", colour = "grey") +
  geom_line(size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Avg. no. cohorts by forest type", y = "no. cohorts",
       subtitle = "no. fires") +
  facet_grid(scenario ~ noFires)

plotData <- summaryBurnCohortDataVegType[, list(BiomassBySpecies = sum(BiomassBySpecies)),
                                         by = .(scenario, year, vegType, firePresAbs)]
plot13 <- ggplot(data = plotData,
                 aes(x = year, y = log(BiomassBySpecies + 0.000001), colour = as.factor(vegType))) +
  # geom_vline(xintercept = fireYears, size = 1, linetype = "dashed", colour = "grey") +
  geom_line(size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Total landscape biomass by forest type", y = expression(paste("log-Biomass"~~"(g/m"^2, ")")),
       subtitle = "presence/absence of fire") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot14 <- ggplot(data = summaryBurnCohortDataVegType,
                 aes(x = year, y = log(BiomassBySpecies + 0.000001), colour = as.factor(vegType))) +
  # geom_vline(xintercept = fireYears, size = 1, linetype = "dashed", colour = "grey") +
  geom_line(size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Total landscape biomass by forest type", y = expression(paste("log-Biomass"~~"(g/m"^2, ")")),
       subtitle = "no. fires") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ noFires)

plotData <- summaryBurnCohortDataVegType[, list(MortalityBySpecies = sum(MortalityBySpecies)),
                                         by = .(scenario, year, vegType, firePresAbs)]
plot13 <- ggplot(data = plotData,
                 aes(x = year, y = log(MortalityBySpecies + 0.000001), colour = as.factor(vegType))) +
  geom_line(size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Total landscape mortality by forest type", y = expression(paste("log-Biomass"~~"(g/m"^2, ")")),
       subtitle = "presence/absence of fire") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot14 <- ggplot(data = summaryBurnCohortDataVegType,
                 aes(x = year, y = log(MortalityBySpecies + 0.000001), colour = as.factor(vegType))) +
  # geom_vline(xintercept = fireYears, size = 1, linetype = "dashed", colour = "grey") +
  geom_line(size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Total landscape mortality by forest type", y = expression(paste("log-Biomass"~~"(g/m"^2, ")")),
       subtitle = "no. fires") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ noFires)


## SAVE PLOTS
amc::.gc()
ggpubr::ggarrange(plot1,
                  plot2 +
                    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
                    labs(title = "", subtitle = ""),
                  widths = c(0.4, 0.6),
                  legend = "bottom", common.legend = TRUE)
ggsave(filename = file.path(figOutputPath, "results_landscapeB.tiff"),
       width = 14, height = 7)

ggpubr::ggarrange(plot13,
                  plot14 +
                    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
                    labs(title = "", subtitle = ""),
                  widths = c(0.4, 0.6),
                  legend = "bottom", common.legend = TRUE)
ggsave(filename = file.path(figOutputPath, "results_landscapeBVegType.tiff"),
       width = 14, height = 7)

ggpubr::ggarrange(plot3,
                  plot4 +
                    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
                    labs(title = "", subtitle = ""),
                  widths = c(0.4, 0.6),
                  legend = "bottom", common.legend = TRUE)
ggsave(filename = file.path(figOutputPath, "results_landscapeMort.tiff"),
       width = 14, height = 7)

ggpubr::ggarrange(plot5,
                  plot6 +
                    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
                    labs(title = "", subtitle = ""),
                  widths = c(0.4, 0.6),
                  legend = "bottom", common.legend = TRUE)
ggsave(filename = file.path(figOutputPath, "results_landscapeAge.tiff"),
       width = 14, height = 7)

ggpubr::ggarrange(plot9.2,
                  plot10.2 +
                    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
                    labs(title = "", subtitle = ""),
                  widths = c(0.4, 0.6),
                  legend = "bottom", common.legend = TRUE)
ggsave(filename = file.path(figOutputPath, "results_landscapeVegTypes.tiff"),
       width = 14, height = 7)

ggpubr::ggarrange(plot7,
                  plot8 +
                    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
                    labs(title = "", subtitle = ""),
                  widths = c(0.4, 0.6),
                  legend = "bottom", common.legend = TRUE)
ggsave(filename = file.path(figOutputPath, "results_noCohorts.tiff"),
       width = 14, height = 7)
ggpubr::ggarrange(plot11,
                  plot12 +
                    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
                    labs(title = "", subtitle = ""),
                  widths = c(0.4, 0.6),
                  legend = "bottom", common.legend = TRUE)
ggsave(filename = file.path(figOutputPath, "results_noCohortsVegType.tiff"),
       width = 14, height = 7)


