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

## path to figure folder
figOutputPath <- "C:/Users/Ceres Barros/Google Drive/Shared/Landscapes In Motion/ModellingTeam/reportFigs/"

## GET SIM LISTS
simList_PM <- qread("R/SpaDES/outputs/PM_newSppParams_fullSA/simList_PM_newSppParams_fullSA.qs")
simList_noPM <- qread("R/SpaDES/outputs/noPM_newSppParams_fullSA/simList_noPM_newSppParams_fullSA.qs")

outputs_PM <- as.data.table(outputs(simList_PM))
outputs_noPM <- as.data.table(outputs(simList_noPM))

## GET FILE NAMES
cohortDataFiles <- c(outputs_noPM[objectName == "cohortData", file],
                     outputs_PM[objectName == "cohortData", file]) %>%
  unique(.)

rstCurrentBurnFiles <- c(outputs_noPM[objectName == "rstCurrentBurn", file],
                         outputs_PM[objectName == "rstCurrentBurn", file]) %>%
  unique(.)
pixelGroupMapFiles <- c(outputs_noPM[objectName == "pixelGroupMap", file],
                        outputs_PM[objectName == "pixelGroupMap", file]) %>%
  unique(.)
vegTypeMapFiles <- c(outputs_noPM[objectName == "vegTypeMap", file],
                     outputs_PM[objectName == "vegTypeMap", file]) %>%
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

names(rstCurrentBurnStk_noPM) <- sub(".*year", "year", sub("\\.rds", "", grep("noPM", rstCurrentBurnFiles, value = TRUE)))
names(rstCurrentBurnStk_PM) <- sub(".*year", "year", sub("\\.rds", "", grep("noPM", rstCurrentBurnFiles, value = TRUE, invert = TRUE)))
names(pixelGroupMapStk_noPM) <- sub(".*year", "year", sub("\\.rds", "", grep("noPM", pixelGroupMapFiles, value = TRUE)))
names(pixelGroupMapStk_PM) <- sub(".*year", "year", sub("\\.rds", "", grep("noPM", pixelGroupMapFiles, value = TRUE, invert = TRUE)))
names(vegTypeMapStk_noPM) <- sub(".*year", "year", sub("\\.rds", "", grep("noPM", vegTypeMapFiles, value = TRUE)))
names(vegTypeMapStk_PM) <- sub(".*year", "year", sub("\\.rds", "", grep("noPM", vegTypeMapFiles, value = TRUE, invert = TRUE)))

## BUILD TABLES OF RESULTS ------------------------
## pixelCohortData tables
files <- grep("noPM", cohortDataFiles, value = TRUE)
pixelCohortData_noPM <- lapply(files, FUN = function(ff, pixelGroupMapStk) {
  cohortData <- readRDS(ff)
  yr <- sub(".*year", "",  sub(".rds", "", ff))
  pixelGroupMap <- pixelGroupMapStk[[paste0("year", yr)]]
  cohortData <- addPixels2CohortData(cohortData, pixelGroupMap)
  cohortData[, year := as.integer(yr)]
  return(cohortData)
}, pixelGroupMapStk = pixelGroupMapStk_noPM) %>%
  rbindlist(fill = TRUE, l = .)

files <- grep("noPM", cohortDataFiles, value = TRUE, invert = TRUE)
pixelCohortData_PM <- lapply(files, FUN = function(ff, pixelGroupMapStk) {
  cohortData <- readRDS(ff)
  yr <- sub(".*year", "",  sub(".rds", "", ff))
  pixelGroupMap <- pixelGroupMapStk[[paste0("year", yr)]]
  cohortData <- addPixels2CohortData(cohortData, pixelGroupMap)
  cohortData[, year := as.integer(yr)]
  return(cohortData)
}, pixelGroupMapStk = pixelGroupMapStk_PM) %>%
  rbindlist(fill = TRUE, l = .)

## vegTypeData tables
vegTypeSubset <- intersect(names(vegTypeMapStk_noPM), names(pixelGroupMapStk_noPM))
vegTypeData_noPM <- lapply(vegTypeSubset, FUN = function(x) {
  data.table(pixelIndex = seq_len(ncell(vegTypeMapStk_noPM[[x]])),
             pixelGroup = getValues(pixelGroupMapStk_noPM[[x]]),
             vegType = vegTypeMapStk_noPM[[x]][],
             year = as.integer(sub("year", "", x)))
}) %>%
  rbindlist(.)
vegTypeData_noPM <- vegTypeData_noPM[!is.na(pixelGroup)]

vegTypeSubset <- intersect(names(vegTypeMapStk_PM), names(pixelGroupMapStk_PM))
vegTypeData_PM <- lapply(vegTypeSubset, FUN = function(x) {
  data.table(pixelIndex = seq_len(ncell(vegTypeMapStk_PM[[x]])),
             pixelGroup = getValues(pixelGroupMapStk_PM[[x]]),
             vegType = vegTypeMapStk_PM[[x]][],
             year  = as.integer(sub("year", "", x)))
})  %>%
  rbindlist(.)
vegTypeData_PM <- vegTypeData_PM[!is.na(pixelGroup)]

## pixelBurnData tables
fireSubset <- intersect(names(rstCurrentBurnStk_noPM), names(pixelGroupMapStk_noPM))
pixelBurnData_noPM <- lapply(fireSubset, FUN = function(x) {
  data.table(pixelIndex = seq_len(ncell(rstCurrentBurnStk_noPM[[x]])),
             pixelGroup = getValues(pixelGroupMapStk_noPM[[x]]),
             burnt = as.integer(!is.na(rstCurrentBurnStk_noPM[[x]][])),
             year = as.integer(sub("year", "", x)))
}) %>%
  rbindlist(.)
pixelBurnData_noPM <- pixelBurnData_noPM[!is.na(pixelGroup)]

fireSubset <- intersect(names(rstCurrentBurnStk_PM), names(pixelGroupMapStk_PM))
pixelBurnData_PM <- lapply(fireSubset, FUN = function(x) {
  data.table(pixelIndex = seq_len(ncell(rstCurrentBurnStk_PM[[x]])),
             pixelGroup = getValues(pixelGroupMapStk_PM[[x]]),
             burnt = as.integer(!is.na(rstCurrentBurnStk_PM[[x]][])),
             year = as.integer(sub("year", "", x)))
})  %>%
  rbindlist(.)
pixelBurnData_PM <- pixelBurnData_PM[!is.na(pixelGroup)]

rm(fireSubset, vegTypeSubset)

## join tables, add scenario col and rbind.
## note that vegTypeData and pixelBurnData have more pixels because PGs of 0 are there, but not on cohortData
## so join by keeping all pixels
pixelCohortData_noPM <- merge(pixelCohortData_noPM, vegTypeData_noPM,
                              by = c("pixelIndex", "pixelGroup", "year"), all = TRUE)
pixelCohortData_noPM <- merge(pixelCohortData_noPM, pixelBurnData_noPM,
                              by = c("pixelIndex", "pixelGroup", "year"), all = TRUE)
pixelCohortData_noPM[is.na(burnt), burnt := 0]
pixelCohortData_noPM[, scenario := "noPM"]

pixelCohortData_PM <- merge(vegTypeData_PM, pixelCohortData_PM,
                            by = c("pixelIndex", "pixelGroup", "year"), all = TRUE)
pixelCohortData_PM <- merge(pixelBurnData_PM, pixelCohortData_PM,
                            by = c("pixelIndex", "pixelGroup", "year"), all = TRUE)
pixelCohortData_PM[is.na(burnt), burnt := 0]
pixelCohortData_PM[, scenario := "PM"]

allPixelCohortData <- rbind(pixelCohortData_noPM, pixelCohortData_PM, use.names = TRUE)
rm(list = grep("^[pixel|vegType].*Data", ls(), value = TRUE))

## ADD MISSING SPECIES IN YEAR/SCENARIO/PIXEL COMBINATION
## cohortData doens't track absent cohorts, so they need to ba added back
## for now pixels from the 0s pixelGroup  have one entry with NAs for speciesCode
## they will be ignored for now and removed later, after adding one species entry for each of these pixels.
## for reporting consistency add to show losses in B
combinations <- unique(allPixelCohortData[, .(scenario, year, pixelIndex, pixelGroup, burnt)])
spp <- as.character(na.omit(unique(allPixelCohortData$speciesCode)))
combinations <- lapply(spp, FUN = function(x) {
  data.table(combinations,
             speciesCode = x)
}) %>%
  rbindlist(., use.names = TRUE)

## join while keeping all combos, NA species will now disappear.
allPixelCohortData <- allPixelCohortData[combinations,
                                         on = .(scenario, year, pixelIndex,
                                                pixelGroup, burnt, speciesCode)]
rm(spp, combinations)
amc::.gc()


## NO. FIRES ---------------------------------------
## add no. fires per pixel
## how many times did each pixel burn?
noFiresPixel <- unique(allPixelCohortData[, .(pixelIndex, scenario, year, burnt)])
noFiresPixel[, noFires := sum(burnt), by = .(pixelIndex, scenario)]
noFiresPixel <- unique(noFiresPixel[, .(pixelIndex, scenario, noFires)])

test <- sapply(unique(noFiresPixel$scenario), FUN = function(x){
  any(duplicated(noFiresPixel[scenario == x, pixelIndex]))
})
if (any(test))
  stop("Each pixel should only have one record of no. fires per scenario")

allPixelCohortData <- noFiresPixel[allPixelCohortData, on = .(pixelIndex, scenario)]

## add presence/absence of fire across simulation per pixel/scenario
amc::.gc()
allPixelCohortData[, firePresAbs := as.integer(any(noFires > 0)), by = .(scenario, pixelIndex)]
rm(noFiresPixel)
amc::.gc()

## SUBSET A FEW YEARS ------------------------------
subsetYrs <- c(seq(range(allPixelCohortData$year)[1], range(allPixelCohortData$year)[2],
                   5), range(allPixelCohortData$year)[2])
allPixelCohortData <- allPixelCohortData[year %in% subsetYrs]
amc::.gc()

test <- length(unique(allPixelCohortData[, length(unique(pixelIndex)), by = .(scenario, year)]$V1)) == 1
test2 <- length(unique(allPixelCohortData[, length(unique(speciesCode)), by = .(scenario, year, pixelIndex)]$V1)) == 1
test3 <- any(is.na(allPixelCohortData$speciesCode))

if (!test | !test2 | test3)
  stop("something's wrong")

## add ecoregion group where it's missing
allPixelCohortData[, ecoregionGroup := unique(na.omit(ecoregionGroup)),
                   by = .(scenario, year, pixelGroup)]
amc::.gc()

## add vegType where it's missing, but it's a pixel with some veg
allPixelCohortData[, vegType := max(vegType, na.rm = TRUE),
                   by = .(scenario, year, pixelGroup)]
amc::.gc()

## replace NAs of cohortData by 0s and add missing spp.
cols <- c("age", "B", "mortality", "aNPPAct", "vegType")
replaceNAs <- function(x, val = 0) {
  x[is.na(x)] <- val
  x
}

allPixelCohortData[, (cols) := lapply(.SD, replaceNAs), .SDcols = cols]
amc::.gc()

## add no. cohorts per pixel(group)
allPixelCohortData[B > 0, noCohorts := length(unique(paste(speciesCode, age))), by = .(scenario, year, pixelGroup)]
allPixelCohortData[, noCohorts := max(noCohorts, na.rm = TRUE), by = .(scenario, year, pixelGroup)]
allPixelCohortData[is.na(noCohorts), noCohorts := 0]
amc::.gc()

## SUMMARY ACROSS LANDSCAPE -----------------------------------
## remember that biomass is multiplied by 100 in *boreal*, this will revert the units to tonnes/ha
summaryBurnCohortData <- allPixelCohortData[, list(BiomassBySpecies = as.numeric(sum((B/100), na.rm = TRUE)),
                                                   MortalityBySpecies = as.numeric(sum((mortality/100), na.rm = TRUE)),
                                                   aNPPBySpecies = as.numeric(sum((aNPPAct/100), na.rm = TRUE)),
                                                   AgeBySppWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                                   sum((B/100), na.rm = TRUE))),
                                            by = .(scenario, year, noFires, speciesCode)]
summaryBurnCohortData[is.nan(AgeBySppWeighted), AgeBySppWeighted := 0]
summaryBurnCohortData[, firePresAbs := as.integer(noFires != 0)]
amc::.gc(); amc::.gc()

## make species labels/colours
speciesLabels <- LandR::equivalentName(value = unique(summaryBurnCohortData$speciesCode), column = "EN_generic_full",
                                       df = simList_noPM$sppEquiv)
names(speciesLabels) <- unique(summaryBurnCohortData$speciesCode)

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


## PLOTS -----------------------------------------------
fireYears <- as.numeric(sub("year", "", names(rstCurrentBurnStk_PM)))

plotData <- summaryBurnCohortData[, list(BiomassBySpecies = sum(BiomassBySpecies)),
                                  by = .(scenario, year, speciesCode, firePresAbs)]
plot1 <- ggplot(data = plotData,
                aes(x = year, y = log(BiomassBySpecies + 0.000001), colour = speciesCode)) +
  # geom_vline(xintercept = fireYears, size = 1, linetype = "dashed", colour = "grey") +
  geom_line(size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours, labels = speciesLabels) +
  labs(title = "Total landscape biomass", y = expression(paste("log-Biomass"~~"(g/m"^2, ")")),
       subtitle = "presence/absence of fire") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot2 <- ggplot(data = summaryBurnCohortData,
                aes(x = year, y = log(BiomassBySpecies + 0.000001), colour = speciesCode)) +
  # geom_vline(xintercept = fireYears, size = 1, linetype = "dashed", colour = "grey") +
  geom_line(size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours, labels = speciesLabels) +
  labs(title = "Total landscape biomass", y = expression(paste("log-Biomass"~~"(g/m"^2, ")")),
       subtitle = "no. fires") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ noFires)

plotData <- summaryBurnCohortData[, list(MortalityBySpecies = sum(MortalityBySpecies)),
                                  by = .(scenario, year, speciesCode, firePresAbs)]
plot3 <- ggplot(data = plotData,
                aes(x = year, y = log(MortalityBySpecies + 0.000001), colour = speciesCode)) +
  geom_line(size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours, labels = speciesLabels) +
  labs(title = "Total landscape mortality", y = expression(paste("log-Biomass"~~"(g/m"^2, ")")),
       subtitle = "presence/absence of fire") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot4 <- ggplot(data = summaryBurnCohortData,
                aes(x = year, y = log(MortalityBySpecies + 0.000001), colour = speciesCode)) +
  # geom_vline(xintercept = fireYears, size = 1, linetype = "dashed", colour = "grey") +
  geom_line(size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours, labels = speciesLabels) +
  labs(title = "Total landscape mortality", y = expression(paste("log-Biomass"~~"(g/m"^2, ")")),
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
  labs(title = "Avg. age across landscape (biomass-weighted)", y = "years",
       subtitle = "presence/absence of fire") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot6 <- ggplot(data = summaryBurnCohortData,
                aes(x = year, y = AgeBySppWeighted, colour = speciesCode)) +
  # geom_vline(xintercept = fireYears, size = 1, linetype = "dashed", colour = "grey") +
  geom_line(size = 1) +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours, labels = speciesLabels) +
  labs(title = "Avg. age across landscape (biomass-weighted)", y = "years",
       subtitle = "no. fires") +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ noFires)

amc::.gc()
plotData <- allPixelCohortData[, list(noPixelsVeg = length(unique(pixelIndex))),
                               by = .(scenario, year, firePresAbs, vegType)]
plot7 <- ggplot(data = plotData,
                aes(x = year, y = noPixelsVeg, fill = as.factor(vegType))) +
  geom_area(stat = "identity", position = "stack") +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Dominant species", y = "no. pixels",
       subtitle = "presence/absence of fire") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

plot7.2 <- ggplot(data = plotData,
                aes(x = year, y = noPixelsVeg, fill = as.factor(vegType))) +
  geom_area(stat = "identity", position = "fill") +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Dominant species", y = "prop. pixels",
       subtitle = "presence/absence of fire") +
  facet_grid(scenario ~ firePresAbs,
             labeller = labeller(firePresAbs = c("0" = "no fire", "1" = "fire")))

amc::.gc()
plotData <- allPixelCohortData[, list(noPixelsVeg = length(unique(pixelIndex))),
                               by = .(scenario, year, noFires, vegType)]
plot8 <- ggplot(data = plotData,
                aes(x = year, y = noPixelsVeg, fill = as.factor(vegType))) +
  # geom_vline(xintercept = fireYears, size = 1, linetype = "dashed", colour = "grey") +
  geom_area(stat = "identity", position = "stack") +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Dominant species", y = "no. pixels",
       subtitle = "no. fires") +
  facet_grid(scenario ~ noFires)

plot8.2 <- ggplot(data = plotData,
                aes(x = year, y = noPixelsVeg, fill = as.factor(vegType))) +
  # geom_vline(xintercept = fireYears, size = 1, linetype = "dashed", colour = "grey") +
  geom_area(stat = "identity", position = "fill") +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Dominant species", y = "prop. pixels",
       subtitle = "no. fires") +
  facet_grid(scenario ~ noFires)

## SAVE PLOTS
ggpubr::ggarrange(plot1,
                  plot2 +
                    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
                    labs(title = "", subtitle = ""),
                  widths = c(0.4, 0.6),
                  legend = "bottom", common.legend = TRUE)
ggsave(filename = file.path(figOutputPath, "results_landscapeB.tiff"),
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

ggpubr::ggarrange(plot7.2,
                  plot8.2 +
                    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
                    labs(title = "", subtitle = ""),
                  widths = c(0.4, 0.6),
                  legend = "bottom", common.legend = TRUE)
ggsave(filename = file.path(figOutputPath, "results_landscapeVegTypes.tiff"),
       width = 14, height = 7)


