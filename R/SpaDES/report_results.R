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
             vegType = as.integer(!is.na(vegTypeMapStk_noPM[[x]][])),
             year = as.integer(sub("year", "", x)))
}) %>%
  rbindlist(.)
vegTypeData_noPM <- vegTypeData_noPM[!is.na(pixelGroup)]

vegTypeSubset <- intersect(names(vegTypeMapStk_PM), names(pixelGroupMapStk_PM))
vegTypeData_PM <- lapply(vegTypeSubset, FUN = function(x) {
  data.table(pixelIndex = seq_len(ncell(vegTypeMapStk_PM[[x]])),
             pixelGroup = getValues(pixelGroupMapStk_PM[[x]]),
             vegType = as.integer(!is.na(vegTypeMapStk_PM[[x]][])),
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

## join tables, add scenario col and rbind
pixelCohortData_noPM <- vegTypeData_noPM[pixelCohortData_noPM, on = .(pixelIndex, pixelGroup, year)]
pixelCohortData_noPM <- pixelBurnData_noPM[pixelCohortData_noPM, on = .(pixelIndex, pixelGroup, year)]
pixelCohortData_noPM[is.na(burnt), burnt := 0]
pixelCohortData_noPM[, scenario := "noPM"]

pixelCohortData_PM <- vegTypeData_PM[pixelCohortData_PM, on = .(pixelIndex, pixelGroup, year)]
pixelCohortData_PM <- pixelBurnData_PM[pixelCohortData_PM, on = .(pixelIndex, pixelGroup, year)]
pixelCohortData_PM[is.na(burnt), burnt := 0]
pixelCohortData_PM[, scenario := "PM"]

allPixelCohortData <- rbind(pixelCohortData_noPM, pixelCohortData_PM, use.names = TRUE)
rm(list = grep("^[pixel|vegType].*Data", ls(), value = TRUE))

## TABLE OF VEG TYPE LABELS -----------------------------------
vegTypeTable <- as.data.table(levels(vegTypeMapStk_noPM[[1]]))
allPixelCohortData <- vegTypeTable[, .(ID, VALUE)][allPixelCohortData, on = "ID==vegType"]
setnames(allPixelCohortData, old = c("ID", "VALUE"),
         new = c("vegTypeID", "vegType"))


## SUMMARY ACROSS LANDSCAPE -----------------------------------
## remember that biomass is multiplied by 100 in *boreal*
## also, B needs to be rescaled, to avoid integer overflow
allPixelCohortData[, noPixels := length(pixelIndex), by = pixelGroup]

summaryBurnCohortData <- allPixelCohortData[, list(BiomassBySpecies = as.numeric(sum((B/100) * noPixels, na.rm = TRUE)),
                                                   MortalityBySpecies = as.numeric(sum((mortality/100) * noPixels, na.rm = TRUE)),
                                                   aNPPBySpecies = as.numeric(sum((aNPPAct/100) * noPixels, na.rm = TRUE)),
                                                   AgeBySppWeighted = as.numeric(sum(age * (B/100) * noPixels, na.rm = TRUE) /
                                                                                   sum((B/100) * noPixels, na.rm = TRUE)),
                                                   noCohorts = as.numeric(length(unique(age)))),
                                            by = .(scenario, year, burnt, speciesCode)]
## missing species in an EXISTING year, scenario, burn combo are of B == 0,
## add them to show losses in B
combinations <- lapply(as.character(unique(summaryBurnCohortData$speciesCode)),
                       FUN = function(x) {
                         combinations <- unique(summaryBurnCohortData[, .(year, scenario, burnt)])
                         combinations[, speciesCode := x]
                         return(combinations)
                       }) %>%
  rbindlist(., use.names = TRUE)

## join while keeping all combos
summaryBurnCohortData <- summaryBurnCohortData[combinations,
                                      on = c("year", "scenario", "burnt", "speciesCode")]
## replace NA's to 0s by converting to matrix
summaryBurnCohortData <- as.matrix(summaryBurnCohortData)
summaryBurnCohortData[is.na(summaryBurnCohortData)] <- 0
summaryBurnCohortData <- as.data.table(summaryBurnCohortData)
cols <- grep("scenario|species", names(summaryBurnCohortData), invert = TRUE, value = TRUE)
summaryBurnCohortData <- summaryBurnCohortData[, (cols) := lapply(.SD, as.numeric),
                                               .SDcols = cols]
summaryBurnCohortData <- unique(summaryBurnCohortData)

## make species labels/colours
speciesLabels <- LandR::equivalentName(value = unique(summaryBurnCohortData$speciesCode), column = "EN_generic_full",
                                       df = simList_noPM$sppEquiv)
names(speciesLabels) <- unique(summaryBurnCohortData$speciesCode)

speciesColours <- levels(vegTypeMapStk_noPM[[1]])[[1]]$colors
names(speciesColours) <- levels(vegTypeMapStk_noPM[[1]])[[1]]$VALUE


fireYears <- as.numeric(sub("year", "", names(rstCurrentBurnStk_PM)))
plot1 <- ggplot(data = summaryBurnCohortData,
                aes(x = year, y = log(BiomassBySpecies/100 + 0.000001), colour = speciesCode)) +
  # geom_vline(xintercept = fireYears, size = 1, linetype = "dashed", colour = "grey") +
  geom_line(size = 1) +
  theme_pubr(base_size = 16, legend = "right") +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours,
                      labels = speciesLabels) +
  labs(title = "Total landscape biomass", y = expression(paste("log-Biomass"~~"(g/m"^2, ")"))) +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  facet_grid(scenario ~ burnt,
             labeller = labeller(burnt = c("0" = "no fire", "1" = "fire")))

plot2 <- ggplot(data = summaryBurnCohortData,
                aes(x = year, y = log(MortalityBySpecies), colour = speciesCode)) +
  geom_vline(xintercept = 5, size = 1.5, linetype = "dashed", colour = "red") +
  geom_line(size = 1.5) +
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours,
                      labels = speciesLabels) +
  facet_grid(burnt ~ scenario,
             labeller = labeller(burnt = c("0" = "no fire", "1" = "fire")))

plot3 <- ggplot(data = allPixelCohortData,
                aes(x = year, fill = vegType)) +
  geom_area(stat = "count", position = "fill") +
  geom_vline(xintercept = 5, size = 1.5, linetype = "dashed", colour = "red") +
  scale_fill_manual(values = speciesColours,
                    labels = speciesLabels) +
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.title = element_blank()) +
  labs(title = "No. pixels per vegetation type", y = "g/m^2") +
  facet_grid(burnt ~ scenario,
             labeller = labeller(burnt = c("0" = "no fire", "1" = "fire")))

plot4 <- ggplot(data = summaryBurnCohortData,
                aes(x = year, y = AgeBySppWeighted, colour = speciesCode)) +
  geom_vline(xintercept = 5, size = 1.5, linetype = "dashed", colour = "red") +
  geom_line(size = 1.5) +
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours,
                      labels = speciesLabels) +
  facet_grid(burnt ~ scenario,
             labeller = labeller(burnt = c("0" = "no fire", "1" = "fire")))

plot5 <- ggplot(data = summaryBurnCohortData,
                aes(x = year, y = noCohorts, colour = speciesCode)) +
  geom_vline(xintercept = 5, size = 1.5, linetype = "dashed", colour = "red") +
  geom_line(size = 1.5) +
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours,
                      labels = speciesLabels) +
  facet_grid(burnt ~ scenario,
             labeller = labeller(burnt = c("0" = "no fire", "1" = "fire")))

ggpubr::ggarrange(plot1, plot2, plot3, nrow = 3,
                  legend = "right", common.legend = TRUE)
ggsave(filename = "R/SpaDES/outputs/blogSep2019_noPM_PM_BMortVegType.tiff",
       width = 10, height = 15)

## GIFS ------------------------------------------------------
## make GIFs of vegetation maps
## individual pics first
speciesLabels <- c("Fir", "Larch", "En. spruce", "Wh. spruce",
                   "Bl. spruce", "Lo. pine", "Aspen","Douglas-fir")
names(speciesLabels) <- levels(vegTypeMapStk_noPM[[1]])[[1]]$ID
speciesColours <- levels(vegTypeMapStk_noPM[[1]])[[1]]$colors
names(speciesColours) <-  levels(vegTypeMapStk_noPM[[1]])[[1]]$ID

foothillsMask <- simList_noPM$rawBiomassMap
foothillsMask[!is.na(foothillsMask)] <- 1
foothillsMaskDF <- as.data.frame(as(foothillsMask, "SpatialPixelsDataFrame"))
names(foothillsMaskDF) <- c("value", "x", "y")

makePNGs <- function(id, rasterStack, filePrefix, gif.dir) {
  suppressWarnings(dir.create(gif.dir, recursive = TRUE))
  cat("Making ", id, "\n")
  rasterVis::gplot(rasterStack[[id]],
                   maxpixels = ncell(rasterStack[[id]])) +
    geom_tile(data = foothillsMaskDF,
              aes(x = x, y = y), fill = "grey95") +
    geom_tile(aes(fill = as.factor(value))) +
    scale_fill_manual(values = speciesColours,
                      labels = speciesLabels) +
    theme_void() + theme(legend.position = "none", text = element_text(size = 20)) +
    labs(title = sub("year", "year ", names(rasterStack)[id])) +
    coord_equal()
  ggsave(file.path(gif.dir, paste0(filePrefix, id, ".png")),
         device = "png", width=5, height=10, dpi = 300, units = "in")
}

map_df(.x = 1:nlayers(vegTypeMapStk_noPM), .f = makePNGs,
       rasterStack = vegTypeMapStk_noPM,
       filePrefix = "vegTypeMapStk_noPM",
       gif.dir = "R/SpaDES/outputs/blogSep2019_noPM_oneFire/gif")
map_df(.x = 1:nlayers(vegTypeMapStk_PM), .f = makePNGs,
       rasterStack = vegTypeMapStk_PM,
       filePrefix = "vegTypeMapStk_PM",
       gif.dir = "R/SpaDES/outputs/blogSep2019_PM_oneFire/gif")

makeGIF <- function(gif.dir, gifPrefix, ...) {
  ## get file list and the file numbers (to sort numerically, rather than alphabetically)
  PNGlist <- list.files(path = gif.dir, pattern = "*.png")
  fileNos <- as.numeric(sub("\\.png", "", sub("^\\D*(\\d)", "\\1", PNGlist)))

  ## make GIF
  file.path(gif.dir, PNGlist[order(fileNos)]) %>%
    map(image_read) %>% # reads each path file
    image_join() %>% # joins image
    image_animate(...) %>% # animates, can opt for number of loops
    image_write(file.path(gif.dir, paste0(gifPrefix, ".gif")))
}

makeGIF(gif.dir = "R/SpaDES/outputs/blogSep2019_noPM_oneFire/gif",
        gifPrefix = "vegTypeMapStk_noPM",
        fps = 2)
makeGIF(gif.dir = "R/SpaDES/outputs/blogSep2019_PM_oneFire/gif",
        gifPrefix = "vegTypeMapStk_PM",
        fps = 2)

## OTHER PLOTS --------------------------------------------------------
## topo and climate examples
slopeRas <- projectRaster(simList_noPM$slopeRas, foothillsMask)
topoPal <- colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "BrBG"))
plot(foothillsMask, col = "grey90", axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE)
plot(slopeRas, col = topoPal(20), axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE, add = TRUE)

temperatureRas <- projectRaster(simList_noPM$temperatureRas, foothillsMask)
tempPal <- colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "RdBu"))
plot(foothillsMask, col = "grey90", axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE)
plot(temperatureRas, col = tempPal(20), axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE, add = TRUE)

precipitationRas <- projectRaster(simList_noPM$precipitationRas, foothillsMask)
plot(foothillsMask, col = "grey90", axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE)
plot(precipitationRas, col = tempPal(20), axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE, add = TRUE)

## spp cover, age and biomass examples
sppPal <- colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "Blues"))
plot(foothillsMask, col = "grey90", axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE)
plot(simList_noPM$speciesLayers[[1]], col = sppPal(20), axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE, add = TRUE)

plot(foothillsMask, col = "grey90", axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE)
plot(simList_noPM$rawBiomassMap, axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE, add = TRUE)

plot(foothillsMask, col = "grey90", axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE)
agePal <- colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "Greens"))
plot(simList_noPM$standAgeMap, col = agePal(20), axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE, add = TRUE)

## ecodistricts
## https://www.statcan.gc.ca/eng/subjects/standard/environment/elc/12-607-x2018001-eng.pdf

ecoDistSF <- sf::st_as_sf(simList_noPM$ecoDistrict)
ecoDistSF$ECODISTRIC <- factor(ecoDistSF$ECODISTRIC,
                               levels = c(798, 800, 801, 799, 793, 750, 631, 1018, 1017, 1019))
canada <- sf::st_as_sf(shapefile("data/CA_admin/gpr_000a11a_e.shp"))
canada <- sf::st_transform(canada, crs = crs(ecoDistSF))
alberta <- canada[canada$PRENAME %in% "Alberta",]

ggplot(ecoDistSF) +
  geom_sf(data = alberta) +
  geom_sf(data = ecoDistSF, aes(fill = as.factor(ECODISTRIC))) +
  scale_fill_brewer(palette = "Paired",
                    labels = c("631" = "W AB upland - Foothills",
                               "750" = "Aspen parkland - Upland",
                               "793" = "Moist mixed grassland - Plain",
                               "798" = "Fescue grassland - Plain",
                               "799" = "Fescue grassland - Upland",
                               "800" = "Fescue grassland - Plain",
                               "801" = "Fescue grassland - Foothills",
                               "1017" = "N Cont. Divide - Mountains",
                               "1018" = "N Cont. Divide - Foothills",
                               "1019" = "N Cont. Divide - Mountains")) +
  theme_void() +
  theme(text = element_text(colour = "white"),
        plot.background = element_rect(fill = "black")) +
  labs(fill = "Ecoregion - ecodistrict") +
  coord_sf()

plot(sf::st_as_sf(simList_noPM$ecoDistrict["ECODISTRIC"]))


save(list = grep("model", ls(), value = TRUE), file = "E:/GitHub/LandscapesInMotion/analyses/modelsGAMLSS_0-3Days_goodSample_Nov1_v2.RData")
