## clean workspace
rm(list=ls()); amc::.gc()

library(data.table)
library(ggplot2)
library(raster)
library(quickPlot)
library(SpaDES)
library(purrr)
library(magick)
simList_PM <- readRDS("R/SpaDES/outputs/blogSep2019_PM_oneFire/simList_fakeRstCurrentBurnblogSep2019_PM_oneFire.rds")
simList_noPM <- readRDS("R/SpaDES/outputs/blogSep2019_noPM_oneFire/simList_fakeRstCurrentBurnblogSep2019_noPM_oneFire.rds")

cohortDataFiles <- c(list.files("R/SpaDES/outputs/blogSep2019_noPM_oneFire", pattern = "cohortData", full.names = TRUE),
                     list.files("R/SpaDES/outputs/blogSep2019_PM_oneFire", pattern = "cohortData", full.names = TRUE))

allCohortData <- rbindlist(fill = TRUE, l = lapply(cohortDataFiles, FUN = function(ff) {
  cohortData <- readRDS(ff)
  yr <- as.integer(sub(".*year", "",  sub(".rds", "", ff)))
  scen <- if (grepl("_noPM", ff)) "noPM" else "PM"
  cohortData[, Year := yr]
  cohortData[, Scenario := scen]
  return(cohortData)
}))

rstCurrentBurnFiles <- c(list.files("R/SpaDES/outputs/blogSep2019_noPM_oneFire/", pattern = "rstCurrentBurn", full.names = TRUE),
                         list.files("R/SpaDES/outputs/blogSep2019_PM_oneFire/", pattern = "rstCurrentBurn", full.names = TRUE))
pixelGroupMapFiles <- c(list.files("R/SpaDES/outputs/blogSep2019_noPM_oneFire/", pattern = "pixelGroupMap", full.names = TRUE),
                        list.files("R/SpaDES/outputs/blogSep2019_PM_oneFire/", pattern = "pixelGroupMap", full.names = TRUE))
vegTypeMapFiles <- c(list.files("R/SpaDES/outputs/blogSep2019_noPM_oneFire/", pattern = "vegTypeMap", full.names = TRUE),
                     list.files("R/SpaDES/outputs/blogSep2019_PM_oneFire/", pattern = "vegTypeMap", full.names = TRUE))

rstCurrentBurnStk_noPM <- stack(lapply(grep("_noPM", rstCurrentBurnFiles, value = TRUE), readRDS))
rstCurrentBurnStk_PM <- stack(lapply(grep("_PM", rstCurrentBurnFiles, value = TRUE), readRDS))
pixelGroupMapStk_noPM <- stack(lapply(grep("_noPM", pixelGroupMapFiles, value = TRUE), readRDS))
pixelGroupMapStk_PM <- stack(lapply(grep("_PM", pixelGroupMapFiles, value = TRUE), readRDS))
vegTypeMapStk_noPM <- stack(lapply(grep("_noPM", vegTypeMapFiles, value = TRUE), readRDS))
vegTypeMapStk_PM <- stack(lapply(grep("_PM", vegTypeMapFiles, value = TRUE), readRDS))

names(rstCurrentBurnStk_noPM) <- sub(".*year", "year", sub("\\.rds", "", grep("_noPM", rstCurrentBurnFiles, value = TRUE)))
names(rstCurrentBurnStk_PM) <- sub(".*year", "year", sub("\\.rds", "", grep("_PM", rstCurrentBurnFiles, value = TRUE)))
names(pixelGroupMapStk_noPM) <- sub(".*year", "year", sub("\\.rds", "", grep("_noPM", pixelGroupMapFiles, value = TRUE)))
names(pixelGroupMapStk_PM) <- sub(".*year", "year", sub("\\.rds", "", grep("_PM", pixelGroupMapFiles, value = TRUE)))
names(vegTypeMapStk_noPM) <- sub(".*year", "year", sub("\\.rds", "", grep("_noPM", vegTypeMapFiles, value = TRUE)))
names(vegTypeMapStk_PM) <- sub(".*year", "year", sub("\\.rds", "", grep("_PM", vegTypeMapFiles, value = TRUE)))

## cheat

if (grepl("_oneFire", rstCurrentBurnFiles[1])) {
  ## compile all rasters except rstCurrentBurn
  pixelBurnCohortData <- rbind(
    rbindlist(lapply(names(vegTypeMapStk_noPM), FUN = function(x) {
      data.table(pixelID = seq_len(ncell(pixelGroupMapStk_noPM[[x]])),
                 pixelGroup = getValues(pixelGroupMapStk_noPM[[x]]),
                 vegType = getValues(vegTypeMapStk_noPM[[x]]),
                 Year = as.integer(sub("year", "", x)),
                 Scenario = "noPM")
    })),
    rbindlist(lapply(names(vegTypeMapStk_PM), FUN = function(x) {
      data.table(pixelID = seq_len(ncell(pixelGroupMapStk_PM[[x]])),
                 pixelGroup = getValues(pixelGroupMapStk_PM[[x]]),
                 vegType = getValues(vegTypeMapStk_PM[[x]]),
                 Year = as.integer(sub("year", "", x)),
                 Scenario = "PM")
    }))
  )

  burnPix <- rbind(
    data.table(pixelID = seq_len(ncell(rstCurrentBurnStk_noPM)),
               burnt = as.integer(!is.na(rstCurrentBurnStk_noPM[[1]][])),
               Year = as.integer(sub("year", "", names(rstCurrentBurnStk_noPM))),
               Scenario = "noPM"),
    data.table(pixelID = seq_len(ncell(rstCurrentBurnStk_PM)),
               burnt = as.integer(!is.na(rstCurrentBurnStk_noPM[[1]][])),
               Year = as.integer(sub("year", "", names(rstCurrentBurnStk_PM))),
               Scenario = "PM")
  )

  ## checks
  if (length(setdiff(burnPix$pixelID, pixelBurnCohortData$pixelID)) |
      length(setdiff(pixelBurnCohortData$pixelID, burnPix$pixelID)))
    stop("pixel IDs differ between rstCurrentBurn and other output rasters")
  pixelBurnCohortData <- burnPix[, .(Scenario, pixelID, burnt)][pixelBurnCohortData, on = c("Scenario", "pixelID")]

  # pixelBurnCohortData <- na.omit(pixelBurnCohortData)

  ## checks
  testList <- split(pixelBurnCohortData[, .(burnt, Year)], by = "Year", keep.by = FALSE)
  testList <- lapply(testList, FUN = function(DT) DT[["burnt"]])

  test <- sapply(seq_along(testList), FUN =  function(n) {
    length(setdiff(testList[[n]], unlist(testList[-n])))
  })
  if (any(test > 0))
    stop("burnt pixel IDs differ between years")

} else {
  pixelBurnCohortData <- rbind(
    rbindlist(lapply(names(rstCurrentBurnStk_noPM), FUN = function(x) {
      data.table(pixelID = seq_len(ncell(rstCurrentBurnStk_noPM[[x]])),
                 pixelGroup = getValues(pixelGroupMapStk_noPM[[x]]),
                 burnt = as.integer(!is.na(rstCurrentBurnStk_PM[[x]][])),
                 vegType = getValues(vegTypeMapStk_noPM[[x]]),
                 Year = as.integer(sub("year", "", x)),
                 Scenario = "noPM")
    })),
    rbindlist(lapply(names(rstCurrentBurnStk_PM), FUN = function(x) {
      data.table(pixelID = seq_len(ncell(rstCurrentBurnStk_PM[[x]])),
                 pixelGroup = getValues(pixelGroupMapStk_PM[[x]]),
                 burnt = as.integer(!is.na(rstCurrentBurnStk_PM[[x]][])),
                 vegType = getValues(vegTypeMapStk_PM[[x]]),
                 Year = as.integer(sub("year", "", x)),
                 Scenario = "PM")
    }))
  )
  # pixelBurnCohortData <- na.omit(pixelBurnCohortData)
}

pixelBurnCohortData <- allCohortData[pixelBurnCohortData, nomatch = 0,
                                     on = c("Scenario", "Year", "pixelGroup")]

vegTypeTable <- as.data.table(levels(vegTypeMapStk_noPM[[1]]))
pixelBurnCohortData <- vegTypeTable[, .(ID, VALUE)][pixelBurnCohortData, on = "ID==vegType"]
setnames(pixelBurnCohortData, old = c("ID", "VALUE"),
         new = c("vegTypeID", "vegType"))

pixelBurnCohortData[, noPixels := length(pixelID), by = pixelGroup]

summaryBurnCohortData <- pixelBurnCohortData[, list(BiomassBySpecies = as.numeric(sum(B * noPixels, na.rm = TRUE)),
                                                    MortalityBySpecies = as.numeric(sum(mortality * noPixels, na.rm = TRUE)),
                                                    aNPPBySpecies = as.numeric(sum(aNPPAct * noPixels, na.rm = TRUE)),
                                                    AgeBySppWeighted = as.numeric(sum(age * B * noPixels, na.rm = TRUE) /
                                                                                    sum(B * noPixels, na.rm = TRUE)),
                                                    noCohorts = as.numeric(length(unique(age)))),
                                             by = .(Scenario, Year, burnt, speciesCode)]
## missing species in a year are of B == 0,
## add them to show losses in B
combinations <- as.data.table(expand.grid(unique(summaryBurnCohortData$Year), summaryBurnCohortData$Scenario,
                                          unique(summaryBurnCohortData$burnt), unique(summaryBurnCohortData$speciesCode)))
setnames(combinations, c("Year", "Scenario", "burnt", "speciesCode"))

summaryBurnCohortData <- summaryBurnCohortData[combinations,
                                               on = c("Year", "Scenario", "burnt", "speciesCode")]
## replace NA's to 0s by converting to matrix
summaryBurnCohortData <- as.matrix(summaryBurnCohortData)
summaryBurnCohortData[is.na(summaryBurnCohortData)] <- 0
summaryBurnCohortData <- as.data.table(summaryBurnCohortData)
cols <- grep("Scenario|species", names(summaryBurnCohortData), invert = TRUE, value = TRUE)
summaryBurnCohortData <- summaryBurnCohortData[, (cols) := lapply(.SD, as.numeric),
                                               .SDcols = cols]
summaryBurnCohortData <- unique(summaryBurnCohortData)

## make species labels
speciesLabels <- c("Abie_sp" = "Fir", "Lari_lar" = "Larch",
                   "Pice_eng" = "En. spruce", "Pice_gla" = "Wh. spruce",
                   "Pice_mar" = "Bl. spruce", "Pinu_sp" = "Lo. pine",
                   "Popu_sp" = "Aspen", "Pseu_men" = "Douglas-fir")
speciesColours <- levels(vegTypeMapStk_noPM[[1]])[[1]]$colors
names(speciesColours) <- levels(vegTypeMapStk_noPM[[1]])[[1]]$VALUE

plot1 <- ggplot(data = summaryBurnCohortData,
                aes(x = Year, y = log(BiomassBySpecies + 0.000001), colour = speciesCode)) +
  geom_vline(xintercept = 5, size = 1.5, linetype = "dashed", colour = "red") +
  geom_line(size = 1.5) +
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours,
                      labels = speciesLabels) +
  labs(title = "Total landscape biomass", y = "log-Biomass (g/m^2)") +
  facet_grid(Scenario ~ burnt,
             labeller = labeller(burnt = c("0" = "no fire", "1" = "fire")))

plot2 <- ggplot(data = summaryBurnCohortData,
                aes(x = Year, y = log(MortalityBySpecies), colour = speciesCode)) +
  geom_vline(xintercept = 5, size = 1.5, linetype = "dashed", colour = "red") +
  geom_line(size = 1.5) +
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours,
                      labels = speciesLabels) +
  facet_grid(burnt ~ Scenario,
             labeller = labeller(burnt = c("0" = "no fire", "1" = "fire")))

plot3 <- ggplot(data = pixelBurnCohortData,
                aes(x = Year, fill = vegType)) +
  geom_area(stat = "count", position = "fill") +
  geom_vline(xintercept = 5, size = 1.5, linetype = "dashed", colour = "red") +
  scale_fill_manual(values = speciesColours,
                      labels = speciesLabels) +
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.title = element_blank()) +
  labs(title = "No. pixels per vegetation type", y = "g/m^2") +
  facet_grid(burnt ~ Scenario,
             labeller = labeller(burnt = c("0" = "no fire", "1" = "fire")))

plot4 <- ggplot(data = summaryBurnCohortData,
                 aes(x = Year, y = AgeBySppWeighted, colour = speciesCode)) +
  geom_vline(xintercept = 5, size = 1.5, linetype = "dashed", colour = "red") +
  geom_line(size = 1.5) +
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours,
                      labels = speciesLabels) +
  facet_grid(burnt ~ Scenario,
             labeller = labeller(burnt = c("0" = "no fire", "1" = "fire")))

plot5 <- ggplot(data = summaryBurnCohortData,
                aes(x = Year, y = noCohorts, colour = speciesCode)) +
  geom_vline(xintercept = 5, size = 1.5, linetype = "dashed", colour = "red") +
  geom_line(size = 1.5) +
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours,
                      labels = speciesLabels) +
  facet_grid(burnt ~ Scenario,
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
    labs(title = sub("year", "Year ", names(rasterStack)[id])) +
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
