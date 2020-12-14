## ----------------------------------------
## PLOTS FOR REPORT
## ----------------------------------------
library(qs)
library(ggplot2)
library(ggpubr)
library(ggspatial)
library(SpaDES)
library(raster)
library(data.table)
library(magick)
library(purrr)

source("R/R_tools/plotFunction.R")

## path to figure folder
figOutputPath <- "C:/Users/Ceres Barros/Google Drive/Shared/McIntire-lab/Manuscripts_inPrep/LIMmodel_paper/figures"

## choose one final run for some plots.
runName <- c("noPM_rep1")

## load simLists
simOutPreSim <- qread(file.path("R/SpaDES/outputs/preSim", "preSimList.qs"))
simOutSim <- qread(file.path("R/SpaDES/outputs", runName, paste0("simList_", runName, ".qs")))

## FIRE LOCATIONS -------------------------
plot1 <- ggplot() +
  layer_spatial(data = simOutPreSim$studyArea, fill = "grey", colour = "black") +
  layer_spatial(data = simOutPreSim$fireLocations, colour = "darkred") +
  annotation_north_arrow(style = north_arrow_minimal,
                         location = "tr", which_north = "true") +
  theme_pubr(margin = FALSE) +
  labs(x = "longitude", y = "latitude")
ggsave(file.path(figOutputPath, "fireLocations.tiff"), plot1,
       width = 4, height = 7, dpi = 300)

## FUEL TYPES -------------------------
rasFuels <- setValues(raster(simOutPreSim$fuelTypesMaps$finalFuelType),
                      getValues(simOutPreSim$fuelTypesMaps$finalFuelType))
FTlabs <- raster::levels(simOutPreSim$fuelTypesMaps$finalFuelType)[[1]][,2]
names(FTlabs) <- raster::levels(simOutPreSim$fuelTypesMaps$finalFuelType)[[1]][,1]

FTlabs <- sub("O1b", "O1b: Grass - Standing", FTlabs)
FTlabs <- sub("NF", "NF: Non-fuel", FTlabs)
FTlabs <- sub("C2", "C2: Boreal spruce", FTlabs)
FTlabs <- sub("C3", "C3: Mature Jack or Lodgepole pine", FTlabs)
FTlabs <- sub("C4", "C4: Immature Jack or Lodgepole pine", FTlabs)
FTlabs <- sub("C7", "C7: Ponderosa pine/Douglas-fir", FTlabs)
FTlabs <- sub("M2", "M2: Boreal mixedwood – green", FTlabs)
FTlabs <- sub("D2", "D2: Green Aspen", FTlabs)

plot2 <- ggplot() +
  layer_spatial(rasFuels, aes(fill = stat(band1))) +
  layer_spatial(data = simOutPreSim$studyArea, fill = "transparent", colour = "black") +
  annotation_north_arrow(style = north_arrow_minimal,
                         location = "tr", which_north = "true") +
  theme_pubr(legend = "right", margin = FALSE) +
  scale_fill_distiller(palette = "Paired", breaks = sort(unique(rasFuels[])),
                       na.value = "transparent", guide = "legend",
                       labels = FTlabs) +
  labs(x = "longitude", y = "latitude", fill = "Fuel type") +
  guides(guide_legend(ncol = 1))
ggsave(file.path(figOutputPath, "fuelsMap.tiff"), plot2,
       width = 7, height = 7, dpi = 300)


## FUEL TYPE COVER
plotList <- lapply(unstack(simOutPreSim$fuelTypesCoverStk), plotFunction,
                   studyArea = simOutPreSim$studyArea, limits = c(0,1))
plot3 <- ggarrange(plotlist = plotList, ncol = 3, nrow = 3,
                   legend = "bottom", common.legend = TRUE)

ggsave(file.path(figOutputPath, "Fueltype_cover.tiff"), plot3,
       width = 9, height = 14, dpi = 300)


## TOPO, VEG, CLIMATE MAPS ----------------
## elevation
## if simLists are from experiment(), they will have no objects.
if (!is.null(simOutSim$DEMRas)) {
  if (!inMemory(simOutSim$DEMRas)) {
    DEMRas <- sub(".*LandscapesInMotion/", "",
                  filename(simOutSim$DEMRas))
    DEMRas <- raster(DEMRas)
  } else
    DEMRas <- simOutSim$DEMRas
} else {
  tempDT <- as.data.table(showCache(cachePath(simOutSim), "DEMRas"))
  DEMRas <- unique(tempDT[createdDate == max(createdDate)]$cacheId) %>%
    reproducible::loadFromCache(cachePath = cachePath(simOutSim), cacheId = .)
}

elevPlot <- ggplot() +
  layer_spatial(DEMRas, aes(fill = stat(band1))) +
  layer_spatial(data = simOutPreSim$studyArea, col = "black",
                fill = "transparent") +
  annotation_north_arrow(style = north_arrow_minimal,
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         location = "tr", which_north = "true") +
  theme_pubr(legend = "right") +
  theme(plot.margin = unit(c(0,0,0,0), units = "mm")) +
  scale_fill_distiller(palette = "YlOrBr", na.value = "transparent",
                       direction = 1) +
  labs(x = "longitude", y = "latitude", fill = "m a.s.l.",
       title = "Elevation")

## stand biomass (raw)
if (!is.null(simOutPreSim$rawBiomassMap)){
  if (!inMemory(simOutPreSim$rawBiomassMap)) {
    biomassRas <- sub(".*LandscapesInMotion/", "",
                      filename(simOutPreSim$rawBiomassMap))
    biomassRas <- raster(biomassRas)
  } else
    biomassRas <- simOutPreSim$rawBiomassMap
} else {
  tempDT <- as.data.table(showCache(cachePath(simOutPreSim), "rawBiomassMap"))
  biomassRas <- unique(tempDT[createdDate == max(createdDate)]$cacheId) %>%
    reproducible::loadFromCache(cachePath = cachePath(simOutPreSim), cacheId = .)
}

biomassPlot <- ggplot() +
  layer_spatial(biomassRas, aes(fill = stat(band1))) +
  layer_spatial(data = simOutPreSim$studyArea, col = "black",
                fill = "transparent") +
  annotation_north_arrow(style = north_arrow_minimal,
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         location = "tr", which_north = "true") +
  theme_pubr(legend = "right") +
  theme(plot.margin = unit(c(0,0,0,0), units = "mm")) +
  scale_fill_distiller(palette = "Greens", na.value = "transparent",
                       direction = 1) +
  labs(x = "longitude", y = "latitude", fill = expression("g/m"^2),
       title = "Stand biomass")

## Stand age (pre-corrections)
if (!is.null(simOutPreSim$standAgeMap)) {
  if (!inMemory(simOutPreSim$standAgeMap)) {
    ageRas <- sub(".*LandscapesInMotion/", "",
                  filename(simOutPreSim$standAgeMap))
    ageRas <- raster(ageRas)
  } else
    ageRas <- simOutPreSim$standAgeMap
} else {
  tempDT <- as.data.table(showCache(cachePath(simOutPreSim), "standAgeMap"))
  ageRas <- unique(tempDT[createdDate == max(createdDate)]$cacheId) %>%
    reproducible::loadFromCache(cachePath = cachePath(simOutPreSim), cacheId = .)
}

agePlot <- ggplot() +
  layer_spatial(ageRas, aes(fill = stat(band1))) +
  layer_spatial(data = simOutPreSim$studyArea, col = "black",
                fill = "transparent") +
  annotation_north_arrow(style = north_arrow_minimal,
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         location = "tr", which_north = "true") +
  theme_pubr(legend = "right") +
  theme(plot.margin = unit(c(0,0,0,0), units = "mm")) +
  scale_fill_distiller(palette = "Greens", na.value = "transparent",
                       direction = 1) +
  labs(x = "longitude", y = "latitude", fill = "years",
       title = "Stand age")

plot4 <- ggarrange(elevPlot, agePlot, biomassPlot, ncol = 3, nrow = 1)
ggsave(file.path(figOutputPath, "elev_StandAge_B.tiff"), plot4,
       width = 12, height = 7, dpi = 300)


## SPECIES COVER MAPS ---------------------
if (!is.null(simOutPreSim$speciesLayers)){
  if (!inMemory(simOutPreSim$speciesLayers)) {
    sppStk <- unstack(simOutPreSim$speciesLayers)
    sppStk <- sub(".*LandscapesInMotion/", "",
                  sapply(sppStk, filename))
    sppStk <- stack(sppStk)
    names(sppStk) <- sub(".*\\.", "", sub(".tif", "", names(sppStk)))
  } else
    sppStk <- simOutPreSim$speciesLayers
}  else {
  tempDT <- as.data.table(showCache(cachePath(simOutPreSim), "speciesLayers"))
  sppStk <- unique(tempDT[createdDate == max(createdDate)]$cacheId) %>%
    reproducible::loadFromCache(cachePath = cachePath(simOutPreSim), cacheId = .)
}


names(sppStk) <- LandR::equivalentName(value = names(sppStk), column = "EN_generic_full",
                                       df = simOutPreSim$sppEquiv)

## mask spp layers to remove 0s around SA
sppStk <- mask(sppStk, simOutPreSim$studyArea)

plotList <- lapply(unstack(sppStk), plotFunction,
                   studyArea = simOutPreSim$studyArea, limits = c(0, 100))
plot5 <- ggarrange(plotlist = plotList, ncol = 3, nrow = 2,
                   legend = "bottom", common.legend = TRUE)

ggsave(file.path(figOutputPath, "Species_cover.tiff"), plot5,
       width = 9, height = 10, dpi = 300)


## ECOREGIONS ---------------------------
## ecological zones
if (!is.null(simOutPreSim$ecoregionLayer)) {
    ecoDistSF <- sf::st_as_sf(simOutPreSim$ecoregionLayer)
} else {
  tempDT <- as.data.table(showCache(cachePath(simOutPreSim), "ecoregionLayer"))
  ecoDistSF <- unique(tempDT[createdDate == max(createdDate)]$cacheId) %>%
    reproducible::loadFromCache(cachePath = cachePath(simOutPreSim), cacheId = .) %>%
    sf::st_as_sf(.)
}

ecoDistSF$ecozoneCode2 <- factor(paste(ecoDistSF$NRNAME, ecoDistSF$NSRNAME, sep = " - "))
canada <- sf::st_as_sf(shapefile("data/CA_admin/gpr_000a11a_e.shp"))
canada <- sf::st_transform(canada, crs = crs(ecoDistSF))
alberta <- canada[canada$PRENAME %in% "Alberta",]

ecoDistPlot <- ggplot(ecoDistSF) +
  # geom_sf(data = alberta) +
  geom_sf(data = ecoDistSF, aes(fill = ecozoneCode2)) +
  layer_spatial(data = simOutPreSim$studyArea, fill = "transparent", colour = "black") +
  scale_fill_brewer(palette = "RdYlGn", direction = 1) +
  annotation_north_arrow(style = north_arrow_minimal,
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         location = "tr", which_north = "true") +
  theme_pubr(legend = "right") +
  theme(plot.margin = unit(c(0,0,0,0), units = "mm")) +
  labs(x = "longitude", y = "latitude", fill = "",
       title = "Natural Regions and Subregions of Alberta")

# ggsave(file.path(figOutputPath, "natRegSubRegAB.tiff"), plot6,
#        width = 6, height = 7, dpi = 300)

## land-cover classes
if (!is.null(simOutPreSim$rstLCCRTM)) {
  if (!inMemory(simOutPreSim$rstLCCRTM)) {
    lccRas <- sub(".*LandscapesInMotion/", "",
                   filename(simOutPreSim$rstLCCRTM))
    lccRas <- raster(lccRas)
  } else
    lccRas <- simOutPreSim$rstLCCRTM
} else {
  tempDT <- as.data.table(showCache(cachePath(simOutPreSim), "rstLCCRTM"))
  lccRas <- unique(tempDT[createdDate == max(createdDate)]$cacheId) %>%
    reproducible::loadFromCache(cachePath = cachePath(simOutPreSim), cacheId = .)
}

## get original lcc, for colortable
tempLCC <- raster("R/SpaDES/inputs/LCC2005_V1_4a.tif")
colortable(lccRas) <- colortable(tempLCC)
rm(tempLCC)

lccLabels <- c("Temperate/subpolar needle-leaved evergreen, closed canopy",
               "Cold deciduous closed canopy",
               "Mixed needle-leaved evergreen > deciduous closed canopy",
               "Mixed needle-leaved evergreen > deciduous closed young canopy",
               "Mixed cold deciduous > needle-leaved evergreen closed canopy",
               "Temperate/subpolar needle-leaved evergreen medium density, moss-shrub",
               "Temperate/subpolar needle-leaved evergreen medium density, lichen-shrub",
               "Temperate/subpolar needle-leaved evergreen low density, moss-shrub",
               "Temperate/subpolar needle-leaved evergreen low density, lichen-rock",
               "Temperate/subpolar needle-leaved evergreen low density, poorly drained",
               "Cold deciduous broadleaved, low-medium density",
               "Cold deciduous broadleaved, medium density, young regenerating",
               "Mixed needle-leaved evergreen > deciduous, low-medium density",
               "Mixed needle-leaved deciduous > evergreen, low-medium density",
               "Low regenerating young mixed cover",
               "High-low shrub dominated",
               "Grassland",
               "Herb-shrub-bare cover",
               "Wetlands",
               "Sparse needle-leaved evergreen, herb-shrub cover",
               "Polar grassland, herb-shrub",
               "Shrub-herb-lichen-bare",
               "Herb-shrub poorly drained",
               "Lichen-shrub-herb-bare soil",
               "Low vegetation cover",
               "Cropland-woodland",
               "High biomass cropland",
               "Medium biomass cropland",
               "Low biomass cropland",
               "Lichen barren",
               "Lichen-sedge-moss-low shrub wetland",
               "Lichen-spruce bog",
               "Rock outcrops",
               "Recent burns",
               "Old burns",
               "Urban and built-up",
               "Water bodies",
               "Mixes of water and land",
               "Snow/ice")
lccLabels <- paste(as.character(seq_along(lccLabels)),
                   lccLabels, sep = " - ")
names(lccLabels) <- as.character(seq_along(lccLabels))

## colour vector
lccColours <- colortable(lccRas)
names(lccColours) <- as.character(0:255)

## plot using DT
rasData <- data.table(coordinates(lccRas), vals = as.factor(getValues(lccRas)))
lccPlot <- ggplot() +
  layer_spatial(data = simOutPreSim$studyArea, col = "black",
                fill = "black", size = 1.5) +
  geom_raster(data = rasData, mapping = aes(x, y, fill = vals)) +
  annotation_north_arrow(style = north_arrow_minimal,
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         location = "tr", which_north = "true") +
  scale_fill_manual(values = lccColours, na.value = "transparent",
                    labels = lccLabels, breaks = as.integer(names(lccLabels))) +
  theme_pubr(legend = "right") +
  theme(plot.margin = unit(c(0,0,0,0), units = "mm"),
        legend.key = element_rect(colour = "black")) +
  labs(x = "longitude", y = "latitude", fill = "",
       title = "Land cover") +
  guides(fill = guide_legend(ncol = 1))

plot6 <- ggarrange(ecoDistPlot, lccPlot, ncol = 2, nrow = 1,
                   widths = c(0.6,1), labels = "auto",
                   font.label = list(size = 20))
ggsave(file.path(figOutputPath, "Ecologicalzones_LCC.tiff"), plot6,
       width = 20, height = 10, dpi = 300)


## STUDY AREA IN ALBERTA ---------------------------
plot7 <- ggplot() +
  geom_sf(data = alberta) +
  layer_spatial(data = simOutPreSim$studyArea, fill = "forestgreen", colour = "black") +
  annotation_north_arrow(style = north_arrow_minimal,
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         location = "tr", which_north = "true") +
  theme_pubr(legend = "right") +
  theme(plot.background = element_rect(fill = "transparent")) +
  labs(x = "longitude", y = "latitude",
       title = "Study area")

ggsave(file.path(figOutputPath, "LIM_studyArea.tiff"), plot7,
       width = 4, height = 7, dpi = 300)


## JULY MONTHLY DROUGHT CODE ---------------------------
plot8 <- ggplot() +
  layer_spatial(simOutPreSim$weatherDataMDCStk[[1]], aes(fill = stat(band1))) +
  layer_spatial(data = simOutPreSim$studyArea, fill = "transparent", colour = "black") +
  annotation_north_arrow(style = north_arrow_minimal,
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         location = "tr", which_north = "true") +
  scale_fill_distiller(palette = "YlOrRd",  direction = 1, na.value = "transparent") +
  theme_pubr(legend = "right") +
  theme(plot.background = element_rect(fill = "transparent")) +
  labs(x = "longitude", y = "latitude", fill = "",
       title = "Drought code")

ggsave(file.path(figOutputPath, "julMDCyr1.tiff"), plot8,
       width = 4, height = 7, dpi = 300)


## GIFS ------------------------------------------------------
## make GIFs of vegetation maps
## individual pics
outputs_PM <- list.dirs("R/SpaDES/outputs", full.names = TRUE, recursive = FALSE) %>%
  grep("/PM_rep1", ., value = TRUE)
outputs_noPM <- list.dirs("R/SpaDES/outputs", full.names = TRUE, recursive = FALSE) %>%
  grep("/noPM_rep1", ., value = TRUE)

vegTypeMapFiles <- c(sapply(outputs_PM, FUN = function(x) list.files(x, pattern = "vegTypeMap", full.names = TRUE)),
                     sapply(outputs_noPM, FUN = function(x) list.files(x, pattern = "vegTypeMap", full.names = TRUE))) %>%
  unique(.)

vegTypeMapStk_noPM <- lapply(grep("noPM", vegTypeMapFiles, value = TRUE), readRDS) %>%
  stack(.)
vegTypeMapStk_PM <- lapply(grep("noPM", vegTypeMapFiles, value = TRUE, invert = TRUE), readRDS) %>%
  stack(.)
names(vegTypeMapStk_noPM) <- sub(".*year", "year", sub("\\.rds", "", grep("noPM", vegTypeMapFiles, value = TRUE)))
names(vegTypeMapStk_PM) <- sub(".*year", "year", sub("\\.rds", "", grep("noPM", vegTypeMapFiles, invert = TRUE, value = TRUE)))


speciesLabels <- LandR::equivalentName(value = levels(vegTypeMapStk_noPM[[1]])[[1]]$VALUE,
                                       column = "EN_generic_full", df = simOutPreSim$sppEquiv)
names(speciesLabels) <- unique(levels(vegTypeMapStk_noPM[[1]])[[1]]$ID)

speciesColours <- levels(vegTypeMapStk_noPM[[1]])[[1]]$colors
names(speciesColours) <- levels(vegTypeMapStk_noPM[[1]])[[1]]$ID

foothillsMask <- simOutPreSim$rasterToMatch
foothillsMask[!is.na(foothillsMask[])] <- 1
foothillsMaskDF <- as.data.frame(as(foothillsMask, "SpatialPixelsDataFrame"))
names(foothillsMaskDF) <- c("value", "x", "y")

makePNGs <- function(id, rasterStack, labPrefix, labSufix, filePrefix, gif.dir) {
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
    labs(title = paste(labPrefix, sub("year", "year ", names(rasterStack)[id]), labSufix, sep = "")) +
    coord_equal()
  ggsave(file.path(gif.dir, paste0(filePrefix, id, ".png")),
         device = "png", width=5, height=10, dpi = 300, units = "in")
}

map_df(.x = 1:nlayers(vegTypeMapStk_noPM), .f = makePNGs,
       rasterStack = vegTypeMapStk_noPM,
       labPrefix = "stand replacement\n",
       labSufix = "",
       filePrefix = "vegTypeMapStk_noPM",
       gif.dir = "R/SpaDES/outputs/noPM_rep1/gif")
map_df(.x = 1:nlayers(vegTypeMapStk_PM), .f = makePNGs,
       rasterStack = vegTypeMapStk_PM,
       labPrefix = "partial mortality\n",
       labSufix = "",
       filePrefix = "vegTypeMapStk_PM",
       gif.dir = "R/SpaDES/outputs/PM_rep1/gif")

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

makeGIF(gif.dir = "R/SpaDES/outputs/noPM_rep1/gif",
        gifPrefix = "vegTypeMapStk_noPM",
        fps = 2)
makeGIF(gif.dir = "R/SpaDES/outputs/PM_rep1/gif",
        gifPrefix = "vegTypeMapStk_PM",
        fps = 2)

