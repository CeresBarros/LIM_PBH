## ----------------------------------------
## PLOTS FOR REPORT
## ----------------------------------------
library(qs)
library(ggplot2)
library(ggpubr)
library(ggspatial)
library(SpaDES)
library(raster)

## path to figure folder
figOutputPath <- "C:/Users/Ceres Barros/Google Drive/Shared/Landscapes In Motion/ModellingTeam/reportFigs/"

## run name(s)
runName <- c("noPM_newSppParams_fullSA", "PM_newSppParams_fullSA")

## load simLists
simOutPreSim <- qread(file.path("R/SpaDES/outputs/preSim/preSimList.qs"))
if (length(runName) > 1) {
  for(run in runName) {
   eval(parse(text =
    paste0("simOut_", run, " <- qread('R/SpaDES/outputs/",
         run, "/", paste0("simList_", run, ".qs"), "')")
   ))
  }
  rm(run)
}

## FIRE LOCATIONS -------------------------
plot1 <- ggplot() +
  layer_spatial(data = simOutPreSim$studyArea, fill = "grey") +
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

plot2 <- ggplot() +
  layer_spatial(data = simOutPreSim$studyArea, fill = "grey") +
  layer_spatial(rasFuels, aes(fill = stat(band1))) +
  annotation_north_arrow(style = north_arrow_minimal,
                         location = "tr", which_north = "true") +
  theme_pubr(legend = "bottom", margin = FALSE) +
  scale_fill_distiller(palette = "Paired", breaks = sort(unique(rasFuels[])),
                       na.value = "transparent", guide = "legend",
                      labels = FTlabs) +
  labs(x = "longitude", y = "latitude", fill = "Fuel type")
ggsave(file.path(figOutputPath, "fuelsMap.tiff"), plot2,
       width = 4, height = 7, dpi = 300)


## FUEL TYPE COVER
plotFunction <- function(ras, limits) {
  ggplot() +
    layer_spatial(ras, aes(fill = stat(band1))) +
    layer_spatial(data = simOutPreSim$studyArea, fill = "transparent", colour = "black") +
    annotation_north_arrow(style = north_arrow_minimal,
                           height = unit(1, "cm"), width = unit(1, "cm"),
                           location = "tr", which_north = "true") +
    theme_pubr(legend = "bottom") +
    theme(plot.margin = unit(c(0,0,0,0), units = "mm")) +
    scale_fill_distiller(palette = "Greys", na.value = "transparent",
                         direction = 1,
                         breaks = seq(limits[1], limits[2], length.out = 6),
                         limits = limits) +
    labs(x = "longitude", y = "latitude", fill = "Cover",
         title = sub("\\.|_", " ", names(ras)))
}
plotList <- lapply(unstack(simOutPreSim$fuelTypesCoverStk), plotFunction,
                   limits = c(0,1))
plot3 <- ggarrange(plotlist = plotList, ncol = 3, nrow = 3,
          legend = "bottom", common.legend = TRUE)

ggsave(file.path(figOutputPath, "Fueltype_cover.tiff"), plot3,
       width = 9, height = 14, dpi = 300)
rm(ras)


## TOPO, VEG, CLIMATE MAPS ----------------
## elevation
if (!inMemory(simOut_noPM_newSppParams_fullSA$DEMRas)) {
  DEMRas <- sub(".*LandscapesInMotion/", "",
                filename(simOut_noPM_newSppParams_fullSA$DEMRas))
  DEMRas <- raster(DEMRas)
} else
  DEMRas <- simOut_noPM_newSppParams_fullSA$DEMRas

elevPlot <- ggplot() +
  layer_spatial(DEMRas, aes(fill = stat(band1))) +
  layer_spatial(data = simOutPreSim$studyArea, col = "black",
                fill = "transparent") +
  annotation_north_arrow(style = north_arrow_minimal,
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         location = "tr", which_north = "true") +
  theme_pubr(legend = "right") +
  theme(plot.margin = unit(c(0,0,0,0), units = "mm")) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent",
                       direction = -1) +
  labs(x = "longitude", y = "latitude", fill = "m a.s.l.",
       title = "Elevation")

## stand biomass (raw)
if (!inMemory(simOut_noPM_newSppParams_fullSA$rawBiomassMap)) {
  biomassRas <- sub(".*LandscapesInMotion/", "",
                filename(simOut_noPM_newSppParams_fullSA$rawBiomassMap))
  biomassRas <- raster(biomassRas)
} else
  biomassRas <- simOut_noPM_newSppParams_fullSA$rawBiomassMap

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
if (!inMemory(simOut_noPM_newSppParams_fullSA$standAgeMap)) {
  ageRas <- sub(".*LandscapesInMotion/", "",
                    filename(simOut_noPM_newSppParams_fullSA$standAgeMap))
  ageRas <- raster(ageRas)
} else
  ageRas <- simOut_noPM_newSppParams_fullSA$standAgeMap

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
if (!inMemory(simOut_noPM_newSppParams_fullSA$speciesLayers)) {
  sppStk <- unstack(simOut_noPM_newSppParams_fullSA$speciesLayers)
  sppStk <- sub(".*LandscapesInMotion/", "",
                sapply(sppStk, filename))
  sppStk <- stack(sppStk)
  names(sppStk) <- sub(".*\\.", "", sub(".tif", "", names(sppStk)))
} else
  sppStk <- simOut_noPM_newSppParams_fullSA$speciesLayers

names(sppStk) <- LandR::equivalentName(value = names(sppStk), column = "EN_generic_full",
                                       df = simOutPreSim$sppEquiv)

plotList <- lapply(unstack(sppStk), plotFunction, limits = c(0, 100))
plot5 <- ggarrange(plotlist = plotList, ncol = 3, nrow = 2,
                     legend = "bottom", common.legend = TRUE)

ggsave(file.path(figOutputPath, "Species_cover.tiff"), plot5,
       width = 9, height = 10, dpi = 300)


## ECODISTRICTS ---------------------------

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
