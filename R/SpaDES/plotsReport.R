## ----------------------------------------
## PLOTS FOR REPORT
## ----------------------------------------

library(qs)
library(ggplot2)
library(ggpubr)
library(ggspatial)
library(SpaDES)

## path to figure folder
figOutputPath <- "C:/Users/Ceres Barros/Google Drive/Shared/Landscapes In Motion/ModellingTeam/reportFigs/"

## run name(s)
runName <- c("noPM_newSppParams_fullSA", "PM_newSppParams_fullSA")

## load simLists
simOutPreSim <- qread(file.path("R/SpaDES/outputs/preSim/preSimList.qs"))
if (length(runName) > 1) {
  for(run in runName) {
   eval(parse = text(
    paste0("simOut_", run, " <- qread(R/SpaDES/outputs/",
         run, "/", paste0("simList_", runName, ".qs"))
   ))
  }
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
plotFunction <- function(ras) {
  ggplot() +
    layer_spatial(data = simOutPreSim$studyArea, fill = "grey") +
    layer_spatial(ras, aes(fill = stat(band1))) +
    annotation_north_arrow(style = north_arrow_minimal,
                           height = unit(1, "cm"), width = unit(1, "cm"),
                           location = "tr", which_north = "true", ) +
    theme_pubr(legend = "bottom") +
    theme(plot.margin = unit(c(0,0,0,0), units = "mm")) +
    scale_fill_distiller(palette = "Greys", na.value = "transparent",
                         direction = 1, breaks = seq(0, 1, 0.2),
                         limits = c(0,1)) +
    labs(x = "longitude", y = "latitude", fill = "Cover",
         title = names(ras))
}
plotList <- lapply(unstack(simOutPreSim$fuelTypesCoverStk), plotFunction)
plot3 <- ggarrange(plotlist = plotList, ncol = 3, nrow = 3,
          legend = "bottom", common.legend = TRUE)

ggsave(file.path(figOutputPath, "Fueltype_cover.tiff"), plot3,
       width = 9, height = 14, dpi = 300)
rm(ras)


## TOPO, VEG, CLIMATE MAPS ----------------
simOutPreSim$

## SPECIES COVER MAPS ---------------------
