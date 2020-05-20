## -------------------------------------------------------------
## HELPER FUNCTION TO PLOT RASTER LAYERS AGAINST A STUDY AREA
## this function plots a raster layer against a study area
## shapefile and is intended to be iterated across many layers.
##
## ras is the raster layer with values to show,
##  should already be masked to studyArea
## studyArea is a polygon of the study area
## limits is a vector of length 2 that forces the limits of the
##   labelling so that the legend scale is uniform across
##   plots when iterating the function through many layers.
##   If ommitted it will take the min and max values of the current ras

plotFunction <- function(ras, studyArea, limits = NULL) {
  if (is.null(limits))
    limits <- range(getValues(ras), na.rm = TRUE)

  ggplot() +
    layer_spatial(ras, aes(fill = stat(band1))) +
    layer_spatial(data = studyArea, fill = "transparent", colour = "black") +
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
