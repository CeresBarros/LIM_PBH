
# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "simplifyLCCVeg",
  description = "Simplification of LCC vegetation classes",
  keywords = c("vegetation", "LCC", "reclassification"),
  authors = person("Ceres", "Barros", email = "cbarros@mail.ubc.ca", role = c("aut", "cre")),
  childModules = character(0),
  version = numeric_version("0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "simplifyLCCVeg.Rmd"),
  reqdPkgs = list("raster"),
  parameters = rbind(
    defineParameter(".plotMaps", "logical", FALSE, NA, NA, "This describes whether maps should be plotted or not"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA, "This describes the simulation time interval between plot events")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "studyArea", objectClass = "SpatialPolygonsDataFrame",
                 desc = "Shapefile of study area. Default is a random polygon somewhere in the landcover map",
                 sourceURL = NA),
    expectsInput(objectName = "vegetationRas", objectClass = "RasterLayer",
                 desc = "Land cover map in study area, default is LCC2010",
                 sourceURL = "http://www.cec.org/sites/default/files/Atlas/Files/Land_Cover_2010/Land_Cover_2010_TIFF.zip")
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = "vegetation", objectClass = "list",
                  desc = "List of vegetation states rasters, with simplified classes")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.simplifyLCCVeg = function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      ## Initialise module
      sim <- fireInit(sim)
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

### module initialization - part of this may pass to another module of data prep
fireInit <- function(sim) {
  ## make raster storage lists
  sim$vegetation <- sim$severity_ras <- sim$spreadRas <- list()

  ## VEGETATION CLASSES ----------------------------------------
  ## Non-burnable:
  ## LCC2005 ---
  ## wetlands (19), wet tundra (23) cropland/woodlands (26, 27,28,29), lichen dominated (30, 31,32)
  ## rock outcrops (33), recent burns (34), cities (36), water (37, 38), snow/ice (39)
  # non_burn <- c(19, 26, 23, 26, 27, 28, 29, 30, 31, 32, 33, 34, 36, 37, 38, 39)

  ## Grasslands/open habs:
  ## open habs (17, 18, 20, 21, 22, 24, 25)
  # grass <- c(17, 18, 20, 21, 22, 24, 25)

  ## Shrublands, recently burnt forest (within last 10y) = old burns (35), shrublands (16)
  # shrubs <- c(16, 35)

  ## Deciduous Forest
  # decid_forst <- c(2, 11, 12)

  ## Mixed forest
  # mixed_forst <- c(3, 4, 5, 13, 14, 15)

  ## Coniferous forest
  # conif_forst <- c(1, 6, 7, 8, 9, 10)

  ## LCC2010 --
  ## wetlands (14), sub-polar grassland/shrubland (11,12) cropland (15), lichen dominated (13)
  ## barren lands (16), urban (17), water (18), snow/ice (19)
  non_burn <- c(11, 12, 13, 14, 15, 16, 17, 18, 19)

  ## Grasslands/open habs:
  grass <- c(9, 10)

  ## Shrublands
  shrubs <- c(7,8)

  ## Deciduous Forest
  decid_forst <- c(3, 4, 5)

  ## Mixed forest
  mixed_forst <- c(6)

  ## Coniferous forest
  conif_forst <- c(1, 2)

  ## RECLASSIFY VEGETATION ------------------------------------
  reclass_mat <- as.matrix(data.frame(old = c(non_burn, grass, shrubs, decid_forst, mixed_forst, conif_forst),
                                      new = c(rep(0, length(non_burn)), rep(1, length(grass)), rep(2, length(shrubs)),
                                              rep(3, length(decid_forst)), rep(4, length(mixed_forst)), rep(5, length(conif_forst)))))
  vegetation_prefire <- reclassify(vegetation_prefire, rcl = reclass_mat)

  sim$vegetation[[start(sim)]] <- vegetation_prefire

  ## clean workspace
  rm(non_burn, grass, shrubs, decid_forst, mixed_forst, conif_forst, reclass_mat, vegetation_prefire)

  return(invisible(sim))
}



## Inputs
.inputObjects <- function(sim) {
  ## PRE-FIRE VEGETATION ---------------------------------------
  ## get first vegetation state from original LCC raster
  if(is.null(sim$vegetation_prefire)){
    sim$vegetationRas <- prepInputs(targetFile = "NA_LandCover_2010_25haMMU.tif",
                                    url = "http://www.cec.org/sites/default/files/Atlas/Files/Land_Cover_2010/Land_Cover_2010_TIFF.zip",
                                    useCache = TRUE, cacheRepo = cachePath(sim))

    if(G(sim)$.useCache) {
      sim$vegetation_prefire <- reproducible::Cache(makePrefireVegetation,
                                                    area = sim$studyArea, vegRas = sim$vegetationRas,
                                                    cacheRepo = paths(sim)$cachePath)
    } else {
      sim$vegetation_prefire <- makePrefireVegetation(area = sim$studyArea, vegRas = sim$vegetationRas)
    }
  }


  if(is.null(sim$studyArea)) {
    sim$studyArea <- randomPolygon(SpatialPoints(cbind(-110, 59)), 1e4)
    sim$studyArea <- sp::spTransform(x = studyArea, CRSobj = crs(sim$vegetationRas))
  }

  return(invisible(sim))
}
