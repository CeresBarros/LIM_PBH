
# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fire_STSM",
  description = "Fire state-and-transition [toy] model",
  keywords = c("fire", "STSM", "state-and-transition"),
  authors = person("Ceres", "Barros", email = "cbarros@mail.ubc.ca", role = c("aut", "cre")),
  childModules = character(0),
  version = numeric_version("0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fire_STSM.Rmd"),
  reqdPkgs = list("raster"),
  parameters = rbind(
    defineParameter("fireSize", "numeric", 1000, NA, NA, desc = "Fire size in pixels"),
    defineParameter("noStartPix", "numeric", 100, NA, NA, desc = "Number of fire events"),
    defineParameter(".plotMaps", "logical", FALSE, NA, NA, "This describes whether maps should be plotted or not"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA, "This describes the simulation time interval between plot events")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "studyArea", objectClass = "SpatialPolygonsDataFrame",
                 desc = paste("Needs to be provided by user."),
                 sourceURL = NA),
    expectsInput(objectName = "vegetationRas", objectClass = "RasterLayer",
                 desc = "2010 land classification map in study area, default is canada national land classification in 2010",
                 sourceURL = "http://www.cec.org/sites/default/files/Atlas/Files/Land_Cover_2010/Land_Cover_2010_TIFF.zip"),
    expectsInput(objectName = "fire_transitprobs", objectClass = "matrix",
                 desc = "Habitat X habitat matrix of transition probabilities (rows = initial, cols = final)", sourceURL = NA)
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = "vegetation", objectClass = "list",
                  desc = "List of vegetation states rasters, with simplified classes"),
    createsOutput(objectName = "burnable_areas", objectClass = "RasterLayer",
                  desc = "Raster of areas that are susceptible to burning"),
    createsOutput(objectName = "ROS_map", objectClass = "RasterLayer",
                  desc = "Raster of rate of spread per pixel"),
    createsOutput(objectName = "spreadRas", objectClass = "list",
                  desc = "List of rasters of fire spread"),
    createsOutput(objectName = "severity_ras", objectClass = "RasterLayer",
                  desc = "Raster of fire severity")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.fire_STSM = function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      ## Initialise module
      sim <- fireInit(sim)

      ## schedule future event(s)
      sim <- scheduleEvent(sim, start(sim) + 1, "fire_STSM", "fireSpread", eventPriority = 1)
      sim <- scheduleEvent(sim, start(sim) + 1, "fire_STSM", "vegetationTransition", eventPriority = 2)
      sim <- scheduleEvent(sim, start(sim) + 1, "fire_STSM", "severityMap", eventPriority = 3)
      if(P(sim)$.plotMaps) {
        sim <- scheduleEvent(sim, start(sim) + 1, "fire_STSM", "plot", eventPriority = 4)
      }
    },

    fireSpread = {
      ## calculate fire spread
      sim <- do.Fire(sim)

      ## schedule future events
      sim <- scheduleEvent(sim, time(sim) + 1, "fire_STSM", "fireSpread", eventPriority = 1)
    },
    vegetationTransition = {
      ## convert vegetation states
      sim <- do.VegetationTransit(sim)

      ## schedule future events
      sim <- scheduleEvent(sim, time(sim) + 1, "fire_STSM", "vegetationTransition", eventPriority = 2)
    },
    severityMap = {
      ## calculate severity
      sim <- do.Severity(sim)

      ## schedule future events
      sim <- scheduleEvent(sim, time(sim) + 1, "fire_STSM", "severityMap", eventPriority = 3)
    },
    plot = {
      ## Plot severity and vegetation changes
      # sim <- fire_STSMPlot(sim = sim ,
      #                      severity = sim$severity_ras[[time(sim)]],
      #                      preVegetation = sim$vegetation[[time(sim) - 1]],
      #                      postVegetation = sim$vegetation[[time(sim)]])

      sim <- fire_STSMPlot(sim)

      ## schedule next plot
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "fire_STSM", "plot", eventPriority = 4)
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

### module initialization - part of this may pass to another module of data prep
fireInit <- function(sim) {
  ## make raster storage lists
  sim$vegetation = sim$severity_ras = sim$spreadRas = list()

  ## PRE-FIRE VEGETATION ---------------------------------------
  ## get first vegetation state from original LCC raster
  if(G(sim)$.useCache) {
    vegetation_prefire <- reproducible::Cache(cacheRepo = paths(sim)$cachePath, makePrefireVegetation,
                                              area = studyArea, vegRas = vegetationRas)
  } else {
    vegetation_prefire <- makePrefireVegetation(area = studyArea, vegRas = vegetationRas)
  }



  ## VEGETATION CLASSES ----------------------------------------
  ## Non-burnable:
  ## wetlands (19), wet tundra (23) cropland/woodlands (26, 27,28,29), lichen dominated (30, 31,32)
  ## rock outcrops (33), recent burns (34), cities (36), water (37, 38), snow/ice (39)
  non_burn <- c(19, 26, 23, 26, 27, 28, 29, 30, 31, 32, 33, 34, 36, 37, 38, 39)

  ## Grasslands/open habs:
  ## open habs (17, 18, 20, 21, 22, 24, 25)
  grass <- c(17, 18, 20, 21, 22, 24, 25)

  ## Shrublands, recently burnt forest (within last 10y) = old burns (35), shrublands (16)
  shrubs <- c(16, 35)

  ## Deciduous Forest
  decid_forst <- c(2, 11, 12)

  ## Mixed forest
  mixed_forst <- c(3, 4, 5, 13, 14, 15)

  ## Coniferous forest
  conif_forst <- c(1, 6, 7, 8, 9, 10)

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

### Simulate fire spread
do.Fire <- function(sim) {
  ## MAKE BURNABLE AREAS RASTER -------------------------------
  sim$burnable_areas <- sim$vegetation[[time(sim) - 1]]   ## using previous year vegetation map
  sim$burnable_areas[sim$vegetation[[time(sim) - 1]][] == 0] <- NA

  ## MAKE RATE OF SPREAD RASTER -------------------------------
  sim$ROS_map <- sim$vegetation[[time(sim) - 1]]
  reclass_mat2 <- matrix(c(0:5, 0, 0.9, 0.7, 0.2, 0.4, 0.5), byrow = FALSE, nrow = 6, ncol = 2)
  sim$ROS_map <- reclassify(sim$ROS_map, rcl = reclass_mat2)
  rm(reclass_mat2)

  ## MAKE RASTER OF FIRE SPREAD -------------------------------
  ## note that this function two random components: selection of starting pixels and fire spread
  if(G(sim)$.useCache){
    sim$spreadRas[[time(sim) - 1]] <- reproducible::Cache(cacheRepo = paths(sim)$cachePath,
                                         FUN = makeFireSpreadRas,
                                         fMask = sim$burnable_areas, pixels = P(sim)$noStartPix,
                                         fROS = sim$ROS_map, fSize = P(sim)$fireSize)
  } else {
    sim$spreadRas[[time(sim) - 1]] <- makeFireSpreadRas(fMask = sim$burnable_areas, pixels = P(sim)$noStartPix,
                                       fROS = sim$ROS_map, fSize = P(sim)$fireSize)
  }

  return(invisible(sim))
}

### Change vegetation states (present year) in function of fire
## to really useCache, I will need to cache the fire rasters - these have a random component now.
do.VegetationTransit <- function(sim) {
  if(G(sim)$.useCache) {
    sim$vegetation[[time(sim)]] <- reproducible::Cache(cacheRepo = paths(sim)$cachePath,
                                                       FUN = calc,
                                                       x = stack(sim$vegetation[[time(sim) - 1]], sim$spreadRas[[time(sim) - 1]]),
                                                       fun = vegTransition)
  } else {
  sim$vegetation[[time(sim)]] <- calc(stack(sim$vegetation[[time(sim) - 1]], sim$spreadRas[[time(sim) - 1]]), fun = vegTransition)
  }

  return(invisible(sim))
}

### Calculate severity based on vegetation state transitions
do.Severity <- function(sim){
  ## convert fire spread raster to a mask
  fireMask <- sim$spreadRas[[time(sim) - 1]]
  fireMask[!is.na(fireMask)] <- 1

  if(G(sim)$.useCache) {
    sim$severity_ras[[time(sim)]] <- reproducible::Cache(cacheRepo = paths(sim)$cachePath,
                                                         FUN = calc,
                                                         x = stack(sim$vegetation[[time(sim) - 1]], sim$vegetation[[time(sim)]], fireMask),
                                                         fun = calculateSeverity)
  } else {
    sim$severity_ras[[time(sim)]] <- calc(stack(sim$vegetation[[time(sim) - 1]], sim$vegetation[[time(sim)]], fireMask), fun = calculateSeverity)
  }

  return(invisible(sim))
}

### Plot fire severity and vegetation
fire_STSMPlot <- function(sim) {
  Plot(sim$severity_ras[[time(sim)]], new = TRUE,
       title = "Fire Severity Map",
       cols = c("grey", "green", "yellow", "red"))

  Plot(sim$vegetation[[time(sim) - 1]], new = TRUE,
       title = "Vegetation pre-fire")

  Plot(sim$vegetation[[time(sim)]], new = TRUE,
       title = "Vegetation post-fire")

  return(invisible(sim))
}


## OTHER INPUTS AND FUNCTIONS --------------------------------
.inputObjects = function(sim) {
  sim$fire_transitprobs <- matrix(c(c(1, 0, 0, 0, 0, 0),
                         rep(c(0, 1, 0, 0, 0, 0), 2),
                         rep(c(0, 0, 0, 1, 0, 0), 2),
                         c(0, 0, 0, 0, 1, 0)),
                       nrow = 6, ncol = 6, byrow = TRUE,
                       dimnames = list(paste0("hab", as.character(0:5)), paste0("hab", as.character(0:5))))
  return(invisible(sim))
}

