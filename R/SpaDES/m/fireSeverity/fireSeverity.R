# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireSeverity",
  description = "Calculate fire severity from the type of vegetation transition occuring after a fire",
  keywords = c("fire", "severity", "state-transition", "vegetation"),
  authors = person("Ceres", "Barros", email = "cbarros@mail.ubc.ca", role = c("aut", "cre")),
  childModules = character(0),
  version = numeric_version("0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSeverity.Rmd"),
  reqdPkgs = list("raster"),
  parameters = rbind(
    defineParameter(name = "fireTimestep", class = "numeric", default = 2L,
                    desc = "The number of time units between successive fire events in a fire module"),
    defineParameter(".plotMaps", "logical", FALSE, NA, NA, "Controls whether maps should be plotted or not"),
    defineParameter(name = ".saveInitialTime", class = "numeric", default = 0,
                    min = NA, max = NA, desc = "This describes the simulation time at which the
                    first save event should occur"),
    defineParameter(".useCache", "logical", "init", NA, NA,
                    desc = "use caching for the spinup simulation?")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "ecoregionMap", objectClass = "RasterLayer",
                 desc = "ecoregion map that has mapcodes match ecoregion table and speciesEcoregion table",
                 sourceURL = ""),
    expectsInput(objectName = "simulatedBiomassMap", objectClass = "RasterLayer",
                 desc = "Biomass map at each succession time step"),
    expectsInput(objectName = "rstCurrentBurn", objectClass = "RasterLayer",
                 desc = "Binary raster of fire spread"),
    expectsInput(objectName = "fireYear", objectClass = "numeric", desc =  "Next fire year"),
    expectsInput("studyArea", "SpatialPolygonsDataFrame",
                 desc = paste("Polygon to use as the study area.",
                              "Defaults to  an area in Southwestern Alberta, Canada."),
                 sourceURL = ""),
    expectsInput(objectName = "studyAreaFBP", objectClass = "SpatialPolygonsDataFrame",
                 desc = "same as studyArea,  but on FBP-compatible projection", sourceURL = "")
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = "biomassMapPreFire", objectClass = "RasterLayer",
                  desc = "Biomass map from before the last fire."),
    createsOutput(objectName = "biomassMapPostFire", objectClass = "RasterLayer",
                  desc = "Biomass map from directly after the last fire."),
    createsOutput(objectName = "severityMap", objectClass = "RasterLayer",
                  desc = "Raster of fire severity")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.fireSeverity = function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      ## Initialise module
      sim <- fireInit(sim)

      stop("add summaryBGM to severity module")
      if (!is.null(sim$rstCurrentBurn)) {
        ## schedule post-fire (pre-growth/mortality/dispersal) summaryBGM
        ## if fires are simulated
        fireTimestep <- grep("fireTimestep", names(unlist(sim@params, recursive = TRUE))) %>%
          unlist(sim@params, recursive = TRUE)[.] %>%
          unique(.) %>%
          as.numeric(.)

        if (length(fireTimestep) > 1)
          stop(paste("Different values for fireTimestep parameter provided.",
                     "Make sure all modules have the same fireTimestep value", sep = "\n"))

        sim <- scheduleEvent(sim, time(sim) + fireTimestep,
                             "LBMR", "summaryBGM", eventPriority = 4)
      }

      ## schedule events
      if(!is.null(sim$rstCurrentBurn)) {   ## only if fire module is "active"
        sim <- scheduleEvent(sim, params(sim)$fireSpread$fireInitialTime, "fireSeverity", "getBiomassPreFire", eventPriority = 2)  ## before regeneration
        sim <- scheduleEvent(sim, params(sim)$fireSpread$fireInitialTime, "fireSeverity", "calcSeverity", eventPriority = 4.25)

        if(P(sim)$.plotMaps)
          sim <- scheduleEvent(sim, params(sim)$fireSpread$fireInitialTime, "fireSeverity", "severityPlot", eventPriority = 8.25)

        if(!any(is.na(P(sim)$.saveInitialTime)))
          sim <- scheduleEvent(sim, params(sim)$fireSpread$fireInitialTime,
                               "fireSeverity", "saveSeverity", eventPriority = 4.50)
      }
    },
    getBiomassPreFire = {
      ## get the pre-fire biomass maps before any fire event
      sim <- doBiomassPreFire(sim)

      ## schedule future events - can't use fireyear because it's not updated at this point
      sim <- scheduleEvent(sim, eventTime = time(sim) + P(sim)$fireTimestep, moduleName = "fireSeverity",
                           eventType = "getBiomassPreFire", eventPriority = 2)
    },
    calcSeverity = {
      if(!all(is.na(sim$rstCurrentBurn[]))) {
        ## calculate severity
        sim <- doSeverity(sim)

        ## schedule future events
        sim <- scheduleEvent(sim, eventTime = sim$fireYear, moduleName = "fireSeverity",
                             eventType = "calcSeverity", eventPriority = 4.25)
      }
    },
    severityPlot = {
      if(!all(is.na(sim$rstCurrentBurn[]))) {
        ## Plot severity and vegetation changes
        sim <- doSeverityPlot(sim)

        ## schedule next plot
        sim <- scheduleEvent(sim, sim$fireYear, "fireSeverity", "severityPlot", eventPriority = 8.25)
      }
    },
    saveSeverity = {
      if(!all(is.na(sim$rstCurrentBurn[]))) {
        ## Plot severity and vegetation changes
        sim <- doSaveSeverity(sim)

        ## schedule next plot
        sim <- scheduleEvent(sim, time(sim) + P(sim)$fireTimestep, "fireSeverity", "saveSeverity", eventPriority = 4.50)
      }
    },

    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

### module initialization - part of this may pass to another module of data prep
fireInit <- function(sim) {
  ## store initial biomassMap as first pre-fire biomass raster
  sim$biomassMapPreFire <- sim$simulatedBiomassMap

  return(invisible(sim))
}

### Calculate severity based on vegetation state transitions
doSeverity <- function(sim) {
  ## if scheduling is correct simulatedBiomassMap should have been updated after fire and before growth/mortality/dispersal
  ## the units for sumB, sumANPP, sumMortality are g/m2, g/m2/year, g/m2/year, respectively.
  biomassMapPostFire <- sim$simulatedBiomassMap

  ## fire severity in changes in biomass
  ## convert fire spread raster to a mask
  fireMask <- sim$rstCurrentBurn
  fireMask[!is.na(fireMask)] <- 1

  severityMap <- ((biomassMapPostFire - sim$biomassMapPreFire)/abs(sim$biomassMapPreFire)) * 100
  severityMap <- setValues(severityMap, values = round(getValues(severityMap)))
  severityMap <- raster::mask(severityMap, mask = fireMask)

  ## export to sim
  sim$biomassMapPostFire <- biomassMapPostFire
  sim$severityMap <- severityMap

  return(invisible(sim))
}

### Plot fire severity and vegetation
doSeverityPlot <- function(sim) {
  Plot(sim$severityMap, new = TRUE,
       title = "Fire severity",
       cols = heat.colors(10))
  # Plot(sim$severityMap, arr = c(2,2), new = TRUE)  ## doesn't re-arrange


  return(invisible(sim))
}

doSaveSeverity <- function(sim) {
  raster::projection(sim$severityMap) <- raster::projection(sim$ecoregionMap)
  writeRaster(sim$severityMap,
              file.path(outputPath(sim), paste("severityMap_Year", round(time(sim)), ".tif",sep="")),
              datatype='INT2S', overwrite = TRUE)
  return(invisible(sim))
}

doBiomassPreFire <- function(sim){
  ## before regeneration occurs, save the pre-fire biomass conditions
  sim$biomassMapPreFire <- sim$simulatedBiomassMap
  return(invisible(sim))
}

## OTHER INPUTS AND FUNCTIONS --------------------------------
.inputObjects = function(sim) {
  dPath <- dataPath(sim)
  cacheTags <- c(currentModule(sim), "function:.inputObjects")

  ## project to Lat/Long (decimal degrees) for compatibility with FBP system
  ## TODO: this results in data loss - but LandR doesn't deal well with lat/long
  ## need to find long term solution
  latLong <- "+proj=longlat +datum=WGS84"

  if (!suppliedElsewhere("studyAreaFBP", sim)) {
    if (!suppliedElsewhere("studyArea", sim)) {
      if (getOption("LandR.verbose", TRUE) > 0)
        message("'studyArea' was not provided by user. Using a polygon (6250000 m^2) in southwestern Alberta, Canada")
      sim$studyArea <- randomStudyArea(seed = 1234, size = (250^2)*100)
    }
    sim$studyAreaFBP <- sim$studyArea
  }

  ## if necessary reproject to lat/long - for compatibility with FBP
  if (!identical(latLong, crs(sim$studyAreaFBP))) {
    sim$studyAreaFBP <- spTransform(sim$studyAreaFBP, latLong) #faster without Cache
  }

  ## load ecoregion map
  ## TODO: make equal to LBMR.
  if (!suppliedElsewhere("ecoregionMap", sim )) {
    ## LANDIS-II demo data:
    ## Don't load ecoregions to sim even if they don't exist
    if (!suppliedElsewhere("ecoregion", sim)) {
      ecoregion <- prepInputsEcoregion(url = extractURL("ecoregion"),
                                       dPath = dPath, cacheTags = cacheTags)
    } else {
      ecoregion <- sim$ecoregion
    }

    ## Dummy version with spatial location in Canada
    ras <- projectExtent(sim$studyAreaFBP, crs = crs(sim$studyAreaFBP))
    res(ras) = 250
    ecoregionMap <- rasterize(sim$studyAreaFBP, ras)
    ecoregionMap[!is.na(getValues(ecoregionMap))][] <- sample(ecoregion$mapcode,
                                                              size = sum(!is.na(getValues(ecoregionMap))),
                                                              replace = TRUE)
    sim$ecoregionMap <- ecoregionMap
  }

  return(invisible(sim))
}
