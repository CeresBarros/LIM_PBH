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
    defineParameter(".plotMaps", "logical", FALSE, NA, NA, "This describes whether maps should be plotted or not")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "biomassMap", objectClass = "RasterLayer",
                 desc = "Biomass map at each succession time step"),
    expectsInput(objectName = "biomassMapPreFire", objectClass = "RasterLayer",
                 desc = "Biomass map from before the last fire."),
    expectsInput(objectName = "fireSpreadRas", objectClass = "list",
                 desc = "List of rasters of fire spread"),
    expectsInput(objectName = "fireYear", objectClass = "list",
                 desc =  "Next fire year")
  ),
  outputObjects = bind_rows(
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
      
      ## schedule events
      if(!is.null(sim$fireSpreadRas)) {   ## only if fire module is "active"
        sim <- scheduleEvent(sim, params(sim)$fireSpread$fireStart, "fireSeverity", "calcSeverity", eventPriority = 8)
        if(P(sim)$.plotMaps) {
          sim <- scheduleEvent(sim, params(sim)$fireSpread$fireStart, "fireSeverity", "severityPlot", eventPriority = 8.5)
        }
        sim <- scheduleEvent(sim, params(sim)$fireSpread$fireStart,
                             "fireSeverity", "saveSeverity", eventPriority = 8.75)
      }
    },
    calcSeverity = {
      if(!all(is.na(sim$fireSpreadRas[[time(sim)]][]))) {
        ## calculate severity
        sim <- doSeverity(sim)
        
        ## schedule future events
        sim <- scheduleEvent(sim, sim$fireYear, "fireSeverity", "calcSeverity", eventPriority = 8)
      }
    },
    severityPlot = {
      if(!all(is.na(sim$fireSpreadRas[[time(sim)]][]))) {
        ## Plot severity and vegetation changes
        sim <- doSeverityPlot(sim)
        
        ## schedule next plot
        sim <- scheduleEvent(sim, sim$fireYear, "fireSeverity", "severityPlot", eventPriority = 8.5)
      }
    },
    saveSeverity = {
      if(!all(is.na(sim$fireSpreadRas[[time(sim)]][]))) {
        ## Plot severity and vegetation changes
        sim <- doSaveSeverity(sim)
        
        ## schedule next plot
        sim <- scheduleEvent(sim, sim$fireYear, "fireSeverity", "saveSeverity", eventPriority = 8.75)
      }
    },
    
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

### module initialization - part of this may pass to another module of data prep
fireInit <- function(sim) {
  return(invisible(sim))
}

### Calculate severity based on vegetation state transitions
doSeverity <- function(sim){
  ## convert fire spread raster to a mask
  fireMask <- sim$fireSpreadRas[[time(sim)]]
  fireMask[!is.na(fireMask)] <- 1
  
  ## fire severity in % mortality ((pre - post)/pre)
  sim$severityMap <- (sim$biomassMapPreFire - sim$biomassMap)/sim$biomassMapPreFire
  sim$severityMap <- setValues(sim$severityMap, values = round(getValues(sim$severityMap)))
  sim$severityMap <- raster::mask(sim$severityMap, mask = fireMask)
  
  return(invisible(sim))
}

### Plot fire severity and vegetation
doSeverityPlot <- function(sim) {  
  browser()
  Plot(sim$severityMap, new = TRUE,
       title = "Fire severity",
       cols = heat.colors(10))
  
  return(invisible(sim))
}

doSaveSeverity = function(sim) {
  raster::projection(sim$severityMap) <- raster::projection(sim$ecoregionMap)
  writeRaster(sim$severityMap,
              file.path(outputPath(sim), paste("severityMap_Year", round(time(sim)), ".tif",sep="")), datatype='INT2S',
              overwrite = TRUE)
  return(invisible(sim))
}

## OTHER INPUTS AND FUNCTIONS --------------------------------
.inputObjects = function(sim) {
  
  dPath <- dataPath(sim)
  cacheTags = c(currentModule(sim), "function:.inputObjects")
  
  if (!suppliedElsewhere("ecoregionMap", sim )) {
    # load ecoregion map
    sim$ecoregionMap <- raster(file.path(dPath, "ecoregions.gis"))
  }
  
  # if(!suppliedElsewhere("biomassMap", sim)) {
  #   sim$biomassMap <- Cache(prepInputs,
  #                           targetFile = biomassMapFilename,
  #                           archive = asPath(c("kNN-StructureBiomass.tar",
  #                                              "NFI_MODIS250m_kNN_Structure_Biomass_TotalLiveAboveGround_v0.zip")),
  #                           url = extractURL("biomassMap", sim), 
  #                           destinationPath = dPath,
  #                           studyArea = sim$shpStudySubRegion,
  #                           useSAcrs = TRUE,
  #                           method = "bilinear",
  #                           datatype = "INT2U",
  #                           filename2 = TRUE,
  #                           userTags = cacheTags)
  # } 
  
  return(invisible(sim))
}