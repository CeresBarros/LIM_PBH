
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
    defineParameter(".plotMaps", "logical", FALSE, NA, NA, "This describes whether maps should be plotted or not"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA, "This describes the simulation time interval between plot events")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "vegetation", objectClass = "list",
                 desc = "List of vegetation states rasters, with simplified classes"),
    expectsInput(objectName = "spreadRas", objectClass = "list",
                 desc = "List of rasters of fire spread")
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = "severity_ras", objectClass = "RasterLayer",
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
      
      ## schedule future event(s)
      sim <- scheduleEvent(sim, start(sim) + 1, "fireSeverity", "severityMap", eventPriority = 3)
      if(P(sim)$.plotMaps) {
        sim <- scheduleEvent(sim, start(sim) + 1, "fireSeverity", "plot", eventPriority = 4)
      }
    },
    severityMap = {
      ## calculate severity
      sim <- do.Severity(sim)
      
      ## schedule future events
      sim <- scheduleEvent(sim, time(sim) + 1, "fireSeverity", "severityMap", eventPriority = 3)
    },
    plot = {
      ## Plot severity and vegetation changes
      sim <- fire_STSMPlot(sim)
      
      ## schedule next plot
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "fireSeverity", "plot", eventPriority = 4)
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

### module initialization - part of this may pass to another module of data prep
fireInit <- function(sim) {
  sim$severity_ras <- list()
  
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


