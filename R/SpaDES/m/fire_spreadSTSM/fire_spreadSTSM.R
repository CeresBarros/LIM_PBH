
# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fire_spreadSTSM",
  description = "Fire state-transition [toy] model using SpaDES::spread()",
  keywords = c("fire", "state-transition model"),
  authors = person("Ceres", "Barros", email = "cbarros@mail.ubc.ca", role = c("aut", "cre")),
  childModules = character(0),
  version = numeric_version("0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fire_spreadSTSM.Rmd"),
  reqdPkgs = list("raster"),
  parameters = rbind(
    defineParameter("fireSize", "numeric", 1000, NA, NA, desc = "Fire size in pixels"),
    defineParameter("noStartPix", "numeric", 100, NA, NA, desc = "Number of fire events"),
    defineParameter(".plotMaps", "logical", FALSE, NA, NA, "This describes whether maps should be plotted or not"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA, "This describes the simulation time interval between plot events")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "vegetation", objectClass = "list",
                 desc = "List of vegetation states rasters, with simplified classes"),
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
                  desc = "List of rasters of fire spread")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.fire_spreadSTSM = function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      ## Initialise module
      sim <- fireInit(sim)
      
      ## schedule future event(s)
      sim <- scheduleEvent(sim, start(sim) + 1, "fire_spreadSTSM", "fireSpread", eventPriority = 1)
      sim <- scheduleEvent(sim, start(sim) + 1, "fire_spreadSTSM", "transitions", eventPriority = 2)
    },
    
    fireSpread = {
      ## calculate fire spread
      sim <- do.Fire(sim)
      
      ## schedule future events
      sim <- scheduleEvent(sim, time(sim) + 1, "fire_spreadSTSM", "fireSpread", eventPriority = 1)
    },
    transitions = {
      ## convert vegetation states
      sim <- do.VegetationTransit(sim)
      
      ## schedule future events
      sim <- scheduleEvent(sim, time(sim) + 1, "fire_spreadSTSM", "transitions", eventPriority = 2)
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

### module initialization
fireInit <- function(sim) {
  ## make raster storage lists 
  sim$spreadRas = list()
  
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
  
  startPix <- sample(which(!is.na(sim$burnable_areas[])), P(sim)$noStartPix)
  
  ## calculate spread - potentially iterate spread2 to allow updating parameters (instead of runnin interactions for fixed argument values)
  # sim$spreadRas[[time(sim) - 1]] <- spread2(landscape = sim$burnable_areas, spreadProbRel = sim$ROS_map,
  #                                           start = startPix, maxSize =  P(sim)$fireSize, plot.it = FALSE)
  
  ## Favier's model:
  sim$spreadRas[[time(sim) - 1]] <- spread2(landscape = sim$burnable_areas, spreadProb = 0.5, 
                                            start = startPix, maxSize =  P(sim)$fireSize, plot.it = FALSE)

  
  ## remove fires outside burnable areas
  sim$spreadRas[[time(sim) - 1]][is.na(sim$burnable_areas)] <- NA
  
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

