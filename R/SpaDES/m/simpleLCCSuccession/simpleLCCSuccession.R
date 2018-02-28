
# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "simpleLCCSuccession",
  description = "VERY simple vegetation succession model based on LCC 2010 classes", #"insert module description here",
  keywords = c("succession", "vegetation", "LCC2010"), # c("insert key words here"),
  authors = person("Ceres", "Barros", email = "cbarros@mail.ubc.ca", role = c("aut", "cre")),
  childModules = character(0),
  version = list(SpaDES.core = "0.1.1.9005", simpleLCCSuccession = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "simpleLCCSuccession.Rmd"),
  reqdPkgs = list("raster"),
  parameters = rbind(
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between plot events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "vegetation", objectClass = "list",
                 desc = "List of vegetation states rasters, with simplified classes"),
    expectsInput(objectName = "transitionMatrix", objectClass = "matrix",
                 desc = "Habitat X habitat matrix of transition probabilities (rows = initial, cols = final)", sourceURL = NA),
    expectsInput(objectName = "spreadRas", objectClass = "list",
                 desc = "List of rasters of fire spread")
  ),
  outputObjects = bind_rows(
    expectsInput(objectName = "vegetation", objectClass = "list",
                 desc = "List of vegetation states rasters, with simplified classes")
  )
))


doEvent.simpleLCCSuccession = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      # schedule future event(s)
      sim <- scheduleEvent(sim, start(sim) + 1, "simpleLCCSuccession", "succession", eventPriority = 2)
    },
    succession = {
      ## convert vegetation states
      sim <- do.VegetationTransit(sim)
      
      ## schedule future events
      sim <- scheduleEvent(sim, time(sim) + 1, "simpleLCCSuccession", "succession", eventPriority = 2)
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

## Initialisation
successionInit <- function(sim) {
  return(invisible(sim))
}

## Change vegetation states (present year) in function of fire
## to really useCache, I will need to cache the fire rasters - these have a random component now.
do.VegetationTransit <- function(sim) {
  if(G(sim)$.useCache) {
    sim$vegetation[[time(sim)]] <- reproducible::Cache(cacheRepo = paths(sim)$cachePath, 
                                                       FUN = calc, 
                                                       x = stack(sim$vegetation[[time(sim) - 1]], sim$spreadRas[[time(sim) - 1]]),
                                                       fun = vegTransition)
  } else {
    ## note: transitionMatrix is being used by vegTransition 
    sim$vegetation[[time(sim)]] <- calc(stack(sim$vegetation[[time(sim) - 1]], sim$spreadRas[[time(sim) - 1]]),
                                        fun = vegTransition)
  }
  
  return(invisible(sim))
}



## OTHER INPUTS AND FUNCTIONS --------------------------------
.inputObjects = function(sim) {
  if(is.null(sim$transitionMatrix)){
    sim$transitionMatrix <- matrix(c(c(1, 0, 0, 0, 0, 0),
                                     rep(c(0, 1, 0, 0, 0, 0), 2),
                                     c(0, 0, 0, 1, 0, 0),
                                     c(0, 0, 1, 0, 0, 0),
                                     c(0, 0, 0, 0, 1, 0)),
                                     nrow = 6, ncol = 6, byrow = TRUE,
                                     dimnames = list(paste0("hab", as.character(0:5)), paste0("hab", as.character(0:5))))
  }
  return(invisible(sim))
}