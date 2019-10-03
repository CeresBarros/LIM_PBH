
# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireSpread",
  description = "Fire spread model using Favier (2004) percolation model, via SpaDES.tools::spread2() ",
  keywords = c("fire", "percolation model", "probability of persistence", "probability of spread"),
  authors = person("Ceres", "Barros", email = "cbarros@mail.ubc.ca", role = c("aut", "cre")),
  childModules = character(0),
  version = numeric_version("0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSpread.Rmd"),
  reqdPkgs = list("raster"),
  parameters = rbind(
    defineParameter("fireSize", "numeric", 1000, NA, NA, desc = "Fire size in pixels"),
    defineParameter("noStartPix", "numeric", 100, NA, NA, desc = "Number of fire events"),
    defineParameter("fireFreq", "numeric", 1, NA, NA, desc = "Fire recurrence in years"),
    defineParameter(".plotMaps", "logical", FALSE, NA, NA, "This describes whether maps should be plotted or not"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA, "This describes the simulation time interval between plot events")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "vegetation", objectClass = "list",
                 desc = "List of vegetation states rasters, with simplified classes"),
    expectsInput(objectName = "spreadProb_mat", objectClass = "matrix",
                 desc = "Two-column matrix of baseline spread probs. per vegetation class. 
                 1st column should contain vegetation class codes, with corresponding probs. in the 2nd column"),
    expectsInput(objectName = "persistProb_mat", objectClass = "matrix",
                 desc = "Two-column matrix of baseline persistence probs. per vegetation class. 
                 1st column should contain vegetation class codes, with corresponding probs. in the 2nd column"),
    expectsInput(objectName = "climate", objectClass = "RasterLayer",
                 desc = "A raster of climate values that will increase fire spread probability and severity; 
                 defaults to the Climate Moisture Index for Canada - Reference Period (1981-2010)",
                 sourceURL = "https://apps-scf-cfs.rncan.gc.ca/opendata/Forest%20Change/Forest%20impact%20maps%20and%20data/Drought/Climate%20Moisture%20Index%20(CMI)/Reference%20Period/CMI%20-%20IHC%20-%201981-2010%20-%20Raster.zip")
    
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = "burnable_areas", objectClass = "RasterLayer",
                  desc = "Raster of areas that are susceptible to burning"),
    createsOutput(objectName = "spreadProb_map", objectClass = "RasterLayer",
                  desc = "Raster of rate of spread per pixel"),
    createsOutput(objectName = "spreadRas", objectClass = "list",
                  desc = "List of rasters of fire spread")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.fireSpread = function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      ## Initialise module
      sim <- fireInit(sim)
      
      ## schedule future event(s)
      sim <- scheduleEvent(sim, start(sim) + 1, "fireSpread", "fireSpread", eventPriority = 1)
    },
    
    fireSpread = {
      ## calculate fire spread
      sim <- do.Fire(sim)
      
      ## schedule future events
      sim <- scheduleEvent(sim, time(sim) + 1, "fireSpread", "fireSpread", eventPriority = 1)
    },
    
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

### module initialization
fireInit <- function(sim) {
  ## make raster storage lists 
  sim$spreadRas <- list()
  
  ## define first fire year
  sim$fireYear <- start(sim) + 1
  
  return(invisible(sim))
}

### Simulate fire spread
do.Fire <- function(sim) {
  if(time(sim) == tail(sim$fireYear, 1)){
    ## MAKE BURNABLE AREAS RASTER -------------------------------
    sim$burnable_areas <- sim$vegetation[[time(sim) - 1]]   ## using previous year vegetation map
    sim$burnable_areas[sim$vegetation[[time(sim) - 1]][] == 0] <- NA
    
    ## MAKE RASTER OF SPREAD PROABILITIES -------------------------------
    sim$spreadProb_map <- sim$vegetation[[time(sim) - 1]]
    sim$spreadProb_map <- reclassify(sim$spreadProb_map, rcl = sim$spreadProb_mat)
    
    ## Climate modulates spread
    sim$spreadProb_map <- sim$spreadProb_map * sim$climate ## drier pixels have "high climate" values/influence
    
    ## MAKE RASTER OF PERSISTENCE PROABILITIES -------------------------------
    sim$persistProb_map <- sim$vegetation[[time(sim) - 1]]
    sim$persistProb_map <- reclassify(sim$persistProb_map, rcl = sim$persistProb_mat)
    
    ## MAKE RASTER OF FIRE SPREAD -------------------------------
    ## note that this function two random components: selection of starting pixels and fire spread
    startPix <- sample(which(!is.na(sim$burnable_areas[])), P(sim)$noStartPix)
    
    ## Favier's model:
    sim$spreadRas[[time(sim) - 1]] <-  spread2(landscape = sim$burnable_areas, spreadProb = sim$spreadProb_map,
                                               persistProb = sim$persistProb_map,
                                               start = startPix, 
                                               maxSize =  P(sim)$fireSize,
                                               plot.it = FALSE)
    
    ## remove fires outside burnable areas
    sim$spreadRas[[time(sim) - 1]][is.na(sim$burnable_areas)] <- NA
    
    ## define next fire year
    sim$fireYear[length(sim$fireYear)+1] <- tail(sim$fireYear, 1) + P(sim)$fireFreq
  } else {
    ## if this is not a fire year make a fire raster with NAs
    sim$spreadRas[[time(sim) - 1]] <- sim$vegetation[[time(sim) - 1]]
    sim$spreadRas[[time(sim) - 1]][] <- NA
  }
  
  
  return(invisible(sim))  
}


## OTHER INPUTS AND FUNCTIONS --------------------------------
.inputObjects = function(sim) {
  if(is.null(sim$spreadProb_mat)){
    sim$spreadProb_mat <- matrix(c(0:5, 0, 0.9, 0.7, 0.2, 0.4, 0.5), byrow = FALSE, nrow = 6, ncol = 2,
                           dimnames = list(paste0("hab", as.character(0:5)), c("hab", "spreadProb")))
    
  }
  
  if(is.null(sim$persistProb_mat)){
    sim$persistProb_mat <- matrix(c(0:5, 0, 0.1, 0.4, 0.3, 0.6, 0.8), byrow = FALSE, nrow = 6, ncol = 2,
                                 dimnames = list(paste0("hab", as.character(0:5)), c("hab", "persistProb")))
    
  }
  
  ## make a  raster template from vegetation - to change later
  if(is.null(sim$climate)) {
    ## default cliamte raster is the Climate Moisture Index 1981-2010
    sim$climate <- prepInputs(targetFile = "cmi_ihc_1981-2010.tif",
                              url = "https://apps-scf-cfs.rncan.gc.ca/opendata/Forest%20Change/Forest%20impact%20maps%20and%20data/Drought/Climate%20Moisture%20Index%20(CMI)/Reference%20Period/CMI%20-%20IHC%20-%201981-2010%20-%20Raster.zip", 
                              destinationPath = dataPath(sim), overwrite = FALSE,
                              fun = "raster::raster", useCache = TRUE, cacheRepo = cachePath(sim))
    
    sim$climate <- if(G(sim)$.useCache) {
      Cache(cropToStudyArea, study.area = sim$vegetation_prefire, tocrop = sim$climate, cacheRepo = cachePath(sim)) 
    } else {
      cropToStudyArea(study.area = sim$vegetation_prefire, tocrop = sim$climate)
    }
    
    ## reverse scale so that drier pixels have higher levels and rescale 0 to 1
    sim$climate[] <- -sim$climate[]
    sim$climate[] <- (sim$climate[] - min(sim$climate[], na.rm = TRUE)) / (max(sim$climate[], na.rm = TRUE) - min(sim$climate[], na.rm = TRUE))
  }
  
  # if(is.null(sim$DEM)) {
  #   sim$DEM <- prepInputs(targetFile = "",
  #                             url = "", 
  #                             destinationPath = dataPath(sim), overwrite = FALSE,
  #                             fun = "raster::raster", useCache = TRUE, cacheRepo = cachePath(sim))
  # }
  
  return(invisible(sim))
}

