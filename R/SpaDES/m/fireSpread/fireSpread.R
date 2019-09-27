# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireSpread",
  description = "Fire spread model using Favier (2004) percolation model, where percolation probabilities are
  conditioned on fire (behaviour) properties calculated using, e.g, the Canadian Forest Fire Behaviour Prediction System",
  keywords = c("fire", "percolation model", "fire-vegetation feedbacks", "fire-climate feedbacks", "FBP system"),
  authors = person("Ceres", "Barros", email = "cbarros@mail.ubc.ca", role = c("aut", "cre")),
  childModules = character(0),
  version = numeric_version("0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSpread.Rmd"),
  reqdPkgs = list("R.utils", "raster", "data.table", "dplyr", "humidity", #"cffdrs",  ## cffdrs is causing installation problems
                  "sf", "PredictiveEcology/SpaDES.core@development",
                  "PredictiveEcology/SpaDES.tools@development",
                  "PredictiveEcology/reproducible@development"),
  parameters = rbind(
    defineParameter("fireSize", "integer", 1000L, NA, NA, desc = "Fire size in pixels"),
    defineParameter("noStartPix", "integer", 100L, NA, NA, desc = "Number of fire events"),
    defineParameter(name = "fireInitialTime", class = "numeric", default = 2L,
                    desc = "The event time that the first fire disturbance event occurs"),
    defineParameter(name = "fireTimestep", class = "numeric", default = 2L,
                    desc = "The number of time units between successive fire events in a fire module"),
    defineParameter(".plotMaps", "logical", FALSE, NA, NA, "This describes whether maps should be plotted or not"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA, "This describes the simulation time interval between plot events"),
    defineParameter(".useCache", "logical", "init", NA, NA,
                    desc = "use caching for the spinup simulation?")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "fireCFBRas", objectClass = "RasterLayer",
                 desc = "Raster of crown fraction burnt"),
    expectsInput(objectName = "fireIntRas", objectClass = "RasterLayer",
                 desc = "Raster of equilibrium head fire intensity [kW/m]"),
    expectsInput(objectName = "fireROSRas", objectClass = "RasterLayer",
                 desc = "Raster of equilibrium rate of spread [m/min]"),
    expectsInput(objectName = "fireRSORas", objectClass = "RasterLayer",
                 desc = "Critical spread rate for crowning [m/min]"),
    expectsInput(objectName = "fireTFCRas", objectClass = "RasterLayer",
                 desc = "Raster of total fuel consumed [kg/m^2]"),
    expectsInput(objectName = "rasterToMatch", "RasterLayer",
                 desc = "a raster of the studyArea in the same resolution and projection as biomassMap ",
                 sourceURL = NA),
    expectsInput(objectName = "simulatedBiomassMap", objectClass = "RasterLayer",
                 desc = "Biomass map at each succession time step. Default is Canada national biomass map",
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-StructureBiomass.tar"),
    expectsInput("studyArea", "SpatialPolygonsDataFrame",
                 desc = paste("Polygon to use as the study area.",
                              "Defaults to  an area in Southwestern Alberta, Canada."),
                 sourceURL = "")
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = "fireCFBRas", objectClass = "RasterLayer",
                  desc = "Raster of crown fraction burnt"),
    createsOutput(objectName = "fireIntRas", objectClass = "RasterLayer",
                  desc = "Raster of equilibrium head fire intensity [kW/m]"),
    createsOutput(objectName = "fireROSRas", objectClass = "RasterLayer",
                  desc = "Raster of equilibrium rate of spread [m/min]"),
    createsOutput(objectName = "fireRSORas", objectClass = "RasterLayer",
                  desc = "Critical spread rate for crowning [m/min]"),
    createsOutput(objectName = "fireTFCRas", objectClass = "RasterLayer",
                  desc = "Raster of total fuel consumed [kg/m^2]"),
    createsOutput(objectName = "fireYear", objectClass = "numeric", desc = "Next fire year"),
    createsOutput(objectName = "pixelGroupMapFBP", objectClass = "RasterLayer",
                  desc = "updated community map at each succession time step, on FBP-compatible projection"),
    createsOutput(objectName = "rstCurrentBurn", objectClass = "RasterLayer",
                  desc = "Binary raster of fire spread"),
    createsOutput(objectName = "startPix", objectClass = "vector",
                  desc = "List of starting fire pixels")
  )
))

doEvent.fireSpread = function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      ## Initialise module
      sim <- fireSpreadInit(sim)

      ## schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$fireInitialTime, "fireSpread",
                           "doFireSpread", eventPriority = 2.5) ## always schedule fire
    },
    doFireSpread = {
      ## calculate fire spread in fire years
      if(time(sim) == sim$fireYear) {
        sim <- doFireSpread(sim)

        ## define next fire year
        sim$fireYear <- time(sim) + P(sim)$fireTimestep
      } else {
        ## No fire
        sim <- doNoFire(sim)
      }

      ## schedule future event(s) - always schedule fire
      sim <- scheduleEvent(sim, time(sim) + 1, "fireSpread",
                           "doFireSpread", eventPriority = 2.5)  ## always schedule fire
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

### module initialization
fireSpreadInit <- function(sim) {
  return(invisible(sim))
}

## Fire spread event in fire years - rasters should be back in LBMR projection
doFireSpread <- function(sim) {
  ## MAKE BURNABLE AREAS RASTER -------------------------------
  ## only areas with biomass can burn
  ## if no simulatedBiomassMap is supplied then generate one from raw data
  ## at the start
  if (time(sim) == P(sim)$fireInitialTime) {
    if (is.null(sim$simulatedBiomassMap)) {
      if (is.null(sim$biomassMap)) {
        cacheTags <- c(currentModule(sim), current(sim))
        dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)

        # If biomassMap is not present either, get rawBiomassMap, but crop it to studyArea/RTM instead of SALarge/RTMLarge
        rawBiomassMapFilename <- file.path(dPath, "NFI_MODIS250m_kNN_Structure_Biomass_TotalLiveAboveGround_v0.tif")
        rawBiomassMapURL <- "http://tree.pfc.forestry.ca/kNN-StructureBiomass.tar"
        burnableAreas <- Cache(prepInputs,
                               targetFile = asPath(basename(rawBiomassMapFilename)),
                               archive = asPath(c("kNN-StructureBiomass.tar",
                                                  "NFI_MODIS250m_kNN_Structure_Biomass_TotalLiveAboveGround_v0.zip")),
                               url = rawBiomassMapURL,
                               destinationPath = dPath,
                               studyArea = sim$studyArea,
                               rasterToMatch = sim$rasterToMatch,
                               maskWithRTM = TRUE,
                               useSAcrs = FALSE,
                               method = "bilinear",
                               datatype = "INT2U",
                               filename2 = TRUE, overwrite = TRUE,
                               userTags = cacheTags,
                               omitArgs = c("destinationPath", "targetFile", cacheTags, "stable"))
        rm(cacheTags)
      } else {
        burnableAreas <- sim$biomassMap
      }
    } else {
      burnableAreas <- sim$simulatedBiomassMap
    }
  } else {
    burnableAreas <- sim$simulatedBiomassMap
  }

  vals <- data.table(B = getValues(burnableAreas))   ## making a mask is probably faster with data.table
  vals <- vals[B > 0, B := 1]
  vals <- vals[B <= 0, B := NA]
  burnableAreas[] <- vals$B

  ## MAKE RASTER OF SPREAD PROABILITIES
  ## spread probability is the combination of ROS and intensity, which have an additive effect
  ## and their sum is scaled to 0-0.23
  ## TODO: the scaling should guarantee an average value of 0.23
  ## TODO: ROS and intensity should be combined differently
  # browser()
  spreadProb_map <- sim$fireROSRas + sim$fireIntRas
  spreadProb_map <- mask(spreadProb_map, burnableAreas)

  vals <- data.table(spreadP = getValues(spreadProb_map))   ## making a mask is probably faster with data.table
  vals[!is.na(spreadP), spreadP := scale(spreadP, scale = FALSE) + 0.20]
  spreadProb_map[] <- vals$spreadP

  ## NAs get 0 probability - not necessary
  # spreadProb_map[is.na(getValues(spreadProb_map))] <- 0

  ## MAKE RASTER OF PERSISTENCE PROABILITIES
  ## persistence probability is the combination of TFC and intensity, which have an additive effect
  ## and their sum is scaled to 0-1
  ## TODO: TFC and intensity should be combined differently
  persistProb_map <- sim$fireTFCRas + sim$fireIntRas
  persistProb_map <- mask(persistProb_map, burnableAreas)

  vals <- data.table(persisP = getValues(persistProb_map))   ## making a mask is probably faster with data.table
  vals[!is.na(persisP), persisP := scales::rescale(persisP, to = c(0,1))]
  persistProb_map[] <- vals$persisP

  ## check if NAs match
  if (any(!is.na(spreadProb_map[is.na(persistProb_map[])])))
    stop("spread and persistence probability rasters have unmatching NAs")
  ## redo burnable areas if missing fire probabilities
  if (any(!is.na(burnableAreas[is.na(spreadProb_map[])])))
    burnableAreas <- mask(burnableAreas, spreadProb_map)


  ## MAKE RASTER OF FIRE SPREAD -------------------------------
  ## note that this function two random components: selection of starting pixels and fire spread
  sim$startPix <- sample(which(!is.na(getValues(burnableAreas))), P(sim)$noStartPix)

  ## Favier's model:
  rstCurrentBurn <- spread2(landscape = burnableAreas,
                            spreadProb = spreadProb_map,
                            persistProb = persistProb_map,
                            start = sim$startPix,
                            maxSize =  P(sim)$fireSize,
                            plot.it = FALSE)

  ## remove fires that spread beyond burnable areas
  if (any(!is.na(rstCurrentBurn[is.na(burnableAreas[])])))
    rstCurrentBurn <- mask(rstCurrentBurn, burnableAreas)

  ## convert to mask
  rstCurrentBurn[!is.na(rstCurrentBurn[])][] <- 1

  ## remove pixels that didn't burn from fire property rasters
  sim$fireCFBRas <- mask(sim$fireCFBRas, rstCurrentBurn)
  sim$fireIntRas <- mask(sim$fireIntRas, rstCurrentBurn)
  sim$fireROSRas <- mask(sim$fireROSRas, rstCurrentBurn)
  sim$fireRSORas <- mask(sim$fireRSORas, rstCurrentBurn)
  sim$fireTFCRas <- mask(sim$fireTFCRas, rstCurrentBurn)

  ## export to sim
  sim$rstCurrentBurn <- rstCurrentBurn
  return(invisible(sim))
}

## What to do in no fire years
doNoFire <- function(sim) {
  sim$rstCurrentBurn <- NULL

  return(invisible(sim))
}

## OTHER INPUTS AND FUNCTIONS --------------------------------
.inputObjects = function(sim) {
  ## TODO: ADD DUMMIES FOR FIRE PROPERTIES
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  cacheTags <- c(currentModule(sim), "function:.inputObjects")

  if (!suppliedElsewhere("studyArea", sim)) {
    message("'studyArea' was not provided by user. Using a polygon (6250000 m^2) in southwestern Alberta, Canada")
    sim$studyArea <- randomStudyArea(seed = 1234, size = (250^2)*100)
  }

  ## DEFAULT RASTER TO MATCH
  needRTM <- FALSE
  if (is.null(sim$rasterToMatch)) {
    if (!suppliedElsewhere("rasterToMatch", sim)) {      ## if one is not provided, re do both (safer?)
      needRTM <- TRUE
      message("There is no rasterToMatch supplied; will attempt to use rawBiomassMap")
    } else {
      stop("rasterToMatch is going to be supplied, but ", currentModule(sim), " requires it ",
           "as part of its .inputObjects. Please make it accessible to ", currentModule(sim),
           " in the .inputObjects by passing it in as an object in simInit(objects = list(rasterToMatch = aRaster)",
           " or in a module that gets loaded prior to ", currentModule(sim))
    }
  }

  if (needRTM) {
    if(!suppliedElsewhere("rawBiomassMap", sim)) {
      rawBiomassMapFilename <- file.path(dPath, "NFI_MODIS250m_kNN_Structure_Biomass_TotalLiveAboveGround_v0.tif")
      rawBiomassMap <- Cache(prepInputs,
                             targetFile = asPath(basename(rawBiomassMapFilename)),
                             archive = asPath(c("kNN-StructureBiomass.tar",
                                                "NFI_MODIS250m_kNN_Structure_Biomass_TotalLiveAboveGround_v0.zip")),
                             url = extractURL("rawBiomassMap"),
                             destinationPath = dPath,
                             studyArea = sim$studyArea,
                             rasterToMatch = if (!needRTM) sim$rasterToMatch else NULL,
                             maskWithRTM = if (!needRTM) TRUE else FALSE,
                             useSAcrs = FALSE,     ## never use SA CRS
                             method = "bilinear",
                             datatype = "INT2U",
                             filename2 = TRUE, overwrite = TRUE, userTags = cacheTags,
                             omitArgs = c("destinationPath", "targetFile", "userTags", "stable"))

    } else {
      rawBiomassMap <- Cache(postProcess,
                             x = sim$rawBiomassMap,
                             studyArea - sim$studyArea,
                             maskWithRTM = FALSE,
                             useSAcrs = FALSE,     ## never use SA CRS
                             method = "bilinear",
                             datatype = "INT2U",
                             filename2 = TRUE, overwrite = TRUE, userTags = cacheTags,
                             omitArgs = c("destinationPath", "userTags"))
    }

    ## if we need rasterToMatch, that means a) we don't have it, but b) we will have rawBiomassMap
    ## even if one of the rasterToMatch is present re-do both.
    ## if we need rasterToMatch, that means a) we don't have it, but b) we will have rawBiomassMap
    sim$rasterToMatch <- rawBiomassMap
    RTMvals <- getValues(sim$rasterToMatch)
    sim$rasterToMatch[!is.na(RTMvals)] <- 1

    sim$rasterToMatch <- Cache(writeOutputs, sim$rasterToMatch,
                               filename2 = file.path(cachePath(sim), "rasters", "rasterToMatch.tif"),
                               datatype = "INT2U", overwrite = TRUE)
  }

  ## DEFAULT FIRE PROPERTIES RASTERS
  if (any(!suppliedElsewhere("fireRSORas", sim),
          !suppliedElsewhere("fireROSRas", sim),
          !suppliedElsewhere("fireIntRas", sim),
          !suppliedElsewhere("fireTFCRas", sim),
          !suppliedElsewhere("fireCFBRas", sim))) {
    message(crayon::red(paste0("fireSpread is missing one/several of the following rasters:\n",
                               "  fireRSORas, fireROSRas, fireIntRas, fireTFCRas and fireCFBRas.\n",
                               "  DUMMY RASTERS will be used - if this is not intended, please \n",
                               "  use a fire module that provides them (e.g. fireSpread)")))
    vals <- getValues(sim$rasterToMatch)
    valsCFB <- valsInt <- valsROS <- valsRSO <- valsTFC <- integer(0)
    valsCFB[!is.na(vals)] <- runif(sum(!is.na(vals)), 0, 1)
    valsROS[!is.na(vals)] <- as.integer(round(runif(sum(!is.na(vals)), 0, 100)))
    valsRSO[!is.na(vals)] <- as.integer(round(runif(sum(!is.na(vals)), 0, 100)))

    ## TODO: decide on best values
    browser()

    valsInt[!is.na(vals)] <- runif(sum(!is.na(vals)), 0, 1)
    valsTFC[!is.na(vals)] <- runif(sum(!is.na(vals)), 0, 1)

    sim$fireCFBRas <- setValues(sim$rasterToMatch, valsCFB)
    sim$fireIntRas <- setValues(sim$rasterToMatch, valsInt)
    sim$fireROSRas <- setValues(sim$rasterToMatch, valsROS)
    sim$fireRSORas <- setValues(sim$rasterToMatch, valsRSO)
    sim$fireTFCRas <- setValues(sim$rasterToMatch, valsTFC)
  }
  return(invisible(sim))
}
