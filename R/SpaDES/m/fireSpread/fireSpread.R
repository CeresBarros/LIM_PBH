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
                  "sf", "scales", "crayon",
                  "PredictiveEcology/SpaDES.core@development",
                  "PredictiveEcology/SpaDES.tools@development",
                  "PredictiveEcology/reproducible@development",
                  "PredictiveEcology/LandR@development (>= 1.0.0.9003)"),
  parameters = rbind(
    defineParameter("noStartPix", "integer", 100L, 0L, NA,
                    desc = paste("Number of fire events. Only used if fireIgnitionProb is not available")),
    defineParameter(name = "fireInitialTime", class = "numeric", default = 2L,
                    desc = "The event time that the first fire disturbance event occurs"),
    defineParameter(name = "fireTimestep", class = "numeric", default = 2L,
                    desc = "The number of time units between successive fire events in a fire module"),
    defineParameter(name = "persistProbRange", class = "numeric", default = c(0, 1), min = 0, max = 1,
                    desc = paste("Range of values to use when rescaling presistence probability, q,",
                                 "calculated from fire properties. Defaults to 0-1, the minimum and",
                                 "maximum values possible. Must be a sorted vector of length 2",
                                 "If min(persistProbRange) != 0, 0s are forced to 0 after scaling")),
    defineParameter(name = "spreadProbRange", class = "numeric", default = c(0.20, 0.25), min = 0, max = 1,
                    desc = paste("Range of values to use when rescaling spread probability, p,",
                                 "calculated from fire properties. Defaults to 0.20-0.25, which",
                                 "produces sensible fire sizes, but can be 0-1. Must be a sorted vector of length 2.",
                                 "If min(spreadProbRange) != 0, 0s in the original scale are forced to 0 after scaling")),
    defineParameter(".plotMaps", "logical", FALSE, NA, NA, "This describes whether maps should be plotted or not"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA, "This describes the simulation time interval between plot events"),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used. If NA, a hash of studyArea will be used."),
    defineParameter(".useCache", "logical", "init", NA, NA,
                    desc = "use caching for the spinup simulation?")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "fireCFBRas", objectClass = "RasterLayer",
                 desc = "Raster of crown fraction burnt"),
    expectsInput(objectName = "fireIntRas", objectClass = "RasterLayer",
                 desc = "Raster of equilibrium head fire intensity [kW/m]"),
    expectsInput(objectName = "fireIgnitionProb", objectClass = "RasterLayer",
                 desc = paste("Raster of expected no. of ignitions (lambda in a poisson ditribution). Optional.",
                              "If not present, will use 'noStartPix' to randomly start a given number of fires",
                              "in the landscape")),
    expectsInput(objectName = "fireROSRas", objectClass = "RasterLayer",
                 desc = "Raster of equilibrium rate of spread [m/min]"),
    expectsInput(objectName = "fireRSORas", objectClass = "RasterLayer",
                 desc = "Critical spread rate for crowning [m/min]"),
    expectsInput(objectName = "fireTFCRas", objectClass = "RasterLayer",
                 desc = "Raster of total fuel consumed [kg/m^2]"),
    expectsInput(objectName = "fireSize", objectClass = "integer",
                 desc = paste("Fire size in pixels. Defaults to maximum fire size obtained from fire perimeter",
                              "records in the study area (using all fire perimeter polygons in the Canadian National Fire Database)."),
                 sourceURL = "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/current_version/NFDB_poly.zip"),
    expectsInput(objectName = "rasterToMatch", "RasterLayer",
                 desc = "a raster of the studyArea in the same resolution and projection as biomassMap ",
                 sourceURL = NA),
    expectsInput(objectName = "simulatedBiomassMap", objectClass = "RasterLayer",
                 desc = paste("Biomass map at each succession time step. If not supplied, will use Canadian Forestry",
                              "Service, National Forest Inventory, kNN-derived total aboveground biomass map",
                              "from 2001. See https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990",
                              "for metadata"),
                 sourceURL = NA),
    expectsInput("studyArea", "SpatialPolygonsDataFrame",
                 desc = paste("Polygon to use as the study area.",
                              "Defaults to  an area in Southwestern Alberta, Canada."),
                 sourceURL = NA)
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
    createsOutput(objectName = "rstCurrentFires", objectClass = "RasterLayer",
                  desc = "Raster of fire spread with fire IDs (not binary)"),
    createsOutput(objectName = "startPix", objectClass = "vector",
                  desc = "List of starting fire pixels")
  )
))

doEvent.fireSpread = function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      ## Initialise module
      sim <- Init(sim)

      ## schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$fireInitialTime, "fireSpread",
                           "doFireSpread", eventPriority = 2.5) ## always schedule fire
    },
    doFireSpread = {
      ## calculate fire spread in fire years
      if (time(sim) == sim$fireYear) {
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
Init <- function(sim) {
  ## checks
  if (start(sim) == P(sim)$fireInitialTime)
    warning(red("start(sim) and P(sim)$fireInitialTime are the same.\nThis may create bad scheduling with init events"))

  if (length(P(sim)$spreadProbRange) != 2) {
    stop("'P(sim)$spreadProbRange' must be a vector of length = 2")
  }

  if (length(P(sim)$persistProbRange) != 2) {
    stop("'P(sim)$persistProbRange' must be a vector of length = 2")
  }

  if (min(P(sim)$spreadProbRange) < 0 | max(P(sim)$spreadProbRange) > 1) {
    stop("'P(sim)$spreadProbRange' must be between [0,1]")
  }

  if (min(P(sim)$persistProbRange) < 0 | max(P(sim)$persistProbRange) > 1) {
    stop("'P(sim)$persistProbRange' must be between [0,1]")
  }

  if (P(sim)$spreadProbRange[1] > P(sim)$spreadProbRange[2]) {
    stop("'P(sim)$spreadProbRange[1]' must be <= 'P(sim)$spreadProbRange[2]'")
  }

  if (P(sim)$persistProbRange[1] > P(sim)$persistProbRange[2]) {
    stop("'P(sim)$persistProbRange[1]' must be <= 'P(sim)$persistProbRange[2]'")
  }

  ## try to make fireIgnitionProb from fireSense_IgnitionPredicted
  ## if fireSense_IgntionPredicted not present, try again later
  mod$useFireSense <- FALSE
  if (is.null(sim$fireIgnitionProb)) {
    if (suppliedElsewhere("fireSense_IgnitionPredicted", sim)) {
      message(blue(paste("'fireIgnitionProb' raster was not supplied, but 'fireSense_IgnitionPredicted' exists in sim.",
                         "Using 'fireSense_IgnitionPredicted' as 'fireIgnitionProb' in", currentModule(sim))))
      sim$fireIgnitionProb <- sim$fireSense_IgnitionPredicted
      mod$useFireSense <- TRUE
    }
  }

  ## check if ignition raster matches RTM
  ## e.g. fireSense_IgnitionPredict projects at lower res. and this may need to be checked at each fireTimeStep
  if (!is.null(sim$fireIgnitionProb)) {
    if (!compareRaster(sim$fireIgnitionProb, sim$rasterToMatch, stopiffalse = FALSE)) {
      message(blue(paste("Properties of 'fireIgnitionProb' and 'rasterToMatch' differ.",
                         "Projecing/masking 'fireIgnitionProb' to 'rasterToMatch")))
      origRes <- res(sim$fireIgnitionProb)[1]
      finalRes <- res(sim$rasterToMatch)[1]
      sim$fireIgnitionProb <- postProcess(sim$fireIgnitionProb,
                                          rasterToMatch = sim$rasterToMatch,
                                          maskWithRTM = TRUE,
                                          method = "bilinear",
                                          filename2 = NULL, ## don't save
                                          useCache = FALSE)  ## don't cache

      if (origRes != finalRes) {
        message(blue(paste("Resolution of 'fireIgnitionProb' and 'rasterToMatch' differed.",
                           "Rescaling values in fireIgnitionProb")))
        sim$fireIgnitionProb[] <- sim$fireIgnitionProb[] * (finalRes/origRes ^ 2)
      }
    }
  }

  return(invisible(sim))
}

## Fire spread event in fire years - rasters should be back in LandR Biomass projection
doFireSpread <- function(sim) {
  ## make or update fireIgnitionProb from fireSense_IgnitionPredicted
  if (time(sim) == P(sim)$fireInitialTime) {
    if (is.null(sim$fireIgnitionProb)) {
      if (suppliedElsewhere("fireSense_IgnitionPredicted", sim)) {
        message(blue(paste("'fireIgnitionProb' raster was not supplied, but 'fireSense_IgnitionPredicted' exists in sim.",
                           "Using 'fireSense_IgnitionPredicted' as 'fireIgnitionProb' in", currentModule(sim))))
        sim$fireIgnitionProb <- sim$fireSense_IgnitionPredicted
        mod$useFireSense <- TRUE
      } else {
        message(blue(paste("'fireIgnitionProb' raster was not supplied. Fires will be ignited",
                           "randomly across the landscape, in number = to 'noStartPix'")))
      }
    }
  } else {
    if (!is.null(sim$fireIgnitionProb)) {
      if (isTRUE(mod$useFireSense)) {
        message(blue("Updating 'fireIgnitionProb' with 'fireSense_IgnitionPredicted'"))
        sim$fireIgnitionProb <- sim$fireSense_IgnitionPredicted
      }
    }
  }


  ## check if ignition raster matches RTM
  ## e.g. fireSense_IgnitionPredict projects at lower res. and this may need to be checked at each fireTimeStep
  ## thise needs to be done each year, in case fireIgnitionsProb is being updated
  if (!is.null(sim$fireIgnitionProb)) {
    if (!compareRaster(sim$fireIgnitionProb, sim$rasterToMatch, stopiffalse = FALSE)) {
      message(blue(paste("Properties of 'fireIgnitionProb' and 'rasterToMatch' differ.",
                         "Projecing/masking 'fireIgnitionProb' to 'rasterToMatch")))
      origRes <- res(sim$fireIgnitionProb)[1]
      finalRes <- res(sim$rasterToMatch)[1]
      sim$fireIgnitionProb <- postProcess(sim$fireIgnitionProb,
                                          rasterToMatch = sim$rasterToMatch,
                                          maskWithRTM = TRUE,
                                          method = "bilinear",
                                          filename2 = NULL, ## don't save
                                          useCache = FALSE)  ## don't cache

      if (origRes != finalRes) {
        message(blue(paste("Resolution of 'fireIgnitionProb' and 'rasterToMatch' differed.",
                           "Rescaling values in fireIgnitionProb")))
        sim$fireIgnitionProb[] <- sim$fireIgnitionProb[] * (finalRes/origRes ^ 2)
      }
    }
  }

  ## MAKE BURNABLE AREAS RASTER -------------------------------
  ## only areas with biomass can burn if no non-forest fire spread is allowed
  ## if no simulatedBiomassMap is supplied then generate one from raw data
  ## at the start
  ## when non-forest fires are allowed, fire can spread across the whole SA
  if (is.null(sim$pixelNonForestFuels)) {
    if (is.null(sim$simulatedBiomassMap)) {
      if (is.null(sim$biomassMap)) {
        cacheTags <- c(currentModule(sim), current(sim))
        dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)

        # If biomassMap is not present either, get rawBiomassMap, but crop it to studyArea/RTM instead of SALarge/RTMLarge
        rawBiomassMapURL <- paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                   "canada-forests-attributes_attributs-forests-canada/",
                                   "2001-attributes_attributs-2001/",
                                   "NFI_MODIS250m_2001_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif")
        rawBiomassMapFilename <- "NFI_MODIS250m_2001_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif"
        rawBiomassMap <- Cache(prepInputs,
                               targetFile = rawBiomassMapFilename,
                               url = rawBiomassMapURL,
                               destinationPath = dPath,
                               studyArea = sim$studyArea,
                               rasterToMatch = sim$rasterToMatch,
                               maskWithRTM = TRUE,
                               useSAcrs = FALSE,
                               method = "bilinear",
                               datatype = "INT2U",
                               filename2 = NULL,
                               userTags = c(cacheTags, "rawBiomassMap"),
                               omitArgs = c("destinationPath", "targetFile", cacheTags, "stable"))

        burnableAreas <- rawBiomassMap
        rm(cacheTags, rawBiomassMap)
      } else {
        burnableAreas <- sim$biomassMap
      }
    } else {
      burnableAreas <- sim$simulatedBiomassMap
    }
  } else {
    burnableAreas <- sim$rasterToMatch
  }

  vals <- data.table(B = getValues(burnableAreas))   ## making a mask is probably faster with data.table
  vals <- vals[B > 0, B := 1]
  vals <- vals[B <= 0, B := NA]
  burnableAreas[] <- vals$B

  ## MAKE RASTER OF SPREAD PROBABILITIES
  ## spread probability is the combination of ROS and intensity, which have an multiplicative effect
  ## TODO: ROS and intensity should be combined differently
  spreadProb_map <- sim$fireROSRas * sim$fireIntRas
  spreadProb_map <- mask(spreadProb_map, burnableAreas)

  vals <- data.table(spreadP = getValues(spreadProb_map))   ## making a mask is probably faster with data.table
  ## before:
  # vals[!is.na(spreadP) & spreadP > 0, spreadPsc := scales::rescale(spreadP, to = P(sim)$spreadProbRange)]

  ## now: there are real 0s, so these shouldn't be excluded but need to stay 0s at the end
  ## this rescaling is not affected by the shape of the distribution, as the relative positions of the values
  ## and dist. shape are maintained in the new scale
  ## if the min(P(sim)$spreadProbRange) != 0, small values will only exist for original values close to 0.
  ## same rationale for persitence prob.
  vals[!is.na(spreadP), spreadPsc := scales::rescale(spreadP, to = P(sim)$spreadProbRange)]
  vals[spreadP == 0, spreadPsc := 0]
  spreadProb_map[] <- vals$spreadPsc

  ## MAKE RASTER OF PERSISTENCE PROBABILITIES
  ## persistence probability is the combination of TFC and intensity, as a ratio
  ## (higher intensity fires should a same amount of biomass for less time tan a low intensity fire)
  ## Note that 0 denominator will create NAs, so these need to be zeroed
  ## and their ratio is scaled to 0-1
  ## TODO: TFC and intensity should be combined differently
  persistProb_map <- sim$fireTFCRas / sim$fireIntRas
  persistProb_map[sim$fireIntRas[] == 0] <- 0
  persistProb_map <- mask(persistProb_map, burnableAreas)

  vals <- data.table(persisP = getValues(persistProb_map))   ## making a mask is probably faster with data.table
  vals[!is.na(persisP), persisPsc := scales::rescale(persisP, to = P(sim)$persistProbRange)]
  vals[persisP == 0, persisPsc := 0]
  persistProb_map[] <- vals$persisPsc

  ## check if NAs match
  if (getOption("LandR.assertions"))
    if (any(!is.na(spreadProb_map[is.na(persistProb_map[])])))
      stop("spread and persistence probability rasters have unmatching NAs")

  ## redo burnable areas if missing fire probabilities
  if (any(!is.na(burnableAreas[is.na(spreadProb_map[])])))
    burnableAreas <- mask(burnableAreas, spreadProb_map)


  ## MAKE RASTER OF FIRE SPREAD -------------------------------
  ## note that this function has two random components: selection of starting pixels and fire spread
  ## Favier's model:

  if (is.null(sim$fireIgnitionProb)) {
    message(blue(paste("'fireIgnitionProb' raster was not supplied. Igniting fires",
                       "randomly across the landscape, in number = to 'noStartPix'")))
    sim$startPix <- sample(which(!is.na(getValues(burnableAreas))), P(sim)$noStartPix)
  } else {
    ## draw prob of having a fire, assess "winners", convert to vector (also export to sim)
    startPix <- mask(sim$fireIgnitionProb, burnableAreas)
    startPix <- rbinom(n = ncell(startPix), size = 1, prob = pmin(startPix[], 1))
    startPix <- which(startPix > 0) ## winners are 0 or larger.
    sim$startPix <- sample(startPix)  ## randomize order so that first fires aren't always at top of landscape
  }

  if (length(sim$startPix)) {
    ## run spread2 for one iteration to simulate escape - a bit like fireSense.R
    escapedFires <- spread2(landscape = burnableAreas,
                            spreadProb = spreadProb_map,
                            start = sim$startPix,
                            iterations = 1,
                            plot.it = FALSE,
                            asRaster = FALSE)

    rstCurrentBurn <- spread2(landscape = burnableAreas,
                              spreadProb = spreadProb_map,
                              persistProb = persistProb_map,
                              start = escapedFires,
                              # start = sim$startPix,
                              maxSize = sim$fireSize,
                              plot.it = FALSE)
  } else {
    rstCurrentBurn <- setValues(burnableAreas, rep(NA, ncell(burnableAreas)))
  }

  ## remove fires that only burned one pixel - these didn't really spread
  if (sum(rstCurrentBurn[], na.rm = TRUE)) {
    noSpreadFires <- as.data.table(table(rstCurrentBurn[]))
    noSpreadFires <- noSpreadFires[N == 1]
    rstCurrentBurn[rstCurrentBurn[] %in% noSpreadFires$V1] <- NA
  }

  ## remove fires that spread beyond burnable areas
  if (any(!is.na(rstCurrentBurn[is.na(burnableAreas[])]))) {
    rstCurrentBurn <- mask(rstCurrentBurn, burnableAreas)
  }

  ## output fire raster with fire IDs
  sim$rstCurrentFires <- rstCurrentBurn
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
    stop("Please provide a 'studyArea' polygon")
    # message("'studyArea' was not provided by user. Using a polygon (6250000 m^2) in southwestern Alberta, Canada")
    # sim$studyArea <- randomStudyArea(seed = 1234, size = (250^2)*100)  # Jan 2021 we agreed to force user to provide a SA/SAL
  }

  if (is.na(P(sim)$.studyAreaName)) {
    params(sim)[[currentModule(sim)]][[".studyAreaName"]] <- reproducible::studyAreaName(sim$studyArea)
    message("The .studyAreaName is not supplied; derived name from sim$studyArea: ",
            params(sim)[[currentModule(sim)]][[".studyAreaName"]])
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
    ## if rawBiomassMap exists, it needs to match SALarge, if it doesn't make it
    if (!suppliedElsewhere("rawBiomassMap", sim) ||
        !compareRaster(sim$rawBiomassMap, sim$studyArea, stopiffalse = FALSE)) {
      rawBiomassMapURL <- paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                 "canada-forests-attributes_attributs-forests-canada/",
                                 "2001-attributes_attributs-2001/",
                                 "NFI_MODIS250m_2001_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif")
      rawBiomassMapFilename <- "NFI_MODIS250m_2001_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif"
      rawBiomassMap <- Cache(prepInputs,
                             targetFile = rawBiomassMapFilename,
                             url = rawBiomassMapURL,
                             destinationPath = dPath,
                             studyArea = sim$studyArea,
                             rasterToMatch = NULL,
                             maskWithRTM = FALSE,
                             useSAcrs = FALSE,     ## never use SA CRS
                             method = "bilinear",
                             datatype = "INT2U",
                             filename2 = NULL,
                             userTags = c(cacheTags, "rawBiomassMap"),
                             omitArgs = c("destinationPath", "targetFile", "userTags", "stable"))
    } else {
      rawBiomassMap <- Cache(postProcess,
                             x = sim$rawBiomassMap,
                             studyArea = sim$studyArea,
                             useSAcrs = FALSE,
                             maskWithRTM = FALSE,   ## mask with SA
                             method = "bilinear",
                             datatype = "INT2U",
                             filename2 = NULL,
                             overwrite = TRUE,
                             userTags = cacheTags,
                             omitArgs = c("destinationPath", "targetFile", "userTags", "stable"))
    }

    ## if we need rasterToMatch, that means a) we don't have it, but b) we will have rawBiomassMap
    ## even if one of the rasterToMatch is present re-do both.
    ## if we need rasterToMatch, that means a) we don't have it, but b) we will have rawBiomassMap
    sim$rasterToMatch <- rawBiomassMap
    RTMvals <- getValues(sim$rasterToMatch)
    sim$rasterToMatch[!is.na(RTMvals)] <- 1

    sim$rasterToMatch <- Cache(writeOutputs, sim$rasterToMatch,
                               filename2 = .suffix(file.path(dPath, "rasterToMatch.tif"),
                                                   paste0("_", P(sim)$.studyAreaName)),
                               datatype = "INT2U", overwrite = TRUE)
  }

  if (!compareCRS(sim$studyArea, sim$rasterToMatch)) {
    warning(paste0("studyArea and rasterToMatch projections differ.\n",
                   "studyArea will be projected to match rasterToMatch"))
    sim$studyArea <- spTransform(sim$studyArea, crs(sim$rasterToMatch))
    sim$studyArea <- fixErrors(sim$studyArea)
  }

  ## DEFAULT FIRE PROPERTIES RASTERS
  if (any(!suppliedElsewhere("fireRSORas", sim),
          !suppliedElsewhere("fireROSRas", sim),
          !suppliedElsewhere("fireIntRas", sim),
          !suppliedElsewhere("fireTFCRas", sim),
          !suppliedElsewhere("fireCFBRas", sim))) {
    message(crayon::red(paste0("fireSpread is missing one/several of the following rasters:\n",
                               "  fireRSORas, fireROSRas, fireIntRas, fireTFCRas, fireCFBRas or fireIgnitionProb.\n",
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
    valvalsIgnitsTFC[!is.na(vals)] <- runif(sum(!is.na(vals)), 0, 1)

    sim$fireCFBRas <- setValues(sim$rasterToMatch, valsCFB)
    sim$fireIntRas <- setValues(sim$rasterToMatch, valsInt)
    sim$fireROSRas <- setValues(sim$rasterToMatch, valsROS)
    sim$fireRSORas <- setValues(sim$rasterToMatch, valsRSO)
    sim$fireTFCRas <- setValues(sim$rasterToMatch, valsTFC)
  }

  ## try to make fireIgnitionProb from fireSense_IgnitionPredicted
  ## if fireSense_IgnitionPredicted not available, try again during init.
  if (!suppliedElsewhere("fireIgnitionProb", sim)) {
    if (suppliedElsewhere("fireSense_IgnitionPredicted", sim)) {
      message(blue(paste("'fireIgnitionProb' raster was not supplied, but 'fireSense_IgnitionPredicted' exists in sim.",
                         "Using 'fireSense_IgnitionPredicted' as 'fireIgnitionProb' in", currentModule(sim))))
      sim$fireIgnitionProb <- sim$fireSense_IgnitionPredicted
    }
  }

  ## try to make fireIgnitionProb from fireSense_IgnitionPredicted
  ## if fireSense_IgnitionPredicted not available, try again during init.
  if (!suppliedElsewhere("fireIgnitionProb", sim)) {
    if (suppliedElsewhere("fireSense_IgnitionPredicted", sim)) {
      message(blue(paste("'fireIgnitionProb' raster was not supplied, but 'fireSense_IgnitionPredicted' exists in sim.",
                         "Using 'fireSense_IgnitionPredicted' as 'fireIgnitionProb' in", currentModule(sim))))
      sim$fireIgnitionProb <- sim$fireSense_IgnitionPredicted
    }
  }

  ## calculate fire sizes from CWFIS polygon data
  if (!suppliedElsewhere("fireSize", sim)) {
    fireRaster <- Cache(prepInputsFireYear,
                        url = extractURL("fireSize"),
                        destinationPath = dPath,
                        rasterToMatch = sim$rasterToMatch,
                        maskWithRTM = TRUE,
                        method = "ngb",
                        datatype = "INT2U",
                        filename2 = NULL,
                        fireField = "CFS_REF_ID",
                        earliestYear = 1,
                        fun = "sf::st_read",
                        userTags = cacheTags,
                        omitArgs = c("destinationPath", "targetFile", "userTags"))

    sim$fireSize <- max(table(fireRaster[]))
  }

  return(invisible(sim))
}
