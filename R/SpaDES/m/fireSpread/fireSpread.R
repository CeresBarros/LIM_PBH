# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireSpread",
  description = "Fire spread model using Favier (2004) percolation model, where percolation probabilities are
  conditioned on vegetation, climate and topography conditions via the Canadian Forest Fire Behaviour System",
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
    defineParameter("vegFeedback", "logical", TRUE, NA, NA, desc = "Should vegetation feedbacks unto fire be simulated? Defaults to TRUE"),
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
    expectsInput(objectName = "aspectRas", objectClass = "RasterLayer",
                 desc = "Raster of aspect values - needs to be previously downloaded at this point"),
    expectsInput(objectName = "FuelTypes", objectClass = "data.table",
                 desc = "Table of Fuel Type parameters, with  base fuel type, species (in LANDIS code), their - or + contribution ('negSwitch'),
                 min and max age for each species"),
    expectsInput(objectName = "FWIinit", objectClass = "data.frame",
                 desc = "Initalisation parameter values for FWI calculations. Defaults to default values in cffdrs::fwi.
                 This table should be updated every year"),
    expectsInput(objectName = "pixelFuelTypes", objectClass = "data.table",
                 desc = "Fuel types per pixel group, calculated from cohort biomasses"),
    expectsInput(objectName = "pixelGroupMap", objectClass = "RasterLayer",
                 desc = "updated community map at each succession time step"),
    expectsInput(objectName = "precipitationRas", objectClass = "RasterLayer",
                 desc = "Raster of precipitation values",
                 sourceURL = "http://biogeo.ucdavis.edu/data/worldclim/v2.0/tif/base/wc2.0_2.5m_prec.zip"),
    expectsInput(objectName = "simulatedBiomassMap", objectClass = "RasterLayer",
                 desc = "Biomass map at each succession time step. Default is Canada national biomass map",
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-StructureBiomass.tar"),
    expectsInput(objectName = "slopeRas", objectClass = "RasterLayer",
                 desc = "Raster of slope values - needs to be previously downloaded at this point"),
    expectsInput("studyArea", "SpatialPolygonsDataFrame",
                 desc = paste("multipolygon to use as the study area,",
                              "with attribute LTHFC describing the fire return interval.",
                              "Defaults to a square shapefile in Southwestern Alberta, Canada."),
                 sourceURL = ""),
    expectsInput(objectName = "studyAreaFBP", objectClass = "SpatialPolygonsDataFrame",
                 desc = "same as studyArea,  but on FBP-compatible projection", sourceURL = ""),
    expectsInput(objectName = "temperatureRas", objectClass = "RasterLayer",
                 desc = "Raster of temperature values",
                 sourceURL = "http://biogeo.ucdavis.edu/data/worldclim/v2.0/tif/base/wc2.0_2.5m_tmax.zip")
  ),
  outputObjects = bind_rows(
    # createsOutput(objectName = "FBPinputs", objectClass = "RasterLayer",
    #               desc = "Fire behaviour prediction system inputs table"),
    # createsOutput(objectName = "FBPoutputs", objectClass = "list",
    #               desc = "Fire weather outputs table"),
    createsOutput(objectName = "aspectRas", objectClass = "RasterLayer",
                  desc = "Raster of aspect values - reprojected/cropped"),
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
    # createsOutput(objectName = "FWIinputs", objectClass = "RasterLayer",
    #               desc = "Fire weather inputs table"),
    # createsOutput(objectName = "FWIoutputs", objectClass = "list",
    #               desc = "Fire weather outputs table"),
    createsOutput(objectName = "pixelGroupMapFBP", objectClass = "RasterLayer",
                  desc = "updated community map at each succession time step, on FBP-compatible projection"),
    createsOutput(objectName = "precipitationRas", objectClass = "RasterLayer",
                  desc = "Raster of precipitation values - reprojected/cropped"),
    createsOutput(objectName = "rstCurrentBurn", objectClass = "RasterLayer",
                  desc = "Binary raster of fire spread"),
    createsOutput(objectName = "slopeRas", objectClass = "RasterLayer",
                  desc = "Raster of slope values - reprojected/cropped"),
    createsOutput(objectName = "startPix", objectClass = "vector",
                  desc = "List of starting fire pixels"),
    createsOutput(objectName = "temperatureRas", objectClass = "RasterLayer",
                  desc = "Raster of temperature values - reprojected/cropped"),
    createsOutput(objectName = "topoClimData", objectClass = "data.table",
                  desc = "Climate data table with temperature, precipitation and relative humidity for each pixelGroup")
  )
))

doEvent.fireSpread = function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      ## Initialise module
      sim <- fireInit(sim)

      ## schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$fireInitialTime, "fireSpread",
                           "fireParams", eventPriority = 2.25) ## always calculate fire parameters before the first fire time
      sim <- scheduleEvent(sim, P(sim)$fireInitialTime, "fireSpread",
                           "fireSpread", eventPriority = 2.5) ## always schedule fire
    },
    fireParams = {
      ## in the first year of fire always calculate parameters
      if(time(sim) == P(sim)$fireInitialTime) {
        ## calculate fire parameters
        sim <- FPBPercParams(sim)
      }

      ## in subsequent years evaluate if parameters are to be calculated again (veg feedbacks = TRUE)
      if(time(sim) == sim$fireYear) {
        if(P(sim)$vegFeedback) {
          ## calculate fire parameters
          sim <- FPBPercParams(sim)

          ## schedule future event(s)
          ## only calculate parameters in fire years.
          sim <- scheduleEvent(sim, time(sim) + P(sim)$fireTimestep, "fireSpread", "fireParams", eventPriority = 2.25)
        }
      }

    },
    fireSpread = {
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
                           "fireSpread", eventPriority = 2.5)
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

### module initialization
fireInit <- function(sim) {
  ## check package is installed
  if (!"cffdrs" %in% installed.packages())
    stop(paste("Please install the cffdrs R package to use",
               currentModule(sim)))

  cacheTags <- c("fireSpread", "fireInit")

  ## define first fire year
  sim$fireYear <- as.integer(P(sim)$fireInitialTime)

  ## project all inputs to Lat/Long (decimal degrees)
  ## for compatibility with FBP system

  ## buffer pixelGroupMap to prevent data loss, using twice the pixel size and distance
  ## then reproject to FBP compatible projection and mask buffer out (later)
  ## note: don't mask to area until the end.

  pixelGroupMapFBP <- buffer(sim$pixelGroupMap,
                             width = res(sim$pixelGroupMap)[1]*2)
  pixelGroupMapFBP <- projectRaster(pixelGroupMapFBP,
                                    crs = crs(sim$studyAreaFBP))

  ## PROJECT CLIMATE/TOPO RASTERS
  message(blue("Processing climate data for fire weather and fuel calculation"))
  sim$temperatureRas <- Cache(postProcess,
                              x = sim$temperatureRas,
                              rasterToMatch = pixelGroupMapFBP,
                              maskWithRTM = TRUE,
                              method = "bilinear",
                              filename2 = NULL,
                              userTags = c(cacheTags, "topoClimRas"),
                              omitArgs = c("userTags"))
  sim$precipitationRas <- Cache(postProcess,
                                x = sim$precipitationRas,
                                rasterToMatch = pixelGroupMapFBP,
                                maskWithRTM = TRUE,
                                method = "bilinear",
                                filename2 = NULL,
                                userTags = c(cacheTags, "topoClimRas"),
                                omitArgs = c("userTags"))
  sim$slopeRas <- Cache(postProcess,
                        x = sim$slopeRas,
                        rasterToMatch = pixelGroupMapFBP,
                        maskWithRTM = TRUE,
                        method = "bilinear",
                        filename2 = NULL,
                        userTags = c(cacheTags, "topoClimRas"),
                        omitArgs = c("userTags"))
  sim$aspectRas <- Cache(postProcess,
                         x = sim$aspectRas,
                         rasterToMatch = pixelGroupMapFBP,
                         maskWithRTM = TRUE,
                         method = "bilinear",
                         filename2 = NULL,
                         userTags = c(cacheTags, "topoClimRas"),
                         omitArgs = c("userTags"))

  ## TOPOCLIMDATA TABLE ----------------------
  topoClimData <- data.table(ID = 1:length(pixelGroupMapFBP),
                             pixelGroup = getValues(pixelGroupMapFBP),
                             temp = getValues(sim$temperatureRas), precip = getValues(sim$precipitationRas),
                             slope = getValues(sim$slopeRas), aspect = getValues(sim$aspectRas),
                             lat = coordinates(pixelGroupMapFBP)[,2],
                             long = coordinates(pixelGroupMapFBP)[,1])
  ## relative humidity
  ## using dew point between -3 and 20%, quarterly seasonal for Jun 2013
  ## https://calgary.weatherstats.ca/metrics/dew_point.html
  topoClimData[, relHum := RH(t = topoClimData$temp, Td = runif(nrow(topoClimData), -3, 20), isK = FALSE)]

  ## export to sim
  sim$pixelGroupMapFBP <- pixelGroupMapFBP
  sim$topoClimData <- topoClimData

  return(invisible(sim))
}

## Derive fire parameters from FBP system - rasters need to be in lat/long
FPBPercParams <- function(sim) {
  cacheTags <- c("fireSpread", "FBPPercParams")

  ## redo only if pixelGroupMap has been updated
  if (time(sim) != start(sim)) {
    pixelGroupMapFBP <- buffer(sim$pixelGroupMap,
                               width = res(sim$pixelGroupMap)[1]*2)
    pixelGroupMapFBP <- projectRaster(pixelGroupMapFBP,
                                      crs = crs(sim$studyAreaFBP))
    ## export to sim and clean ws
    sim$pixelGroupMapFBP <- pixelGroupMapFBP
    rm(pixelGroupMapFBP); amc::.gc()
  }

  ## FUEL TYPES ------------------------------
  ## rasterize fuel types table
  fuelTypesMaps <- rasterizeReduced(sim$pixelFuelTypes, sim$pixelGroupMap,
                                    newRasterCols = c("finalFuelType" , "coniferDom"),
                                    mapcode = "pixelGroup")

  ## now reproject to FBP-compatible crs
  fuelTypeRas <- Cache(postProcess,
                       x = fuelTypesMaps$finalFuelType,
                       rasterToMatch = sim$pixelGroupMapFBP,
                       maskWithRTM = TRUE,
                       method = "ngb",
                       filename2 = NULL,
                       userTags = c(cacheTags, "fuelTypeRas"),
                       omitArgs = c("userTags"))

  coniferDomRas <- Cache(postProcess,
                         x = fuelTypesMaps$coniferDom,
                         rasterToMatch = sim$pixelGroupMapFBP,
                         maskWithRTM = TRUE,
                         method = "bilinear",
                         filename2 = NULL,
                         userTags = c(cacheTags, "coniferDomRas"),
                         omitArgs = c("userTags"))
  ## make table of final fuel types
  FTs <- data.table(ID = 1:length(sim$pixelGroupMapFBP),
                    FuelType = getValues(fuelTypeRas),
                    coniferDom = getValues(coniferDomRas))
  ## add FBP fuel type names
  FTs <- unique(sim$FuelTypes[, .(FuelTypeFBP, FuelType)])[FTs, on = "FuelType"]
  FTs <- FTs[!is.na(FuelType)]

  ## check for duplicates (there shouldn't be any)
  if (any(duplicated(FTs))) {
    stop("Duplicated pixels found in fuel types table.",
         " Please debug fireSpread FPBPercParams() event function")
  }

  ## FWI ------------------------------
  ## make/update table of FWI inputs
  FWIinputs <- data.frame(id = sim$topoClimData$ID,
                          lat = sim$topoClimData$lat,
                          long = sim$topoClimData$long,
                          mon = 7,
                          temp = sim$topoClimData$temp,
                          rh = sim$topoClimData$relHum,
                          ws = 0,
                          prec = sim$topoClimData$precip)

  ## calculate FW indices
  FWIoutputs <- suppressWarnings({
    cffdrs::fwi(input = na.omit(FWIinputs),
                init = na.omit(sim$FWIinit),
                batch = FALSE, lat.adjust = TRUE)
  })
  FWIoutputs <- data.table(FWIoutputs)

  ## FBP -----------------------------
  ## make inputs dataframe for FBI
  ## add fuel types and conifer dominance to FWIOutputs
  ## note that because climate/topo data is "larger" there are pixels that have no fuels - these are removed.
  FWIoutputs <- FTs[FWIoutputs, on = "ID", nomatch = 0]

  ## add slope and aspect
  ## again, only keep pixels that have fuels
  FWIoutputs <- sim$topoClimData[, .(ID, slope, aspect)][FWIoutputs, on = "ID", nomatch = 0]

  FBPinputs <- data.frame(id = FWIoutputs$ID,
                          FuelType = FWIoutputs$FuelTypeFBP,
                          LAT = FWIoutputs$LAT,
                          LONG = FWIoutputs$LONG,
                          FFMC = FWIoutputs$FFMC,
                          BUI = FWIoutputs$BUI,
                          WS = FWIoutputs$WS,
                          GS = FWIoutputs$slope,
                          Dj = rep(180, nrow(FWIoutputs)),
                          Aspect = FWIoutputs$aspect,
                          PC = FWIoutputs$coniferDom)

  FBPoutputs <- suppressWarnings({
    cffdrs::fbp(input = na.omit(FBPinputs), output = "All")
  })
  FBPoutputs <- data.table(FBPoutputs)

  ## FBP OUTPUTS TO SPATIALPOINTS
  FBPOutputsPts <- FBPoutputs[, .(ID, CFB, ROS, RSO, HFI, TFC)]
  FBPOutputsPts <- FBPOutputsPts[FWIoutputs[, .(ID, LAT, LONG)],
                                 on = "ID", nomatch = 0][, ID := NULL]

  FBPOutputsSf <- st_as_sf(FBPOutputsPts, coords = c("LONG", "LAT"),
                           crs =  as.character(crs(sim$studyAreaFBP)),
                           agr = "constant")

  ## REPROJECT TO ORIGINAL CRS without data loss and convert to raster
  ## note that after rasterizing it is safer to mask, in case some of
  ## the buffer artefacts come through
  FBPOutputsSf <- st_transform(FBPOutputsSf,
                               crs = as.character(crs(sim$pixelGroupMap)))
  FBPOutputsPoly <- as_Spatial(FBPOutputsSf)

  ## Crown fraction burnt
  sim$fireCFBRas <- rasterize(FBPOutputsPoly, sim$pixelGroupMap,
                              field = "CFB", fun = function(x, na.rm = TRUE) max(x))
  sim$fireCFBRas <- mask(sim$fireCFBRas, sim$pixelGroupMap)

  ## Head fire intensity
  sim$fireIntRas <- rasterize(FBPOutputsPoly, sim$pixelGroupMap,
                              field = "HFI", fun = function(x, na.rm = TRUE) max(x))
  sim$fireIntRas <- mask(sim$fireIntRas, sim$pixelGroupMap)

  ## Rate of spread
  sim$fireROSRas <- rasterize(FBPOutputsPoly, sim$pixelGroupMap,
                              field = "ROS", fun = function(x, na.rm = TRUE) max(x))
  sim$fireROSRas <- mask(sim$fireROSRas, sim$pixelGroupMap)

  ## Critical spread rate for crowning
  sim$fireRSORas <- rasterize(FBPOutputsPoly, sim$pixelGroupMap,
                              field = "RSO", fun = function(x, na.rm = TRUE) max(x))
  sim$fireRSORas <- mask(sim$fireRSORas, sim$pixelGroupMap)

  ## Total fuel consumption
  sim$fireTFCRas <- rasterize(FBPOutputsPoly, sim$pixelGroupMap,
                              field = "TFC", fun = function(x, na.rm = TRUE) max(x))
  sim$fireTFCRas <- mask(sim$fireTFCRas, sim$pixelGroupMap)

  ## export to sim
  # sim$FWIinputs <- FWIinputs
  # sim$FWIoutputs <- FWIoutputs
  # sim$FBPinputs <- FBPinputs
  # sim$FBPoutputs <- FBPoutputs

  return(invisible(sim))
}

## Fire spread event in fire years - rasters should be back in LBMR projection
doFireSpread <- function(sim) {

  ## if no simulatedBiomassMap is supplied then generate one from raw data
  ## at the start
  if (time(sim) == P(sim)$fireInitialTime) {
    if (!suppliedElsewhere("simulatedBiomassMap", sim)) {
      if (!suppliedElsewhere("rawBiomassMap", sim)) {
        # Filenames
        rawBiomassMapFilename <- file.path(dPath, "NFI_MODIS250m_kNN_Structure_Biomass_TotalLiveAboveGround_v0.tif")
        rawBiomassMapURL <- "http://tree.pfc.forestry.ca/kNN-StructureBiomass.tar"
        sim$simulatedBiomassMap <- Cache(prepInputs,
                                         targetFile = asPath(basename(rawBiomassMapFilename)),
                                         archive = asPath(c("kNN-StructureBiomass.tar",
                                                            "NFI_MODIS250m_kNN_Structure_Biomass_TotalLiveAboveGround_v0.zip")),
                                         url = rawBiomassMapURL,
                                         destinationPath = dPath,
                                         studyArea = sim$studyArea,
                                         rasterToMatch = NULL,
                                         maskWithRTM = FALSE,
                                         useSAcrs = TRUE,
                                         method = "bilinear",
                                         datatype = "INT2U",
                                         filename2 = TRUE, overwrite = TRUE,
                                         omitArgs = c("destinationPath", "targetFile", cacheTags, "stable"))
      } else {
        sim$simulatedBiomassMap <- sim$rawBiomassMap
      }
    }
  }

  ## MAKE BURNABLE AREAS RASTER -------------------------------
  ## only areas with biomass can burn
  burnableAreas <- sim$simulatedBiomassMap
  vals <- data.table(B = getValues(sim$simulatedBiomassMap))   ## making a mask is probably faster with data.table
  vals <- vals[B > 0, B := 1]
  vals <- vals[B <= 0, B := NA]
  burnableAreas[] <- vals$B

  ## MAKE RASTER OF SPREAD PROABILITIES
  ## spread probability is the combination of ROS and intensity, which have an additive effect
  ## and their sum is scaled to 0-0.23
  ## TODO: the scaling should guarantee an average value of 0.23
  ## TODO: ROS and intensity should be combined differently
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

  ## remove pixels that didn't burn from ROS, Int and TFC rasters
  sim$fireCFBRas <- mask(sim$fireIntRas, rstCurrentBurn)
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
  dPath <- dataPath(sim)
  cacheTags = c(currentModule(sim), "function:.inputObjects")

  ## make raster storage lists
  sim$rstCurrentBurn <- list()

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

  ## DEFAULT TOPO, TEMPERATURE AND PRECIPITATION
  ## these defaults are only necessary if the rasters are not supplied by another module
  ## climate defaults to http://worldclim.org/version2 (temperarute historical ensemble)
  ## slope/aspect defaults obtained from Environment Canada using foothills shp
  if(!suppliedElsewhere("temperatureRas", sim)){
    ## get default temperature values
    fnames <- paste0("wc2.0_2.5m_tmax_0", 5:9, ".tif")
    temperatureStk <- list()
    for(f in fnames) {
      temperatureStk[[f]] <- Cache(prepInputs, targetFile = f,
                                   url = extractURL("temperatureRas", sim),
                                   archive = "wc2.0_2.5m_tmax.zip",
                                   alsoExtract = NA,
                                   destinationPath = getPaths()$inputPath,
                                   datatype = "FLT4S",
                                   filename2 = FALSE,
                                   userTags = cacheTags)
    }
    temperatureStk <- stack(temperatureStk)
    temperatureRas <- raster::mean(temperatureStk)
    sim$temperatureRas <- temperatureRas
  }

  if(!suppliedElsewhere("precipitationRas", sim)){
    ## get default precipitation values
    fnames <- paste0("wc2.0_2.5m_prec_0", 5:9, ".tif")
    precipitationStk <- list()
    for(f in fnames) {
      precipitationStk[[f]] <- Cache(prepInputs, targetFile = f,
                                     url = extractURL("precipitationRas", sim),
                                     archive = "wc2.0_2.5m_prec.zip",
                                     alsoExtract = NA,
                                     destinationPath = getPaths()$inputPath,
                                     datatype = "FLT4S",
                                     filename2 = FALSE,
                                     userTags = cacheTags)
    }
    precipitationStk <- stack(precipitationStk)
    precipitationRas <- raster::mean(precipitationStk)
    sim$precipitationRas <- precipitationRas
  }

  if(!suppliedElsewhere("slopeRas", sim)){
    ## TODO defaults of slope/aspect should cover whole of Canada
    ## get default slope values
    slopeRas <- Cache(prepInputs, targetFile = "dataset/SLOPE.tif",
                      archive = "DEM_Foothills_study_area.zip",
                      alsoExtract = NA,
                      destinationPath = getPaths()$inputPath,
                      datatype = "FLT4S",
                      filename2 = FALSE,
                      userTags = cacheTags)
    sim$slopeRas <- slopeRas

  }

  if(!suppliedElsewhere("aspectRas", sim)){
    ## get default aspect values
    aspectRas <- Cache(prepInputs, targetFile = "dataset/ASPECT.tif",
                       archive = "DEM_Foothills_study_area.zip",
                       alsoExtract = NA,
                       destinationPath = getPaths()$inputPath,
                       datatype = "FLT4S",
                       filename2 = FALSE,
                       userTags = cacheTags)
    sim$aspectRas <- aspectRas
  }

  ## FWI INITIALISATION DATAFRAME
  ## TODO:FWIinit should be updated every year from previous year's/days/months results
  if(!suppliedElsewhere("FWIinit", sim)) {
    sim$FWIinit = data.frame(ffmc = 85,
                             dmc = 6,
                             dc = 15)
  }

  return(invisible(sim))
}
