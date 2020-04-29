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
  reqdPkgs = list("raster", "data.table", "dplyr", "humidity", #"cffdrs",  ## cffdrs is causing installation problems
                  "PredictiveEcology/SpaDES.core@development",
                  "PredictiveEcology/SpaDES.tools@development",
                  "CeresBarros/reproducible@development"),
  parameters = rbind(
    defineParameter(".crsUsed", "character", "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0",
                    NA, NA, "CRS to be used. Defaults to the biomassMap projection"),
    defineParameter("fireSize", "integer", 1000L, NA, NA, desc = "Fire size in pixels"),
    defineParameter("vegFeedback", "logical", TRUE, NA, NA, desc = "Should vegetation feedbacks unto fire be simulated? Defaults to TRUE"),
    defineParameter("noStartPix", "integer", 100L, NA, NA, desc = "Number of fire events"),
    defineParameter("fireStart", "integer", 2L, NA, NA, desc = "First fire year. Defaults to the 2nd year of the simulation"),
    defineParameter("fireFreq", "integer", 1L, NA, NA, desc = "Fire recurrence in years"),
    defineParameter(".plotMaps", "logical", FALSE, NA, NA, "This describes whether maps should be plotted or not"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA, "This describes the simulation time interval between plot events")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "shpStudySubRegionFBP", objectClass = "SpatialPolygonsDataFrame",
                 desc = "this shape file contains two informaton: Sub study area with fire return interval attribute.
                 Defaults to a square shapefile in Southwestern Alberta, Canada", sourceURL = ""),
    expectsInput(objectName = "shpStudyRegionFullFBP", objectClass = "SpatialPolygonsDataFrame",
                 desc = "this shape file contains two informaton: Full study area with fire return interval attribute.
                 Defaults to a square shapefile in Southwestern Alberta, Canada", sourceURL = ""),
    expectsInput(objectName = "biomassMap", objectClass = "RasterLayer",
                 desc = paste("Biomass map at each succession time step. Defaults to the Canadian Forestry",
                              "Service, National Forest Inventory, kNN-derived total aboveground biomass map",
                              "from 2001. See https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990",
                              "for metadata"),
                 sourceURL = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                    "canada-forests-attributes_attributs-forests-canada/2001-attributes_attributs-2001/")),
    expectsInput(objectName = "topoClimData", objectClass = "data.table",
                 desc = "Climate data table with temperature, precipitation and relative humidity for each pixelGroup"),
    expectsInput(objectName = "FWIinit", objectClass = "data.frame",
                 desc = "Initalisation parameter values for FWI calculations. Defaults to default values in cffdrs::fwi.
                 This table should be updated every year"),
    expectsInput(objectName = "pixelFuelTypes", objectClass = "data.table",
                 desc = "Fuel types per pixel group, calculated from cohort biomasses"),
    expectsInput(objectName = "pixelGroupMapFBP", objectClass = "RasterLayer",
                 desc = "updated community map at each succession time step, on FBP-compatible projection"),
    expectsInput(objectName = "FuelTypes", objectClass = "data.table",
                 desc = "Table of Fuel Type parameters, with  base fuel type, species (in LANDIS code), their - or + contribution ('negSwitch'),
                 min and max age for each species"),
    # expectsInput(objectName = "FWIinit", objectClass = "data.table",
    #              desc = "Table of Fire Weather Index initalisation parameters, defaults to default values
    #              available in the cffrs::fwi documentation"),
    expectsInput(objectName = "pixelGroupMap", objectClass = "RasterLayer",
                 desc = "updated community map at each succession time step")
    ),
  outputObjects = bind_rows(
    createsOutput(objectName = "biomassMapFBP", objectClass = "RasterLayer",
                 desc = "Biomass map at each succession time step, on FBP-compatible projection.
                 Default is Canada national biomass map"),
    createsOutput(objectName = "topoClimData", objectClass = "data.table",
                  desc = "Climate data table with temperature, precipitation and relative humidity for each pixelGroup"),
    createsOutput(objectName = "fireYear", objectClass = "numeric", desc = "Next fire year"),
    createsOutput(objectName = "pixelGroupMapFBP", objectClass = "RasterLayer",
                 desc = "updated community map at each succession time step, on FBP-compatible projection"),
    createsOutput(objectName = "pixelGroupMap", objectClass = "RasterLayer",
                  desc = "updated community map at each succession time step"),
    createsOutput(objectName = "FWIinputs", objectClass = "RasterLayer",
                  desc = "Fire weather inputs table"),
    createsOutput(objectName = "FWIinit", objectClass = "RasterLayer",
                  desc = "Fire weather initialisation table, updated for the following fire year"),
    createsOutput(objectName = "FWIoutputs", objectClass = "list",
                  desc = "Fire weather outputs table"),
    createsOutput(objectName = "FBPinputs", objectClass = "RasterLayer",
                  desc = "Fire behaviour prediction system inputs table"),
    createsOutput(objectName = "FBPoutputs", objectClass = "list",
                  desc = "Fire weather outputs table"),
    createsOutput(objectName = "fireROSRas", objectClass = "RasterLayer",
                  desc = "Raster of equilibrium rate of spread"),
    createsOutput(objectName = "fireIntRas", objectClass = "RasterLayer",
                  desc = "Raster of equilibrium head fire intensity"),
    createsOutput(objectName = "fireTFCRas", objectClass = "RasterLayer",
                  desc = "Raster of total fuel consumed"),
    createsOutput(objectName = "startPix", objectClass = "vector",
                  desc = "List of starting fire pixels"),
    createsOutput(objectName = "fireSpreadRas", objectClass = "list",
                  desc = "List of rasters of fire spread")
  )
))

doEvent.fireSpread = function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      ## Initialise module
      sim <- fireInit(sim)

      ## schedule future event(s)
      sim <- scheduleEvent(sim, time(sim) + 1, "fireSpread", "Fire", eventPriority = 2) ## always schedule fire
    },
    Fire = {
      ## in the first year of fire calculate parameters and do fire
      if(time(sim) == P(sim)$fireStart) {
        ## calculate fire parameters
        sim <- FPBPercParams(sim)
        ## calculate fire spread
        sim <- doFireSpread(sim)

      } else{
        ## in subsequent years evaluate if parameters are to be aclcualted again (veg feedbacks)
        if(time(sim) == sim$fireYear) {
          if(P(sim)$vegFeedback) {
            ## calculate fire parameters
            sim <- FPBPercParams(sim)
            ## calculate fire spread
            sim <- doFireSpread(sim)
          } else {
            ## calculate fire spread
            sim <- doFireSpread(sim)
          }
        } else {
          ## No fire
          sim <- doNoFire(sim)
        }
      }

      ## define next fire year
      sim$fireYear <- time(sim) + P(sim)$fireFreq

      ## schedule future event(s) - always schedule fire
      sim <- scheduleEvent(sim, time(sim) + 1, "fireSpread", "Fire", eventPriority = 2)
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

### module initialization
fireInit <- function(sim) {
  ## project to Lat/Long (decimal degrees) for compatibility with FBP system
  ## TODO: this results in data loss - but LandR doesn't deal well with lat/long
  ## need to find long term solution
  latLong = "+proj=longlat +datum=WGS84"

  ## define first fire year
  sim$fireYear <- as.integer(P(sim)$fireStart)

  ## get pixelGroupMapFBP and reproject
  pixelGroupMapFBP <- sim$pixelGroupMap
  pixelGroupMapFBP <- postProcess(pixelGroupMapFBP, rasterToMatch = sim$biomassMapFBP,
                                  method = "ngb", filename2 = NULL)

  ## get biomassMap and reproject to previous biomassMapFBP
  biomassMapFBP <- sim$biomassMap
  biomassMapFBP <- postProcess(biomassMapFBP, rasterToMatch = sim$biomassMapFBP,
                               method = "bilinear", filename2 = NULL)


  ## TOPOCLIMDATA TABLE - only needs to be created if not supplied by another module
  ## can't be put in .inputObjects without a default picelGroupMap
  if(!suppliedElsewhere('topoClimData', sim)) {
    topoClimData <- data.table(pixID = 1:length(pixelGroupMapFBP),
                               pixelGroup = getValues(pixelGroupMapFBP),
                               temp = getValues(sim$temperatureRas), precip = getValues(sim$precipitationRas),
                               slope = getValues(sim$slopeRas), aspect = getValues(sim$aspectRas),
                               lat =  coordinates(sim$temperatureRas)[,2],
                               long = coordinates(sim$temperatureRas)[,1])
    ## relative humidity
    ## using dew point between -3 and 20%, quarterly seasonal for Jun 2013
    ## https://calgary.weatherstats.ca/metrics/dew_point.html
    topoClimData[, relHum := RH(t = topoClimData$temp, Td = runif(nrow(topoClimData), -3, 20), isK = FALSE)]
    sim$topoClimData <- topoClimData
  }

  ## export to sim
  sim$pixelGroupMapFBP <- pixelGroupMapFBP
  sim$biomassMapFBP <- biomassMapFBP

  return(invisible(sim))
}

## Derive fire parameters from FBP system - rasters need to be in lat/long
FPBPercParams <- function(sim) {
  require(cffdrs)   ## cffdrs was causing issues when included in metadata

  ## Update pixelGroupMap and biomassMap if not init
  if(time(sim) != start(sim)) {
    ## get pixelGroupMapFBP and reproject
    pixelGroupMapFBP <- sim$pixelGroupMap
    pixelGroupMapFBP <- postProcess(pixelGroupMapFBP, rasterToMatch = sim$biomassMapFBP,
                                    method = "ngb", filename2 = NULL)

    ## get biomassMap and reproject to previous biomassMapFBP
    biomassMapFBP <- sim$biomassMap
    biomassMapFBP <- postProcess(biomassMapFBP, rasterToMatch = sim$biomassMapFBP,
                                 method = "bilinear", filename2 = NULL)

    ## export to sim and clean ws
    sim$pixelGroupMapFBP <- pixelGroupMapFBP
    sim$biomassMapFBP <- biomassMapFBP
    rm(pixelGroupMapFBP, biomassMapFBP)
  }

  ## make/update table of FWI inputs
  FWIinputs <- data.frame(id = sim$topoClimData$pixID,
                          lat = sim$topoClimData$lat,
                          long = sim$topoClimData$long,
                          mon = 7,
                          temp = sim$topoClimData$temp,
                          rh = sim$topoClimData$relHum,
                          ws = 0,
                          prec = sim$topoClimData$precip)

  ## calculate FW indices
  FWIoutputs <- fwi(input = na.omit(FWIinputs),
                    init = na.omit(sim$FWIinit),
                    batch = FALSE, lat.adjust = TRUE) %>%
    data.table

  ## add pixelGroup
  FWIoutputs  <- sim$topoClimData[,.(pixID, pixelGroup)] %>%   ## this is a right join (right table being FWIoutputs)
    .[FWIoutputs , , on = "pixID==ID"]


  ## make table of final fuel types
  FTs <- data.table(pixelGroup = sim$pixelFuelTypes$pixelGroup,
                    FuelType = sim$pixelFuelTypes$finalFuelType,
                    coniferDom = sim$pixelFuelTypes$coniferDom)

  ## merge with FBP fuel type names
  FTs <- sim$FuelTypes[, .(FuelTypeFBP,FuelType)] %>%
    .[!duplicated(.)] %>%
    .[FTs, on = "FuelType"]

  ## add fuel types to FWIOutputs
  FWIoutputs <- FTs[!duplicated(FTs)] %>%
    .[FWIoutputs, on = "pixelGroup"]

  ## add slope and aspect
  FWIoutputs <- sim$topoClimData[, .(pixID, slope, aspect)] %>%
    .[FWIoutputs, on = "pixID"]

  ## make inputs dataframe for FBI
  FBPinputs <- data.frame(id = FWIoutputs$pixID,
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

  FBPoutputs <- fbp(input = na.omit(FBPinputs)) %>%
    data.table

  ## add pixelGroup
  FBPoutputs <- FWIoutputs[,.(pixID, pixelGroup)] %>%   ## this is a right join (right table being FWIoutputs)
    .[FBPoutputs, , on = "pixID==ID"]

  ## FBP OUTPUTS TO RASTERS
  sim$fireROSRas = sim$fireIntRas =  sim$fireTFCRas = sim$pixelGroupMapFBP

  ## Rate of spread
  ROSvals <- as.matrix(data.frame(is = FBPoutputs$pixelGroup, becomes = FBPoutputs$ROS))
  sim$fireROSRas[!getValues(sim$fireROSRas) %in% FBPoutputs$pixelGroup] <- NA
  sim$fireROSRas <- raster::reclassify(sim$fireROSRas, rcl = ROSvals)

  # Head fire intensity
  HFIvals <- as.matrix(data.frame(is = FBPoutputs$pixelGroup, becomes = FBPoutputs$HFI))
  sim$fireIntRas[!getValues(sim$fireIntRas) %in% FBPoutputs$pixelGroup] <- NA
  sim$fireIntRas <- raster::reclassify(sim$fireIntRas, rcl = HFIvals)

  ## Total fuel consumption
  TFCvals <- as.matrix(data.frame(is = FBPoutputs$pixelGroup, becomes = FBPoutputs$TFC))
  sim$fireTFCRas[!getValues(sim$fireTFCRas) %in% FBPoutputs$pixelGroup] <- NA
  sim$fireTFCRas <- raster::reclassify(sim$fireTFCRas, rcl = TFCvals)

  ## transform pixel IDs in tables to Biomass_core compatible
  tempRas <- sim$pixelGroupMapFBP
  tempRas[] <- 1:ncell(tempRas)
  tempRas <- postProcess(tempRas, rasterToMatch = sim$pixelGroupMap,
                         res = res(sim$pixelGroupMap), method = "ngb", filename2 = NULL,
                         useCache = FALSE)
  pixCorresp <- data.table(pixFBP = getValues(tempRas), pixelIndex = 1:ncell(tempRas))

  FWIinputs  <- merge(FWIinputs, pixCorresp, by.x = "id", by.y = "pixFBP", all.y = TRUE)
  FWIoutputs <- merge(FWIoutputs, pixCorresp, by.x = "pixID", by.y = "pixFBP", all.y = TRUE)
  FBPinputs <- merge(FBPinputs, pixCorresp, by.x = "id", by.y = "pixFBP", all.y = TRUE)
  FBPoutputs <- merge(FBPoutputs, pixCorresp, by.x = "pixID", by.y = "pixFBP", all.y = TRUE)

  FWIinputs$id <- NULL; FBPinputs$id <- NULL
  FWIoutputs[, pixID := NULL]; FBPoutputs[, pixID := NULL]

  ## reproject rasters and maps for Biomass_core compatibility
  sim$fireROSRas <- postProcess(sim$fireROSRas, rasterToMatch = sim$biomassMap, method = "bilinear", filename2 = NULL)
  sim$fireIntRas <- postProcess(sim$fireIntRas, rasterToMatch = sim$biomassMap, method = "bilinear", filename2 = NULL)
  sim$fireTFCRas <- postProcess(sim$fireTFCRas, rasterToMatch = sim$biomassMap, method = "bilinear", filename2 = NULL)

  ## export to sim
  sim$FWIinputs <- FWIinputs
  sim$FWIoutputs <- FWIoutputs
  sim$FBPinputs <- FBPinputs
  sim$FBPoutputs <- FBPoutputs

  return(invisible(sim))
}

## Fire spread event in fire years - rasters should be back in Biomass_core projection
doFireSpread <- function(sim) {
  ## MAKE BURNABLE AREAS RASTER -------------------------------
  ## only areas with biomass can burn
  burnableAreas <- sim$biomassMap
  vals <- data.table(B = getValues(sim$biomassMap))   ## making a mask is probably faster with data.table
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
  vals[!is.na(spreadP), spreadP := scales::rescale(spreadP, to = c(0.1,0.3))]
  spreadProb_map[] <- vals$spreadP

  ## MAKE RASTER OF PERSISTENCE PROABILITIES
  ## persistence probability is the combination of TFC and intensity, which have an additive effect
  ## and their sum is scaled to 0-1
  ## TODO: TFC and intensity should be combined differently
  persistProb_map <- sim$fireTFCRas + sim$fireIntRas
  persistProb_map <- mask(persistProb_map , burnableAreas)

  vals <- data.table(persisP = getValues(persistProb_map))   ## making a mask is probably faster with data.table
  vals[!is.na(persisP), persisP := scales::rescale(persisP, to = c(0,1))]
  persistProb_map[] <- vals$persisP

  ## MAKE RASTER OF FIRE SPREAD -------------------------------
  ## note that this function two random components: selection of starting pixels and fire spread
  sim$startPix <- sample(which(!is.na(getValues(burnableAreas))), P(sim)$noStartPix)

  ## Favier's model:
  x <- try(expr = {
    fireSpreadRas <- spread2(landscape = burnableAreas,
                             spreadProb = spreadProb_map,
                             persistProb = persistProb_map,
                             start = sim$startPix,
                             maxSize =  P(sim)$fireSize,
                             plot.it = FALSE)
  }, silent = TRUE)

  while (class(x) == "try-error") {
    x <- try(expr = {
      fireSpreadRas <- spread2(landscape = burnableAreas,
                               spreadProb = spreadProb_map,
                               persistProb = persistProb_map,
                               start = startPix,
                               maxSize =  P(sim)$fireSize,
                               plot.it = FALSE)
    }, silent = TRUE)
  }



  ## remove fires that spread beyond burnable areas
  fireSpreadRas <- mask(fireSpreadRas, burnableAreas)

  ## export to sim
  sim$burnableAreas
  sim$fireSpreadRas[[time(sim)]] <- fireSpreadRas
  return(invisible(sim))
}

## What to do in no fire years
doNoFire <- function(sim) {
  sim$fireSpreadRas[[time(sim)]] <- NA

  return(invisible(sim))
}

## OTHER INPUTS AND FUNCTIONS --------------------------------
.inputObjects = function(sim) {

  dPath <- dataPath(sim)
  cacheTags <- c(currentModule(sim), "function:.inputObjects")

  ## make raster storage lists
  sim$fireSpreadRas <- list()

  ## project to Lat/Long (decimal degrees) for compatibility with FBP system
  ## TODO: this results in data loss - but LandR doesn't deal well with lat/long
  ## need to find long term solution
  latLong <- "+proj=longlat +datum=WGS84"

  if(!suppliedElsewhere("shpStudyRegionFullFBP", sim)) {
    if (!suppliedElsewhere("shpStudyRegionFull", sim)) {
      message("'shpStudyRegionFull' was not provided by user. Using a polygon in Southwestern Alberta, Canada")

      canadaMap <- Cache(getData, 'GADM', country = 'CAN', level = 1, path = asPath(dPath),
                         cacheRepo = getPaths()$cachePath, quick = FALSE)
      smallPolygonCoords = list(coords = data.frame(x = c(-115.9022,-114.9815,-114.3677,-113.4470,-113.5084,-114.4291,-115.3498,-116.4547,-117.1298,-117.3140),
                                                    y = c(50.45516,50.45516,50.51654,50.51654,51.62139,52.72624,52.54210,52.48072,52.11243,51.25310)))

      sim$shpStudyRegionFull <- SpatialPolygons(list(Polygons(list(Polygon(smallPolygonCoords$coords)), ID = "swAB_polygon")),
                                                proj4string = crs(canadaMap))

      ## use CRS of biomassMap
      sim$shpStudyRegionFull <- spTransform(sim$shpStudyRegionFull, CRSobj = P(sim)$.crsUsed)

    }
    sim$shpStudyRegionFullFBP <- sim$shpStudyRegionFull
  }

  if (!suppliedElsewhere("shpStudySubRegionFBP", sim)) {
    if (!suppliedElsewhere("shpStudySubRegion", sim)) {
      message("'shpStudySubRegion' was not provided by user. Using the same as 'shpStudyRegionFull'")
      sim$shpStudySubRegion <- sim$shpStudyRegionFull
    }
    sim$shpStudySubRegionFBP <- sim$shpStudySubRegion
  }

  ## if necessary reproject to lat/long - for compatibility with FBP
  if (!compareCRS(latLong, crs(sim$shpStudyRegionFullFBP))) {
    sim$shpStudyRegionFullFBP <- spTransform(sim$shpStudyRegionFullFBP, latLong) #faster without Cache
  }

  if (!compareCRS(latLong, crs(sim$shpStudySubRegionFBP))) {
    sim$shpStudySubRegionFBP <- spTransform(sim$shpStudySubRegionFBP, latLong) #faster without Cache
  }

  if(!suppliedElsewhere("biomassMap", sim)) {
    sim$biomassMap <- Cache(prepInputs,
                            targetFile = biomassMapFilename,
                            archive = asPath(c("kNN-StructureBiomass.tar",
                                               "NFI_MODIS250m_kNN_Structure_Biomass_TotalLiveAboveGround_v0.zip")),
                            url = extractURL("biomassMap", sim),
                            destinationPath = dPath,
                            studyArea = sim$shpStudySubRegion,
                            useSAcrs = TRUE,
                            method = "bilinear",
                            datatype = "INT2U",
                            filename2 = TRUE,
                            userTags = cacheTags)
  }

  ## project biomassMap for next prepinput calls.
  biomassMapFBP <- sim$biomassMap
  biomassMapFBP <- postProcess(biomassMapFBP, rasterToMatch = sim$biomassMapFBP,
                               method = "bilinear", filename2 = NULL)

  ## DEFAULT TOPO, TEMPERATURE AND PRECIPITATION
  ## these defaults are only necessary if the topoClimData is not supplied by another module
  ## climate defaults to http://worldclim.org/version2 (temperarute historical ensemble)
  ## slope/aspect defaults obtained from Environment Canada using foothills shp
  if(!suppliedElsewhere("topoClimData", sim)){
    ## get default temperature values
    fnames <- paste0("wc2.0_2.5m_tmax_0", 5:9, ".tif")
    temperatureStk <- list()
    for(f in fnames) {
      temperatureStk[[f]] <- Cache(prepInputs, targetFile = f,
                                   url = "http://biogeo.ucdavis.edu/data/worldclim/v2.0/tif/base/wc2.0_2.5m_tmax.zip",
                                   archive = "wc2.0_2.5m_tmax.zip",
                                   alsoExtract = NA,
                                   destinationPath = getPaths()$inputPath,
                                   studyArea = sim$shpStudySubRegionFBP,
                                   rasterToMatch = biomassMapFBP,
                                   method = "bilinear",
                                   datatype = "FLT4S",
                                   filename2 = FALSE,
                                   userTags = cacheTags)
    }
    temperatureStk <- stack(temperatureStk)
    temperatureRas <- raster::mean(temperatureStk)

    ## get default precipitation values
    fnames <- paste0("wc2.0_2.5m_prec_0", 5:9, ".tif")
    precipitationStk <- list()
    for(f in fnames) {
      precipitationStk[[f]] <- Cache(prepInputs, targetFile = f,
                                     url = "http://biogeo.ucdavis.edu/data/worldclim/v2.0/tif/base/wc2.0_2.5m_prec.zip",
                                     archive = "wc2.0_2.5m_prec.zip",
                                     alsoExtract = NA,
                                     destinationPath = getPaths()$inputPath,
                                     studyArea = sim$shpStudySubRegionFBP,
                                     rasterToMatch = biomassMapFBP,
                                     method = "bilinear",
                                     datatype = "FLT4S",
                                     filename2 = FALSE,
                                     userTags = cacheTags)
    }
    precipitationStk <- stack(precipitationStk)
    precipitationRas <- raster::mean(precipitationStk)

    ## TODO defaults of slope/aspect should cover whole of Canada
    ## get default slope values
    slopeRas <- Cache(prepInputs, targetFile = "dataset/SLOPE.tif",
                      archive = "DEM_Foothills_study_area.zip",
                      alsoExtract = NA,
                      destinationPath = getPaths()$inputPath,
                      studyArea = sim$shpStudySubRegionFBP,
                      rasterToMatch =biomassMapFBP,
                      method = "bilinear",
                      datatype = "FLT4S",
                      filename2 = FALSE,
                      userTags = cacheTags)

    ## get default aspect values
    aspectRas <- Cache(prepInputs, targetFile = "dataset/ASPECT.tif",
                       archive = "DEM_Foothills_study_area.zip",
                       alsoExtract = NA,
                       destinationPath = getPaths()$inputPath,
                       studyArea = sim$shpStudySubRegionFBP,
                       rasterToMatch = biomassMapFBP,
                       method = "bilinear",
                       datatype = "FLT4S",
                       filename2 = FALSE,
                       userTags = cacheTags)

    ## make a copy to reproject for FBP and base file loads and export to sim
    sim$temperatureRas <- temperatureRas
    sim$precipitationRas <- precipitationRas
    sim$slopeRas <- slopeRas
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
