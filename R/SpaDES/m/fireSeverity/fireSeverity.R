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
    defineParameter(".crsUsed", "character", "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0",
                    NA, NA, "CRS to be used. Defaults to the biomassMap projection"),
    defineParameter(".plotMaps", "logical", FALSE, NA, NA, "This describes whether maps should be plotted or not"),
    defineParameter(name = ".saveInitialTime", class = "numeric", default = 0,
                    min = NA, max = NA, desc = "This describes the simulation time at which the
                    first save event should occur"),
    defineParameter(name = "fireTimestep", class = "numeric", default = 2L,
                    desc = "The number of time units between successive fire events in a fire module")
    ),
  inputObjects = bind_rows(
    expectsInput(objectName = "pixelGroupMap", objectClass = "RasterLayer",
                  desc = "updated community map at each succession time step"),
    expectsInput(objectName = "ecoregionMap", objectClass = "RasterLayer",
                 desc = "ecoregion map that has mapcodes match ecoregion table and speciesEcoregion table",
                 sourceURL = ""),
    expectsInput(objectName = "cohortData", objectClass = "data.table",
                  desc = "age cohort-biomass table hooked to pixel group map by pixelGroupIndex at
                  succession time step"),
    expectsInput(objectName = "biomassMap", objectClass = "RasterLayer",
                 desc = "Biomass map at each succession time step"),
    expectsInput(objectName = "biomassMapPreFire", objectClass = "RasterLayer",
                 desc = "Biomass map from before the last fire."),
    expectsInput(objectName = "rstCurrentBurn", objectClass = "RasterLayer",
                 desc = "Binary raster of fire spread"),
    expectsInput(objectName = "fireYear", objectClass = "numeric", desc =  "Next fire year"),
    expectsInput(objectName = "shpStudySubRegion", objectClass = "SpatialPolygonsDataFrame",
                 desc = "this shape file contains two informaton: Sub study area with fire return interval attribute. 
                 Defaults to a shapefile in Southwestern Alberta, Canada", sourceURL = ""),
    expectsInput(objectName = "shpStudyRegionFull", objectClass = "SpatialPolygonsDataFrame",
                 desc = "this shape file contains two informaton: Full study area with fire return interval attribute.
                 Defaults to a shapefile in Southwestern Alberta, Canada", sourceURL = "")
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
      if(!is.null(sim$rstCurrentBurn)) {   ## only if fire module is "active"
        sim <- scheduleEvent(sim, params(sim)$fireSpread$fireInitialTime, "fireSeverity", "calcSeverity", eventPriority = 4)
        
        if(P(sim)$.plotMaps)
          sim <- scheduleEvent(sim, params(sim)$fireSpread$fireInitialTime, "fireSeverity", "severityPlot", eventPriority = 4.5)
        
        if(!any(is.na(P(sim)$.saveInitialTime)))
          sim <- scheduleEvent(sim, params(sim)$fireSpread$fireInitialTime,
                               "fireSeverity", "saveSeverity", eventPriority = 4.75)
      }
    },
    calcSeverity = {
      if(!all(is.na(sim$rstCurrentBurn[]))) {
        ## calculate severity
        sim <- doSeverity(sim)
        
        ## schedule future events
        sim <- scheduleEvent(sim, eventTime = sim$fireYear, moduleName = "fireSeverity", 
                             eventType = "calcSeverity", eventPriority = 4)
      }
    },
    severityPlot = {
      if(!all(is.na(sim$rstCurrentBurn[]))) {
        ## Plot severity and vegetation changes
        sim <- doSeverityPlot(sim)
        
        ## schedule next plot
        sim <- scheduleEvent(sim, sim$fireYear, "fireSeverity", "severityPlot", eventPriority = 4.5)
      }
    },
    saveSeverity = {
      if(!all(is.na(sim$rstCurrentBurn[]))) {
        ## Plot severity and vegetation changes
        sim <- doSaveSeverity(sim)
        
        ## schedule next plot
        sim <- scheduleEvent(sim, time(sim) + P(sim)$fireTimestep, "fireSeverity", "saveSeverity", eventPriority = 4.75)
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
  fireMask <- sim$rstCurrentBurn
  fireMask[!is.na(fireMask)] <- 1
  
  ## Make raster of post-fire biomass  
  ## note that for this, this event must come before dispersal/regeneration events in LBMR
  cohortData <- sim$cohortData
  pixelAll <- cohortData[,.(uniqueSumB = as.integer(sum(B, na.rm=TRUE))), by=pixelGroup]
  biomassMapPostFire <- rasterizeReduced(pixelAll, sim$pixelGroupMap, "uniqueSumB")
  raster::projection(biomassMapPostFire) <- raster::projection(sim$biomassMap)
  
  ## fire severity in % mortality ((pre - post)/pre) -- this nedes to be changed as post-fire biomass cannot include dispersal/regeneration
  severityMap <- ((sim$biomassMapPreFire - biomassMapPostFire)/sim$biomassMapPreFire) * 100
  severityMap <- setValues(severityMap, values = round(getValues(severityMap)))
  severityMap <- raster::mask(severityMap, mask = fireMask)

  ## fire severity in % TFC
  # sim$severityMap <- sim$fireTFCRas
  # sim$severityMap <- setValues(sim$severityMap, values = scales::rescale(getValues(sim$severityMap), to = c(0,1))*100)
  # sim$severityMap <- raster::mask(sim$severityMap, mask = fireMask)

  ## export to sim
  sim$severityMap <- severityMap
  
  return(invisible(sim))
}

### Plot fire severity and vegetation
doSeverityPlot <- function(sim) {  
  Plot(sim$severityMap, new = TRUE,
       title = "Fire severity",
       cols = heat.colors(10))
  
  return(invisible(sim))
}

doSaveSeverity = function(sim) {
  raster::projection(sim$severityMap) <- raster::projection(sim$ecoregionMap)
  writeRaster(sim$severityMap,
              file.path(outputPath(sim), paste("severityMap_Year", round(time(sim)), ".tif",sep="")), 
              datatype='INT2S', overwrite = TRUE)
  return(invisible(sim))
}

## OTHER INPUTS AND FUNCTIONS --------------------------------
.inputObjects = function(sim) {
  
  dPath <- dataPath(sim)
  cacheTags = c(currentModule(sim), "function:.inputObjects")
  
  if (!suppliedElsewhere("shpStudyRegionFull", sim)) {
    message("'shpStudyRegionFull' was not provided by user. Using a polygon in Southwestern Alberta, Canada")
    
    canadaMap <- Cache(getData, 'GADM', country = 'CAN', level = 1, path = asPath(dPath),
                       cacheRepo = getPaths()$cachePath, quick = FALSE) 
    smallPolygonCoords = list(coords = data.frame(x = c(-115.9022,-114.9815,-114.3677,-113.4470,-113.5084,-114.4291,-115.3498,-116.4547,-117.1298,-117.3140), 
                                                  y = c(50.45516,50.45516,50.51654,50.51654,51.62139,52.72624,52.54210,52.48072,52.11243,51.25310)))
    
    sim$shpStudyRegionFull <- SpatialPolygons(list(Polygons(list(Polygon(smallPolygonCoords$coords)), ID = "swAB_polygon")),
                                              proj4string = crs(canadaMap))
    
    ## use CRS of biomassMap
    sim$shpStudyRegionFull <- spTransform(sim$shpStudyRegionFull,
                                          CRSobj = P(sim)$.crsUsed)
    
  }
  
  if (!suppliedElsewhere("shpStudySubRegion", sim)) {
    message("'shpStudySubRegion' was not provided by user. Using the same as 'shpStudyRegionFull'")
    sim$shpStudySubRegion <- sim$shpStudyRegionFull
  }
  
  if (!identical(P(sim)$.crsUsed, crs(sim$shpStudyRegionFull))) {
    sim$shpStudyRegionFull <- spTransform(sim$shpStudyRegionFull, P(sim)$.crsUsed) #faster without Cache
  }
  
  if (!identical(P(sim)$.crsUsed, crs(sim$shpStudySubRegion))) {
    sim$shpStudySubRegion <- spTransform(sim$shpStudySubRegion, P(sim)$.crsUsed) #faster without Cache
  }
  
  ## load ecoregion map
  if (!suppliedElsewhere("ecoregionMap", sim )) {
    ## LANDIS-II demo data:
    
    # sim$ecoregionMap <- Cache(prepInputs,
    #                           url = extractURL("ecoregionMap"),
    #                           destinationPath = dPath,
    #                           targetFile = "ecoregions.gis",
    #                           fun = "raster::raster")
    
    ## Dummy version with spatial location in Canada
    ras <- projectExtent(sim$shpStudySubRegion, crs = sim$shpStudySubRegion)
    res(ras) = 250
    ecoregionMap <- rasterize(sim$shpStudySubRegion, ras)
    
    ecoregionMap[!is.na(getValues(ecoregionMap))][] <- sample(ecoregion$mapcode, 
                                                              size = sum(!is.na(getValues(ecoregionMap))), 
                                                              replace = TRUE) 
    sim$ecoregionMap <- ecoregionMap
  }
  
  if(!suppliedElsewhere("biomassMap", sim)) {
    biomassMapFilename <- file.path(dPath, "NFI_MODIS250m_kNN_Structure_Biomass_TotalLiveAboveGround_v0.tif")
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
  
  return(invisible(sim))
}