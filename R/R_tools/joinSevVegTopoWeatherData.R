## ----------------------------------------
## JOIN FIRE SEVERITY, VEGETATION, TOPO AND WEATHER DATA
##
## Ceres Dec 11th 2018
## ----------------------------------------

## this script should be sourced
## sevDataSf - an sf object of severity polygons
## vegDataSf - an sf object of vegetation polygons
## topoDataSf - an sf object of topography data polygons that match the fire severity polygons (IDs).
## weatherDataDt - a dataframe of weather data per fire ID which matches the severity fire IDs.
## resolution is the resolution of output raster in meters. Defaults to 100.
## doAll = controls whether all data joins must be re-done/saved (joins are made by fire)
## saveDir is the directory where joint data tables are saved for each fire

## FUNCTION WRAPPER JOIN VEGETATION, TOPOGRAPHY AND WEATHER DATASETS, AFTER THEY'VE BEEN PREP'ED
joinSevVegTopoWeatherData <- function(sevDataSf, vegDataSf, topoDataSf, weatherDataDt,
                                      resolution = 100, doAll = FALSE, saveDir) {
  ## make a template raster that will be used to extract pixIDs
  ## note that IDs run 1:ncell(templateRas)
  ## not enought memory to rasterize the full extent to 10x10res
  rasterToMatch <- raster(sevDataSf, resolution = resolution,
                          crs = crs(sevDataSf))
  if(!compareCRS(crs(sevDataSf), crs(rasterToMatch))) {
    crs(rasterToMatch) <- crs(sevDataSf)
    rasterToMatch <- projectRaster(rasterToMatch, crs = crs(sevDataSf))
  }
  rasterToMatch <- fasterize(sf = sevDataSf, raster = rasterToMatch)
  rasterToMatch <- setValues(rasterToMatch, values = 1:ncell(rasterToMatch))

  ## do computations per fire
  fireNames <- as.character(unique(sevDataSf$FIRE_NAME))
  setnames(weatherDataDt, "fireName", "FIRE_NAME")

  ## save data.tables in a temporary folder
  dir.create(saveDir)

  if(!doAll) {
    firesDone <- sub(".RData", "", sub("dataTable_", "",
                                       list.files(saveDir, pattern = "dataTable")))
    firesDone <- gsub("_", " ", firesDone)
    fireNames <- fireNames[!gsub("_", " ", fireNames) %in% firesDone]
  }

  if(length(fireNames))
    message(paste0(length(fireNames), " fires to do.")) else
      message("All fires have been joined.")
  ## do joins and save output tables by fire
  invisible(utils::memory.limit(64000))
  for(x in fireNames) {
    joinPerFire(smallSevDataSf = sevDataSf[sevDataSf$FIRE_NAME == x,],
                vegDataSf = vegDataSf, topoDataSf = topoDataSf,
                weatherDataDt = weatherDataDt,
                rasterToMatch = rasterToMatch, saveDir = saveDir)
  }

  message("Binding tables...")

  allDataDT <- rbindlist(lapply(list.files(saveDir, pattern = "dataTable", full.names = TRUE),
                                FUN = function(x) {
                                  load(x)
                                  return(dataDT)
                                }))
  message("... done!")
  return(allDataDT)
}

## FUNCTION TO JOIN VEGETATION, TOPOGRAPHY AND WEATHER DATASETS, AFTER THEY'VE BEEN PREP'ED
## BY FIRE - INTERNAL
joinPerFire <- function(smallSevDataSf, vegDataSf, topoDataSf, weatherDataDt,
                        rasterToMatch, saveDir) {
  amc::.gc()
  message(paste0("Joining data for: ", as.character(unique(smallSevDataSf$FIRE_NAME))))
  ## do data "joins" by intersecting SF objects with a raster and making a data.table with pixIDs
  ## severity polygons to then extract veg, topo, weather data

  ## To make a templateraster with RTM pix IDs need to:
  ## 1) extract pix IDs in RTM from an extent, 2) crop RTM to extent,
  ## 3) rasterize, 4) subset pix IDs where polygons exist 1s
  ## the is because direct rasterization to RTM doesn't work on small fires
  ## and cropping doesn't keep RTM cell values as it should, creating duplicates
  pixID <- extract(rasterToMatch, extent(as_Spatial(smallSevDataSf)))

  templateRas <- crop(rasterToMatch, as_Spatial(smallSevDataSf))
  templateRas <- fasterize(sf = smallSevDataSf,
                           raster = templateRas)
  pixID <- pixID[!is.na(getValues(templateRas))]

  coords <- xyFromCell(templateRas, cell = which(!is.na(getValues(templateRas))))
  dataDT <- data.table(pixID = pixID,
                       Long = coords[, 1], Lat = coords[, 2])
  setkey(dataDT, pixID)

  ## convert templateRas to points, provide CRS
  templatePoints <- st_as_sf(dataDT, coords = c("Long", "Lat"))

  if(!compareCRS(crs(templatePoints), crs(templateRas))) {
    st_crs(templatePoints) <- as.character(crs(templateRas))
    templatePoints <- st_transform(templatePoints, crs =  as.character(crs(templateRas)))
  }

  ## GET SEVERITY DATA
  message("... joining fire severity data")
  ## simplify data to polygons - should be faster
  smallSevDataSf <- st_cast(smallSevDataSf, "POLYGON")

  ## make sure CRS is the same
  if(!compareCRS(crs(templatePoints), crs(smallSevDataSf)))
    smallSevDataSf <- st_transform(smallSevDataSf, crs = st_crs(templatePoints))

  ## extract polygon info per point
  sevDataPoints <- st_join(x = templatePoints,
                           y = smallSevDataSf)  ## >10xs faster than st_intersection

  ## make datatable and join to master table
  sevDataPoints <- data.table(st_set_geometry(sevDataPoints, NULL))
  setkey(sevDataPoints, pixID)
  dataDT <- sevDataPoints[dataDT, nomatch = 0]

  ## JOIN TOPOGRAPHY DATA
  message("... joining topographic data")
  if(!compareCRS(crs(topoDataSf), crs(smallSevDataSf)))
    topoDataSf <- st_transform(topoDataSf, crs = st_crs(smallSevDataSf))

  ## use prioritizr::fast_extract to extract raster IDs per point
  templateRasV <- templateRas
  templateRasV[!is.na(getValues(templateRasV))] <- pixID

  ## if extents to not intersect (some fires are beyong the topo SF extent)
  ## then make NA IDs
  topoDataPixID <- if (prioritizr:::intersecting_extents(templateRasV, topoDataSf)) {
    fast_extract(templateRasV, topoDataSf)   ## faster than raster::extract()
  } else {
    rep(NA, dim(topoDataSf)[1])
  }

  topoDataPoints <- data.table(st_set_geometry(topoDataSf, NULL))
  topoDataPoints$pixID <- topoDataPixID
  setkey(topoDataPoints, pixID)

  dataDT <- topoDataPoints[dataDT, nomatch = 0]

  ## clean-up to free memory:
  rm(topoDataPoints, topoDataPixID,
     templateRasV, templateRas,
     sevDataPoints, pixID, coords)
  amc::.gc()

  ## JOIN VEGETATION DATA
  message("... joining pre-fire vegetation data")

  if(!compareCRS(crs(templatePoints), crs(vegDataSf)))
    vegDataSf <- st_transform(vegDataSf, crs = st_crs(templatePoints))

  vegDataSf <- st_cast(vegDataSf, "POLYGON")

  vegDataPoints <- st_join(x = templatePoints,
                           y = vegDataSf)

  vegDataPoints <- data.table(st_set_geometry(vegDataPoints, NULL))
  setkey(vegDataPoints, pixID)
  setkey(dataDT, pixID)
  dataDT <- vegDataPoints[dataDT, nomatch = 0]

  ## JOIN WEATHER DATA
  message("... joining weather data")
  dataDT[FIRE_NAME == "Alfred", FIRE_NAME := "Alfred Lake"]
  dataDT[, FIRE_NAME := toupper(FIRE_NAME)]

  setkey(weatherDataDt, FIRE_NAME)
  setkey(dataDT, FIRE_NAME)
  dataDT <- weatherDataDt[dataDT, allow.cartesian = TRUE, nomatch = 0]

  message("... saving temp file")
  tempFile <- paste0("dataTable_",
                     sub(" ", "_", as.character(unique(smallSevDataSf$FIRE_NAME))),
                     ".RData")
  save("dataDT", file = file.path(saveDir, tempFile))
  message("... done!")
}
