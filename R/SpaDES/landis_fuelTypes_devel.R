library(cffdrs)
library(humidity)

## TEMPERATURE (MONTHLY MAX)
## using http://worldclim.org/version2 (temperarute historical ensemble)
fnames <- paste0("wc2.0_2.5m_tmax_0", 5:9, ".tif")
temperatureStk <- list()
for(f in fnames) {
  temperatureStk[[f]] <- Cache(prepInputs, targetFile = f,
                               url = "http://biogeo.ucdavis.edu/data/worldclim/v2.0/tif/base/wc2.0_2.5m_tmax.zip",
                               archive = "wc2.0_2.5m_tmax.zip", 
                               alsoExtract = NA,
                               destinationPath = getPaths()$inputPath, 
                               studyArea = LBMR_testSimout$shpStudySubRegion,
                               rasterToMatch = LBMR_testSimout$biomassMap,
                               method = "bilinear",
                               datatype = "FLT4S",
                               postProcessedFilename = FALSE)
}
temperatureStk <- stack(temperatureStk)
temperatureRas <- mean(temperatureStk)

## PRECIPITATION (MONTHLY CUMMULATIVE)
fnames <- paste0("wc2.0_2.5m_prec_0", 5:9, ".tif")
precipitationStk <- list()
for(f in fnames) {
  precipitationStk[[f]] <- Cache(prepInputs, targetFile = f,
                               url = "http://biogeo.ucdavis.edu/data/worldclim/v2.0/tif/base/wc2.0_2.5m_prec.zip",
                               archive = "wc2.0_2.5m_prec.zip", 
                               alsoExtract = NA,
                               destinationPath = getPaths()$inputPath, 
                               studyArea = LBMR_testSimout$shpStudySubRegion,
                               rasterToMatch = LBMR_testSimout$biomassMap,
                               method = "bilinear",
                               datatype = "FLT4S",
                               postProcessedFilename = FALSE)
}
precipitationStk <- stack(precipitationStk)
precipitationRas <- mean(precipitationStk)

## project to Lat/Long (decimal degrees) for compatibility with FBP system
latLong = "+proj=longlat +datum=WGS84"
temperatureRas <- projectRaster(temperatureRas, crs = latLong)
precipitationRas <- projectRaster(precipitationRas, crs = latLong)

## make data.tables of temp and precip per pixel group
climData <- data.table(temp = temperatureRas[], precip = precipitationRas[],
                       pixelGroup = LBMR_testSimout@.envir$pixelGroupMap[],
                       lat =  coordinates(temperatureRas)[,2],
                       long = coordinates(temperatureRas)[,1])
## relative humidity
## using dew point between -3 and 20%, quarterly seasonal for Jun 2013
## https://calgary.weatherstats.ca/metrics/dew_point.html
climData[, relHum := RH(t = climData$temp, Td = runif(nrow(climData), -3, 20), isK = FALSE)]   

FWIinputs <- data.frame(id = climData$pixelGroup,
                        lat = climData$lat,
                        long = climData$long,
                        mon = 7,
                        temp = climData$temp,
                        rh = climData$relHum,
                        ws = 0,
                        prec = climData$precip)

## using defaults
FWIinit = data.frame(ffmc = 85,
                      dmc = 6,
                      dc = 15)

## CALCULATE FIRE WEATHER INDICES
FWIoutputs <- fwi(input = na.omit(FWIinputs), init = na.omit(FWIinit), batch = FALSE, lat.adjust = TRUE) %>%
  data.table

## merge FWI outputs with calculated fuel types
FTs <- data.table(pixelGroup = LBMR_testSimout$pixelFuelTypes$pixelGroup,
                  FuelType = LBMR_testSimout$pixelFuelTypes$finalFuelType)
FTs <- LBMR_testSimout$FuelTypes[, .(FuelTypeFBP,FuelType)] %>%
  .[!duplicated(.)] %>%
  FTs[., on = "FuelType", nomatch = 0]

FWIoutputs <- FWIoutputs[FTs[!duplicated(FTs)], on = "ID==pixelGroup", nomatch = 0]

## dataframe of minimum FBP inputs set to their defaults. See ?fbp for info
FBPinputs <- data.frame(id = FWIoutputs$ID,
                        FuelType = FWIoutputs$FuelTypeFBP,
                        LAT = FWIoutputs$LAT,
                        LONG = FWIoutputs$LONG, 
                        FFMC = FWIoutputs$FFMC,
                        BUI = FWIoutputs$BUI,
                        WS = FWIoutputs$WS,
                        GS = rep(0, nrow(FWIoutputs)),
                        Dj = rep(180, nrow(FWIoutputs)),
                        Aspect = rep(0, nrow(FWIoutputs)))
FBPoutputs <- fbp(input = FBPinputs) %>%
  data.table

## FBP OUTPUTS TO RASTERS
ras <- LBMR_testSimout@.envir$pixelGroupMap
## Rate of spread
ROSras <- reclassify(ras, as.matrix(data.frame(is = FBPoutputs$ID, becomes = FBPoutputs$ROS)))
# Head fire intensity
IntRas <- reclassify(ras, as.matrix(data.frame(is = FBPoutputs$ID, becomes = FBPoutputs$HFI)))
## Total fuel consumption
TFCRas <- reclassify(ras, as.matrix(data.frame(is = FBPoutputs$ID, becomes = FBPoutputs$TFC)))

dev();
clearPlot()
Plot(ROSras); Plot(IntRas); Plot(TFCRas)
