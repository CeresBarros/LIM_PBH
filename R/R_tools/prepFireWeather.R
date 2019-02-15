## ----------------------------------------
## FIRE WEATHER PREP
## 
## Ceres Nov 26th 2018
## Adapted from Colin Ferster Oct 6, 2014
## ----------------------------------------

## this script should be sourced

## FUNCTION TO COMPILE FIRE WEATHER DATA - adapted from Colin Ferter's script (processFireWeather.R)
## folder is the directory path where the converted shapefiles will be saved

prepFireWeather <- function(folder) {
  ## LOAD DATA -----------------------------------------
  ## files to load
  files <- c("interp_50s.tab", "interp_60s.tab","interp_70s.tab","interp_80s.tab", "interp_90s.tab")
  weatherDataFolder <- file.path(folder, "20DayInterpolations/extracted") 
  weatherCodeFolder <- file.path(folder, "../fireWeatherCode")  
  
  ## metadata downloaded from FTP Site follows: (the 21 days of fire weather seems to fit)
  ## Fire Description line from LFDB:
  ## num, fireyr, firemon, fireday, prov, lat, lon, cause, size(ha), ecoregion, ecozone
  ## 21 days of weather lines:
  ## julDay, temp, rh, ws, rain, ffmc, dmc, dc, isi, bui, fwi
  
  weatherLinesNames <- c("julDay", "temp", "rh", "ws", "rain", "ffmc", "dmc", 
                         "dc", "isi", "bui", "fwi")
  fireDescriptionsNames <- c("num", "fireyr", "firemon", "fireday", "prov", "lat", "lon", 
                             "cause", "size", "ecoregion", "ecozone", "fireNum")
  
  provinceCodes <- c("AB","BC","MB","NB","NL","NS","NT","NU","ON","PE","QC","SK","YT","MN","NWT","PQ","WBNP")
  provinceInterestCodes <- c("AB","SK")
  
  fireWeatherTables <- lapply(files, FUN = function(x) {
    ## this file initially contained both fire description and weather data.
    tempWeatherLines <- as.data.table(read.delim2(file.path(weatherDataFolder, x), 
                                                  sep = "", header = FALSE, stringsAsFactors = FALSE))
    setnames(tempWeatherLines, weatherLinesNames)  ## set names, even if data.table is still a mix between description and weather

    ## add fire number (julDay) - needs to be added line by line
    tempWeatherLines[, fireNum := as.numeric(NA)]
    fireNo <- NA
    for (i in 1:nrow(tempWeatherLines)) {
      if (tempWeatherLines[i, rain] %in% provinceCodes) {
        fireNo <- tempWeatherLines[i, julDay]
      }
      set(tempWeatherLines, i = i, j = "fireNum", value = fireNo)
    }
    
    ## separate fire descriptions, from the weather data
    tempFireDescription <- copy(tempWeatherLines[rain %in% provinceCodes])
    setnames(tempFireDescription, fireDescriptionsNames[])
    tempWeatherLines <- copy(tempWeatherLines[!rain %in% provinceCodes])
    
    ## subset by provinces of interest:
    tempFireDescription <- tempFireDescription[prov %in% provinceInterestCodes]
    tempWeatherLines <- tempWeatherLines[fireNum %in% unique(tempFireDescription$fireNum)]
    
    return(list(fireDescriptions = tempFireDescription, weatherLines = tempWeatherLines))
  })
 
  ## bind all tables
  fireDescriptions <- rbindlist(lapply(fireWeatherTables, FUN = function(x) x$fireDescriptions))
  weatherLines <- rbindlist(lapply(fireWeatherTables, FUN = function(x) x$weatherLines))
  
  ## convert to appropriate data types
  cols <- c("fireyr", "firemon", "fireday", "lat", "lon", "size")
  fireDescriptions[, (cols) := lapply(.SD, function(x) as.numeric(as.character(x))), .SDcols = cols]
  
  cols <- c("temp", "rh", "ws", "rain", "ffmc", "dmc", "dc", "isi", "bui", "fwi")
  weatherLines[, (cols) := lapply(.SD, function(x) as.numeric(as.character(x))), .SDcols = cols]
  
  ## convert to full year notation
  fireDescriptions[, fireyr := fireyr + 1900]
  
  
  ## ADD OTHER FIRE DETAILS ----------------
  ## Find matches identified using GIS, for non-matches Knn will be used
  fireWeatherMatches <- as.data.table(read.dbf(file.path(folder, "fire_weather.dbf")))
  firesOverview <- as.data.table(read.csv(file.path(weatherCodeFolder, "Overview.csv"), stringsAsFactors = FALSE))
  setnames(firesOverview, paste0("overview_", names(firesOverview)))
  setnames(firesOverview, "overview_Fire.Name", "FireName")
  firesOverview[, FireName := toupper(FireName)]
  
  ## Link the Algar and the Algarita up to the Baseline fire, since they happened at the same time and in the same area
  ## duplicate the baseline x 2, rename to algar and agarita
  baselineRow <- fireWeatherMatches[FireName %in% "BASELINE"]
  baselineRow[, FireName := "ALGAR"]
  fireWeatherMatches <- rbind(fireWeatherMatches, baselineRow)
  baselineRow[, FireName := "ALGARITA"]
  fireWeatherMatches <- rbind(fireWeatherMatches, baselineRow)
  rm(baselineRow)
  
  exactMatchesFires <- fireWeatherMatches[firesOverview, on = "FireName", nomatch = NA]
  setdiff(fireWeatherMatches$FireName, exactMatchesFires$FireName) ## two fires had no match in GIS
  
  ## make an ID column
  exactMatchesFires[, fireWeatherID := fireNum]
  
  
  ## find large enough fires (> 175 ha) between 1959 - 2001 with no match
  if (nrow(exactMatchesFires[is.na(RowNum) & 
                             overview_Area > 175 & 
                             overview_fyear >= 1959 & 
                             overview_fyear <= 2001]))  {
    warning(paste("No weather data for these recent (1959-2001) and large (>175 ha) fires:\n",
                  paste(exactMatchesFires[is.na(RowNum) & 
                                            overview_Area > 175 & 
                                            overview_fyear >= 1959 & 
                                            overview_fyear <= 2001, FireName], collapse = " "), 
                  "\nPlease check."))
  }

  ## these values do not make sense - change to NA
  firesOverview[overview_fday == 0, overview_fday := NA]
  firesOverview[overview_fmonth == 0, overview_fday := NA]
  
  ## add fire names to fire weather table and clean up 
  weatherLines <- fireWeatherMatches[, .(fireNum, FireName)][weatherLines, on = "fireNum",
                                                             nomatch = 0]
  weatherLines <- weatherLines[!is.na(FireName)]
  
  ## Get fire duration and correct one fire name
  fireDuration <- as.data.table(read.csv(file.path(weatherCodeFolder, "fireDuration.csv"), stringsAsFactors = FALSE))
  fireDuration[fireName == "ALFRED", fireName := "ALFRED LAKE"]
 
  ## join fields of interest with to fire weather table
  weatherLines <- fireDuration[, .(fireName, Duration, Area, Pre.)][weatherLines, on = "fireName==FireName", nomatch = 0]
  
  return(list(fireDescriptions = fireDescriptions, fireWeather = weatherLines))
}