## -----------------------------------
## LOAD/MAKE FIRE WEATHER
## -----------------------------------

## this script makes a pre-simulation object that gets and filters fire weather data
## by running fireWeather This is module involves huge data objects and
## unless the study area or the data change, it should only
## be run once (even if other things change, like the simulation rep,
## or other modules). That's why caching is kept separate from the rest
## of the simulation

fireWeatherPaths <-list(cachePath = file.path("R/SpaDES/cache/LIM_tests", "fireWeather"),
                        modulePath = file.path("R/SpaDES/m"),
                        inputPath = file.path("R/SpaDES/inputs"))

fireWeatherParameters <- list(
  fireWeather = list(
    ".useCache" = eventCaching
  )
)

simOutFireWeather <- Cache(simInitAndSpades
                           , times = list(start = 0, end = 1)
                           , params = fireWeatherParameters
                           , modules = "fireWeather"
                           , paths = fireWeatherPaths
                           , debug = TRUE
                           , .plotInitialTime = NA
                           # , useCache = "overwrite"
                           , cacheRepo = fireWeatherPaths$cachePath
                           , userTags = "simInitFireWeather"
                           , omitArgs = c("userTags"))
