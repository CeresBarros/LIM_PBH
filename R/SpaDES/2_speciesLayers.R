## -----------------------------------
## LOAD/MAKE SPECIES LAYERS
## -----------------------------------

## this script makes a pre-simulation object that makes species layers
## by running Biomass_speciesData. This is the longest module to run and,
## unless the study area or the species needed change, it whould only
## be run once (even if other things change, like the simulation rep,
## or other modules). That's why caching is kept separate from the rest
## of the simulation

speciesPaths <- list(cachePath = file.path(simPaths$cachePath, "speciesLayers"),
                     modulePath = simPaths$modulePath,
                     inputPath = simPaths$inputPath,
                     outputPath = file.path(simPaths$outputPath, "speciesLayers"))

speciesParameters <- list(
  Biomass_speciesData = list(
    "dataYear" = 2011L
    , "types" = c("KNN", "CASFRI", "ForestInventory")   ## Pickell has no data here having errors with the other two. extents and NA data.
    , "sppEquivCol" = sppEquivCol
    , ".useCache" = eventCaching
  )
)

speciesObjects <- list(
  "sppEquiv" = sppEquivalencies_CA
  # , "studyAreaLarge" = foothillsMED
  , "studyAreaLarge" = foothills
  , "rasterToMatchLarge" = rasterToMatchLarge
)

simOutSpeciesLayers <- Cache(simInitAndSpades
                             , times = list(start = 0, end = 1)
                             , params = speciesParameters
                             , modules = "Biomass_speciesData"
                             , objects = speciesObjects
                             , paths = speciesPaths
                             , debug = TRUE
                             , .plotInitialTime = NA
                             , cacheRepo = speciesPaths$cachePath
                             , userTags = "simInitSpeciesLayers"
                             , omitArgs = c("userTags", ".plotInitialTime", "debug"))
