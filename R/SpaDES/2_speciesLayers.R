## -----------------------------------
## LOAD/MAKE SPECIES LAYERS
## -----------------------------------

## this script makes a pre-simulation object that makes species layers
## by running Biomass_speciesData. This is the longest module to run and,
## unless the study area or the species needed change, it whould only
## be run once (even if other things change, like the simulation rep,
## or other modules). That's why caching is kept separate from the rest
## of the simulation

speciesPaths <-list(cachePath = file.path("R/SpaDES/cache/AI_report", "speciesLayers"),
                    modulePath = file.path("R/SpaDES/m"),
                    inputPath = file.path("R/SpaDES/inputs"),
                    outputPath = file.path("R/SpaDES/outputs", "speciesLayers"))

speciesParameters <- list(
  Biomass_speciesData = list(
    "types" = c("KNN", "CASFRI", "Pickell", "ForestInventory")
    , "sppEquivCol" = sppEquivCol
    , ".useCache" = c(".inputObjects", "init")
  )
)

speciesObjects <- list(
  "sppEquiv" = sppEquivalencies_CA
  , "sppColorVect" = sppColorVect
  # , "studyAreaLarge" = foothillsMED
  , "studyAreaLarge" = foothills
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
                             , omitArgs = c("userTags"))
