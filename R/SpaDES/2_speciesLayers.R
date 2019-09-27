## -----------------------------------
## LOAD/MAKE SPECIES LAYERS
## -----------------------------------

## this script makes a pre-simulation object that makes species layers
## by running BiomassSpeciesData. This is the longest module to run and,
## unless the study area or the species needed change, it whould only
## be run once (even if other things change, like the simulation rep,
## or other modules). That's why caching is kept separate from the rest
## of the simulation

speciesPaths <-list(cachePath = file.path("R/SpaDES/cache/LIM_tests", "speciesLayers"),
                    modulePath = file.path("R/SpaDES/m"),
                    inputPath = file.path("R/SpaDES/inputs"),
                    outputPath = file.path("R/SpaDES/outputs", runName))

speciesParameters <- list(
  BiomassSpeciesData = list(
    "types" = c("KNN", "CASFRI", "Pickell", "ForestInventory")
    , "sppEquivCol" = sppEquivCol
    , ".useCache" = TRUE
  )
)

speciesObjects <- list(
  "sppEquiv" = sppEquivalencies_CA
  , "sppColorVect" = sppColorVect
  , "studyAreaLarge" = foothillsMED
)

simOutSpeciesLayers <- Cache(simInitAndSpades
                             , times = list(start = 0, end = 1)
                             , params = speciesParameters
                             , modules = "BiomassSpeciesData"
                             , objects = speciesObjects
                             , paths = speciesPaths
                             , debug = TRUE
                             , .plotInitialTime = NA)
