## ------------------------------------------------------
## FIRE MODELLING WITH SpaDES
## TESTS
##
## Ceres: Nov 2017
## ------------------------------------------------------

## clean workspace
rm(list=ls()); amc::.gc()

## requires as of June 18th 2019
# loading reproducible     0.2.10.9000
# loading quickPlot        0.1.6.9000
# loading SpaDES.core      0.2.6.9000
# loading SpaDES.tools     0.3.2.9002
# loading SpaDES.addins    0.1.2

# devtools::install_github("PredictiveEcology/reproducible@development", upgrade = "always")
# devtools::install_github("achubaty/amc@development", upgrade = "always")
# devtools::install_github("PredictiveEcology/pemisc@development", upgrade = "always")
# devtools::install_github("PredictiveEcology/map@development", upgrade = "always")
# devtools::install_github("PredictiveEcology/LandR@development", upgrade = "always")
# devtools::install_github("PredictiveEcology/quickPlot@development", upgrade = "always")
# devtools::install_github("PredictiveEcology/SpaDES.tools@development", upgrade = "always")
# devtools::install_github("PredictiveEcology/SpaDES.core@development", upgrade = "always")
library(SpaDES)
library(LandR)

## testing packages
# try(detach("package:LandR", unload = TRUE))
# try(detach("package:SpaDES.core", unload = TRUE))
# try(detach("package:SpaDES.tools", unload = TRUE))
# try(detach("package:reproducible", unload = TRUE))
# devtools::load_all("../reproducible")
# devtools::load_all("../SpaDES.tools")
# devtools::load_all("../SpaDES.core")
# devtools::load_all("../LandR")

source("R/R_tools/Useful_functions.R")

## define paths
setPaths(modulePath = file.path("R/SpaDES/m"),
         inputPath = file.path("R/SpaDES/inputs"),
         outputPath = file.path("R/SpaDES/outputs"))


## STUDY AREA(S) ---------------------------------------

## Foothills and a smaller region for testing
## prepInputs doens't work with kmz, so download and unzipping need to be done externally.
foothills <- Cache(prepKMZ2shapefile,
                   url = "https://drive.google.com/open?id=1OCqRRIjRNFi6LmxY6m8QH4gMBOLTNeDs",
                   archive = "Foothills_study_area.zip",
                   destinationPath = "data/maps",
                   cacheRepo = "data/cache")
foothills <- spTransform(foothills,
                         "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
foothillsSMALL <- raster::buffer(foothills, width = -30000)

## -----------------------------------------------
## SIMULATION SETUP
## -----------------------------------------------

## Set up sppEquiv  ---------------------------
data("sppEquivalencies_CA", package = "LandR")
sppEquivalencies_CA[grep("Pin", LandR), `:=`(EN_generic_short = "Pine",
                                             EN_generic_full = "Pine",
                                             Leading = "Pine leading")]

## Make LIM spp equivalencies column
sppEquivalencies_CA[, LIM := c(Abie_bal = "Abie_sp", Abie_las = "Abie_sp", Abie_sp = "Abie_sp",
                               Lari_lar = "Lari_lar",
                               Pice_mar = "Pice_mar", Pice_gla = "Pice_gla", Pice_eng = "Pice_eng",
                               Pinu_con = "Pinu_sp",
                               Popu_tre = "Popu_sp", Betu_pap = "Popu_sp", Popu_bal = "Popu_sp",
                               Pseu_men = "Pseu_men")[LandR]]

sppEquivalencies_CA[, EN_generic_short := c(Abie_sp = "Fir",
                                            Lari_lar = "Larch",
                                            Pice_mar = "Bl spruce", Pice_gla = "Wh spruce", Pice_eng = "Eng spruce",
                                            Pinu_sp = "Pine",
                                            Popu_sp = "Decid",
                                            Pseu_men = "Doug-fir")[LIM]]
sppEquivalencies_CA[, EN_generic_full := c(Abie_sp = "Fir",
                                           Lari_lar = "Larch",
                                           Pice_mar = "Black spruce", Pice_gla = "White spruce", Pice_eng = "Engelmann spruce",
                                           Pinu_sp = "Pine",
                                           Popu_sp = "Deciduous",
                                           Pseu_men = "Doug-fir")[LIM]]

sppEquivalencies_CA[, FI_layers := c(Abie_sp = "Fir",
                                     Lari_lar = "",
                                     Pice_mar = "Black.Spruce", Pice_gla = "White.Spruce", Pice_eng = "",
                                     Pinu_sp = "Pine",
                                     Popu_sp = "Deciduous",
                                     Pseu_men = "")[LIM]]

sppEquivalencies_CA[LIM == "Abie_sp", Leading := "Fir leading"]
sppEquivalencies_CA[LIM == "Popu_sp", Leading := "Deciduous leading"]
sppEquivalencies_CA[LIM == "Pinu_sp", Leading := "Pine leading"]

## define spp column to use for model
sppEquivCol <- "LIM"
sppEquivalencies_CA <- na.omit(sppEquivalencies_CA, cols = sppEquivCol)

## create color palette for species used in model
sppColorVect <- sppColors(sppEquivalencies_CA, sppEquivCol,
                          newVals = "Mixed", palette = "Accent")

## Set up modelling parameters  ---------------------------
options('reproducible.useNewDigestAlgorithm' = TRUE)
runName <- "testFeb2019"
eventCaching <- c(".inputObjects", "init")
useParallel <- FALSE

## paths
pathsSim <- getPaths()
# pathsSim$outputPath <- "R/SpaDES/outputs/vegFB_0"
# pathsSim$outputPath <- file.path(pathsSim$outputPath, "vegFB_1/tests/Regen")
# pathsSim$cachePath <- file.path("R/SpaDES/cache/LIM_tests/Regen")
# pathsSim$outputPath <- file.path(pathsSim$outputPath, "allSPP")
# pathsSim$cachePath <- file.path("R/SpaDES/cache/LIM_tests/allSPP")
pathsSim$outputPath <- file.path(pathsSim$outputPath, runName)
pathsSim$cachePath <- file.path("R/SpaDES/cache/LIM_tests", runName)

## simulation params
timesSim <- list(start = 0, end = 15)
eventCaching <- c(".inputObjects", "init")


vegLeadingProportion <- 0 # indicates what proportion the stand must be in one species group for it to be leading.
# If all are below this, then it is a "mixed" stand
fireTimestep <- 2L
successionTimestep <- 1L

modulesSim <- list("BiomassSpeciesData"
                   , "Boreal_LBMRDataPrep"
                   , "LBMR"
                   , "Biomass_fuels"
                   , "Biomass_regenerationPM"
                   # , "Biomass_regeneration"
                   , "fireSpread"
                   , "fireSeverity"
)

objectsSim <- list("studyArea" = foothillsSMALL
                   , "sppEquiv" = sppEquivalencies_CA
                   , "sppColorVect" = sppColorVect
)

# outputs <- data.frame(expand.grid(objectName = c("cohortData"),
#                                   saveTime = seq(2, 50, by = 5),
#                                   stringsAsFactors = FALSE))
# outputs <- rbind(outputs, data.frame(objectName = "rstCurrentBurn",
#                                      saveTime = tail(seq(2, 50, by = 5), 1)))

paramsSim <- list(
  Boreal_LBMRDataPrep = list(
    "sppEquivCol" = sppEquivCol
    , "forestedLCCClasses" = c(1:15, 34:36)
    # next two are used when assigning pixelGroup membership; what resolution for
    #   age and biomass
    , "pixelGroupAgeClass" = successionTimestep * 10L
    , "pixelGroupBiomassClass" = 100
    , "runName" = runName
    , "useCloudCacheForStats" = FALSE
    , "cloudFolderID" = NA
    , ".useCache" = eventCaching
  )
  , LBMR = list(
    "calcSummaryBGM" = c("start")
    , "initialBiomassSource" = "cohortData" # can be 'biomassMap' or "spinup" too
    , ".plotInitialTime" = timesSim$start
    , "seedingAlgorithm" = "wardDispersal"
    , "sppEquivCol" = sppEquivCol
    , "successionTimestep" = successionTimestep * 10L
    , "vegLeadingProportion" = vegLeadingProportion
    , ".plotInterval" = 1
    , ".plotMaps" = FALSE
    , ".saveInitialTime" = NA
    , ".useCache" = eventCaching[eventCaching] # seems slower to use Cache for both
    , ".useParallel" = useParallel
  )
  , BiomassSpeciesData = list(
    "types" = c("KNN", "CASFRI", "Pickell", "ForestInventory")
    , "sppEquivCol" = sppEquivCol
    , ".useCache" = TRUE
  )
  , Biomass_fuels = list(
      "fireInitialTime" = fireTimestep
    , "fireTimestep" = fireTimestep
    , "sppEquivCol" = sppEquivCol
    , ".useCache" = eventCaching
  )
  , Biomass_regenerationPM = list(
    "fireInitialTime" = fireTimestep
    , "fireTimestep" = fireTimestep
    , "successionTimestep" = successionTimestep
  )
  # , Biomass_regeneration = list(
  #   "fireInitialTime" = fireTimestep
  #   , "fireTimestep" = fireTimestep
  #   , "successionTimestep" = successionTimestep
  # )
  , fireSpread = list(
    "fireInitialTime" = fireTimestep
    , "fireTimestep" = fireTimestep
    , "fireSize" = 1000L
    , "noStartPix" = 10
    , "vegFeedback" = TRUE
    , ".useCache" = eventCaching
  ),
  fireSeverity = list(
    "fireTimestep" = fireTimestep
    , ".plotMaps" = TRUE
    , ".saveInitialTime" = 1
    , ".useCache" = eventCaching
  )
)

# showCache(pathsSim$cachePath, after = "2018-09-26 00:00:00")
# reproducible::clearCache(pathsSim$cachePath, userTags = c("prepInputsLCC2005_rtm", "Boreal_LBMRDataPrep"))

## TODO CHANGE FIRE MODULES TO USE COHORT DATA RATHER THAN SUMMARY BMG OUTPUTS, LIKE BIOMASSMAP
## TODO: Biomass_fuels doFuelTypes is very slow. check.
options(spades.moduleCodeChecks = TRUE)
graphics.off()

## LBMR - only
# pathsSim$cachePath <- "R/SpaDES/cache/LBMRonly/testJun2019"
# pathsSim$outputPath <- "R/SpaDES/outputs/LBMRonly/testJun2019"
# options("reproducible.overwrite" = TRUE)
reproducible::clearCache(pathsSim$cachePath, userTags = c("^LBMR$", "init"), ask = FALSE)
set.seed(524326)
LBMR_testSim <- simInitAndSpades(times = timesSim
                                 , params = paramsSim
                                 , modules = modulesSim[c(3, 1, 2, 5, 4, 6, 7)]
                                 , objects = objectsSim
                                 , paths = pathsSim
                                 , debug = "TRUE"
                                 # , .plotInitialTime = NA
)


## TEST WITH FAKE FIRE MAP
## make fake fire map
rstCurrentBurn <- LBMR_testSimout@.envir$pixelGroupMap
rstCurrentBurn[rstCurrentBurn[]>0] <- 1
rstCurrentBurn[rstCurrentBurn[] <= 0] <- NA
IDs <- which(rstCurrentBurn[] == 1)
rstCurrentBurn[IDs[1:round(length(IDs)/2)]] <- NA

objectsSim["rstCurrentBurn"] <- rstCurrentBurn

modulesSim <- list("BiomassSpeciesData", "Boreal_LBMRDataPrep",   ## biomassSpeciesData needs a data prep -can't cope with LBMR defaults
                   "LBMR", "LandR_BiomassGMOrig",
                   "Biomass_fuels",
                   "LandR_BiomassRegen", "fireSpread",
                   "fireSeverity")

LBMR_testSim <- simInit(times = timesSim, params = paramsSim, modules = modulesSim,
                        objects = objectsSim, outputs = outputs, paths = pathsSim)

graphics.off()
dev()
clearPlot()
LBMR_testSimout <- spades(LBMR_testSim, cache = TRUE, debug = TRUE)   ## debug = TRUE activates automatic browsing when errors occur
completed(LBMR_testSimout)


