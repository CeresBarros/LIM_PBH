## ------------------------------------------------------
## FIRE MODELLING WITH SpaDES
## TESTS
##
## Ceres: Nov 2017
## ------------------------------------------------------

## clean workspace
rm(list=ls()); amc::.gc()

## requires as of Feb 27nd 2019
# loading reproducible     0.2.6.9003
# loading quickPlot        0.1.6
# loading SpaDES.core      0.2.5
# loading SpaDES.tools     0.3.1.9000
# loading SpaDES.addins    0.1.2
# devtools::install_github("PredictiveEcology/reproducible@development")
# devtools::install_github("achubaty/amc@development")
# devtools::install_github("PredictiveEcology/pemisc@development")
# devtools::install_github("PredictiveEcology/map@development")
# devtools::install_github("CeresBarros/LandR@development")
# devtools::install_github("PredictiveEcology/quickPlot@development")
# devtools::install_github("PredictiveEcology/SpaDES.tools@development")
# devtools::install_github("PredictiveEcology/SpaDES.core@development")
library(SpaDES)
library(LandR)

## testing packages
# try(detach("package:LandR", unload = TRUE))
# try(detach("package:SpaDES.core", unload = TRUE))
# try(detach("package:SpaDES.tools", unload = TRUE))
# try(detach("package:reproducible", unload = TRUE))
# devtools::load_all("C:/Ceres/GitHub/reproducible")
# devtools::load_all("C:/Ceres/GitHub/SpaDES.tools")
# devtools::load_all("C:/Ceres/GitHub/SpaDES.core")
# devtools::load_all("C:/Ceres/GitHub/LandR")

source("R/R_tools/Useful_functions.R")

## define paths
setPaths(modulePath = file.path("R/SpaDES/m"),
         inputPath = file.path("R/SpaDES/inputs"),
         outputPath = file.path("R/SpaDES/outputs"))


## STUDY AREA(S) ---------------------------------------

## Foothills and a smaller region for testing
foothills <- raster::shapefile("data/maps/Foothills_study_area.shp")
foothills <- spTransform(foothills,
                         "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

foothillsSMALL <- raster::buffer(foothills, width = -20000)

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
                               Pinu_con = "Pinu_sp", Pinu_ban = "Pinu_sp",
                               Popu_tre = "Popu_sp", Betu_pap = "Popu_sp", Popu_bal = "Popu_sp",
                               Pseu_men = "Pseu_men")[LandR]]

sppEquivalencies_CA[LIM == "Abie_sp", EN_generic_full := "Fir"]
sppEquivalencies_CA[LIM == "Abie_sp", EN_generic_short := "Fir"]
sppEquivalencies_CA[LIM == "Abie_sp", Leading := "Fir leading"]

sppEquivalencies_CA[LIM == "Popu_sp", EN_generic_full := "Deciduous"]
sppEquivalencies_CA[LIM == "Popu_sp", EN_generic_short := "Decid"]
sppEquivalencies_CA[LIM == "Popu_sp", Leading := "Deciduous leading"]

sppEquivalencies_CA <- sppEquivalencies_CA[!is.na(LIM),]


## create color palette for species used in model
sppEquivCol <- "LIM"
sppColors <- sppColors(sppEquivalencies_CA, sppEquivCol,
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
timesSim <- list(start = 1, end = 10)
.plotInitialTime = timesSim$start
# eventCaching <- c(".inputObjects")


vegLeadingProportion <- 0.8 # indicates what proportion the stand must be in one species group for it to be leading.
# If all are below this, then it is a "mixed" stand
fireTimestep <- 2L
successionTimestep <- 1L

modulesSim <- list("BiomassSpeciesData"
                   , "Biomass_regeneration"
                   , "Boreal_LBMRDataPrep"   ## biomassSpeciesData needs a data prep -can't cope with LBMR defaults
                   , "LandR_BiomassGMOrig"
                   , "LandR_BiomassFuels"
                   , "LBMR"
                   , "fireSpread"
                   , "fireSeverity"
)



objectsSim <- list("studyArea" = foothillsSMALL,
                   "sppEquiv" = sppEquivalencies_CA,
                   "sppColors" = sppColors)

# outputs <- data.frame(expand.grid(objectName = c("cohortData"),
#                                   saveTime = seq(2, 50, by = 5),
#                                   stringsAsFactors = FALSE))
# outputs <- rbind(outputs, data.frame(objectName = "rstCurrentBurn",
#                                      saveTime = tail(seq(2, 50, by = 5), 1)))

paramsSim <- list(
  Boreal_LBMRDataPrep = list(
    "sppEquivCol" = sppEquivCol
    # next two are used when assigning pixelGroup membership; what resolution for
    #   age and biomass
    , "pixelGroupAgeClass" = successionTimestep
    , "pixelGroupBiomassClass" = 100
    # , "establishProbAdjFacResprout" = 1 ## TODO: check defaults - not necessary with Boreal_.*@rmEstablishProbAdj
    # , "establishProbAdjFacNonResprout" = 1 ## TODO: check defaults
    , "runName" = runName
    , "useCloudCacheForStats" = FALSE
    , "cloudFolderID" = NA
    , ".useCache" = eventCaching[1]
  )
  , LBMR = list(
    "initialBiomassSource" = "cohortData" # can be 'biomassMap' or "spinup" too
    , "seedingAlgorithm" = "wardDispersal"
    , "sppEquivCol" = sppEquivCol
    , "successionTimestep" = successionTimestep
    , ".plotInitialTime" = 1
    , ".useCache" = FALSE # seems slower to use Cache for both
    , ".useParallel" = useParallel
  )
  , BiomassSpeciesData = list(
    "types" = c("KNN", "CASFRI", "Pickell", "ForestInventory")
    , "sppEquivCol" = sppEquivCol
    , "omitNonTreePixels" = FALSE
    , .useCache = TRUE
    )
  , LandR_BiomassGMOrig = list(
    ".useParallel" = useParallel
  )
  , LandR_BiomassFuels = list(
    "sppEquivCol" = sppEquivCol
  )
  , Biomass_regeneration = list(
    "fireInitialTime" = fireTimestep
    , "fireTimestep" = fireTimestep
    , "successionTimestep" = successionTimestep
    )
  , fireSpread = list(
    "fireSize" = 1000L
    , "noStartPix" = 10
    , "fireTimestep" = fireTimestep
    , "vegFeedback" = TRUE
    # , ".useCache" = eventCaching[1]
  ),
  fireSeverity = list(
    "fireTimestep" = fireTimestep
    , ".plotMaps" = TRUE
    , ".saveInitialTime" = 1)
)


# showCache(pathsSim$cachePath, after = "2018-09-26 00:00:00")
# reproducible::clearCache(pathsSim$cachePath, userTags = ".inputObjects")

## TODO: LandR_BiomassFuels doFuelTypes is very slow. check.
LBMR_testSim <- simInit(times = timesSim, params = paramsSim, modules = modulesSim,
                        objects = objectsSim, paths = pathsSim)

graphics.off()
dev()
clearPlot()
LBMR_testSimout <- spades(LBMR_testSim, cache = TRUE, debug = TRUE)   ## debug = TRUE activates automatic browsing when errors occur


## NOTE WITH "ALL" SPP, LAST YEAR IS THROWING AN ERROR.

dev(3)
clearPlot(dev = 3)

biomStk <- stack(list.files(path = "R/SpaDES/outputs/vegFB_1/tests/Regen/", pattern = "biomass", full.names = TRUE))
sevStk <- stack(list.files(path = "R/SpaDES/outputs/vegFB_1/tests/Regen/", pattern = "severity", full.names = TRUE))

plot(biomStk)
plot(sevStk)

# fireMap <- sevStk[[1]]
# fireMap[!is.na(fireMap[])] <- 1
# fireMap <- rasterToPolygons(fireMap, dissolve = TRUE)
# clearPlot()
# par(mfrow = c(1,2))
# plot(biomStk[[9]], main = "Biomass - pre-fire"); plot(fireMap, add = TRUE, border = "black")
# plot(biomStk[[1]], main = "Biomass - post-fire"); plot(fireMap, add = TRUE, border = "black")




# events(LBMR_testSimout)
# year <- 5L
# LBMR_testSimout <- scheduleEvent(LBMR_testSimout, eventTime = year, eventType = "save", moduleName = "LBMR")
# events(LBMR_testSimout)
#
# year <- 5.0
# LBMR_testSimout <- scheduleEvent(LBMR_testSimout, eventTime = year, eventType = "save", moduleName = "LBMR")
# events(LBMR_testSimout)
#
# year <- time(LBMR_testSimout) + 2
# LBMR_testSimout <- scheduleEvent(LBMR_testSimout, eventTime = year, eventType = "save", moduleName = "LBMR")
# events(LBMR_testSimout)


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
                   "LandR_BiomassFuels",
                   "LandR_BiomassRegen", "fireSpread",
                   "fireSeverity")

LBMR_testSim <- simInit(times = timesSim, params = paramsSim, modules = modulesSim,
                        objects = objectsSim, outputs = outputs, paths = pathsSim)

graphics.off()
dev()
clearPlot()
LBMR_testSimout <- spades(LBMR_testSim, cache = TRUE, debug = TRUE)   ## debug = TRUE activates automatic browsing when errors occur
completed(LBMR_testSimout)


