## ------------------------------------------------------
## FIRE MODELLING WITH SpaDES
## TESTS
##
## Ceres: Nov 2017
## ------------------------------------------------------

## clean workspace
rm(list=ls()); gc(reset = TRUE)

## requires as of july 23rd 2018
# loading reproducible     0.2.2
# loading quickPlot        0.1.4
# loading SpaDES.core      0.2.1
# loading SpaDES.tools     0.3.0
# loading SpaDES.addins    0.1.1
# devtools::install_github("PredictiveEcology/reproducible@development")
# devtools::install_github("PredictiveEcology/quickPlot@development")
# devtools::install_github("PredictiveEcology/SpaDES.core@development")
# devtools::install_github("Predictive/SpaDES.tools@development")


library(SpaDES)
source("R/R_tools/Useful_functions.R")

## define paths
setPaths(cachePath = file.path("R/SpaDES/cache"),
         modulePath = file.path("R/SpaDES/m"),
         inputPath = file.path("R/SpaDES/inputs"),
         outputPath = file.path("R/SpaDES/outputs"))


## STUDY AREA(S) ---------------------------------------

## Foothills and a smaller region for testing
foothills <- prepInputs(targetFile = "data/maps/Foothills_study_area.shp",
                        url = "https://drive.google.com/file/d/10vFcsyMu_-UF3PEcDngKsU72gH7jciGk/view?usp=sharing",
                        destinationPath = "data/maps",
                        fun = "raster::shapefile")

foothills <-  raster::shapefile("data/maps/Foothills_study_area.shp")

foothillsSMALL <- raster::buffer(foothills, width = -0.3)

## Dave's AB and SK fire datasets - get from GDrive - NOT WORKING YET
# firesABSK <- prepInputs(targetFile = "data/fires_Dave/albertafires1_postfire.shp",
#                         alsoExtract = c("data/fires_Dave/albertafires1_postfire.",
#                                         "data/fires_Dave/albertafires2_prefire.",
#                                         "data/fires_Dave/albertafires2_postfire.",
#                                         "data/fires_Dave/saskatchewanfires_prefire.",
#                                         "data/fires_Dave/saskatchewanfires_postfire."),
#                         url = "https://drive.google.com/file/d/1wGCqI_X1t-PDM5eO6JWQW9hpNv_8zlum/view?usp=sharing",
#                         destinationPath = "data/fires_Dave",
#                         fun = "st_read", pkg = "sf")
# 
# firesABSK <- Cache(loadBindSpatialObjs, 
#                    files = c("albertafires1_postfire.shp", "albertafires2_postfire.shp", "saskatchewanfires_postfire.shp"),
#                    destinationPath = "data/fires_Dave", 
#                    urls = "https://drive.google.com/file/d/1wGCqI_X1t-PDM5eO6JWQW9hpNv_8zlum/view?usp=sharing",
#                    cacheRepo = getPaths()$cachePath)
# 
# ## buffer around the polygons will be the study area
# firesABSK.buf <- Cache(outerBuffer, 
#                        x = firesABSK, 
#                        cacheRepo = getPaths()$cachePath)
# 
# raster::crs(firesABSK.buf) = raster::crs(firesABSK)
# dev()
# sp::plot(firesABSK) ; sp::plot(firesABSK.buf, add =TRUE)

## HERE ##
## TO DO: adapt study areas to fire datasets.
## FULL REGION - From LandWeb 
# source("R/R_tools/inputMaps.R")   ## functions for input maps
# # studyArea = firesABSK.buf
# 
# # These are used in inputTables.R for filling the tables of parameters in
# subStudyRegionName <- "RIA"
# landisInputs <- readRDS("../LandWeb/inputs/landisInputs.rds")
# spEcoReg <- readRDS("../LandWeb/inputs/SpEcoReg.rds")
# 
# # The CRS for the Study -- spTransform converts this first one to the second one, they are identical geographically
# crsStudyRegion <- sp::CRS(paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0",
#                                 "+ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
# 
# 

# studyRegionFilePath <- {
#   studyRegionFilename <- if ("RIA" %in% subStudyRegionName) {
#     "RIA_SE_ResourceDistricts_Clip.shp"
#   } else {
#     "studyarea-correct.shp"
#   }
#   file.path(getPaths()$inputPath, studyRegionFilename)
# }
# 
# 
# studyRegionsShps <- Cache(loadStudyRegions, shpStudyRegionCreateFn = shpStudyRegionCreate,
#                           asPath(studyRegionFilePath),
#                           fireReturnIntervalMap = asPath(file.path(getPaths()$inputPath, "ltfcmap correct.shp")),
#                           subStudyRegionName = subStudyRegionName,
#                           crsStudyRegion = crsStudyRegion, cacheRepo = getPaths()$cachePath)
# list2env(studyRegionsShps, envir = environment()) # shpStudyRegion & shpStudyRegion


## Species list
# speciesList <- as.matrix(data.frame(speciesnamesRaw = c("Abie_Las", "Abie_Bal", "Betu_Pap", "Betu_Spp", "Lari_Lya", "Lari_Occ", "Pice_Eng", 
#                                                         "Pice_Gla", "Pice_Mar", "Pinus_Alb","Pinu_Ban", "Pinu_Con", "Pinu_Fle", "Pinu_Pon",
#                                                         "Popu_Del", "Popu_Tre", "Pseu_Men"),
#                      speciesNamesEnd =  c("Abie_sp", "Abie_sp", "Betu_pap", "Betu_sp", "Lari_sp", "Lari_sp", "Pice_eng",
#                                           "Pice_gla", "Pice_mar", "Pinu_sp", "Pinu_ban", "Pinu_con", "Pinu_sp", "Pinu_pon", 
#                                           "Popu_tre", "Popu_tre", "Pseu_men")))

## Betu_pap was causing spades to hang on writing CASFRI raster to disk
speciesList <- as.matrix(read.csv(file.path(getPaths()$inputPath, "speciesList.csv"), header = TRUE))

## SIMULATION SETUP ------------------------------

## simulation parameters

pathsSim <- getPaths()
timesSim <- list(start = 1, end = 3)

modulesSimtoy <- list("simplifyLCCVeg", "simpleLCCSuccession", "fireSpread", "fireSeverity")
objectsSimtoy <-  list(studyArea = foothillsSMALL)
paramsSimtoy <- list(
  .globals = list(.useCache = TRUE),
  fireSpread = list(fireSize = 1000, noStartPix = 100, fireFreq = 3),
  fireSeverity = list(.plotMaps = TRUE, .plotInterval = 1),
  fireStats = list(.plotStats = TRUE)
)


# eventCaching <- c(".inputObjects")
modulesSimLBMR <- list("BiomassSpeciesData", "Boreal_LBMRDataPrep",
                       "LBMR", "LandR_BiomassFuels")
                       
objectsSimLBMR <- list("shpStudyRegionFull" = foothillsSMALL,
                       "shpStudySubRegion" = foothillsSMALL,
                       "speciesList" = speciesList)

paramsSimLBMR <- list(
  # Boreal_LBMRDataPrep = list(.useCache = eventCaching),
  LBMR = list(successionTimestep = 1,
              .plotInitialTime = timesSim$start,
              .saveInitialTime = timesSim$start, 
              seedingAlgorithm = "wardDispersal"#,
              # seedingAlgorithm = "universalDispersal",
              # .useCache = eventCaching
              )
)


showCache(getPaths()$cachePath)
# clearCache(getPaths()$cachePath#, userTags = "Boreal_LBMRDataPrep"
           # )

LBMR_testSim <- simInit(times = timesSim, params = paramsSimLBMR, modules = modulesSimLBMR,
                        objects = objectsSimLBMR, paths = pathsSim)

# moduleDiagram(LBMR_testSim)
# objectDiagram(LBMR_testSim)
# moduleDiagram(LBMR_testSim)
# events(LBMR_testSim)

dev()
clearPlot()
LBMR_testSimout <- spades(LBMR_testSim, cache = TRUE, debug = TRUE)   ## debug = TRUE activates automatic browsing when errors occur
# events(LBMR_testSim)





