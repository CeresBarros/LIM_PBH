## -----------------------------------
##  DATA PREP FOR ANALYSES OF RESULTS
## -----------------------------------

if (!exists("yearSubset")) {
  stop("please provide yearSubset vector")
}

## load one of the preSim lists to get years and later the ecolocations map
preSimList <- loadSimList(file.path(simPaths$outputPath, "noPM", "LIM_simInit_noPM"))

## LOAD DATA (RESULTS)  -------------------
## Given the size of the data put together in a pixel-based format, results were sampled every 10 years (instead of the 5-year interval used for saving),
paramsResults <- list("LIM_resultsDataPrep" = list("endYear" = as.integer(end(preSimList)),
                                                   "parallel" = FALSE,
                                                   "reps" = 1L:10L,
                                                   "startYear" = start(preSimList),
                                                   "yearSubset" = as.integer(unique(c(seq(2011, 2111, 5), 2111))),
                                                   ".useCache" = c(".inputObjects", "init")))

objectsResults <- list("ecoregionLayer" = preSimList$ecoregionLayer,
                       "rasterToMatch" = preSimList$rasterToMatch,
                       "sppEquiv" = preSimList$sppEquiv)

outputsResults <- data.frame(expand.grid(objectName = c("allPixelBurnData"),
                                         saveTime = 1,
                                         eventPriority = 10,
                                         stringsAsFactors = FALSE))
outputsResults <- rbind(outputsResults, data.frame(objectName = "allPixelCohortData",
                                                   saveTime = 1,
                                                   eventPriority = 10))
outputsResults <- rbind(outputsResults, data.frame(objectName = "allPixelCohortDataMnt",
                                                   saveTime = 1,
                                                   eventPriority = 10))
# options("LandR.assertions" = FALSE)
# simOut <- Cache(simInitAndSpades,
#                 times = list(start = 1, end = 1),
#                 params = paramsResults,
#                 modules = "LIM_resultsDataPrep",
#                 objects = objectsResults,
#                 outputs = outputsResults,
#                 paths = simPaths,
#                 cacheRepo = simPaths$cachePath,
#                 userTags = c("simInitAndSpades", "LIM_resultsDataPrep"),
#                 omitArgs = "userTags")
# rm(simOut)   ## cache not working. maybe because the data's too big.

## alternatively:
allPixelBurnData <- readRDS(list.files(simPaths$outputPath, "allPixelBurnData", full.names = TRUE))
allPixelCohortData <- readRDS(list.files(simPaths$outputPath, "allPixelCohortData_", full.names = TRUE))
allPixelCohortDataMnt <- readRDS(list.files(simPaths$outputPath, "allPixelCohortDataMnt", full.names = TRUE))

amc::.gc()

## GET CAMERON'S AGE DATA AND STAND VEG TYPES -----------------------
ageDataCN <- fread("data/CameronsAgeData/treelist_outputs_for Ceres.csv")
patchVegTypeCN <- fread("data/CameronsAgeData/patch outputs_for Ceres.csv")

ageDataCN <- patchVegTypeCN[ageDataCN, on = .(Patch.ID)]
ageDataCN$Cover.dendro <- sub("Mixedwood", "mixedwood", ageDataCN$Cover.dendro)
ageDataCN$Cover.dendro <- sub("Broadleaf", "broadleaf", ageDataCN$Cover.dendro)
ageDataCN$Cover.dendro <- sub("-", "", ageDataCN$Cover.dendro)
rm(patchVegTypeCN)

## remove a record that seems funky (maybe it's a new cohort?)
ageDataCN <- ageDataCN[Reconstructed.age != 2018]

## calculate no. cohorts:
ageDataCN[, noCohorts := length(unique(Est.bin)) , by = .(Patch.ID)]

ageDataCN[, firePresAbs := ifelse(is.na(mean.FI), 0, 1)]
