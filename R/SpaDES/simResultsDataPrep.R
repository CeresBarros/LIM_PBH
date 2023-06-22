## -----------------------------------
##  DATA PREP FOR ANALYSES OF RESULTS
## -----------------------------------
if (!exists("yearSubset")) {
  stop("please provide yearSubset vector")
}

if (!exists("runPrepResultsModule")) {
  stop("please provide runPrepResultsModule TRUE/FALSE")
}


## load one of the preSim lists to get years and later the ecolocations map
preSimList <- loadSimList(file.path(simPaths$outputPath, "noPM", "LIM_simInit_noPM.qs"))

## LOAD DATA (RESULTS)  -------------------
## Given the size of the data put together in a pixel-based format, results were sampled every 10 years (instead of the 5-year interval used for saving),
if (runPrepResultsModule) {
  paramsResults <- list("LIM_resultsDataPrep" = list("startYear" = as.integer(min(yearSubset)),
                                                     "endYear" = as.integer(max(yearSubset)),
                                                     "parallel" = FALSE,
                                                     "reps" = 1L:5L,
                                                     "yearSubset" = as.integer(yearSubset),
                                                     ".useCache" = c("init", "loadSimulationData", "joinSimulationData",
                                                                     "addVegTypesCN", "averageAcrossYears"))) ## parsing error when caching

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

  options("LandR.assertions" = FALSE)
  options("spades.useRequire" = FALSE)
  options("spades.moduleCodeChecks" = FALSE)
  simOut <- Cache(simInitAndSpades,
                  times = list(start = 1, end = 1),
                  params = paramsResults,
                  modules = "LIM_resultsDataPrep",
                  objects = objectsResults,
                  outputs = outputsResults,
                  paths = simPaths,
                  useCache = TRUE,
                  cacheRepo = simPaths$cachePath,
                  userTags = c("simInitAndSpades", "LIM_resultsDataPrep"),
                  omitArgs = "userTags")

  allPixelBurnData <- simOut$allPixelBurnData
  allPixelCohortData <- simOut$allPixelCohortData
  allPixelCohortDataMnt <- simOut$allPixelCohortDataMnt

  rm(simOut)   ## cache not working. maybe because the data's too big.
  gc(reset = TRUE)
} else {
  ## alternatively:
  allPixelBurnData <- readRDS(list.files(simPaths$outputPath, "allPixelBurnData", full.names = TRUE))
  allPixelCohortData <- readRDS(list.files(simPaths$outputPath, "allPixelCohortData_", full.names = TRUE))
  allPixelCohortDataMnt <- readRDS(list.files(simPaths$outputPath, "allPixelCohortDataMnt", full.names = TRUE))
}

## not sure why the veg type isn't there..
if (!"vegType" %in% colnames(allPixelCohortDataMnt)) {
  allPixelCohortDataMnt <- unique(allPixelCohortData[, .(scenario, rep, year, pixelIndex, vegType)])[allPixelCohortDataMnt, on = .(scenario, rep, year, pixelIndex)]
}

gc(reset = TRUE)

