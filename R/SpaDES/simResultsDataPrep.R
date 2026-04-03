## -----------------------------------
##  DATA PREP FOR ANALYSES OF RESULTS
## -----------------------------------
if (!exists("yearSubset")) {
  stop("please provide yearSubset vector")
}

if (!exists("runPrepResultsModule")) {
  stop("please provide runPrepResultsModule TRUE/FALSE")
}

## load one of the preSim lists to get maps, species traits and end(sim)
if (paste0(version$major, ".", strsplit(version$minor, "[.]")[[1]][1]) != "4.1") {
  warning("Simulations were run and saved with R v4.1 -- attempting to load simList with qread")
  preSimList <- qs::qread(file.path(simPaths$outputPath, "LIM_preSimulation.qs"))
} else
preSimList <- loadSimList(file.path(simPaths$outputPath, "LIM_preSimulation.qs"))

speciesTraits <- preSimList$species

## LOAD DATA (RESULTS)  -------------------
## Given the size of the data put together in a pixel-based format, results were sampled every 10 years (instead of the 5-year interval used for saving),
if (runPrepResultsModule) {
  preSimList$rasterToMatch@file@name <- sub(".*/LandscapesInMotion", getwd(), filename(preSimList$rasterToMatch))  ## paths change between machines.

  paramsResults <- list("LIM_resultsDataPrep" = list("startYear" = as.integer(min(yearSubset)),
                                                     "endYear" = as.integer(max(yearSubset)),
                                                     "parallel" = FALSE,
                                                     "reps" = 1L:5L,
                                                     "yearSubset" = as.integer(yearSubset),
                                                     ".useCache" = c("init")))

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

  options("LandR.assertions" = TRUE)
  options("spades.useRequire" = FALSE)
  options("spades.moduleCodeChecks" = FALSE)

  dtthreads <- data.table::getDTthreads()
  data.table::setDTthreads(1)
  simOut <- Cache(simInitAndSpades,
                  times = list(start = 1, end = 1),
                  params = paramsResults,
                  modules = "LIM_resultsDataPrep",
                  objects = objectsResults,
                  outputs = outputsResults,
                  paths = simPaths,
                  # useCache = FALSE,
                  cacheRepo = simPaths$cachePath,
                  # cacheId = "912a9bdc8ce27d3a",
                  showSimilar = TRUE,
                  userTags = c("simInitAndSpades", "LIM_resultsDataPrep"),
                  omitArgs = "userTags")

  allPixelBurnData <- simOut$allPixelBurnData
  allPixelCohortData <- simOut$allPixelCohortData
  allPixelCohortDataMnt <- simOut$allPixelCohortDataMnt

  rm(simOut)   ## cache not working. maybe because the data's too big.
  gc(reset = TRUE)
  data.table::setDTthreads(dtthreads)
} else {
  ## alternatively:
  allPixelBurnData <- readRDS(file.path(simPaths$outputPath, "allPixelBurnData_year1.rds"))
  allPixelCohortData <- readRDS(file.path(simPaths$outputPath, "allPixelCohortData_year1.rds"))
  allPixelCohortDataMnt <- readRDS(file.path(simPaths$outputPath, "allPixelCohortDataMnt_year1.rds"))
}

## not sure why the veg type isn't there..
if (!"vegType" %in% colnames(allPixelCohortDataMnt)) {
  allPixelCohortDataMnt <- unique(allPixelCohortData[, .(scenario, rep, year, pixelIndex, vegType)])[allPixelCohortDataMnt, on = .(scenario, rep, year, pixelIndex)]
}

for (i in 1:3) gc(reset = TRUE)

