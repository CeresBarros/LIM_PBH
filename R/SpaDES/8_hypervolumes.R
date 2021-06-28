## -----------------------------------
##  HYPERVOLUMES OF PYRODIVERISTY
## -----------------------------------

library(SpaDES)
library(ToolsCB)
library(data.table)
library(future.apply)

source("R/R_tools/convertToCNVegType.R")
source("R/R_tools/Useful_functions.R")
source("R/R_tools/hypervolumesWrapper.R")
options("reproducible.useNewDigestAlgorithm" = 2)
options("spades.moduleCodeChecks" = FALSE)
options("reproducible.useCache" = TRUE)
options("reproducible.destinationPath" = normPath("R/SpaDES/inputs"))
options("reproducible.useGDAL" = FALSE)


## general paths
simDirName <- "jun2021Runs"
simPaths <- list(cachePath = file.path("R/SpaDES/cache", simDirName, "postSimAnalyses")
                 , modulePath = file.path("R/SpaDES/m")
                 , inputPath = file.path("R/SpaDES/inputs")
                 , outputPath = file.path("R/SpaDES/outputs", simDirName))

## path to figure folder and cache folder
figOutputPath <- file.path(simPaths$outputPath, "figuresAnalysis")
dir.create(figOutputPath)
HVoutputPath <- file.path(simPaths$outputPath, "hypervolumes")
# bw.outputPath <- file.path(HVoutputPath, "bwTest")

## LOAD DATA (RESULTS)  -------------------
preSimList <- loadSimList(file.path(simPaths$outputPath, "noPM", "LIM_simInit_noPM"))

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

## get rid of simOut
# allPixelBurnData <- simOut$allPixelBurnData
# allPixelCohortDataMnt <- simOut$allPixelCohortDataMnt
# rm(simOut)

## alternatively:
allPixelBurnData <- readRDS(list.files(simPaths$outputPath, "allPixelBurnData", full.names = TRUE))
allPixelCohortDataMnt <- readRDS(list.files(simPaths$outputPath, "allPixelCohortDataMnt", full.names = TRUE))
amc::.gc()

## FIRE ATTRIBUTES HYPERVOLUMES -----------
## Fire properties (fire patch size in pixels, fire frequency, fire severity as biomass loss)

## summarize data first - as in Steel et al 2021, fire properties are summarized across time, but by pixel
## (and by scenario/rep)
## only look at pixels with vegetation dynamics so that we can compare with biodiv. HVs
summaryFireAttributes <- allPixelBurnData[!is.na(pixelGroup), list(meanFreq = mean(fireFreq),
                                                                   meanSev = mean(severity),
                                                                   meanSevB = mean(severityB),
                                                                   meanPatchS = mean(patchSize)),
                                          by = .(scenario, rep, pixelIndex)]
## add vegType per pixel at the first year of fire,
## and add pixels that had no fires
firstFireYr <- P(preSimList)$fireSpread$fireInitialTime
cols <- c("pixelIndex", "vegTypeCN", "scenario", "rep")
summaryFireAttributes <- summaryFireAttributes[unique(allPixelCohortDataMnt[year == firstFireYr, ..cols]),
                                               on = c("scenario", "rep", "pixelIndex")]
## checks
test1 <- any(is.na(summaryFireAttributes$vegTypeCN))
if (test1) {
  stop("NA vegTypeCNs where pixelGroup - i.e. vegetation - exists")
}

test2 <- allPixelBurnData[!is.na(pixelGroup)][summaryFireAttributes[is.na(meanFreq), .(scenario, rep, pixelIndex)],
                                              on = .(scenario, rep, pixelIndex), nomatch = 0]
if (dim(test2)[1]) {
  stop("pixels that had fire and veg data in allPixelBurnData were ",
       "accidentally dropped when adding vegTypeCN")
}

## set fire attributes to 0 in pixels that had no fires
cols <- c("meanFreq", "meanSev", "meanSevB", "meanPatchS")
summaryFireAttributes[, (cols) := lapply(.SD, replaceNAs), .SDcols = cols]
amc::.gc()

## HYPERVOLUMES BY VEGETATION TYPE ----------
## only montane belt
amc::.gc()
ncores <- 5
if (.Platform$OS.type == "windows") {
  plan("multisession", workers = ncores)
} else {
  plan("multicore", workers = ncores)
}
## gaussian HVs were extremely slow
future_lapply(split(summaryFireAttributes, by = c("rep", "vegTypeCN")), FUN = function(allData, bwVal1, bwVal2) {
  r <- unique(allData$rep)
  veg <- unique(allData$vegTypeCN)
  file.suffix <- paste0("fireHVs_", veg, "_rep", r)

  noAxes <- 4
  cols <- c("meanSev", "meanFreq", "meanSevB", "meanPatchS")

  fireHVWrapper(allData, noAxes, cols, file.suffix)
})
future:::ClusterRegistry("stop")


## HYPERVOLUMES ACROSS THE LANDSCAPE ----------------
## only montane belt
amc::.gc()
if (.Platform$OS.type == "windows") {
  plan("multisession", workers = ncores)
} else {
  plan("multicore", workers = ncores)
}

future_lapply(split(summaryFireAttributes, by = c("rep")), FUN = function(allData, bwVal1, bwVal2) {
  r <- unique(allData$rep)
  file.suffix <- paste0("fireHVs_landscape_rep", r)

  noAxes <- 4
  cols <- c("meanSev", "meanFreq", "meanSevB", "meanPatchS")

  fireHVWrapper(allData, noAxes, cols, file.suffix)

})
future:::ClusterRegistry("stop")

## VEGETATION ATTRIBUTES HYPERVOLUMES -----------
## HYPERVOLUMES BY VEGETATION TYPE --------------
## only montane belt

## Use the first fire year to identify the pixels we want to follow in time
## we follow the same pixels that were used to make fire attributes HVs, only now
## we select the start year and end years of the simulation
## the join shouldn't acutally change anything because we already subset the pixels with veg
## in the montane belt (regardless of fire)
pixelIndices <- unique(summaryFireAttributes[,.(scenario, rep, pixelIndex)])
vegDataForHVs <- allPixelCohortDataMnt[year %in% c(start(preSimList), end(preSimList))]

if (getOption("LandR.assertions")) {
  temp <- vegDataForHVs[pixelIndices, on = .(scenario, rep, pixelIndex), nomatch = 0]
  setkey(temp, scenario, rep, pixelIndex)
  setkey(vegDataForHVs, scenario, rep, pixelIndex)

  if (isFALSE(identical(temp, vegDataForHVs))) {
    stop("Something is wrong. summaryFireAttributes should have the same pixelIndex/scenario/rep\n",
         "Combinations as allPixelCohortDataMnt")
  }
}

## HV comparisons per year, between scenarios --------------
## gaussian HVs were extremely slow
amc::.gc()
ncores <- 3
if (.Platform$OS.type == "windows") {
  plan("multisession", workers = ncores)
} else {
  plan("multicore", workers = ncores)
}
future_lapply(split(vegDataForHVs, by = c("rep", "vegTypeCN", "year")), FUN = function(allData, HVoutputPath) {
  r <- unique(allData$rep)
  yr <- unique(allData$year)
  veg <- unique(allData$vegTypeCN)
  file.suffix <- paste0("vegHVs_", veg, "_yr", yr, "_rep", r)
  IDcols <- c("scenario", "rep", "pixelIndex")
  print(file.suffix)

  vegHVWrapper(allData,
               IDcols,
               "scenario",
               file.suffix,
               noAxes = 4,
               HVmethod = "svm",
               no.runs = 3,
               svm.gamma = 0.01,
               do.scale = TRUE,
               outputs.dir = HVoutputPath,
               saveOrdi = TRUE,
               plotOrdi = TRUE,
               plotHV = TRUE,
               verbose = FALSE)
}, HVoutputPath = HVoutputPath)
future:::ClusterRegistry("stop")


## HV comparisons per scenario, between years --------------
## gaussian HVs were extremely slow
amc::.gc()
ncores <- 3
if (.Platform$OS.type == "windows") {
  plan("multisession", workers = ncores)
} else {
  plan("multicore", workers = ncores)
}
future_lapply(split(vegDataForHVs, by = c("rep", "vegTypeCN", "scenario")), FUN = function(allData, HVoutputPath) {
  r <- unique(allData$rep)
  scen <- unique(allData$scenario)
  veg <- unique(allData$vegTypeCN)
  file.suffix <- paste0("vegHVs_", veg, "_", scen, "_rep", r)
  IDcols <- c("year", "rep", "pixelIndex")
  print(file.suffix)

  vegHVWrapper(allData,
               IDcols,
               "year",
               file.suffix,
               noAxes = 4,
               HVmethod = "svm",
               no.runs = 3,
               svm.gamma = 0.01,
               do.scale = TRUE,
               outputs.dir = HVoutputPath,
               saveOrdi = TRUE,
               plotOrdi = TRUE,
               plotHV = TRUE,
               verbose = FALSE)
}, HVoutputPath = HVoutputPath)
future:::ClusterRegistry("stop")

