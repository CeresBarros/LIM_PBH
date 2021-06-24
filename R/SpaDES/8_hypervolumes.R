## -----------------------------------
##  HYPERVOLUMES OF PYRODIVERISTY
## -----------------------------------

library(SpaDES)
library(ToolsCB)

source("R/R_tools/convertToCNVegType.R")
source("R/R_tools/Useful_functions.R")
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
bw.outputPath <- file.path(HVoutputPath, "bwTest")

## LOAD DATA (RESULTS)  -------------------
preSimList <- loadSimList(file.path(simPaths$outputPath, "noPM", "LIM_simInit_noPM"))

## Given the size of the data put together in a pixel-based format, results were sampled every 10 years (instead of the 5-year interval used for saving),
paramsResults <- list("LIM_resultsDataPrep" = list("endYear" = end(preSimList),
                                                   "reps" = 1L:10L,
                                                   "startYear" = start(preSimList),
                                                   "yearSubset" = unique(c(seq(2011L, 2111L, 5), 2111L)),
                                                   ".useCache" = c(".inputObjects", "init",
                                                                   "loadSimulationData", "joinSimulationData",
                                                                   "addVegTypesCN"))
)

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
simOut <- Cache(simInitAndSpades,
                times = list(start = 1, end = 1),
                params = paramsResults,
                modules = "LIM_resultsDataPrep",
                outputs = outputsResults,
                objects = objectsResults,
                paths = simPaths,
                cacheRepo = simPaths$cachePath,
                userTags = c("simInitAndSpades", "LIM_resultsDataPrep"),
                omitArgs = "userTags")

# source("R/SpaDES/6_resultsDataPrep.R")


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
summaryFireAttributes <- summaryFireAttributes[allPixelCohortDataMnt[year == firstFireYr, ..cols],
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

## HYPERVOLUMES BY VEGETATION TYPE -----------------------------------
## BANDWITH ESTIMATES ----
amc::.gc()
parallel_wrapper <- function(ncores, summaryFireAttributes, byVars, bw.outputPath) {
  if (!dir.exists(bw.outputPath)) dir.create(bw.outputPath)

  if (.Platform$OS.type == "windows") {
    plan("multisession", workers = ncores)
  } else {
    plan("multicore", workers = ncores)
  }
  bw_estimates <- future_lapply(X = split(summaryFireAttributes, by = byVars),
                                FUN = function(DT, bw.outputPath) {
                                  r <- unique(DT$rep)
                                  veg <- unique(DT$vegTypeCN)
                                  message(paste("Calculating PCAs and estimating BWs for: rep", r, "and", veg))
                                  file.suffix <- paste0("fireHVs_freeBW_", veg, "_rep", r)

                                  init.vars <- grep("mean", names(DT))

                                  DT <- ToolsCB:::.scaleVars(DT, init.vars)

                                  out <- estimateBW_wrapper(as.data.frame(DT),
                                                            init.vars = init.vars,
                                                            HVidvar = which(names(DT) == "scenario"),
                                                            noAxes = 4,
                                                            ordination = "PCA",
                                                            file.suffix = file.suffix,
                                                            outputs.dir = bw.outputPath)
                                  return(out)
                                },
                                bw.outputPath = bw.outputPath)
  future:::ClusterRegistry("stop")
  return(bw_estimates)
}

bw_estimates <- Cache(parallel_wrapper,
                      summaryFireAttributes = summaryFireAttributes,
                      ncores = 10,
                      byVars = c("rep", "vegTypeCN"),
                      bw.outputPath = bw.outputPath,
                      cacheRepo = simPaths$cachePath,
                      userTags = c("bw_estimates", "hypervolumes", "vegTypeCN"),
                      omitArgs = c("userTags", "ncores", "bw.outputPath"))

bw_estimates <- do.call(rbind.data.frame, bw_estimates)
bw_estimates$vegType <- sub("[[:digit:]]*\\.", "", sub("\\.PC.*", "", row.names(bw_estimates)))
bw_estimates$rep <- sub("\\..*", "", sub("\\.PC.*", "", row.names(bw_estimates)))
bw_estimates <- as.data.table(bw_estimates)
saveRDS(bw_estimates, file.path(bw.outputPath, "BW_estimates_vegType.rds"))

summaryBW <- bw_estimates[, list(SilvermanMean = mean(c(SilvBW_HV1, SilvBW_HV2)),
                                 StDevMean = mean(c(stdev_HV1, stdev_HV2)),
                                 SilvermanMax = max(c(SilvBW_HV1, SilvBW_HV2)),
                                 StDevMax = max(c(stdev_HV1, stdev_HV2))),
                          by = "PC"]

saveRDS(summaryBW, file.path(bw.outputPath, "BW_MeanMax_vegType.rds"))


## fix bandwidth to max of estimated BW
summaryBW <- readRDS(file.path(bw.outputPath, "BW_MeanMax_vegType.rds"))
bwHV <- summaryBW[["StDevMax"]]

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
  hypervolumes(HVdata1 = as.data.frame(allData[scenario == "noPM"]),
               HVdata2 = as.data.frame(allData[scenario == "PM"]),
               HVidvar = which(names(allData) == "scenario"),
               init.vars = which(names(allData) %in% c("meanSev", "meanFreq", "meanSevB", "meanPatchS")),
               HVmethod = "box", no.runs = 3,
               freeBW = FALSE, bwHV1 = bwHV, bwHV2 = bwHV,
               do.scale = TRUE,
               noAxes = 4, outputs.dir = HVoutputPath,
               file.suffix = paste0("fireHVs_freeBW_", veg, "_rep", r),
               saveOrdi = TRUE, plotOrdi = TRUE, plotHV = TRUE)
})
future:::ClusterRegistry("stop")


## HYPERVOLUMES ACROSS THE LANDSCAPE - only montane belt ----------------
## BANDWITH ESTIMATES ----
amc::.gc()
parallel_wrapper <- function(ncores, summaryFireAttributes, byVars, bw.outputPath) {
  if (!dir.exists(bw.outputPath)) dir.create(bw.outputPath)

  if (.Platform$OS.type == "windows") {
    plan("multisession", workers = ncores)
  } else {
    plan("multicore", workers = ncores)
  }
  bw_estimates <- future_lapply(X = split(summaryFireAttributes, by = byVars),
                                FUN = function(DT, bw.outputPath) {
                                  r <- unique(DT$rep)
                                  message(paste("Calculating PCAs and estimating BWs for: rep", r))
                                  file.suffix <- paste0("fireHVs_freeBW_landscape_rep", r)

                                  init.vars <- grep("mean", names(DT))

                                  DT <- ToolsCB:::.scaleVars(DT, init.vars)

                                  out <- estimateBW_wrapper(as.data.frame(DT),
                                                            init.vars = init.vars,
                                                            HVidvar = which(names(DT) == "scenario"),
                                                            noAxes = 4,
                                                            ordination = "PCA",
                                                            file.suffix = file.suffix,
                                                            outputs.dir = bw.outputPath)
                                  return(out)
                                },
                                bw.outputPath = bw.outputPath)
  future:::ClusterRegistry("stop")
  return(bw_estimates)
}
bw_estimates <- Cache(parallel_wrapper,
                      ncores = 10,
                      summaryFireAttributes = summaryFireAttributes,
                      byVars = c("rep"),
                      bw.outputPath = bw.outputPath,
                      cacheRepo = simPaths$cachePath,
                      userTags = c("bw_estimates", "hypervolumes", "landscape"),
                      omitArgs = c("userTags", "ncores", "bw.outputPath"))

bw_estimates <- do.call(rbind.data.frame, bw_estimates)
bw_estimates$rep <- sub("\\..*", "", sub("\\.PC.*", "", row.names(bw_estimates)))
bw_estimates <- as.data.table(bw_estimates)
saveRDS(bw_estimates, file.path(bw.outputPath, "BW_estimates_landscape.rds"))

summaryBW <- bw_estimates[, list(SilvermanMean = mean(c(SilvBW_HV1, SilvBW_HV2)),
                                 StDevMean = mean(c(stdev_HV1, stdev_HV2)),
                                 SilvermanMax = max(c(SilvBW_HV1, SilvBW_HV2)),
                                 StDevMax = max(c(stdev_HV1, stdev_HV2))),
                          by = "PC"]

saveRDS(summaryBW, file.path(bw.outputPath, "BW_MeanMax_landscape.rds"))

## fix bandwidth to max of estimated BW
summaryBW <- readRDS(file.path(bw.outputPath, "BW_MeanMax_landscape.rds"))
bwHV <- summaryBW[["StDevMax"]]

amc::.gc()
if (.Platform$OS.type == "windows") {
  plan("multisession", workers = ncores)
} else {
  plan("multicore", workers = ncores)
}
future_lapply(split(summaryFireAttributes, by = c("rep"))[9], FUN = function(allData, bwVal1, bwVal2) {
  r <- unique(allData$rep)
  hypervolumes(HVdata1 = as.data.frame(allData[scenario == "noPM"]),
               HVdata2 = as.data.frame(allData[scenario == "PM"]),
               HVidvar = which(names(allData) == "scenario"),
               init.vars = which(names(allData) %in% c("meanSev", "meanFreq", "meanSevB", "meanPatchS")),
               HVmethod = "box", no.runs = 3,
               freeBW = FALSE, bwHV1 = bwHV, bwHV2 = bwHV,
               do.scale = TRUE,
               noAxes = 4, outputs.dir = HVoutputPath,
               file.suffix = paste0("fireHVs_freeBW_landscape_rep", r),
               saveOrdi = TRUE, plotOrdi = TRUE, plotHV = TRUE)
})
future:::ClusterRegistry("stop")




