## --------------------------------------------------
##  HYPERVOLUMES OF PYRODIVERSITY AND BIODIVERSITY
## --------------------------------------------------

library(SpaDES)
library(ToolsCB)
library(data.table)

source("R/R_tools/convertToCNVegType.R")
source("R/R_tools/Useful_functions.R")
source("R/R_tools/hypervolumesHelpers.R")
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
HVoutputPath <- file.path(simPaths$outputPath, "hypervolumes")
# bw.outputPath <- file.path(HVoutputPath, "bwTest")

## LOAD DATA (RESULTS)  -------------------
yearSubset <- c(seq(2011, 2111, 5), 2111)
source("R/SpaDES/6_resultsDataPrep.R")

## MERGE DOUGLAS-FIR/DRY-CONIFER STANDS
mergeDMCPSME <- TRUE

if (mergeDMCPSME) {
  HVoutputPath <- file.path(simPaths$outputPath, "hypervolumes/mergeDMCPSME")
  allPixelCohortDataMnt[vegTypeCN %in% c("DMCPSME", "PSME", "dryPSME"), vegTypeCN := "DMCPSME"]
}

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

## Hypervolumes by vegetation type ----------
## only montane belt
lapply(split(summaryFireAttributes, by = c("rep", "vegTypeCN")), FUN = function(allData, HVoutputPath) {
  r <- unique(allData$rep)
  veg <- unique(allData$vegTypeCN)
  file.suffix <- paste0("fireHVs_", veg, "_rep", r)

  noAxes <- 3
  cols <- c("meanFreq", "meanSevB", "meanPatchS")

  fireHVWrapper(allData, cols, file.suffix,
                noAxes = noAxes,
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


## Hypervolumes across the landscape ----------------
## only montane belt
lapply(split(summaryFireAttributes, by = c("rep")), FUN = function(allData, HVoutputPath) {
  r <- unique(allData$rep)
  file.suffix <- paste0("fireHVs_landscape_rep", r)

  noAxes <- 3
  cols <- c("meanFreq", "meanSevB", "meanPatchS")

  fireHVWrapper(allData, cols, file.suffix,
                noAxes = noAxes,
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


## VEGETATION ATTRIBUTES HYPERVOLUMES -----------
## Hypervolumes by vegetation type --------------
## only montane belt

## Use the first fire year to identify the pixels we want to follow in time
## we follow the same pixels that were used to make fire attributes HVs, only now
## we select the start year and end years of the simulation
## the join shouldn't actually change anything because we already subset the pixels with veg
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

  temp <- vegDataForHVs[, length(unique(pixelIndex)), by = .(scenario, rep, year)]
  if (unique(temp$V1) > 1) {
    stop("There should be the same number of pixels every year.")
  }

  temp <- split(vegDataForHVs[, .(pixelIndex, scenario, rep, year)],
                by = c("scenario", "rep", "year"), keep.by = FALSE)
  temp <- lapply(temp, FUN = function(x) unique(x[["pixelIndex"]]))
  test <- lapply(1:length(temp), function(n) setdiff(temp[[n]], unlist(temp[-n])))

  test <- sapply(test, length)
  if (any(test))
    stop("Different pixelIndex between scenario/rep/year combinations")

  temp <- split(vegDataForHVs[year == start(preSimList), .(vegTypeCN, pixelIndex, scenario, rep)],
                by = c("scenario", "rep"), keep.by = FALSE)
  temp <- lapply(temp, FUN = function(x) setkey(x, vegTypeCN, pixelIndex))
  temppix <- lapply(temp, FUN = function(x) x[["pixelIndex"]])
  tempveg <- lapply(temp, FUN = function(x) x[["vegTypeCN"]])

  test <- lapply(1:length(temppix), function(n) setdiff(temppix[[n]], unlist(temppix[-n])))
  test2 <- lapply(1:length(tempveg), function(n) setdiff(tempveg[[n]], unlist(tempveg[-n])))
  test <- sapply(test, length)
  test2 <- sapply(test2, length)

  if (any(test) | any(test2))
    stop("Difference pixelIndex/vegTypeCN combinations between scenario/reps in the first year")
}

## HV comparisons per year, between scenarios --------------
## note that splitting by veg type has to be done on the first
## year as vegTypes can change. Splitting is done by rep only
## as vegTypeCN/pixelIndex combos for the first year have to be
## identical between scenarios (tested above)
## gaussian HVs were extremely slow

pixelIndexList <- split(vegDataForHVs[year == start(preSimList), .(rep, vegTypeCN, pixelIndex)],
                        by = c("rep", "vegTypeCN"))

lapply(pixelIndexList, FUN = function(pixelIndexDT, vegDataForHVs, HVoutputPath) {
  r <- unique(pixelIndexDT$rep)
  veg <- unique(pixelIndexDT$vegTypeCN)

  ## filter data to appropriate pixels, note that vegType may change in the second year
  allData <- vegDataForHVs[pixelIndexDT[, .(rep, pixelIndex)], on = .(rep, pixelIndex)]

  ## now split by year to calculate and compare hypervolumes between
  ## scenarios for each year
  lapply(split(allData, by = "year"), FUN = function(allData, HVoutputPath, r, veg) {
    yr <- unique(allData$year)
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
  }, HVoutputPath = HVoutputPath, r = r, veg = veg)
}, vegDataForHVs = vegDataForHVs, HVoutputPath = HVoutputPath)


## HV comparisons per scenario, between years --------------
## gaussian HVs were extremely slow
lapply(pixelIndexList, FUN = function(pixelIndexDT, vegDataForHVs, HVoutputPath) {
  r <- unique(pixelIndexDT$rep)
  veg <- unique(pixelIndexDT$vegTypeCN)

  ## filter data to appropriate pixels, note that vegType may change in the second year
  allData <- vegDataForHVs[pixelIndexDT[, .(rep, pixelIndex)], on = .(rep, pixelIndex)]

  ## now split by scenario to calculate and compare hypervolumes between
  ## scenarios for each scenario
  lapply(split(allData, by = "scenario"), FUN = function(allData, HVoutputPath, r, veg) {
    scen <- unique(allData$scenario)
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
  }, HVoutputPath = HVoutputPath, r = r, veg = veg)
}, vegDataForHVs = vegDataForHVs, HVoutputPath = HVoutputPath)



## Hypervolumes across the landscape ----------------
## Hypervolumes by vegetation type --------------
## only montane belt

## Use the first fire year to identify the pixels we want to follow in time
## we follow the same pixels that were used to make fire attributes HVs, only now
## we select the start year and end years of the simulation
## the join shouldn't acutally change anything because we already subset the pixels with veg
## in the montane belt (regardless of fire)
vegDataForHVs <- allPixelCohortDataMnt[year %in% c(start(preSimList), end(preSimList))]

if (getOption("LandR.assertions")) {
  pixelIndices <- unique(summaryFireAttributes[,.(scenario, rep, pixelIndex)])
  temp <- vegDataForHVs[pixelIndices, on = .(scenario, rep, pixelIndex), nomatch = 0]
  setkey(temp, scenario, rep, pixelIndex)
  setkey(vegDataForHVs, scenario, rep, pixelIndex)

  if (isFALSE(identical(temp, vegDataForHVs))) {
    stop("Something is wrong. summaryFireAttributes should have the same pixelIndex/scenario/rep\n",
         "Combinations as allPixelCohortDataMnt")
  }

  temp <- vegDataForHVs[, length(unique(pixelIndex)), by = .(scenario, rep, year)]
  if (unique(temp$V1) > 1) {
    stop("There should be the same number of pixels every year.")
  }

  temp <- split(vegDataForHVs[, .(pixelIndex, scenario, rep, year)],
                by = c("scenario", "rep", "year"), keep.by = FALSE)
  temp <- lapply(temp, FUN = function(x) unique(x[["pixelIndex"]]))
  test <- lapply(1:length(temp), function(n) setdiff(temp[[n]], unlist(temp[-n])))

  test <- sapply(test, length)
  if (any(test))
    stop("Different pixelIndex between scenario/rep/year combinations")
}

## HV comparisons per year, between scenarios --------------
## split by year and rep to calculate and compare hypervolumes between
## scenarios for each year
lapply(split(vegDataForHVs, by = c("rep", "year")), FUN = function(allData, HVoutputPath) {
  r <- unique(allData$rep)
  yr <- unique(allData$year)
  file.suffix <- paste0("vegHVs_landscape", "_yr", yr, "_rep", r)
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


## HV comparisons per scenario, between years --------------
## gaussian HVs were extremely slow
## now split by scenario and rep to calculate and compare hypervolumes between
## years for each scenario
lapply(split(vegDataForHVs, by = c("rep","scenario")), FUN = function(allData, HVoutputPath) {
  r <- unique(allData$rep)
  scen <- unique(allData$scenario)
  file.suffix <- paste0("vegHVs_landscape_", scen, "_rep", r)
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
