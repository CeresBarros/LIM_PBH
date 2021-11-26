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

## MERGE MIXED CONIFER AND DOUGLAS-FIR/DRY-CONIFER STANDS? OR JUST DOUGLAS-FIR/DRY-CONIFER STANDS?
mergeDMCPSME <- FALSE  ## merge DMCPSME PSME dryPSME
mergePSME <- TRUE ## merge PSME dryPSME
options("LandR.assertions" = FALSE)
if (mergeDMCPSME) {
  HVoutputPath <- file.path(simPaths$outputPath, "hypervolumes/mergeDMCPSME")
  allPixelCohortDataMnt[vegTypeCN %in% c("DMCPSME", "PSME", "dryPSME"), vegTypeCN := "DMCPSME"]
}

if (mergePSME) {
  HVoutputPath <- file.path(simPaths$outputPath, "hypervolumes/mergePSME")
  allPixelCohortDataMnt[vegTypeCN %in% c("PSME", "dryPSME"), vegTypeCN := "PSME"]
}

## FIRE ATTRIBUTES HYPERVOLUMES -----------
## Fire properties (fire patch size in pixels, fire frequency, fire severity as biomass loss)
## FIRE DATA SUMMARY FOR HVs -----------------------
source("R/R_tools/prepFireData4HVs.R")

## Global pyrodiversity PCA ----------
## a large PCA on the pooled dataset is needed to ensure that
## hypervolume sizes can be compared across repetitions and forest types.
summaryFireAttributes[, dummyHVid := paste(scenario, rep, vegTypeCN, sep = "_")]
cols <- c("meanFreq", "meanSevB", "meanPatchS", "dummyHVid")
firePCA <- summaryFireAttributes[, ..cols] %>%
  ToolsCB:::.scaleVars(., init.vars = c(1:3)) %>%
  Cache(HVordination,
        datatable = .,
        HVidvar = 4,
        noAxes = 3,
        plotOrdi = FALSE,
        saveOrdi = TRUE,
        saveOrdiSumm = TRUE,
        file.suffix = "fireHVs_FULLPCA",
        outputs.dir = file.path(simPaths$outputPath, "hypervolumes"),
        cacheRepo = simPaths$cachePath,
        omitArgs = c("plotOrdi", "saveOrdi"),
        userTags = c("hypervolumes", "pyrodivPCA"))

fireHVdata <- as.data.table(firePCA$HVpoints)
fireHVdata <- cbind(fireHVdata, summaryFireAttributes[, .(scenario, rep, pixelIndex, vegTypeCN)])

## Hypervolumes by vegetation type ----------
## only montane belt
if (mergeDMCPSME) {
  fireHVdata <- fireHVdata[vegTypeCN == "DMCPSME"]
}

if (mergePSME) {
  fireHVdata <- fireHVdata[vegTypeCN == "PSME"]
}

doAll <- FALSE
lapply(split(fireHVdata, by = c("rep", "vegTypeCN"), drop = TRUE),
       FUN = function(allData, HVoutputPath, doAll) {
         r <- unique(allData$rep)
         veg <- unique(allData$vegTypeCN)
         file.suffix <- paste0("fireHVs_", veg, "_rep", r)
         cols <- c("PC1", "PC2", "PC3")
         no.runs <- 3

         skip <- FALSE
         if (!doAll) {
           ## check if all HV intersections were computed already
           computedHVpairs <- grep(paste0(pattern = file.suffix, "_Intersection.*(1|2|3).rds$"),
                                   list.files(HVoutputPath), value = TRUE)
           if (length(computedHVpairs) == no.runs) {
             skip <- TRUE
           }
         }

         if (skip) {
           message(crayon::blue(file.suffix, ": all", no.runs, "done already.",
                                "'doAll' is", doAll,"... Skipping"))
         } else {
           fireHVWrapper(allData, cols, file.suffix,
                         # noAxes = 3,
                         ordination = "none",
                         HVmethod = "svm",
                         no.runs = no.runs,
                         svm.gamma = 0.01,
                         outputs.dir = HVoutputPath,
                         do.scale = FALSE,
                         # do.scale = TRUE,
                         # saveOrdi = TRUE,
                         # plotOrdi = TRUE,
                         plotHV = TRUE,
                         verbose = FALSE,
                         addNoise = TRUE)
         }
       }, HVoutputPath = HVoutputPath, doAll = doAll)   ## if doAll == FALSE, only missing HV intersection pairs will be computed


## Hypervolumes across the landscape ----------------
## only montane belt
doAll <- FALSE
lapply(split(fireHVdata, by = c("rep"), drop = TRUE),
       FUN = function(allData, HVoutputPath, doAll) {
         r <- unique(allData$rep)
         file.suffix <- paste0("fireHVs_landscape_rep", r)
         cols <- c("PC1", "PC2", "PC3")
         no.runs <- 3

         skip <- FALSE
         if (!doAll) {
           ## check if all HV intersections were computed already
           computedHVpairs <- grep(paste0(pattern = file.suffix, "_Intersection.*(1|2|3).rds$"),
                                   list.files(HVoutputPath), value = TRUE)
           if (length(computedHVpairs) == no.runs) {
             skip <- TRUE
           }
         }

         if (skip) {
           message(crayon::blue(file.suffix, ": all", no.runs, "done already.",
                                "'doAll' is", doAll,"... Skipping"))
         } else {
           fireHVWrapper(allData, cols, file.suffix,
                         # noAxes = 3,
                         ordination = "none",
                         HVmethod = "svm",
                         no.runs = no.runs,
                         svm.gamma = 0.01,
                         outputs.dir = HVoutputPath,
                         do.scale = FALSE,
                         # do.scale = TRUE,
                         # saveOrdi = TRUE,
                         # plotOrdi = TRUE,
                         plotHV = TRUE,
                         verbose = FALSE,
                         addNoise = TRUE)
         }
       }, HVoutputPath = HVoutputPath, doAll = doAll)


## VEGETATION ATTRIBUTES HYPERVOLUMES -----------
## Hypervolumes by vegetation type --------------
## only montane belt

## Use the first fire year to identify the pixels we want to follow in time
## we follow the same pixels that were used to make fire attributes HVs, only now
## we select the start year and end years of the simulation
## the join shouldn't actually change anything because we already subset the pixels with veg
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


  ## checks at landscape scale:
  pixelIndices <- unique(summaryFireAttributes[,.(scenario, rep, pixelIndex)])
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

## prep data for hypervolumes
## calculate stand-age as the mean biomass-weighted age
vegDataForHVs[, `:=`(meanStandAge = mean(as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                      sum((B/100), na.rm = TRUE))),
                     sdStandAge = sd(as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                  sum((B/100), na.rm = TRUE)))),
              by = c("scenario", "rep", "year", "pixelIndex", "vegTypeCN")]
vegDataForHVs[is.na(meanStandAge), meanStandAge := 0]  ## NAs come from stands with 0 B and 0 age
vegDataForHVs[is.na(sdStandAge), sdStandAge := 0]  ## NAs come from stands with 0 B and 0 age or just one value of standAge

## calculate relative species B (across speciesCohorts)
vegDataForHVs[, standB := sum(B, na.rm = TRUE), by =  c("scenario", "rep", "year", "pixelIndex", "vegTypeCN")]
vegDataForHVs[, relB := sum(B) / standB, by = c("scenario", "rep", "year", "pixelIndex", "vegTypeCN", "speciesCode")]
vegDataForHVs[relB == "NaN", relB := 0]

## expand data
cols <- c("meanStandAge", "sdStandAge", "relB", "speciesCode", "scenario", "rep", "year", "pixelIndex", "vegTypeCN")   ## keep rep for wrapper.
vegDataForHVs <- unique(vegDataForHVs[, ..cols])
vegDataForHVs <- dcast.data.table(vegDataForHVs, as.formula("... ~ speciesCode"),
                                  value.var = "relB")

## Global biodiversity PCA ----------
## a large PCA on the pooled dataset is needed to ensure that
## hypervolume sizes can be compared across repetitions and forest types.
vegDataForHVs[, dummyHVid := paste(scenario, rep, year, vegTypeCN, sep = "_")]
cols <- setdiff(names(vegDataForHVs),
                c("scenario", "rep", "year", "pixelIndex", "vegTypeCN"))
vegPCA <- vegDataForHVs[, ..cols] %>%
  ToolsCB:::.scaleVars(., init.vars = c(1:3)) %>%
  Cache(HVordination,
        datatable = .,
        HVidvar = 9,
        noAxes = 4,
        plotOrdi = TRUE,
        saveOrdiSumm = TRUE,
        saveOrdi = TRUE,   ## save actual PCA
        file.suffix = "vegHVs_FULLPCA",
        outputs.dir = file.path(simPaths$outputPath, "hypervolumes"),
        cacheRepo = simPaths$cachePath,
        omitArgs = c("plotOrdi", "saveOrdi"),
        userTags = c("hypervolumes", "biodivPCA"))
vegHVdata <- as.data.table(vegPCA$HVpoints)
vegHVdata <- cbind(vegHVdata, vegDataForHVs[, .(scenario, rep, year, pixelIndex, vegTypeCN)])

if (mergeDMCPSME) {  ## only do the the merged vegTypes
  vegHVdata <- vegHVdata[vegTypeCN == "DMCPSME" & year %in% c(start(preSimList), end(preSimList))]
} else {
  if (mergePSME) {  ## only do the the merged vegTypes
    vegHVdata <- vegHVdata[vegTypeCN == "PSME" & year %in% c(start(preSimList), end(preSimList))]
  } else {
    vegHVdata <- vegHVdata[year %in% c(start(preSimList), end(preSimList))]
  }
}


## subset first 4 axes:
cols <- c(grep("PC(1|2|3|4)", names(vegHVdata), value = TRUE),
          grep("^PC", names(vegHVdata), value = TRUE, invert = TRUE))
vegHVdata <- vegHVdata[, ..cols]

## HV comparisons per year, between scenarios --------------
## note that splitting by veg type has to be done on the first
## year as vegTypes can change (cannot use first fire year, because cohortData will have
## been impacted by fire already). Splitting is done by rep only
## as vegTypeCN/pixelIndex combos for the first year have to be
## identical between scenarios (tested above)
## gaussian HVs were extremely slow
pixelIndexList <- split(vegHVdata[year == start(preSimList), .(rep, vegTypeCN, pixelIndex)],
                        by = c("rep", "vegTypeCN"), drop = TRUE)
doAll <- FALSE
lapply(pixelIndexList,
       FUN = function(pixelIndexDT, vegHVdata, HVoutputPath, doAll) {
         r <- unique(pixelIndexDT$rep)
         veg <- unique(pixelIndexDT$vegTypeCN)

         ## filter data to appropriate pixels, note that vegType may change in the second year
         allData <- vegHVdata[pixelIndexDT[, .(rep, pixelIndex)], on = .(rep, pixelIndex)]

         ## now split by year to calculate and compare hypervolumes between
         ## scenarios for each year
         lapply(split(allData, by = "year"),
                FUN = function(allData, HVoutputPath, r, veg, doAll) {
                  yr <- unique(allData$year)
                  file.suffix <- paste0("vegHVs_", veg, "_yr", yr, "_rep", r)
                  IDcols <- c("scenario", "rep", "pixelIndex", "year", "vegTypeCN")

                  no.runs <- 3
                  skip <- FALSE
                  if (!doAll) {
                    ## check if all HV intersections were computed already
                    computedHVpairs <- grep(paste0(pattern = file.suffix, "_Intersection.*(1|2|3).rds$"),
                                            list.files(HVoutputPath), value = TRUE)
                    if (length(computedHVpairs) == no.runs) {
                      skip <- TRUE
                    }
                  }

                  if (skip) {
                    message(crayon::blue(file.suffix, ": all", no.runs, "done already.",
                                         "'doAll' is", doAll,"... Skipping"))
                  } else {
                    print(file.suffix)
                    vegHVWrapper(allData,
                                 IDcols,
                                 HVIDcol = "scenario",
                                 file.suffix,
                                 # noAxes = 4,
                                 ordination = "none",
                                 HVmethod = "svm",
                                 no.runs = no.runs,
                                 svm.gamma = 0.01,
                                 outputs.dir = HVoutputPath,
                                 do.scale = FALSE,
                                 addNoise = TRUE,
                                 # do.scale = TRUE,
                                 # saveOrdi = TRUE,
                                 # plotOrdi = TRUE,
                                 plotHV = TRUE,
                                 verbose = FALSE)
                  }
                }, HVoutputPath = HVoutputPath, r = r, veg = veg, doAll = doAll)
       }, vegHVdata = vegHVdata, HVoutputPath = HVoutputPath, doAll = doAll)


## HV comparisons per scenario, between years --------------
## gaussian HVs were extremely slow
doAll <- FALSE
lapply(pixelIndexList,
       FUN = function(pixelIndexDT, vegHVdata, HVoutputPath, doAll) {
         r <- unique(pixelIndexDT$rep)
         veg <- unique(pixelIndexDT$vegTypeCN)

         ## filter data to appropriate pixels, note that vegType may change in the second year
         allData <- vegHVdata[pixelIndexDT[, .(rep, pixelIndex)], on = .(rep, pixelIndex)]

         ## now split by scenario to calculate and compare hypervolumes between
         ## scenarios for each scenario
         lapply(split(allData, by = "scenario"),
                FUN = function(allData, HVoutputPath, r, veg, doAll) {
                  scen <- unique(allData$scenario)
                  file.suffix <- paste0("vegHVs_", veg, "_", scen, "_rep", r)
                  IDcols <- c("year", "rep", "pixelIndex", "scenario", "vegTypeCN")

                  no.runs <- 3
                  skip <- FALSE
                  if (!doAll) {
                    ## check if all HV intersections were computed already
                    computedHVpairs <- grep(paste0(pattern = file.suffix, "_Intersection.*(1|2|3).rds$"),
                                            list.files(HVoutputPath), value = TRUE)
                    if (length(computedHVpairs) == no.runs) {
                      skip <- TRUE
                    }
                  }

                  if (skip) {
                    message(crayon::blue(file.suffix, ": all", no.runs, "done already.",
                                         "'doAll' is", doAll,"... Skipping"))
                  } else {
                    print(file.suffix)
                    vegHVWrapper(allData,
                                 IDcols,
                                 HVIDcol = "year",
                                 file.suffix,
                                 # noAxes = 4,
                                 ordination = "none",
                                 HVmethod = "svm",
                                 no.runs = no.runs,
                                 svm.gamma = 0.01,
                                 outputs.dir = HVoutputPath,
                                 do.scale = FALSE,
                                 addNoise = TRUE,
                                 # do.scale = TRUE,
                                 # saveOrdi = TRUE,
                                 # plotOrdi = TRUE,
                                 plotHV = TRUE,
                                 verbose = FALSE)
                  }
                }, HVoutputPath = HVoutputPath, r = r, veg = veg, doAll = doAll)
       }, vegHVdata = vegHVdata, HVoutputPath = HVoutputPath, doAll = doAll)


## Hypervolumes across the landscape ----------------
## only montane belt
## Now we follow all pixels, so there is no need to subset pixels by veg type in
## in the first year

## HV comparisons per year, between scenarios --------------
## split by year and rep to calculate and compare hypervolumes between
## scenarios for each year
doAll <- FALSE
lapply(split(vegHVdata, by = c("rep", "year")),
       FUN = function(allData, HVoutputPath, doAll) {
         r <- unique(allData$rep)
         yr <- unique(allData$year)
         file.suffix <- paste0("vegHVs_landscape", "_yr", yr, "_rep", r)
         IDcols <- c("scenario", "rep", "pixelIndex", "year", "vegTypeCN")
         no.runs <- 3
         skip <- FALSE
         if (!doAll) {
           ## check if all HV intersections were computed already
           computedHVpairs <- grep(paste0(pattern = file.suffix, "_Intersection.*(1|2|3).rds$"),
                                   list.files(HVoutputPath), value = TRUE)
           if (length(computedHVpairs) == no.runs) {
             skip <- TRUE
           }
         }

         if (skip) {
           message(crayon::blue(file.suffix, ": all", no.runs, "done already.",
                                "'doAll' is", doAll,"... Skipping"))
         } else {
           print(file.suffix)

           vegHVWrapper(allData,
                        IDcols,
                        HVIDcol = "scenario",
                        file.suffix,
                        # noAxes = 4,
                        ordination = "none",
                        HVmethod = "svm",
                        no.runs = no.runs,
                        svm.gamma = 0.01,
                        outputs.dir = HVoutputPath,
                        do.scale = FALSE,
                        addNoise = TRUE,
                        # do.scale = TRUE,
                        # saveOrdi = TRUE,
                        # plotOrdi = TRUE,
                        plotHV = TRUE,
                        verbose = FALSE)
         }
       }, HVoutputPath = HVoutputPath, doAll = doAll)


## HV comparisons per scenario, between years --------------
## gaussian HVs were extremely slow
## now split by scenario and rep to calculate and compare hypervolumes between
## years for each scenario
doAll <- FALSE
lapply(split(vegHVdata, by = c("rep","scenario")),
       FUN = function(allData, HVoutputPath, doAll) {
         r <- unique(allData$rep)
         scen <- unique(allData$scenario)
         file.suffix <- paste0("vegHVs_landscape_", scen, "_rep", r)
         IDcols <- c("year", "rep", "pixelIndex", "scenario", "vegTypeCN")
         no.runs <- 3
         skip <- FALSE
         if (!doAll) {
           ## check if all HV intersections were computed already
           computedHVpairs <- grep(paste0(pattern = file.suffix, "_Intersection.*(1|2|3).rds$"),
                                   list.files(HVoutputPath), value = TRUE)
           if (length(computedHVpairs) == no.runs) {
             skip <- TRUE
           }
         }

         if (skip) {
           message(crayon::blue(file.suffix, ": all", no.runs, "done already.",
                                "'doAll' is", doAll,"... Skipping"))
         } else {
           print(file.suffix)

           vegHVWrapper(allData,
                        IDcols,
                        HVIDcol = "year",
                        file.suffix,
                        # noAxes = 4,
                        ordination = "none",
                        HVmethod = "svm",
                        no.runs = no.runs,
                        svm.gamma = 0.01,
                        outputs.dir = HVoutputPath,
                        do.scale = FALSE,
                        addNoise = TRUE,
                        # do.scale = TRUE,
                        # saveOrdi = TRUE,
                        # plotOrdi = TRUE,
                        plotHV = TRUE,
                        verbose = FALSE)
         }
       }, HVoutputPath = HVoutputPath, doAll = doAll)
