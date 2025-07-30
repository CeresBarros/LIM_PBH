## --------------------------------------------------
##  HYPERVOLUMES OF PYRODIVERSITY AND BIODIVERSITY
## --------------------------------------------------

if (!exists("pkgDir")) {
  rver <- paste0(version$major, ".", strsplit(version$minor, "[.]")[[1]][1])
  pkgDir <- file.path(
    if (Sys.info()[["sysname"]] == "Linux" && rver == "4.1") "packages_docker" else "packages",
    version$platform,
    rver
  )

  if (!dir.exists(pkgDir)) {
    dir.create(pkgDir, recursive = TRUE)
  }
}

.libPaths(pkgDir)

library(SpaDES.core)
library(ToolsCB)
library(data.table)

source("R/R_tools/convertToCNVegType.R")
source("R/R_tools/Useful_functions.R")
source("R/R_tools/hypervolumesHelpers.R")
options("reproducible.useNewDigestAlgorithm" = 2)
options("reproducible.useCache" = TRUE)
options("reproducible.destinationPath" = normPath("R/SpaDES/inputs"))
options("reproducible.useGDAL" = FALSE)

## general paths
# simDirName <- "jun2021Runs"
simDirName <- "mar2022Runs"
if (Sys.info()["nodename"] == "W-VIC-A127584") {
  simPaths <- list(cachePath = file.path("F:", basename(getwd()), "R/SpaDES/cache", simDirName, "postSimAnalyses")
                   , modulePath = file.path("R/SpaDES/m")
                   , inputPath = file.path("F:", basename(getwd()), "R/SpaDES/inputs", simDirName)
                   , outputPath = file.path("F:", basename(getwd()), "R/SpaDES/outputs", simDirName))
} else if (grepl("for-cast|ad4d65210e84", Sys.info()["nodename"])) {
  ## settings for for-cast and coco machines
  simPaths <- list(cachePath = file.path("R/SpaDES/cache", simDirName, "postSimAnalyses")
                   , modulePath = file.path("R/SpaDES/m")
                   , inputPath = file.path("R/SpaDES/inputs")
                   , outputPath = file.path("R/SpaDES/outputs", simDirName)
                   , rasterPath = file.path("/mnt/scratch/cbarros", basename(getwd()), "R/SpaDES/scratch/raster")
                   , scratchPath = file.path("/mnt/scratch/cbarros", basename(getwd()), "R/SpaDES/scratch"))
} else {
  simPaths <- list(cachePath = file.path("R/SpaDES/cache", simDirName, "postSimAnalyses")
                   , modulePath = file.path("R/SpaDES/m")
                   , inputPath = file.path("R/SpaDES/inputs")
                   , outputPath = file.path("R/SpaDES/outputs", simDirName))
}

if (grepl("for-cast", Sys.info()["nodename"]) ||
    grepl("4458e1a42ddc", Sys.info()["nodename"])) {
  ## settings for for-cast and coco machines
  data.table::setDTthreads(5)
  options(bitmapType="cairo")
}

## path to figure folder and cache folder
figOutputPath <- file.path(simPaths$outputPath, "figuresAnalysis")
HVoutputPath <- file.path(simPaths$outputPath, "hypervolumes")
# bw.outputPath <- file.path(HVoutputPath, "bwTest")

## LOAD DATA (RESULTS)  -------------------
yearSubset <- unique(as.integer(c(seq(3511, 4011, 5), 4011)))
runPrepResultsModule <- TRUE
source("R/SpaDES/simResultsDataPrep.R")

rm(allPixelCohortData)
gc()

## MERGE MIXED CONIFER AND DOUGLAS-FIR/DRY-CONIFER STANDS? OR JUST DOUGLAS-FIR/DRY-CONIFER STANDS?
mergeDMCPSME <- FALSE  ## merge DMCPSME PSME dryPSME
mergePSME <- TRUE ## merge PSME dryPSME
doMergedOnly <- FALSE ## should HVs only be calculated for the merged veg types?
options("LandR.assertions" = TRUE)
if (mergeDMCPSME) {
  HVoutputPath <- file.path(simPaths$outputPath, "hypervolumes/mergeDMCPSME")
  allPixelCohortDataMnt[vegTypeCN %in% c("DMCPSME", "PSME", "dryPSME"), vegTypeCN := "DMCPSME"]
}

if (mergePSME) {
  HVoutputPath <- file.path(simPaths$outputPath, "hypervolumes/mergePSME")
  allPixelCohortDataMnt[vegTypeCN %in% c("PSME", "dryPSME"), vegTypeCN := "PSME"]
}

dir.create(HVoutputPath, recursive = TRUE)

## FIRE ATTRIBUTES HYPERVOLUMES -----------
## Fire properties (fire patch size in pixels, fire frequency, fire severity as biomass loss)
## FIRE DATA SUMMARY FOR HVs -----------------------
opts <- options("LandR.assertions", FALSE)
source("R/R_tools/prepFireData4HVs.R")
options(opts)

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
        plotOrdi = TRUE,
        saveOrdi = TRUE,
        saveOrdiSumm = TRUE,
        file.suffix = "fireHVs_FULLPCA",
        outputs.dir = file.path(simPaths$outputPath, "hypervolumes"),
        cacheRepo = simPaths$cachePath,
        omitArgs = c("plotOrdi", "saveOrdi", "saveOrdiSumm"),
        userTags = c("hypervolumes", "pyrodivPCA"))

fireHVdata <- as.data.table(firePCA$HVpoints)
fireHVdata <- cbind(fireHVdata, summaryFireAttributes[, .(scenario, rep, pixelIndex, vegTypeCN)])

## Hypervolumes by vegetation type ----------
## only montane belt
if (mergeDMCPSME & doMergedOnly) {
  fireHVdata <- fireHVdata[vegTypeCN == "DMCPSME"]
}

if (mergePSME & doMergedOnly) {
  fireHVdata <- fireHVdata[vegTypeCN == "PSME"]
}

doAll <- FALSE ## if doAll == FALSE, only missing HV intersection pairs will be computed
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
doAll <- FALSE ## if doAll == FALSE, only missing HV intersection pairs will be computed
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

## for each rep, we will randomly draw 5 years from each 100yrs window (the same years are used across scenarios).
## all years enter the same hypervolume, so that it integrates the last 500yrs of simulation
yearSamples <- sample5SimYears(allPixelCohortDataMnt[, .(year, rep)])   ## seed ensures same years are drawn

useFirstLastYear <- FALSE ## use only first/last year of yearSubset, or all years?
source("R/R_tools/prepVegData4HVs.R")

rm(allPixelBurnData, allPixelCohortDataMnt)
gc(reset = TRUE)

## Global biodiversity PCA ----------
## a large PCA on the pooled dataset is needed to ensure that
## hypervolume sizes can be compared across repetitions and forest types.

## with 2000 yrs sims the sampled years from the last 500 years are integrated in a single
## HV (per scenario, rep, vegtype)
vegDataForHVs[, dummyHVid := ifelse(useFirstLastYear, paste(scenario, rep, year, vegTypeCN, sep = "_"),
                                    paste(scenario, rep, vegTypeCN, sep = "_"))]

cols <- setdiff(names(vegDataForHVs),
                c("scenario", "rep", "year", "pixelIndex", "vegTypeCN"))
vegPCA <- vegDataForHVs[, ..cols] %>%
  ToolsCB:::.scaleVars(., init.vars = c(1:3)) %>%
  Cache(HVordination,
        datatable = .,
        HVidvar = 9,
        noAxes = 4,
        plotOrdi = FALSE,
        saveOrdiSumm = TRUE,
        saveOrdi = TRUE,   ## save actual PCA
        file.suffix = "vegHVs_FULLPCA",
        outputs.dir = file.path(simPaths$outputPath, "hypervolumes"),
        cacheRepo = simPaths$cachePath,
        omitArgs = c("plotOrdi", "saveOrdi", "saveOrdiSumm"),
        userTags = c("hypervolumes", "biodivPCA"))
vegHVdata <- as.data.table(vegPCA$HVpoints)
vegHVdata <- cbind(vegHVdata, vegDataForHVs[, .(scenario, rep, year, pixelIndex, vegTypeCN)])

## subset to first/last years
if (useFirstLastYear) {
  vegHVdata <- vegHVdata[year %in% c(min(yearSubset), max(yearSubset))]
}

if (mergeDMCPSME & doMergedOnly) {  ## only do the the merged vegTypes
  vegHVdata <- vegHVdata[vegTypeCN == "DMCPSME"]
}
if (mergePSME & doMergedOnly) {  ## only do the the merged vegTypes
  vegHVdata <- vegHVdata[vegTypeCN == "PSME"]
}


## subset first 4 axes:
cols <- c(grep("PC(1|2|3|4)", names(vegHVdata), value = TRUE),
          grep("^PC", names(vegHVdata), value = TRUE, invert = TRUE))
vegHVdata <- vegHVdata[, ..cols]

## Hypervolumes by vegetation type --------------
## only montane belt

## HV comparisons per year, between scenarios --------------
## note that splitting by veg type has to be done on the last
## year as vegTypes can change (cannot use first fire year, because cohortData will have
## been impacted by fire already). Splitting is done by rep only
## as vegTypeCN/pixelIndex combos for the first year have to be
## identical between scenarios (tested above)
## gaussian HVs were extremely slow

if (useFirstLastYear) {
  pixelIndexList <- split(unique(vegHVdata[year == max(yearSubset), .(rep, vegTypeCN, pixelIndex)]),
                          by = c("rep", "vegTypeCN"), drop = TRUE)
  ## check that there is only one vegType at the start year across scenarios
  ## not relevant when using last 500yrs of a 2000yrs simulation
  test <- sapply(pixelIndexList,
                 FUN = function(pixelIndexDT, vegHVdata) {
                   ## filter data to appropriate pixels, note that vegType may change in the second year
                   allData <- vegHVdata[pixelIndexDT[, .(rep, pixelIndex)], on = .(rep, pixelIndex)]
                   unique(allData$vegTypeCN) == unique(pixelIndexDT$vegTypeCN)

                 }, vegHVdata = vegHVdata[year == min(yearSubset)])
  if (any(!test)) stop("different veg types at start year fround across scenarios.")
} else {
  ## the last year may no longer be 4011
  pixelIndexList <- vegHVdata[, list(year = max(year)), .(rep)]
  pixelIndexList <- vegHVdata[pixelIndexList, on = .(rep, year)]
  pixelIndexList <- as.data.frame(pixelIndexList)  ## for some reason the next line is crashing R with a memory access issue
  pixelIndexList <- unique(pixelIndexList[, c("rep", "vegTypeCN", "pixelIndex")])
  pixelIndexList <- as.data.table(pixelIndexList)
  pixelIndexList <- split(pixelIndexList, by = c("rep", "vegTypeCN"), drop = TRUE)
}

doAll <- FALSE
lapply(pixelIndexList,
       FUN = function(pixelIndexDT, vegHVdata, HVoutputPath, doAll, useFirstLastYear) {
         r <- unique(pixelIndexDT$rep)
         veg <- unique(pixelIndexDT$vegTypeCN)

         ## filter data to appropriate pixels, note that vegType may change in the second year
         allData <- vegHVdata[pixelIndexDT[, .(rep, pixelIndex)], on = .(rep, pixelIndex)]

         if (!useFirstLastYear) {
           allData[, year := NA_integer_] ## don't need year anymore, all will be integrated
         }

         # if necessary split by year to calculate and compare hypervolumes between
         # scenarios for each year
         # not relevant when using last 500yrs of a 2000yrs simulation
         lapply(split(allData, by = "year"),
                FUN = function(allData, HVoutputPath, r, veg, doAll) {
                  yr <- unique(allData$year)
                  if (!is.na(yr)) {
                    file.suffix <- paste0("vegHVs_", veg, "_yr", yr, "_rep", r)
                  } else {
                    file.suffix <- paste0("vegHVs_", veg, "_rep", r)
                  }

                  IDcols <- c("scenario", "rep", "pixelIndex")

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
                                 HVcols = c("PC1", "PC2", "PC3", "PC4"),
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
       },
       vegHVdata = vegHVdata, HVoutputPath = HVoutputPath, doAll = doAll, useFirstLastYear)


## HV comparisons per scenario, between years --------------
## gaussian HVs were extremely slow

## not relevant when using last 500yrs of a 2000yrs simulation
if (useFirstLastYear) {
  doAll <- FALSE
  lapply(pixelIndexList,
         FUN = function(pixelIndexDT, vegHVdata, HVoutputPath, doAll) {
           r <- unique(pixelIndexDT$rep)
           veg <- unique(pixelIndexDT$vegTypeCN)

           ## filter data to appropriate pixels, note that vegType may change in the second year
           allData <- vegHVdata[pixelIndexDT[, .(rep, pixelIndex)], on = .(rep, pixelIndex)]

           ## now split by scenario to calculate and compare hypervolumes between
           ## years for each scenario
           lapply(split(allData, by = "scenario"),
                  FUN = function(allData, HVoutputPath, r, veg, doAll) {
                    scen <- unique(allData$scenario)
                    file.suffix <- paste0("vegHVs_", veg, "_", scen, "_rep", r)
                    IDcols <- c("year", "rep", "pixelIndex")

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
                                   HVcols = c("PC1", "PC2", "PC3", "PC4"),
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
}

## Hypervolumes across the landscape ----------------
## only montane belt
## Now we follow all pixels, so there is no need to subset pixels by veg type in
## in the first year

## HV comparisons per year, between scenarios --------------
## split by year and rep to calculate and compare hypervolumes between
## scenarios for each year

## year split not relevant when using last 500yrs of a 2000yrs simulation

doAll <- FALSE
lapply(ifelse(useFirstLastYear,
              split(vegHVdata, by = c("rep", "year")),
              split(vegHVdata, by = c("rep"))),
       FUN = function(allData, HVoutputPath, doAll, useFirstLastYear) {
         r <- unique(allData$rep)

         if (!useFirstLastYear) {
           allData[, year := NA_integer_]
         }

         yr <- unique(allData$year)
         if (!is.na(yr)) {
           file.suffix <- paste0("vegHVs_landscape", "_yr", yr, "_rep", r)
         } else {
           file.suffix <- paste0("vegHVs_landscape", "_rep", r)
         }

         IDcols <- c("scenario", "rep", "pixelIndex")
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
                        HVcols = c("PC1", "PC2", "PC3", "PC4"),
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
       }, HVoutputPath = HVoutputPath, doAll = doAll, useFirstLastYear)


## HV comparisons per scenario, between years --------------
## gaussian HVs were extremely slow
## now split by scenario and rep to calculate and compare hypervolumes between
## years for each scenario

## not relevant when using last 500yrs of a 2000yrs simulation

if (useFirstLastYear) {
  doAll <- FALSE
  lapply(split(vegHVdata, by = c("rep","scenario")),
         FUN = function(allData, HVoutputPath, doAll) {
           r <- unique(allData$rep)
           scen <- unique(allData$scenario)
           file.suffix <- paste0("vegHVs_landscape_", scen, "_rep", r)
           IDcols <- c("year", "rep", "pixelIndex")
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
                          HVcols = c("PC1", "PC2", "PC3", "PC4"),
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

}
