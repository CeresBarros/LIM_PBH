## --------------------------------------------------
##  HYPERVOLUMES OF PYRODIVERSITY AND BIODIVERSITY
##  ANALYSING PCAs
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

## LOAD DATA (RESULTS)  ---------------------
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

## PREP FIRE AND VEG DATA -------------------
  HVoutputPath <- file.path(simPaths$outputPath, "hypervolumes/mergePSME")
  allPixelCohortDataMnt[vegTypeCN %in% c("PSME", "dryPSME"), vegTypeCN := "PSME"]
}

## PREP FIRE AND VEG DATA -----------
## Fire properties (fire patch size in pixels, fire frequency, fire severity as biomass loss)
source("R/R_tools/prepFireData4HVs.R")
source("R/R_tools/prepVegData4HVs.R")



## LOAD PCAs ---------------------------
fireHVPCA <- grep("fireHVs", list.files(HVoutputPath, "OrdinationObj", full.names = TRUE), value = TRUE)
fireHVPCA <- readRDS(fireHVPCA)
vegHVPCA <- grep("vegHVs", list.files(HVoutputPath, "OrdinationObj", full.names = TRUE), value = TRUE)
vegHVPCA <- readRDS(vegHVPCA)
