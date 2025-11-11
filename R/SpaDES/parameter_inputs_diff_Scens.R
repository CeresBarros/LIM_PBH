## --------------------------------------------------
##  FINDING DIFFERENCES IN INPUTS AND PARAMETERS
##  BETWEEN SCENARIOS
## --------------------------------------------------

library(SpaDES)

options("reproducible.useNewDigestAlgorithm" = 2)
options("spades.moduleCodeChecks" = FALSE)
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
} else if (grepl("for-cast", Sys.info()["nodename"])) {
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
    grepl("45eafed436c8", Sys.info()["nodename"])) {
  data.table::setDTthreads(25)
  options(bitmapType="cairo")
}

preSimListnoPM <- loadSimList(file.path(simPaths$outputPath, "noPM", "LIM_simInit_noPM.qs"))
preSimListPM <- loadSimList(file.path(simPaths$outputPath, "PM", "LIM_simInit_PM.qs"))

paramsnoPM <- params(preSimList)$Biomass_regeneration
paramsPM <- params(preSimList)$Biomass_regenerationPM

inputsnoPM <- depends(preSimListnoPM)@dependencies$Biomass_regeneration@inputObjects
inputsPM <- depends(preSimListPM)@dependencies$Biomass_regenerationPM@inputObjects

diffInputsnoPM <- setdiff(inputsnoPM$objectName, inputsPM$objectName)
diffInputsPM <- setdiff(inputsPM$objectName, inputsnoPM$objectName)
