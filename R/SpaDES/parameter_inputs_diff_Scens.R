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
simPaths <- list(cachePath = file.path("R/SpaDES/cache", simDirName, "postSimAnalyses")
                 , modulePath = file.path("R/SpaDES/m")
                 , inputPath = file.path("R/SpaDES/inputs")
                 , outputPath = file.path("R/SpaDES/outputs", simDirName))

preSimListnoPM <- loadSimList(file.path(simPaths$outputPath, "noPM", "LIM_simInit_noPM.qs"))
preSimListPM <- loadSimList(file.path(simPaths$outputPath, "PM", "LIM_simInit_PM.qs"))

paramsnoPM <- params(preSimList)$Biomass_regeneration
paramsPM <- params(preSimList)$Biomass_regenerationPM

inputsnoPM <- depends(preSimListnoPM)@dependencies$Biomass_regeneration@inputObjects
inputsPM <- depends(preSimListPM)@dependencies$Biomass_regenerationPM@inputObjects

diffInputsnoPM <- setdiff(inputsnoPM$objectName, inputsPM$objectName)
diffInputsPM <- setdiff(inputsPM$objectName, inputsnoPM$objectName)
