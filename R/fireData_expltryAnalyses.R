## ---------------------------------
## DAVID ANDISON FIRE DATA
##
## Exploratory analysis
## ---------------------------------

## clean workspace
rm(list=ls()); gc(reset = TRUE)

## requires
library(SpaDES)
source("R/R_tools/Useful_functions.R")

## define paths
setPaths(cachePath = file.path("R/SpaDES/cache"),
         modulePath = file.path("R/SpaDES/m"),
         inputPath = file.path("R/SpaDES/inputs"),
         outputPath = file.path("R/SpaDES/outputs"))


## LOAD DATA ---------------------------------------
files = c("albertafires1_postfire", "albertafires2_postfire", "saskatchewanfires_postfire")
folder = "data/fires_Dave/Projected_renamed"

for(x in files) {
  eval(parse(text = paste0(
    x, " <- sf::st_read(file.path(folder", ", paste0('", x,"', '.shp')))"
    )))
}

head(albertafires1_postfire)
head(albertafires2_postfire)
head(saskatchewanfires_postfire)

## Saskatchewan attributes table has some repeated columns, check if they are duplicated before deleting
if(all(saskatchewanfires_postfire$FIRE_NAME_ == saskatchewanfires_postfire$FIRE_NAME)) {
  saskatchewanfires_postfire$FIRE_NAME_ <- NULL
}
if(all(saskatchewanfires_postfire$AREA_ == saskatchewanfires_postfire$AREA)) {
  saskatchewanfires_postfire$AREA_ <- NULL
}
if(all(saskatchewanfires_postfire$PERIMETER_ == saskatchewanfires_postfire$PERIMETER)) {
  saskatchewanfires_postfire$PERIMETER_ <- NULL
}

## Uniformize column names
names(albertafires2_postfire)[grep("YEAR", names(albertafires2_postfire))] = 
  names(saskatchewanfires_postfire)[grep("EAR", names(saskatchewanfires_postfire))] = 
  grep("YEAR", names(albertafires1_postfire), value = TRUE)

names(albertafires2_postfire)[grep("NAME", names(albertafires2_postfire))] = 
  names(saskatchewanfires_postfire)[grep("NAME", names(saskatchewanfires_postfire))] = 
  grep("NAME", names(albertafires1_postfire), value = TRUE)

names(albertafires2_postfire)[grep("FIRE_CODE", names(albertafires2_postfire))] = 
  names(albertafires1_postfire)[grep("FIRE_NUM", names(albertafires1_postfire))] = 
  grep("FIRE_ID", names(saskatchewanfires_postfire), value = TRUE)


# firesABSK <- Cache(loadBindSpatialObjs, 
#                    files = c("albertafires1_postfire", "albertafires2_postfire", "saskatchewanfires_postfire"),
#                    folder = "data/fires_Dave/Projected_renamed", 
#                    cacheRepo = getPaths()$cachePath, userTags = "fireData")

## Use Alberta 1 post fire data only for now, as severity classes on other datasets and not yet comparable.
AB1_fireEvents <- Cache(defineFireEvents, 
                    sf.obj = albertafires1_postfire, fireNAMES = "FIRE_NAME", buff.dist = 200L, 
                    PLOT = FALSE, SAVE = FALSE, outputDIR = "analyses/FireEvents", 
                    fileNAME = "Andison_AB1_fireEvents", overwrite = FALSE,
                    cacheRepo = getPaths()$cachePath, userTags = "dataTreat_fireEvents")
