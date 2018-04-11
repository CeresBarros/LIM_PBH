## ---------------------------------
## DAVID ANDISON FIRE DATA
##
## Dataprep for Exploratory analysis
## ---------------------------------


## requires
library(SpaDES); library(sf); library(ggplot2); library(data.table); library(dplyr)
source("R/R_tools/Useful_functions.R")

## define paths
setPaths(cachePath = file.path("R/SpaDES/cache"),
         modulePath = file.path("R/SpaDES/m"),
         inputPath = file.path("R/SpaDES/inputs"),
         outputPath = file.path("R/SpaDES/outputs"))


## LOAD DATA ---------------------------------------

## POST FIRE DATA
files = c("albertafires1_postfire", "albertafires2_postfire", "saskatchewanfires_postfire")
folder = "data/fires_Dave/Projected_renamed"

for(x in files) {
  eval(parse(text = paste0(
    x, " <- st_read(file.path(folder", ", paste0('", x,"', '.shp')))"
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

names(albertafires2_postfire)[grep("ELEMENT", names(albertafires2_postfire))] = 
  names(albertafires1_postfire)[grep("SURVIVAL", names(albertafires1_postfire))] = 
  names(saskatchewanfires_postfire)[grep("CLASS", names(saskatchewanfires_postfire))] = "SEV_CLAS"

## PRE FIRE DATA
files = c("albertafires1_prefire", "albertafires2_prefire", "saskatchewanfires_prefire")
folder = "data/fires_Dave/Projected_renamed"

for(x in files) {
  eval(parse(text = paste0(
    x, " <- st_read(file.path(folder", ", paste0('", x,"', '.shp')))"
  )))
}

## WATER DATA
# files = c("water-abta", "water-sask")
# folder = "data/fires_Dave/Projected_renamed"
# 
# for(x in files) {
#   eval(parse(text = paste0(
#     sub("-", "_", x), " <- st_read(file.path(folder", ", paste0('", x,"', '.shp')))"
#   )))
# }


## DEFINE FIRE EVENTS ----------------------------------------

# firesABSK <- Cache(loadBindSpatialObjs, 
#                    files = c("albertafires1_postfire", "albertafires2_postfire", "saskatchewanfires_postfire"),
#                    folder = "data/fires_Dave/Projected_renamed", 
#                    cacheRepo = getPaths()$cachePath, userTags = "fireData")

## Use Alberta 1 post fire data only for now, as severity classes on other datasets and not yet comparable.
rm(albertafires2_postfire, albertafires2_prefire, saskatchewanfires_postfire, saskatchewanfires_prefire)
gc(reset = TRUE)
AB1_fireEvents <- Cache(defineFireEvents, 
                    sf.obj = albertafires1_postfire, fireNAMES = "FIRE_NAME",
                    # fireVARS = c("FIRE_ID", "FIRE_YEAR", "SEV_CLAS"),   ## this makes the output object huge
                    buff.dist = 200L, 
                    PLOT = FALSE, SAVE = FALSE, outputDIR = "analyses/FireEvents", 
                    fileNAME = "Andison_AB1_fireEvents", overwrite = TRUE,
                    cacheRepo = getPaths()$cachePath, userTags = "dataTreat_fireEvents",
                    omitArgs = c("PLOT", "SAVE", "outputDIR", "fileNAME", "overwrite"))

## add severity (doing it in defineFireEvents seems to produce an overly large polygon)
AB1_distPatchSev <- Cache(st_intersection, 
                          x = AB1_fireEvents[AB1_fireEvents$PatchType == "disturbedPatch",], 
                          y = albertafires1_postfire[, "SEV_CLAS"], userTags = "dataTreat_fireEvents")
AB1_fireEvents[, setdiff(names(AB1_distPatchSev), names(AB1_fireEvents))] <- NA ## make NA columns for binding
AB1_fireEventsSev <- rbind(AB1_fireEvents[AB1_fireEvents$PatchType != "disturbedPatch",], AB1_distPatchSev)

## JOIN WATER, VEGETATION AND FIRE EVENTS --------------------
AB1_vegFireEvents <- Cache(st_intersection, 
                           x = albertafires1_prefire, y = AB1_fireEventsSev,
                           userTags = "dataTreat_fireEvents_wVeg")

## save - not working, R thinks this is sfc instead of sf.
# st_write(st_as_sf(AB1_vegFireEvents), 
#          dsn = "analyses/FireEvents/Andison_AB1_fireEventsVegetation.shp",
#          delete_layer = TRUE)  ##  not working
# raster::shapefile(as_Spatial(st_as_sf(AB1_vegFireEvents)),
#                   filename = "analyses/FireEvents/Andison_AB1_fireEventsVegetation.shp",
#                   overwrite = TRUE)

## not sure what to do about water yet... 
## perhaps just remove these areas with st difference before intersecting with veg?
# AB1_watVegFireInters <- Cache(st_intersection,
#                               x = AB1_vegFireEvents, y = water_abta[, "FEATURE_TY"],
#                               userTags = "dataTreat_fireEvents_wVeg_water")
# names(AB1_watVegFireEvents)[which(names(AB1_watVegFireEvents) == "FEATURE_TY")] <- "WATER_TY"

## extract dataframe only
AB1_VegFireEvents.dt <- as.data.table(AB1_vegFireEvents[,, drop = TRUE])   ## drops geometries

## remove columns with NAs only
NAcols <- sapply(AB1_VegFireEvents.dt, FUN = function(x) {
  return(any(sum(is.na(x)) == length(x)))
})

AB1_VegFireEvents.dt <- AB1_VegFireEvents.dt[, !NAcols, with = FALSE]


## CALCULATE RELATIVE OCCURRENCES OF VEGETATION ATTRIBUTES PER PATCH/FIRE EVENT
sums <- AB1_VegFireEvents.dt %>% group_by(FIRE_NAME, PatchType, SP1) %>%
summarise(Count = n())
totals <- AB1_VegFireEvents.dt %>% group_by(FIRE_NAME, PatchType) %>%
  summarise(Count = n())

lapply(unique(sums$PatchType))




